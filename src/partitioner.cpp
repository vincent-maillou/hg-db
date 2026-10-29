// partitioner.cpp - MT-KaHyPar integration implementation

// Must precede every include: sched_getaffinity / CPU_COUNT are GNU
// extensions and this build compiles with -std=c++20 (CMAKE_CXX_EXTENSIONS
// OFF), so _GNU_SOURCE is not predefined by the toolchain.
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include "hypergraph_reorder/partitioner.hpp"

#include <mtkahypar.h>

#include <algorithm>
#include <iostream>
#include <thread>
#include <unordered_set>

#if defined(__linux__)
#include <sched.h>  // sched_getaffinity / CPU_COUNT (affinity-aware autodetect)
#endif

namespace hypergraph_reorder {

namespace {
// Convert our preset enum to MT-KaHyPar preset type
mt_kahypar_preset_type_t convert_preset(MtKahyparPreset preset) {
  switch (preset) {
    case MtKahyparPreset::QUALITY:
      return QUALITY;
    case MtKahyparPreset::DETERMINISTIC:
      return DETERMINISTIC;
    case MtKahyparPreset::LARGE_K:
      return LARGE_K;
    case MtKahyparPreset::DEFAULT:
    default:
      return DEFAULT;
  }
}

// Thread pool initialization flag
bool s_initialized = false;
size_t s_initialized_num_threads = 0;

// CPUs available to this process (affinity-aware). Falls back to
// hardware_concurrency() off-Linux or on failure (floored at 1).
// Respects SLURM/taskset/cgroup cpuset allocations. Spec WP1 D2.
size_t available_cpus() {
#if defined(__linux__)
  cpu_set_t set;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    const int count = CPU_COUNT(&set);
    if (count > 0) return static_cast<size_t>(count);
  }
#endif
  const unsigned hw = std::thread::hardware_concurrency();
  return hw > 0 ? static_cast<size_t>(hw) : 1;
}
}  // namespace

size_t mtkahypar_pool_threads() {
  return s_initialized ? s_initialized_num_threads : 0;
}

HypergraphPartitioner::HypergraphPartitioner(const PartitionOptions& opts,
                                             int num_threads,
                                             bool suppress_output)
    : opts_(opts),
      num_threads_(num_threads),
      suppress_output_(suppress_output),
      context_(nullptr) {
  init_context();
}

HypergraphPartitioner::~HypergraphPartitioner() { cleanup_context(); }

void HypergraphPartitioner::init_context() {
  // Initialize MT-KaHyPar thread pool (once globally)
  if (!s_initialized) {
    size_t num_threads = num_threads_ > 0
                             ? static_cast<size_t>(num_threads_)
                             : available_cpus();
    mt_kahypar_initialize(num_threads, true);
    s_initialized = true;
    s_initialized_num_threads = num_threads;
  } else if (num_threads_ > 0 &&
             static_cast<size_t>(num_threads_) != s_initialized_num_threads) {
    // MT-KaHyPar's thread pool initializes once per process; a second
    // reorderer asking for a different thread count keeps the first
    // setting. Warn once (spec §4.6).
    std::cerr << "Warning: MT-KaHyPar thread pool already initialized with "
              << s_initialized_num_threads
              << " threads; ignoring requested " << num_threads_ << "."
              << std::endl;
  }

  // Create context from preset
  context_ = mt_kahypar_context_from_preset(convert_preset(opts_.preset));

  // Map the configured objective to the MT-KaHyPar enum. Strict: an
  // invalid objective must fail loudly, not silently fall back (an
  // objective that differs from the requested one would corrupt
  // scientific A/B comparisons).
  mt_kahypar_objective_t objective = KM1;
  switch (opts_.objective) {
    case MtKahyparObjective::kKm1:
      objective = KM1;
      break;
    case MtKahyparObjective::kCut:
      objective = CUT;
      break;
    default:
      throw HypergraphReorderError("Invalid MT-KaHyPar partitioning objective");
  }

  // Configure partitioning parameters
  mt_kahypar_set_partitioning_parameters(
      static_cast<mt_kahypar_context_t*>(context_),
      static_cast<mt_kahypar_partition_id_t>(opts_.n_parts), opts_.imbalance,
      objective
  );

  // Set seed if specified (global setting)
  if (opts_.seed >= 0) {
    mt_kahypar_set_seed(static_cast<size_t>(opts_.seed));
  }

  // Suppress output if requested
  if (suppress_output_) {
    mt_kahypar_error_t error{nullptr, 0, SUCCESS};
    mt_kahypar_set_context_parameter(
        static_cast<mt_kahypar_context_t*>(context_), VERBOSE, "0", &error);
    if (error.status != SUCCESS) {
      mt_kahypar_free_error_content(&error);
    }
  }
}

void HypergraphPartitioner::cleanup_context() {
  if (context_) {
    mt_kahypar_free_context(static_cast<mt_kahypar_context_t*>(context_));
    context_ = nullptr;
  }
}

HypergraphPartition HypergraphPartitioner::partition(const Hypergraph& hg) {
  if (hg.n_nets() == 0) {
    if (!suppress_output_)
      std::cout << "Warning: Empty hypergraph, creating trivial partition"
                << std::endl;

    HypergraphPartition result;
    result.n_parts = opts_.n_parts;
    result.node_partition.resize(hg.n_nodes());
    for (index_t i = 0; i < hg.n_nodes(); ++i) {
      result.node_partition[i] = i % opts_.n_parts;
    }
    result.km1_objective = 0;
    result.cut_objective = 0;
    result.effective_imbalance = opts_.imbalance;  // no relax on this path
    result.part_sizes.resize(opts_.n_parts, 0);
    for (auto p : result.node_partition) {
      result.part_sizes[p]++;
    }
    return result;
  }

  if (!suppress_output_)
    std::cout << "Partitioning CNH: " << hg.n_nodes() << " nodes, " << hg.n_nets()
              << " nets into " << opts_.n_parts << " parts" << std::endl;

  // Prepare hypergraph data for MT-KaHyPar
  const mt_kahypar_hypernode_id_t num_nodes =
      static_cast<mt_kahypar_hypernode_id_t>(hg.n_nodes());
  const mt_kahypar_hyperedge_id_t num_nets =
      static_cast<mt_kahypar_hyperedge_id_t>(hg.n_nets());

  // Convert edge indices to size_t
  std::vector<size_t> hyperedge_indices(hg.edge_indices().begin(),
                                        hg.edge_indices().end());

  // Convert edges (pins) to the expected type
  std::vector<mt_kahypar_hyperedge_id_t> hyperedges(hg.edges().size());
  for (size_t i = 0; i < hg.edges().size(); ++i) {
    hyperedges[i] = static_cast<mt_kahypar_hyperedge_id_t>(hg.edges()[i]);
  }

  // Net weights are unit by construction (WP2 D5): null is passed
  // unconditionally — MT-KaHyPar's documented unit-weight path.

  // Convert node weights (can be null if all weights are 1)
  const mt_kahypar_hypernode_weight_t* node_weights = nullptr;
  std::vector<mt_kahypar_hypernode_weight_t> node_weights_conv;
  if (!hg.node_weights().empty()) {
    node_weights_conv.resize(hg.node_weights().size());
    for (size_t i = 0; i < hg.node_weights().size(); ++i) {
      node_weights_conv[i] =
          static_cast<mt_kahypar_hypernode_weight_t>(hg.node_weights()[i]);
    }
    node_weights = node_weights_conv.data();
  }

  // ---- Epsilon auto-relax (v4, spec D3) --------------------------------
  // MT-KaHyPar's balance constraint is sum_{i in block} w_i <= (1+eps)*W/k.
  // A single hypernode heavier than (1+eps)*W/k makes the instance
  // infeasible; rather than failing, relax eps to the minimum feasible
  // value and warn (unconditional stderr, like the thread-pool warning).
  // eps_eff is reported back through HypergraphPartition.
  const double W = static_cast<double>(hg.total_node_weight());
  double w_max = 0.0;
  for (const auto& w : hg.node_weights()) {
    w_max = std::max(w_max, static_cast<double>(w));
  }
  const double k = static_cast<double>(opts_.n_parts);
  double eps_eff = opts_.imbalance;
  if (W > 0.0) {
    const double min_feasible = k * w_max / W - 1.0;
    if (min_feasible > eps_eff) {
      eps_eff = min_feasible;
      std::cerr << "Warning: node weight " << static_cast<long long>(w_max)
                << " exceeds the block budget at imbalance "
                << opts_.imbalance << "; relaxed imbalance " << opts_.imbalance
                << " -> " << eps_eff << " (" << static_cast<long long>(W)
                << " total node weight, k=" << opts_.n_parts << ")"
                << std::endl;
    }
  }
  // Re-issue the partitioning parameters with the effective epsilon.
  // Always re-issued (not only on relax) so there is exactly one code
  // path that sets epsilon for the partition call; init_context() keeps
  // responsibility for the initial value.
  //
  // Caveat: the LARGE_K preset may apply different block-weight semantics
  // internally (perfectly balanced blocks); the relaxation remains a
  // necessary-condition guard in that case.
  {
    mt_kahypar_objective_t objective = KM1;
    switch (opts_.objective) {
      case MtKahyparObjective::kKm1:
        objective = KM1;
        break;
      case MtKahyparObjective::kCut:
        objective = CUT;
        break;
      default:
        throw HypergraphReorderError(
            "Invalid MT-KaHyPar partitioning objective");
    }
    mt_kahypar_set_partitioning_parameters(
        static_cast<mt_kahypar_context_t*>(context_),
        static_cast<mt_kahypar_partition_id_t>(opts_.n_parts), eps_eff,
        objective);
  }

  // Create MT-KaHyPar hypergraph
  mt_kahypar_error_t error{nullptr, 0, SUCCESS};
  mt_kahypar_hypergraph_t mt_hg = mt_kahypar_create_hypergraph(
      static_cast<mt_kahypar_context_t*>(context_), num_nodes, num_nets,
      hyperedge_indices.data(), hyperedges.data(),
      nullptr,  // net weights: unit by construction (WP2 D5)
      node_weights, &error);

  if (error.status != SUCCESS) {
    std::string err_msg = error.msg ? error.msg : "Unknown error";
    mt_kahypar_free_error_content(&error);
    throw HypergraphReorderError("Failed to create MT-KaHyPar hypergraph: " +
                                 err_msg);
  }

  // Partition the hypergraph
  mt_kahypar_partitioned_hypergraph_t partitioned_hg = mt_kahypar_partition(
      mt_hg, static_cast<mt_kahypar_context_t*>(context_), &error);

  if (error.status != SUCCESS) {
    std::string err_msg = error.msg ? error.msg : "Unknown error";
    mt_kahypar_free_error_content(&error);
    mt_kahypar_free_hypergraph(mt_hg);
    throw HypergraphReorderError("Failed to partition hypergraph: " + err_msg);
  }

  // Extract partition assignment
  std::vector<mt_kahypar_partition_id_t> partition(num_nodes);
  mt_kahypar_get_partition(partitioned_hg, partition.data());

  // Get both metrics of the achieved partition. They are reported for
  // every run regardless of the optimized objective so km1-runs and
  // cut-runs stay comparable.
  const mt_kahypar_hyperedge_weight_t km1 = mt_kahypar_km1(partitioned_hg);
  const mt_kahypar_hyperedge_weight_t cut = mt_kahypar_cut(partitioned_hg);

  // Convert to our result format
  HypergraphPartition result;
  result.n_parts = opts_.n_parts;
  result.node_partition.resize(hg.n_nodes());
  result.part_sizes.resize(opts_.n_parts, 0);

  for (index_t i = 0; i < hg.n_nodes(); ++i) {
    result.node_partition[i] = partition[i];
    result.part_sizes[partition[i]]++;
  }

  result.km1_objective = km1;
  result.cut_objective = cut;
  result.effective_imbalance = eps_eff;

  if (!suppress_output_) {
    std::cout << "Partition objectives: km1=" << km1 << " cut=" << cut
              << std::endl;
    std::cout << "Part sizes: ";
    for (auto sz : result.part_sizes) {
      std::cout << sz << " ";
    }
    std::cout << std::endl;
  }

  // Cleanup MT-KaHyPar structures
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  mt_kahypar_free_hypergraph(mt_hg);

  return result;
}

VertexPartition HypergraphPartitioner::create_vertex_partition(
    const HypergraphPartition& cnh_partition, const CliqueCover& cover,
    index_t n_vertices) {
  if (!suppress_output_)
    std::cout << "Creating vertex separator from CNH partition..." << std::endl;

  VertexPartition result;
  result.n_parts = cnh_partition.n_parts;
  result.parts.resize(result.n_parts);

  // For each vertex, determine if it belongs to a single part or separator
  for (index_t v = 0; v < n_vertices; ++v) {
    auto cliques = cover.get_vertex_cliques(v);

    if (cliques.empty()) {
      // Isolated vertex -> separator
      result.separator.push_back(v);
      continue;
    }

    // Get unique partitions this vertex belongs to
    std::unordered_set<index_t> vertex_parts;
    for (auto clique_id : cliques) {
      vertex_parts.insert(cnh_partition.node_partition[clique_id]);
    }

    if (vertex_parts.size() == 1) {
      // Vertex belongs to single partition
      index_t part_id = *vertex_parts.begin();
      result.parts[part_id].push_back(v);
    } else {
      // Vertex spans multiple partitions -> separator
      result.separator.push_back(v);
    }
  }

  // Sort for consistency
  for (auto& part : result.parts) {
    std::sort(part.begin(), part.end());
  }
  std::sort(result.separator.begin(), result.separator.end());

  if (!suppress_output_) {
    std::cout << "Vertex partition completed:" << std::endl;
    std::cout << "  Part sizes: ";
    for (const auto& part : result.parts) {
      std::cout << part.size() << " ";
    }
    std::cout << std::endl;
    std::cout << "  Separator size: " << result.separator.size() << " ("
              << (100.0 * result.separator.size() / n_vertices) << "%)"
              << std::endl;
  }

  return result;
}

}  // namespace hypergraph_reorder
