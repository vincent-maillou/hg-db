// partitioner.cpp - MT-KaHyPar integration implementation
#include "hypergraph_reorder/partitioner.hpp"

#include <mtkahypar.h>

#include <algorithm>
#include <iostream>
#include <thread>
#include <unordered_set>

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
}  // namespace

HypergraphPartitioner::HypergraphPartitioner(const Options& opts)
    : opts_(opts), context_(nullptr) {
  init_context();
}

HypergraphPartitioner::~HypergraphPartitioner() { cleanup_context(); }

void HypergraphPartitioner::init_context() {
  // Initialize MT-KaHyPar thread pool (once globally)
  if (!s_initialized) {
    size_t num_threads =
        opts_.num_threads > 0
            ? static_cast<size_t>(opts_.num_threads)
            : std::max(1u, std::thread::hardware_concurrency());
    mt_kahypar_initialize(num_threads, true);
    s_initialized = true;
  }

  // Create context from preset
  context_ = mt_kahypar_context_from_preset(convert_preset(opts_.preset));

  // Configure partitioning parameters
  mt_kahypar_set_partitioning_parameters(
      static_cast<mt_kahypar_context_t*>(context_),
      static_cast<mt_kahypar_partition_id_t>(opts_.n_parts), opts_.imbalance,
      KM1  // Use connectivity (km1) metric
  );

  // Set seed if specified (global setting)
  if (opts_.seed >= 0) {
    mt_kahypar_set_seed(static_cast<size_t>(opts_.seed));
  }

  // Suppress output if requested
  if (opts_.suppress_output) {
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
    std::cout << "Warning: Empty hypergraph, creating trivial partition"
              << std::endl;

    HypergraphPartition result;
    result.n_parts = opts_.n_parts;
    result.node_partition.resize(hg.n_nodes());
    for (index_t i = 0; i < hg.n_nodes(); ++i) {
      result.node_partition[i] = i % opts_.n_parts;
    }
    result.objective = 0;
    result.part_sizes.resize(opts_.n_parts, 0);
    for (auto p : result.node_partition) {
      result.part_sizes[p]++;
    }
    return result;
  }

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

  // Convert edge weights (can be null if all weights are 1)
  const mt_kahypar_hyperedge_weight_t* edge_weights = nullptr;
  std::vector<mt_kahypar_hyperedge_weight_t> edge_weights_conv;
  if (!hg.edge_weights().empty()) {
    edge_weights_conv.resize(hg.edge_weights().size());
    for (size_t i = 0; i < hg.edge_weights().size(); ++i) {
      edge_weights_conv[i] =
          static_cast<mt_kahypar_hyperedge_weight_t>(hg.edge_weights()[i]);
    }
    edge_weights = edge_weights_conv.data();
  }

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

  // Create MT-KaHyPar hypergraph
  mt_kahypar_error_t error{nullptr, 0, SUCCESS};
  mt_kahypar_hypergraph_t mt_hg = mt_kahypar_create_hypergraph(
      static_cast<mt_kahypar_context_t*>(context_), num_nodes, num_nets,
      hyperedge_indices.data(), hyperedges.data(), edge_weights, node_weights,
      &error);

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

  // Get objective value (connectivity metric)
  mt_kahypar_hyperedge_weight_t objective = mt_kahypar_km1(partitioned_hg);

  // Convert to our result format
  HypergraphPartition result;
  result.n_parts = opts_.n_parts;
  result.node_partition.resize(hg.n_nodes());
  result.part_sizes.resize(opts_.n_parts, 0);

  for (index_t i = 0; i < hg.n_nodes(); ++i) {
    result.node_partition[i] = partition[i];
    result.part_sizes[partition[i]]++;
  }

  result.objective = objective;

  std::cout << "Partition objective (km1): " << objective << std::endl;
  std::cout << "Part sizes: ";
  for (auto sz : result.part_sizes) {
    std::cout << sz << " ";
  }
  std::cout << std::endl;

  // Cleanup MT-KaHyPar structures
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  mt_kahypar_free_hypergraph(mt_hg);

  return result;
}

VertexPartition HypergraphPartitioner::create_vertex_partition(
    const HypergraphPartition& cnh_partition, const CliqueCover& cover,
    index_t n_vertices) {
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

  std::cout << "Vertex partition completed:" << std::endl;
  std::cout << "  Part sizes: ";
  for (const auto& part : result.parts) {
    std::cout << part.size() << " ";
  }
  std::cout << std::endl;
  std::cout << "  Separator size: " << result.separator.size() << " ("
            << (100.0 * result.separator.size() / n_vertices) << "%)"
            << std::endl;

  return result;
}

}  // namespace hypergraph_reorder
