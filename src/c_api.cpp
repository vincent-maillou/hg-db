#include <cstring>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "hypergraph_reorder.h"
#include "hypergraph_reorder/reorderer.hpp"

/* Fallbacks when the build system does not inject the project version */
#ifndef HGR_VERSION_MAJOR
#define HGR_VERSION_MAJOR 5
#endif
#ifndef HGR_VERSION_MINOR
#define HGR_VERSION_MINOR 2
#endif
#ifndef HGR_VERSION_PATCH
#define HGR_VERSION_PATCH 0
#endif

using namespace hypergraph_reorder;

// Opaque structures
struct hgr_reorderer_s {
  std::unique_ptr<SymmetricDBReorderer> reorderer;
};

struct hgr_result_s {
  SymmetricDBReorderer::Result result;
  // Flattened clique-order histogram, owned by the result so the pointers
  // handed out by hgr_get_clique_order_counts stay valid until
  // hgr_result_free.
  std::vector<int64_t> clique_orders;
  std::vector<int64_t> clique_counts;
  // Flattened maximal-clique-enumeration histogram (v5, BK only), same
  // ownership contract for hgr_get_maximal_clique_enumeration_histogram.
  std::vector<int64_t> enum_clique_orders;
  std::vector<int64_t> enum_clique_counts;
};

// ---- Error channel (spec §4.6) ------------------------------------------
namespace {
thread_local std::string g_last_error;

void set_last_error_current_exception() {
  try {
    throw;  // re-throw the active exception to capture its message
  } catch (const std::exception& e) {
    g_last_error = e.what();
  } catch (...) {
    g_last_error = "unknown error";
  }
}
}  // namespace

// ---- CSR input validation (v5, D9) ---------------------------------------
// Content-only validation of the caller's CSR arrays before any copy or
// construction (fail fast through the error channel instead of undefined
// behavior downstream). Array *lengths* remain trusted by contract: the
// caller states n_rows/nnz and C cannot verify them — reading row_ptr[0],
// row_ptr[n_rows] and the first nnz col_idx entries is within that contract.
namespace {

// Returns an empty string on success, else a precise error message naming
// the offending position and value. Checks, in order (first failure wins):
//   1. n_rows >= 0, nnz >= 0
//   2. row_ptr[0] == 0
//   3. row_ptr non-decreasing over [0, n_rows]
//   4. row_ptr[n_rows] == nnz
//   5. every col_idx[k] in [0, n_rows)
std::string validate_csr_input(int64_t n_rows, int64_t nnz,
                               const int64_t* row_ptr, const int64_t* col_idx) {
  if (n_rows < 0) {
    return "hgr_reorder_csr: n_rows = " + std::to_string(n_rows) +
           " must be >= 0";
  }
  if (nnz < 0) {
    return "hgr_reorder_csr: nnz = " + std::to_string(nnz) + " must be >= 0";
  }
  if (row_ptr[0] != 0) {
    return "hgr_reorder_csr: row_ptr[0] = " + std::to_string(row_ptr[0]) +
           " must be 0";
  }
  for (int64_t i = 0; i < n_rows; ++i) {
    if (row_ptr[i + 1] < row_ptr[i]) {
      return "hgr_reorder_csr: row_ptr[" + std::to_string(i + 1) + "] = " +
             std::to_string(row_ptr[i + 1]) + " < row_ptr[" +
             std::to_string(i) + "] = " + std::to_string(row_ptr[i]) +
             " (must be non-decreasing)";
    }
  }
  if (row_ptr[n_rows] != nnz) {
    return "hgr_reorder_csr: row_ptr[" + std::to_string(n_rows) + "] = " +
           std::to_string(row_ptr[n_rows]) + " != nnz = " +
           std::to_string(nnz);
  }
  for (int64_t k = 0; k < nnz; ++k) {
    if (col_idx[k] < 0 || col_idx[k] >= n_rows) {
      return "hgr_reorder_csr: col_idx[" + std::to_string(k) + "] = " +
             std::to_string(col_idx[k]) + " out of range [0, " +
             std::to_string(n_rows) + ")";
    }
  }
  return "";
}

}  // namespace

extern "C" const char* hgr_get_last_error(void) {
  return g_last_error.c_str();  // never NULL; "" when clean
}

// Create reorderer
extern "C" hgr_reorderer_t* hgr_create(const hgr_options_t* opts) {
  try {
    SymmetricDBReorderer::Options cpp_opts;

    if (opts) {
      // ECC options
      switch (opts->ecc.method) {
        case HGR_ECC_EDGE_COVER:
          cpp_opts.ecc.method = EccMethod::kEdgeCover;
          break;
        case HGR_ECC_TRI_ENUM:
          cpp_opts.ecc.method = EccMethod::kTriEnum;
          break;
        case HGR_ECC_BK:
        default:
          cpp_opts.ecc.method = EccMethod::kBk;
          break;
      }
      cpp_opts.ecc.use_parallel = opts->ecc.use_parallel != 0;

      // Partition options
      cpp_opts.partition.n_parts = opts->partition.n_parts;
      cpp_opts.partition.imbalance = opts->partition.imbalance;
      cpp_opts.partition.seed = opts->partition.seed;
      switch (opts->partition.preset) {
        case HGR_PRESET_QUALITY:
          cpp_opts.partition.preset = MtKahyparPreset::QUALITY;
          break;
        case HGR_PRESET_DETERMINISTIC:
          cpp_opts.partition.preset = MtKahyparPreset::DETERMINISTIC;
          break;
        case HGR_PRESET_LARGE_K:
          cpp_opts.partition.preset = MtKahyparPreset::LARGE_K;
          break;
        case HGR_PRESET_DEFAULT:
        default:
          cpp_opts.partition.preset = MtKahyparPreset::DEFAULT;
          break;
      }
      // Objective: strict mapping — unlike the preset (lenient, historical),
      // an invalid objective must fail loudly. A silent fallback would
      // corrupt A/B benchmark comparisons without any signal.
      switch (opts->partition.objective) {
        case HGR_OBJECTIVE_KM1:
          cpp_opts.partition.objective = MtKahyparObjective::kKm1;
          break;
        case HGR_OBJECTIVE_CUT:
          cpp_opts.partition.objective = MtKahyparObjective::kCut;
          break;
        default:
          throw std::invalid_argument(
              "hgr_create: invalid partition.objective (expected "
              "HGR_OBJECTIVE_KM1 or HGR_OBJECTIVE_CUT)");
      }
      // Node weighting: strict mapping, same rationale as the objective —
      // a silent fallback would corrupt A/B benchmark comparisons.
      switch (opts->partition.node_weighting) {
        case HGR_NODE_WEIGHTING_UNIT:
          cpp_opts.partition.node_weighting = NodeWeighting::kUnit;
          break;
        case HGR_NODE_WEIGHTING_CLIQUE_SIZE:
          cpp_opts.partition.node_weighting = NodeWeighting::kCliqueSize;
          break;
        case HGR_NODE_WEIGHTING_VERTEX_SHARE:
          cpp_opts.partition.node_weighting = NodeWeighting::kVertexShare;
          break;
        case HGR_NODE_WEIGHTING_EDGE_FILL_MASS:
          cpp_opts.partition.node_weighting = NodeWeighting::kEdgeFillMass;
          break;
        default:
          throw std::invalid_argument(
              "hgr_create: invalid partition.node_weighting (expected "
              "HGR_NODE_WEIGHTING_UNIT, HGR_NODE_WEIGHTING_CLIQUE_SIZE, "
              "HGR_NODE_WEIGHTING_VERTEX_SHARE, or "
              "HGR_NODE_WEIGHTING_EDGE_FILL_MASS)");
      }

      // Runtime options
      cpp_opts.runtime.num_threads = opts->runtime.num_threads;
      cpp_opts.runtime.suppress_output = opts->runtime.suppress_output != 0;
      cpp_opts.ecc.num_threads = opts->runtime.num_threads;
      cpp_opts.ecc.suppress_output = opts->runtime.suppress_output != 0;
    }

    auto* handle = new hgr_reorderer_t();
    handle->reorderer = std::make_unique<SymmetricDBReorderer>(cpp_opts);
    return handle;
  } catch (...) {
    set_last_error_current_exception();
    return nullptr;
  }
}

extern "C" void hgr_free(hgr_reorderer_t* reorderer) { delete reorderer; }

extern "C" void hgr_version(int version[3]) {
  version[0] = HGR_VERSION_MAJOR;
  version[1] = HGR_VERSION_MINOR;
  version[2] = HGR_VERSION_PATCH;
}

extern "C" void hgr_default_options(hgr_options_t* opts) {
  if (!opts) return;

  // Single source of defaults: the C++ default-constructed Options.
  const SymmetricDBReorderer::Options cpp_opts{};

  opts->ecc.method = cpp_opts.ecc.method == EccMethod::kEdgeCover
                         ? HGR_ECC_EDGE_COVER
                     : cpp_opts.ecc.method == EccMethod::kTriEnum
                         ? HGR_ECC_TRI_ENUM
                         : HGR_ECC_BK;
  opts->ecc.use_parallel = cpp_opts.ecc.use_parallel ? 1 : 0;

  opts->partition.n_parts = static_cast<int64_t>(cpp_opts.partition.n_parts);
  opts->partition.imbalance = cpp_opts.partition.imbalance;
  switch (cpp_opts.partition.preset) {
    case MtKahyparPreset::QUALITY:
      opts->partition.preset = HGR_PRESET_QUALITY;
      break;
    case MtKahyparPreset::DETERMINISTIC:
      opts->partition.preset = HGR_PRESET_DETERMINISTIC;
      break;
    case MtKahyparPreset::LARGE_K:
      opts->partition.preset = HGR_PRESET_LARGE_K;
      break;
    case MtKahyparPreset::DEFAULT:
    default:
      opts->partition.preset = HGR_PRESET_DEFAULT;
      break;
  }
  opts->partition.objective =
      cpp_opts.partition.objective == MtKahyparObjective::kCut
          ? HGR_OBJECTIVE_CUT
          : HGR_OBJECTIVE_KM1;
  switch (cpp_opts.partition.node_weighting) {
    case NodeWeighting::kUnit:
      opts->partition.node_weighting = HGR_NODE_WEIGHTING_UNIT;
      break;
    case NodeWeighting::kCliqueSize:
      opts->partition.node_weighting = HGR_NODE_WEIGHTING_CLIQUE_SIZE;
      break;
    case NodeWeighting::kVertexShare:
      opts->partition.node_weighting = HGR_NODE_WEIGHTING_VERTEX_SHARE;
      break;
    case NodeWeighting::kEdgeFillMass:
    default:
      opts->partition.node_weighting = HGR_NODE_WEIGHTING_EDGE_FILL_MASS;
      break;
  }
  opts->partition.seed = cpp_opts.partition.seed;

  opts->runtime.num_threads = cpp_opts.runtime.num_threads;
  opts->runtime.suppress_output = cpp_opts.runtime.suppress_output ? 1 : 0;
}

extern "C" int hgr_reorder_csr(hgr_reorderer_t* reorderer, int64_t n_rows,
                               int64_t nnz, const int64_t* row_ptr,
                               const int64_t* col_idx, const double* values,
                               hgr_result_t** result) {
  // Total timer (v5 D7): starts at function entry — before the D9 content
  // validation — so time_total_ms is the wall clock of the whole
  // hgr_reorder_csr call, validation included, as documented in the header.
  const Timer total_timer;

  if (!reorderer || !row_ptr || !col_idx || !result) return -1;

  // Content validation (v5 D9): before any copy or construction, so a
  // malformed CSR fails cleanly through the error channel instead of
  // corrupting memory downstream.
  {
    const std::string validation_error =
        validate_csr_input(n_rows, nnz, row_ptr, col_idx);
    if (!validation_error.empty()) {
      g_last_error = validation_error;
      return -1;
    }
  }

  // Ingestion timer (v5 D7): covers the CSR array copies + CSRMatrix
  // construction here, on top of the graph build timed inside reorder().
  Timer ingest_timer;

  try {
    // Create CSR matrix
    std::vector<index_t> cpp_row_ptr(row_ptr, row_ptr + n_rows + 1);
    std::vector<index_t> cpp_col_idx(col_idx, col_idx + nnz);
    std::vector<value_t> cpp_values;
    if (values) {
      cpp_values.assign(values, values + nnz);
    }

    CSRMatrix matrix(n_rows, n_rows, std::move(cpp_row_ptr),
                     std::move(cpp_col_idx), std::move(cpp_values),
                     true);  // Assume symmetric

    // Stop the ingestion clock after the copies + construction (D7): it
    // must not include the pipeline itself.
    const double ingest_ms = ingest_timer.elapsed_ms();

    // Run reordering
    auto cpp_result = reorderer->reorderer->reorder(matrix);

    // Ingestion = C-side copy + CSRMatrix + graph build (D7); total = the
    // whole call wall clock, ingestion included.
    cpp_result.stats.time_ingestion_ms += ingest_ms;
    cpp_result.stats.time_total_ms = total_timer.elapsed_ms();

    // Create result handle
    auto* res = new hgr_result_t();
    res->result = std::move(cpp_result);
    *result = res;

    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" int hgr_get_permutation(hgr_result_t* result,
                                   const int64_t** permutation, int64_t* n) {
  if (!result || !permutation || !n) return -1;

  try {
    *permutation = result->result.permutation.data();
    *n = static_cast<int64_t>(result->result.permutation.size());
    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" int hgr_get_partition(hgr_result_t* result, int64_t* n_parts,
                                 int64_t* separator_size,
                                 const int64_t** part_sizes) {
  if (!result || !n_parts || !separator_size || !part_sizes) return -1;

  try {
    *n_parts = result->result.stats.n_parts;
    *separator_size = result->result.stats.separator_size;
    *part_sizes = result->result.stats.part_sizes.data();
    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" int hgr_get_statistics(hgr_result_t* result,
                                  hgr_statistics_t* stats) {
  if (!result || !stats) return -1;

  try {
    const Statistics& s = result->result.stats;
    stats->n_rows = s.n_rows;
    stats->n_cols = s.n_cols;
    stats->nnz = s.nnz;
    stats->n_vertices = s.n_vertices;
    stats->n_edges = s.n_edges;
    stats->degeneracy = s.degeneracy;
    stats->ecc_method_used =
        s.ecc_method_used == EccMethod::kEdgeCover
            ? HGR_ECC_EDGE_COVER
            : s.ecc_method_used == EccMethod::kTriEnum ? HGR_ECC_TRI_ENUM
                                                       : HGR_ECC_BK;
    stats->n_cliques = s.n_cliques;
    stats->max_clique_size = s.max_clique_size;
    stats->avg_clique_size = s.avg_clique_size;
    stats->n_maximal_cliques = s.n_maximal_cliques;
    stats->n_hypernodes = s.n_hypernodes;
    stats->n_hyperedges = s.n_hyperedges;
    stats->total_pins = s.total_pins;
    stats->n_parts = s.n_parts;
    stats->separator_size = s.separator_size;
    stats->separator_ratio = s.separator_ratio;
    stats->km1_objective = s.km1_objective;
    stats->cut_objective = s.cut_objective;
    switch (s.node_weighting_used) {
      case NodeWeighting::kUnit:
        stats->node_weighting_used = HGR_NODE_WEIGHTING_UNIT;
        break;
      case NodeWeighting::kCliqueSize:
        stats->node_weighting_used = HGR_NODE_WEIGHTING_CLIQUE_SIZE;
        break;
      case NodeWeighting::kVertexShare:
        stats->node_weighting_used = HGR_NODE_WEIGHTING_VERTEX_SHARE;
        break;
      case NodeWeighting::kEdgeFillMass:
      default:
        stats->node_weighting_used = HGR_NODE_WEIGHTING_EDGE_FILL_MASS;
        break;
    }
    stats->total_node_weight = s.total_node_weight;
    stats->effective_imbalance = s.effective_imbalance;
    stats->achieved_block_imbalance = s.achieved_block_imbalance;
    stats->time_ingestion_ms = s.time_ingestion_ms;
    stats->time_degeneracy_ms = s.time_degeneracy_ms;
    stats->time_enumeration_ms = s.time_enumeration_ms;
    stats->time_greedy_ms = s.time_greedy_ms;
    stats->time_cnh_ms = s.time_cnh_ms;
    stats->time_partition_ms = s.time_partition_ms;
    stats->time_separator_ms = s.time_separator_ms;
    stats->time_permutation_ms = s.time_permutation_ms;
    stats->time_total_ms = s.time_total_ms;
    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" int hgr_get_clique_order_counts(hgr_result_t* result, int64_t* n,
                                           const int64_t** orders,
                                           const int64_t** counts) {
  if (!result || !n || !orders || !counts) return -1;

  try {
    // Flatten the map (ascending order, as std::map iterates sorted).
    result->clique_orders.clear();
    result->clique_counts.clear();
    for (const auto& [order, count] :
         result->result.stats.clique_order_counts) {
      result->clique_orders.push_back(static_cast<int64_t>(order));
      result->clique_counts.push_back(static_cast<int64_t>(count));
    }
    *n = static_cast<int64_t>(result->clique_orders.size());
    *orders = result->clique_orders.data();
    *counts = result->clique_counts.data();
    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" int hgr_get_maximal_clique_enumeration_histogram(
    hgr_result_t* result, int64_t* n, const int64_t** orders,
    const int64_t** counts) {
  if (!result || !n || !orders || !counts) return -1;

  try {
    // Flatten the map (ascending order, as std::map iterates sorted).
    // Empty for non-BK methods (*n == 0; arrays may be null/empty then —
    // consumers must check n before dereferencing).
    result->enum_clique_orders.clear();
    result->enum_clique_counts.clear();
    for (const auto& [order, count] :
         result->result.stats.maximal_clique_enumeration_histogram) {
      result->enum_clique_orders.push_back(static_cast<int64_t>(order));
      result->enum_clique_counts.push_back(static_cast<int64_t>(count));
    }
    *n = static_cast<int64_t>(result->enum_clique_orders.size());
    *orders = result->enum_clique_orders.data();
    *counts = result->enum_clique_counts.data();
    return 0;
  } catch (...) {
    set_last_error_current_exception();
    return -1;
  }
}

extern "C" size_t hgr_get_mtkahypar_pool_threads(void) {
  return mtkahypar_pool_threads();
}

extern "C" void hgr_result_free(hgr_result_t* result) { delete result; }
