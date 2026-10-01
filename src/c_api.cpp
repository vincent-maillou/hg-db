#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "hypergraph_reorder.h"
#include "hypergraph_reorder/reorderer.hpp"

/* Fallbacks when the build system does not inject the project version */
#ifndef HGR_VERSION_MAJOR
#define HGR_VERSION_MAJOR 2
#endif
#ifndef HGR_VERSION_MINOR
#define HGR_VERSION_MINOR 0
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
  opts->partition.seed = cpp_opts.partition.seed;

  opts->runtime.num_threads = cpp_opts.runtime.num_threads;
  opts->runtime.suppress_output = cpp_opts.runtime.suppress_output ? 1 : 0;
}

extern "C" int hgr_reorder_csr(hgr_reorderer_t* reorderer, int64_t n_rows,
                               int64_t nnz, const int64_t* row_ptr,
                               const int64_t* col_idx, const double* values,
                               hgr_result_t** result) {
  if (!reorderer || !row_ptr || !col_idx || !result) return -1;

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

    // Run reordering
    auto cpp_result = reorderer->reorderer->reorder(matrix);

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
    stats->ecc_method_used =
        s.ecc_method_used == EccMethod::kEdgeCover
            ? HGR_ECC_EDGE_COVER
            : s.ecc_method_used == EccMethod::kTriEnum ? HGR_ECC_TRI_ENUM
                                                       : HGR_ECC_BK;
    stats->n_cliques = s.n_cliques;
    stats->max_clique_size = s.max_clique_size;
    stats->avg_clique_size = s.avg_clique_size;
    stats->n_hypernodes = s.n_hypernodes;
    stats->n_hyperedges = s.n_hyperedges;
    stats->total_pins = s.total_pins;
    stats->n_parts = s.n_parts;
    stats->separator_size = s.separator_size;
    stats->separator_ratio = s.separator_ratio;
    stats->km1_objective = s.km1_objective;
    stats->time_graph_ms = s.time_graph_ms;
    stats->time_ecc_ms = s.time_ecc_ms;
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

extern "C" void hgr_result_free(hgr_result_t* result) { delete result; }
