#include <cstring>
#include <memory>

#include "hypergraph_reorder.h"
#include "hypergraph_reorder/reorderer.hpp"

using namespace hypergraph_reorder;

// Opaque structures
struct hgr_reorderer_s {
  std::unique_ptr<SymmetricDBReorderer> reorderer;
};

struct hgr_result_s {
  SymmetricDBReorderer::Result result;
};

// Create reorderer
extern "C" hgr_reorderer_t* hgr_create(const hgr_options_t* opts) {
  try {
    SymmetricDBReorderer::Options cpp_opts;

    if (opts) {
      cpp_opts.n_parts = opts->n_parts;
      cpp_opts.imbalance = opts->imbalance;
      cpp_opts.seed = opts->seed;
      cpp_opts.use_openmp = opts->use_openmp != 0;
      cpp_opts.num_threads = opts->num_threads;
      cpp_opts.suppress_partitioner_output =
          opts->suppress_partitioner_output != 0;
      cpp_opts.suppress_output = opts->suppress_output != 0;

      // Convert preset
      switch (opts->preset) {
        case HGR_PRESET_QUALITY:
          cpp_opts.preset = MtKahyparPreset::QUALITY;
          break;
        case HGR_PRESET_DETERMINISTIC:
          cpp_opts.preset = MtKahyparPreset::DETERMINISTIC;
          break;
        case HGR_PRESET_LARGE_K:
          cpp_opts.preset = MtKahyparPreset::LARGE_K;
          break;
        case HGR_PRESET_DEFAULT:
        default:
          cpp_opts.preset = MtKahyparPreset::DEFAULT;
          break;
      }
    }

    auto* handle = new hgr_reorderer_t();
    handle->reorderer = std::make_unique<SymmetricDBReorderer>(cpp_opts);
    return handle;
  } catch (...) {
    return nullptr;
  }
}

extern "C" void hgr_free(hgr_reorderer_t* reorderer) { delete reorderer; }

extern "C" void hgr_default_options(hgr_options_t* opts) {
  if (!opts) return;

  opts->n_parts = 4;
  opts->imbalance = 0.03;
  opts->preset = HGR_PRESET_DEFAULT;
  opts->seed = -1;
  opts->use_openmp = 1;
  opts->num_threads = 0;  // 0 = auto-detect
  opts->suppress_partitioner_output = 0;
  opts->suppress_output = 0;
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
    return -1;
  }
}

extern "C" int hgr_get_permutation(hgr_result_t* result,
                                   const int64_t** permutation, int64_t* n) {
  if (!result || !permutation || !n) return -1;

  try {
    *permutation = result->result.permutation.data();
    *n = result->result.permutation.size();
    return 0;
  } catch (...) {
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
    return -1;
  }
}

extern "C" int hgr_get_statistics(hgr_result_t* result,
                                  hgr_statistics_t* stats) {
  if (!result || !stats) return -1;

  try {
    stats->n_rows = result->result.stats.n_rows;
    stats->n_cols = result->result.stats.n_cols;
    stats->nnz = result->result.stats.nnz;
    stats->n_parts = result->result.stats.n_parts;
    stats->separator_size = result->result.stats.separator_size;
    stats->separator_ratio = result->result.stats.separator_ratio;
    stats->time_total_ms = result->result.stats.time_total_ms;
    stats->time_clique_cover_ms = result->result.stats.time_clique_cover_ms;
    stats->time_partitioning_ms = result->result.stats.time_partitioning_ms;
    return 0;
  } catch (...) {
    return -1;
  }
}

extern "C" void hgr_result_free(hgr_result_t* result) { delete result; }
