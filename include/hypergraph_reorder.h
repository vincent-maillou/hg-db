/* hypergraph_reorder.h - C API for hypergraph-based matrix reordering */
#ifndef HYPERGRAPH_REORDER_H
#define HYPERGRAPH_REORDER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle types */
typedef struct hgr_reorderer_s hgr_reorderer_t;
typedef struct hgr_result_s hgr_result_t;

/* MT-KaHyPar preset enum */
typedef enum {
  HGR_PRESET_DEFAULT = 0,       /* Fast, good quality */
  HGR_PRESET_QUALITY = 1,       /* Higher quality, slower */
  HGR_PRESET_DETERMINISTIC = 2, /* Deterministic partitioning */
  HGR_PRESET_LARGE_K = 3        /* Optimized for large number of parts */
} hgr_preset_t;

/* Options structure */
typedef struct {
  int64_t n_parts;
  double imbalance;
  hgr_preset_t preset; /* MT-KaHyPar preset */
  int seed;
  int use_openmp;
  int num_threads; /* 0 = auto-detect */
  int suppress_partitioner_output;
  int suppress_output; /* Suppress all reorderer progress output */
} hgr_options_t;

/* Statistics structure */
typedef struct {
  int64_t n_rows;
  int64_t n_cols;
  int64_t nnz;
  int64_t n_parts;
  int64_t separator_size;
  double separator_ratio;
  double time_total_ms;
  double time_clique_cover_ms;
  double time_partitioning_ms;
} hgr_statistics_t;

/* Create/destroy reorderer */
hgr_reorderer_t* hgr_create(const hgr_options_t* opts);
void hgr_free(hgr_reorderer_t* reorderer);

/* Set default options */
void hgr_default_options(hgr_options_t* opts);

/* Reorder from CSR matrix in memory */
int hgr_reorder_csr(hgr_reorderer_t* reorderer, int64_t n_rows, int64_t nnz,
                    const int64_t* row_ptr, const int64_t* col_idx,
                    const double* values, /* can be NULL */
                    hgr_result_t** result);

/* Extract results */
int hgr_get_permutation(hgr_result_t* result, const int64_t** permutation,
                        int64_t* n);

/* Get partition information (diagonal block sizes + separator size) */
int hgr_get_partition(hgr_result_t* result, int64_t* n_parts,
                      int64_t* separator_size, const int64_t** part_sizes);

/* Get statistics */
int hgr_get_statistics(hgr_result_t* result, hgr_statistics_t* stats);

/* Free result */
void hgr_result_free(hgr_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* HYPERGRAPH_REORDER_H */
