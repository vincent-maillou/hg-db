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

/* Ordering method enum */
typedef enum {
  HGR_ORDERING_AMD = 0,
  HGR_ORDERING_CAMD = 1,
  HGR_ORDERING_METIS = 2,
  HGR_ORDERING_NESDIS = 3,
  HGR_ORDERING_COLAMD = 4,
  HGR_ORDERING_NATURAL = 5,
  HGR_ORDERING_NONE = 6
} hgr_ordering_method_t;

/* MT-KaHyPar preset enum */
typedef enum {
  HGR_PRESET_DEFAULT = 0,       /* Fast, good quality */
  HGR_PRESET_QUALITY = 1,       /* Higher quality, slower */
  HGR_PRESET_DETERMINISTIC = 2, /* Deterministic partitioning */
  HGR_PRESET_LARGE_K = 3        /* Optimized for large number of parts */
} hgr_preset_t;

/* Matrix format enum */
typedef enum {
  HGR_FORMAT_MTX = 0,
  HGR_FORMAT_METIS = 1,
  HGR_FORMAT_CSR_BINARY = 2,
  HGR_FORMAT_AUTO = 3
} hgr_matrix_format_t;

/* Options structure */
typedef struct {
  int64_t n_parts;
  double imbalance;
  hgr_preset_t preset; /* MT-KaHyPar preset */
  int seed;
  hgr_ordering_method_t ordering_method;
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
  double time_block_ordering_ms;
} hgr_statistics_t;

/* Create/destroy reorderer */
hgr_reorderer_t* hgr_create(const hgr_options_t* opts);
void hgr_free(hgr_reorderer_t* reorderer);

/* Set default options */
void hgr_default_options(hgr_options_t* opts);

/* Main reordering function */
int hgr_reorder_file(hgr_reorderer_t* reorderer, const char* input_path,
                     hgr_matrix_format_t format, hgr_result_t** result);

/* Reorder from CSR matrix in memory */
int hgr_reorder_csr(hgr_reorderer_t* reorderer, int64_t n_rows, int64_t nnz,
                    const int64_t* row_ptr, const int64_t* col_idx,
                    const double* values, /* can be NULL */
                    hgr_result_t** result);

/* Extract results */
int hgr_get_permutation(hgr_result_t* result, const int64_t** permutation,
                        int64_t* n);

int hgr_get_reordered_matrix(hgr_result_t* result, int64_t* n_rows,
                             int64_t* nnz, const int64_t** row_ptr,
                             const int64_t** col_idx, const double** values);

/* Save result to file */
int hgr_save_matrix(hgr_result_t* result, const char* output_path,
                    hgr_matrix_format_t format);

/* Get statistics */
int hgr_get_statistics(hgr_result_t* result, hgr_statistics_t* stats);

/* Get block sizes after partitioning.
   Sets *part_sizes to an array of n_parts sizes (vertices per block),
   and *n to the number of parts. The pointer is valid until hgr_result_free. */
int hgr_get_part_sizes(hgr_result_t* result, const int64_t** part_sizes,
                       int64_t* n, int64_t* separator)

    /* Free result */
    void hgr_result_free(hgr_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* HYPERGRAPH_REORDER_H */
