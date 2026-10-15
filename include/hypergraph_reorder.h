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

/* Edge-clique cover (ECC) method */
typedef enum {
  HGR_ECC_EDGE_COVER = 0, /* Every edge becomes a 2-clique (ablation baseline) */
  HGR_ECC_TRI_ENUM = 1,   /* Triangle candidates + greedy cover + remainder */
  HGR_ECC_BK = 2          /* Maximal cliques (Bron-Kerbosch) + greedy cover */
} hgr_ecc_method_t;

/* MT-KaHyPar partitioning objective for the CNH stage.
 * Governs what the partition minimizes under the imbalance constraint;
 * separator extraction is objective-agnostic.
 * HGR_OBJECTIVE_KM1 (connectivity metric) is the default and the
 * historical behavior. Value 2 is reserved for a possible future SOED
 * (SOED requires a MT-KaHyPar build with KAHYPAR_ENABLE_SOED_METRIC). */
typedef enum {
  HGR_OBJECTIVE_KM1 = 0, /* sum over nets of (connectivity - 1) */
  HGR_OBJECTIVE_CUT = 1  /* number of nets spanning >= 2 parts */
} hgr_objective_t;

/* Hypernode (clique) weighting for the CNH partition stage: what
 * MT-KaHyPar balances. The imbalance bound applies to part weights under
 * the selected scheme (see hgr_statistics_t::node_weighting_used and the
 * per-scheme W invariants documented in ECC.md "Hypernode weighting"). */
typedef enum {
  HGR_NODE_WEIGHTING_UNIT = 0,          /* w_i = 1 (ablation baseline, v3 behavior) */
  HGR_NODE_WEIGHTING_CLIQUE_SIZE = 1,   /* w_i = |C_i| (default, v4) */
  HGR_NODE_WEIGHTING_VERTEX_SHARE = 2,  /* per-vertex unit split over cliques */
  HGR_NODE_WEIGHTING_EDGE_FILL_MASS = 3 /* w_i = #edges first-covered by C_i */
} hgr_node_weighting_t;

/* ECC options */
typedef struct {
  hgr_ecc_method_t method; /* default HGR_ECC_BK */
  int use_parallel;        /* 1 = OpenMP for enumeration + selection */
} hgr_ecc_options_t;

/* Partitioning options */
typedef struct {
  int64_t n_parts;
  /* Imbalance bound for part weights under the selected node weighting
   * scheme. May be auto-relaxed when a single hypernode exceeds a block's
   * weight budget — the value actually used is reported as
   * hgr_statistics_t::effective_imbalance. */
  double imbalance;
  hgr_preset_t preset;
  int seed;                             /* default 42; -1 = random */
  hgr_objective_t objective;            /* v3 */
  hgr_node_weighting_t node_weighting;  /* v4: what the partitioner balances */
} hgr_partition_options_t;

/* Runtime options */
typedef struct {
  int num_threads;     /* 0 = auto-detect */
  int suppress_output; /* 1 = suppress ALL progress output (including
                          MT-KaHyPar verbose) */
} hgr_runtime_options_t;

/* Options structure (v4: partition.node_weighting added; ABI break vs v3) */
typedef struct {
  hgr_ecc_options_t ecc;
  hgr_partition_options_t partition;
  hgr_runtime_options_t runtime;
} hgr_options_t;

/* Statistics structure (v4: node weighting observables added; ABI break
 * vs v3) */
typedef struct {
  /* Matrix */
  int64_t n_rows, n_cols, nnz;
  /* Graph */
  int64_t n_vertices, n_edges;
  /* ECC */
  hgr_ecc_method_t ecc_method_used;
  int64_t n_cliques, max_clique_size;
  double avg_clique_size;
  /* CNH */
  int64_t n_hypernodes, n_hyperedges, total_pins;
  /* Partition */
  int64_t n_parts, separator_size;
  double separator_ratio;
  /* km1_objective: connectivity-(km1) cut metric of the CNH partition as
   * reported by MT-KaHyPar — the sum over nets of (connectivity - 1).
   * Every net (original vertex) whose cliques span lambda >= 2 parts
   * contributes lambda - 1. This is the quantity MT-KaHyPar minimizes
   * under the imbalance constraint; it is NOT the separator size (a net
   * split across parts puts its vertex in the separator, and vertices in
   * no net — isolated ones — land in the separator without affecting
   * km1). Main configuration levers: n_parts (more parts -> larger km1
   * and typically a larger separator), imbalance (tighter -> larger),
   * preset (quality/deterministic trade km1 vs time/reproducibility),
   * seed (run-to-run variance at fixed config), objective (km1 vs cut). */
  int64_t km1_objective;
  /* Cut-net metric of the CNH partition: number of nets (original
   * vertices) whose cliques span >= 2 parts, as reported by MT-KaHyPar.
   * Reported for every run regardless of the optimized objective.
   * With unit net weights km1_objective >= cut_objective always. */
  int64_t cut_objective;
  /* Node weighting (v4). node_weighting_used echoes the scheme the
   * partitioner balanced under; total_node_weight is W = sum_i w_i under
   * that scheme (unit: n_cliques, clique_size: total_pins, vertex_share:
   * SCALE*n_nonisolated, edge_fill_mass: n_edges); effective_imbalance is
   * the imbalance actually used after auto-relax (>= the requested one);
   * achieved_block_imbalance is max(part sizes) / (n_interior / k) — the
   * vertex-level balance the doubly-bordered output cares about (0.0
   * when there are no interior vertices). */
  hgr_node_weighting_t node_weighting_used;
  int64_t total_node_weight;
  double effective_imbalance;
  double achieved_block_imbalance;
  /* Timings (ms) */
  double time_graph_ms, time_ecc_ms, time_cnh_ms, time_partition_ms,
      time_separator_ms, time_permutation_ms, time_total_ms;
} hgr_statistics_t;

/* Library version: fills version[0..2] = {major, minor, patch} */
void hgr_version(int version[3]);

/* Create/destroy reorderer.
 * hgr_create maps ALL nested option fields; on failure sets the
 * thread-local last-error (see hgr_get_last_error) and returns NULL. */
hgr_reorderer_t* hgr_create(const hgr_options_t* opts);
void hgr_free(hgr_reorderer_t* reorderer);

/* Set default options. Single source of defaults: constructs the C++
 * SymmetricDBReorderer::Options() and transcribes it, so the defaults can
 * never diverge between layers. */
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

/* Clique-order histogram of the final ECC.
 * On success: *orders and *counts point to internal arrays of length *n,
 *   orders[i] = clique order (size), counts[i] = number of cliques of that
 *   order in the final cover. Arrays are sorted by ascending order and are
 *   owned by the result (valid until hgr_result_free). */
int hgr_get_clique_order_counts(hgr_result_t* result, int64_t* n,
                                const int64_t** orders,
                                const int64_t** counts);

/* Returns the message of the last error that occurred on the calling
 * thread, or "" when no error occurred. Thread-local; valid until the next
 * hgr_* call on the same thread. Never returns NULL. */
const char* hgr_get_last_error(void);

/* Free result */
void hgr_result_free(hgr_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* HYPERGRAPH_REORDER_H */
