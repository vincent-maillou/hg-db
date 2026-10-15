// types.hpp - Core type definitions for hypergraph reordering library
#ifndef HYPERGRAPH_REORDER_TYPES_HPP
#define HYPERGRAPH_REORDER_TYPES_HPP

#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace hypergraph_reorder {

// Core index type - int64_t for large matrices (millions of rows)
using index_t = int64_t;

// Value type for matrix entries
using value_t = double;

// Edge ID type
using edge_id_t = int64_t;

// Partition ID type
using part_id_t = int32_t;

// Constants
constexpr index_t INVALID_INDEX = -1;
constexpr part_id_t SEPARATOR_PART = -1;

// Error handling
class HypergraphReorderError : public std::runtime_error {
 public:
  explicit HypergraphReorderError(const std::string& msg)
      : std::runtime_error(msg) {}
};

// Edge-clique cover (ECC) method selection.
// Defined here (not in ecc.hpp) so that Statistics can reference it without
// an include cycle (ecc.hpp -> clique_cover.hpp -> types.hpp); ecc.hpp
// re-exports it via its include of types.hpp.
//  - kEdgeCover: every edge becomes a 2-clique (ablation baseline)
//  - kTriEnum:   triangle candidates + greedy cover + 2-clique remainder
//  - kBk:        maximal-clique candidates (Bron-Kerbosch) + greedy cover +
//                2-clique remainder
enum class EccMethod { kEdgeCover, kTriEnum, kBk };

// MT-KaHyPar preset configuration.
// Lives in types.hpp (dependency-light, shared base) so that both
// PartitionOptions and the C API layer can reference it.
enum class MtKahyparPreset {
  DEFAULT,        // Fast, good quality (default preset)
  QUALITY,        // Higher quality, slower
  DETERMINISTIC,  // Deterministic partitioning
  LARGE_K         // Optimized for large number of parts
};

// MT-KaHyPar partitioning objective. The objective governs what the CNH
// partition minimizes under the imbalance constraint; separator extraction
// and everything downstream are objective-agnostic.
//  - kKm1: connectivity metric, sum over nets of (connectivity - 1).
//    Default (historical behavior).
//  - kCut: cut-net metric, number of nets spanning >= 2 parts.
// SOED is not exposed: the Parallax MT-KaHyPar build does not enable
// KAHYPAR_ENABLE_SOED_METRIC. Value 2 is reserved for a future addition.
enum class MtKahyparObjective { kKm1, kCut };

// Hypernode (clique) weighting for the CNH partition stage. Governs what
// MT-KaHyPar balances; the imbalance parameter bounds part weights under
// the selected scheme (see ECC.md "Hypernode weighting").
enum class NodeWeighting {
  kUnit,          // w_i = 1                     (ablation baseline)
  kCliqueSize,    // w_i = |C_i|                 (default)
  kVertexShare,   // per-vertex unit split over cliques (largest remainder)
  kEdgeFillMass   // w_i = #edges first-covered by C_i
};

// Vertex partition (result of converting a CNH partition to a vertex
// separator). Lives in types.hpp alongside the other shared value types.
struct VertexPartition {
  index_t n_parts;
  std::vector<std::vector<index_t>> parts;  // parts[i] = vertices in part i
  std::vector<index_t> separator;           // Separator vertices

  index_t separator_size() const { return separator.size(); }
  double separator_ratio(index_t n_total) const {
    return static_cast<double>(separator.size()) / n_total;
  }
};

// Timing utilities
struct Timer {
  double start_time;

  Timer();
  void reset();
  double elapsed_ms() const;

 private:
  static double get_time();
};

// Statistics structure: every field is populated by reorder() — no
// placeholders, no dead fields.
struct Statistics {
  // Matrix statistics
  index_t n_rows = 0;
  index_t n_cols = 0;
  index_t nnz = 0;

  // Graph statistics
  index_t n_vertices = 0;
  index_t n_edges = 0;

  // Edge-clique cover statistics
  EccMethod ecc_method_used = EccMethod::kBk;
  index_t n_cliques = 0;
  index_t max_clique_size = 0;
  double avg_clique_size = 0.0;
  // Clique-order histogram of the *final* cover (after greedy selection and
  // the 2-clique remainder): clique order (size) -> number of cliques of
  // that order. Not part of the plain-C statistics struct; exposed through
  // a dedicated C accessor.
  std::map<index_t, index_t> clique_order_counts;

  // Clique-node hypergraph statistics
  index_t n_hypernodes = 0;
  index_t n_hyperedges = 0;
  index_t total_pins = 0;

  // Partition statistics
  index_t n_parts = 0;
  std::vector<index_t> part_sizes;
  index_t separator_size = 0;
  double separator_ratio = 0.0;
  // Connectivity-(km1) cut metric of the CNH partition as reported by
  // MT-KaHyPar: sum over nets of (connectivity - 1). This is the quantity
  // MT-KaHyPar minimizes; it is NOT the separator size.
  index_t km1_objective = 0;
  // Cut-net metric of the CNH partition as reported by MT-KaHyPar:
  // number of nets (original vertices) whose cliques span >= 2 parts.
  // Reported for every run regardless of the optimized objective, so
  // km1-runs and cut-runs stay comparable. With unit net weights
  // km1_objective >= cut_objective always (equality iff no net spans
  // more than 2 parts).
  index_t cut_objective = 0;

  // Hypernode weighting (v4). node_weighting_used echoes the scheme the
  // partitioner balanced under; total_node_weight is W = sum_i w_i under
  // that scheme; effective_imbalance is eps after auto-relax (>= the
  // requested imbalance; see partitioner); achieved_block_imbalance is
  // max(part sizes) / (n_interior / k) — the vertex-level balance the DB
  // output actually cares about (0.0 when there are no interior vertices).
  NodeWeighting node_weighting_used = NodeWeighting::kCliqueSize;
  index_t total_node_weight = 0;         // W under the used scheme
  double effective_imbalance = 0.0;      // eps after auto-relax
  double achieved_block_imbalance = 0.0; // max part / (n_interior / k)

  // Timing breakdown (milliseconds) — seven stages including the total,
  // which is the wall clock around the whole reorder() call.
  double time_graph_ms = 0.0;
  double time_ecc_ms = 0.0;
  double time_cnh_ms = 0.0;
  double time_partition_ms = 0.0;
  double time_separator_ms = 0.0;
  double time_permutation_ms = 0.0;
  double time_total_ms = 0.0;
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_TYPES_HPP
