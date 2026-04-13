// types.hpp - Core type definitions for hypergraph reordering library
#ifndef HYPERGRAPH_REORDER_TYPES_HPP
#define HYPERGRAPH_REORDER_TYPES_HPP

#include <cstdint>
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

// Timing utilities
struct Timer {
  double start_time;

  Timer();
  void reset();
  double elapsed_ms() const;

 private:
  static double get_time();
};

// Statistics structure
struct Statistics {
  // Matrix statistics
  index_t n_rows = 0;
  index_t n_cols = 0;
  index_t nnz = 0;

  // Graph statistics
  index_t n_vertices = 0;
  index_t n_edges = 0;

  // Clique cover statistics
  index_t n_cliques = 0;
  index_t max_clique_size = 0;
  double avg_clique_size = 0.0;

  // Hypergraph statistics
  index_t n_hypernodes = 0;
  index_t n_hyperedges = 0;
  index_t total_pins = 0;

  // Partition statistics
  index_t n_parts = 0;
  std::vector<index_t> part_sizes;
  index_t separator_size = 0;
  double separator_ratio = 0.0;

  // Ordering statistics
  index_t blocks_ordered = 0;
  index_t blocks_failed = 0;

  // Timing breakdown (in milliseconds)
  double time_load_ms = 0.0;
  double time_graph_construction_ms = 0.0;
  double time_clique_cover_ms = 0.0;
  double time_hypergraph_construction_ms = 0.0;
  double time_partitioning_ms = 0.0;
  double time_separator_construction_ms = 0.0;
  double time_block_ordering_ms = 0.0;
  double time_permutation_ms = 0.0;
  double time_total_ms = 0.0;
};

// Ordering method enum
enum class OrderingMethod {
  AMD,      // Approximate Minimum Degree
  CAMD,     // Constrained AMD
  METIS,    // METIS nested dissection
  NESDIS,   // Nested dissection
  COLAMD,   // Column AMD
  NATURAL,  // Natural ordering (no reordering)
  NONE      // Skip block ordering
};

// Convert ordering method to string
inline const char* ordering_method_to_string(OrderingMethod method) {
  switch (method) {
    case OrderingMethod::AMD:
      return "AMD";
    case OrderingMethod::CAMD:
      return "CAMD";
    case OrderingMethod::METIS:
      return "METIS";
    case OrderingMethod::NESDIS:
      return "NESDIS";
    case OrderingMethod::COLAMD:
      return "COLAMD";
    case OrderingMethod::NATURAL:
      return "NATURAL";
    case OrderingMethod::NONE:
      return "NONE";
    default:
      return "UNKNOWN";
  }
}

// Parse ordering method from string
inline OrderingMethod parse_ordering_method(const std::string& str) {
  if (str == "amd" || str == "AMD") return OrderingMethod::AMD;
  if (str == "camd" || str == "CAMD") return OrderingMethod::CAMD;
  if (str == "metis" || str == "METIS") return OrderingMethod::METIS;
  if (str == "nesdis" || str == "NESDIS") return OrderingMethod::NESDIS;
  if (str == "colamd" || str == "COLAMD") return OrderingMethod::COLAMD;
  if (str == "natural" || str == "NATURAL") return OrderingMethod::NATURAL;
  if (str == "none" || str == "NONE") return OrderingMethod::NONE;
  throw HypergraphReorderError("Unknown ordering method: " + str);
}

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_TYPES_HPP
