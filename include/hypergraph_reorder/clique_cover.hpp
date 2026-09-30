// clique_cover.hpp - Clique-cover data structure (dual CSR representation)
#ifndef HYPERGRAPH_REORDER_CLIQUE_COVER_HPP
#define HYPERGRAPH_REORDER_CLIQUE_COVER_HPP

#include <map>
#include <span>
#include <unordered_map>
#include <vector>

#include "graph.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Clique cover data structure with dual representation
class CliqueCover {
 public:
  CliqueCover();

  // Construct from cliques
  CliqueCover(index_t n_vertices, std::vector<std::vector<index_t>> cliques);

  // Accessors
  index_t n_cliques() const { return n_cliques_; }
  index_t n_vertices() const { return n_vertices_; }
  index_t max_clique_size() const { return max_clique_size_; }
  double avg_clique_size() const;

  // Clique-order histogram of this cover: clique order (size) -> number of
  // cliques of that order. O(n_cliques). Used by reorder() to populate
  // Statistics::clique_order_counts.
  std::map<index_t, index_t> clique_order_counts() const;

  // Get members of a clique (zero-copy span)
  std::span<const index_t> get_clique(index_t clique_id) const;

  // Get cliques containing a vertex (zero-copy span)
  std::span<const index_t> get_vertex_cliques(index_t vertex_id) const;

  // Access raw data
  const std::vector<index_t>& clique_ptr() const { return clique_ptr_; }
  const std::vector<index_t>& clique_members() const { return clique_members_; }
  const std::vector<index_t>& vertex_cliques_ptr() const {
    return vertex_cliques_ptr_;
  }
  const std::vector<index_t>& vertex_cliques() const { return vertex_cliques_; }

 private:
  index_t n_cliques_;
  index_t n_vertices_;
  index_t max_clique_size_;

  // Clique-centric representation: clique_id -> members
  std::vector<index_t> clique_ptr_;      // Size: n_cliques + 1
  std::vector<index_t> clique_members_;  // Size: sum of clique sizes

  // Vertex-centric representation: vertex_id -> cliques
  std::vector<index_t> vertex_cliques_ptr_;  // Size: n_vertices + 1
  std::vector<index_t> vertex_cliques_;      // Size: total incidences
};

// The ECC *solvers* live in ecc.hpp / ecc_*.cpp (EdgeCover, TriEnum, BK).

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_CLIQUE_COVER_HPP
