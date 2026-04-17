// clique_cover.hpp - Edge-clique cover algorithm (performance critical!)
#ifndef HYPERGRAPH_REORDER_CLIQUE_COVER_HPP
#define HYPERGRAPH_REORDER_CLIQUE_COVER_HPP

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

// Clique cover solver
class CliqueCoverSolver {
 public:
  struct Options {
    bool use_maximal_cliques;
    bool use_parallel;
    index_t max_clique_enum_vertices;
    index_t max_clique_enum_edges;
    bool greedy_only;
    int num_threads;

    Options()
        : use_maximal_cliques(true),
          use_parallel(true),
          max_clique_enum_vertices(5000),
          max_clique_enum_edges(100000),
          greedy_only(false),
          num_threads(-1) {}
  };

  explicit CliqueCoverSolver(const Options& opts = Options());

  // Main algorithm: compute edge-clique cover
  CliqueCover solve(const Graph& graph);

  // Get solver statistics
  const Statistics& get_stats() const { return stats_; }

 private:
  Options opts_;
  Statistics stats_;

  // Phase 1a: Fast triangle enumeration (for large graphs)
  std::vector<std::vector<index_t>> enumerate_triangles(const Graph& graph);

  // Phase 1b: Find maximal cliques using parallel Bron-Kerbosch
  std::vector<std::vector<index_t>> find_maximal_cliques(const Graph& graph);

  // Bron-Kerbosch algorithm with pivoting
  void bron_kerbosch_pivot(const Graph& graph,
                           std::vector<index_t>& R,  // Current clique
                           std::vector<index_t>& P,  // Candidates
                           std::vector<index_t>& X,  // Already processed
                           std::vector<std::vector<index_t>>& cliques);

  // Phase 2: Greedy clique selection to cover edges
  std::vector<std::vector<index_t>> select_covering_cliques(
      const Graph& graph,
      const std::vector<std::vector<index_t>>& maximal_cliques);

  // Phase 3: Cover remaining edges with 2-cliques
  void cover_remaining_edges(const Graph& graph,
                             const std::vector<bool>& covered_edges,
                             std::vector<std::vector<index_t>>& cliques);

  // Compute edge ID for bitset
  inline index_t edge_id(index_t u, index_t v, index_t n) const {
    if (u > v) std::swap(u, v);
    return u * n + v;
  }
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_CLIQUE_COVER_HPP
