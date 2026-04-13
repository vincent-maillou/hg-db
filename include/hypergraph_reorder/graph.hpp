#ifndef HYPERGRAPH_REORDER_GRAPH_HPP
#define HYPERGRAPH_REORDER_GRAPH_HPP

#include <span>
#include <vector>

#include "sparse_matrix.hpp"
#include "types.hpp"

namespace hypergraph_reorder {
class Graph {
 public:
  Graph();

  explicit Graph(index_t n_vertices);

  Graph(index_t n_vertices, index_t n_edges, std::vector<index_t> adj_ptr,
        std::vector<index_t> adj_list);

  index_t n_vertices() const { return n_vertices_; }
  index_t n_edges() const { return n_edges_; }

  std::span<const index_t> neighbors(index_t v) const;
  index_t degree(index_t v) const;
  bool has_edge(index_t u, index_t v) const;

  const std::vector<index_t> &adj_ptr() const { return adj_ptr_; }
  const std::vector<index_t> &adj_list() const { return adj_list_; }

  static Graph from_symmetric_matrix(const CSRMatrix &matrix);

  bool validate() const;
  std::vector<index_t> compute_degeneracy_ordering() const;
  Graph induced_subgraph(const std::vector<index_t> &vertices) const;

 private:
  index_t n_vertices_;
  index_t n_edges_;

  std::vector<index_t> adj_ptr_;
  std::vector<index_t> adj_list_;
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_GRAPH_HPP
