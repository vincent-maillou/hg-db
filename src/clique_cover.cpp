#include "hypergraph_reorder/clique_cover.hpp"

#include <algorithm>


namespace hypergraph_reorder {

// ===== CliqueCover Implementation =====

CliqueCover::CliqueCover()
    : n_cliques_(0), n_vertices_(0), max_clique_size_(0) {}

CliqueCover::CliqueCover(index_t n_vertices,
                         std::vector<std::vector<index_t>> cliques)
    : n_cliques_(cliques.size()), n_vertices_(n_vertices), max_clique_size_(0) {
  // Build clique-centric representation
  clique_ptr_.resize(n_cliques_ + 1);
  clique_ptr_[0] = 0;

  index_t total_clique_size = 0;
  for (size_t i = 0; i < cliques.size(); ++i) {
    total_clique_size += cliques[i].size();
    clique_ptr_[i + 1] = total_clique_size;
    max_clique_size_ =
        std::max(max_clique_size_, static_cast<index_t>(cliques[i].size()));
  }

  clique_members_.reserve(total_clique_size);
  for (const auto &clique : cliques) {
    clique_members_.insert(clique_members_.end(), clique.begin(), clique.end());
  }

  // Build vertex-centric representation (vertex -> cliques) directly as
  // CSR: count, prefix-sum, fill. Pass 2 iterates cliques in ascending
  // cid order and each clique contributes each member exactly once, so
  // every vertex's entries end up ascending AND deduplicated by
  // construction (WP2 D2; the invariant D4 asserts on the CNH side).

  // Pass 1: incidence counts per vertex.
  std::vector<index_t> count(n_vertices + 1, 0);
  for (index_t cid = 0; cid < n_cliques_; ++cid) {
    for (auto v : get_clique(cid)) {
      ++count[v + 1];
    }
  }
  // Exclusive prefix sum -> vertex_cliques_ptr_.
  vertex_cliques_ptr_.resize(n_vertices + 1);
  vertex_cliques_ptr_[0] = 0;
  for (index_t v = 0; v < n_vertices; ++v) {
    vertex_cliques_ptr_[v + 1] = vertex_cliques_ptr_[v] + count[v + 1];
  }
  // Pass 2: fill member array; cursor[v] is the next write slot of v.
  vertex_cliques_.resize(vertex_cliques_ptr_[n_vertices]);
  std::vector<index_t> cursor(vertex_cliques_ptr_.begin(),
                              vertex_cliques_ptr_.end() - 1);
  for (index_t cid = 0; cid < n_cliques_; ++cid) {
    for (auto v : get_clique(cid)) {
      vertex_cliques_[cursor[v]++] = cid;
    }
  }
}

double CliqueCover::avg_clique_size() const {
  if (n_cliques_ == 0) return 0.0;
  return static_cast<double>(clique_members_.size()) / n_cliques_;
}

std::map<index_t, index_t> CliqueCover::clique_order_counts() const {
  std::map<index_t, index_t> counts;
  for (index_t cid = 0; cid < n_cliques_; ++cid) {
    index_t order = clique_ptr_[cid + 1] - clique_ptr_[cid];
    ++counts[order];
  }
  return counts;
}

std::span<const index_t> CliqueCover::get_clique(index_t clique_id) const {
  if (clique_id < 0 || clique_id >= n_cliques_) {
    throw HypergraphReorderError("Clique ID out of bounds");
  }
  index_t start = clique_ptr_[clique_id];
  index_t end = clique_ptr_[clique_id + 1];
  return std::span<const index_t>(clique_members_.data() + start, end - start);
}

std::span<const index_t> CliqueCover::get_vertex_cliques(
    index_t vertex_id) const {
  if (vertex_id < 0 || vertex_id >= n_vertices_) {
    throw HypergraphReorderError("Vertex ID out of bounds");
  }
  index_t start = vertex_cliques_ptr_[vertex_id];
  index_t end = vertex_cliques_ptr_[vertex_id + 1];
  return std::span<const index_t>(vertex_cliques_.data() + start, end - start);
}

}  // namespace hypergraph_reorder
