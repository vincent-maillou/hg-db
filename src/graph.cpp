#include "hypergraph_reorder/graph.hpp"

#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace hypergraph_reorder {

Graph::Graph() : n_vertices_(0), n_edges_(0) {}

Graph::Graph(index_t n_vertices) : n_vertices_(n_vertices), n_edges_(0) {
  adj_ptr_.resize(n_vertices + 1, 0);
}

Graph::Graph(index_t n_vertices, index_t n_edges, std::vector<index_t> adj_ptr,
             std::vector<index_t> adj_list)
    : n_vertices_(n_vertices),
      n_edges_(n_edges),
      adj_ptr_(std::move(adj_ptr)),
      adj_list_(std::move(adj_list)) {
  if (adj_ptr_.size() != static_cast<size_t>(n_vertices + 1)) {
    throw HypergraphReorderError("Invalid adj_ptr size");
  }
  if (adj_ptr_[n_vertices] != static_cast<index_t>(adj_list_.size())) {
    throw HypergraphReorderError("adj_list size doesn't match adj_ptr");
  }
}

std::span<const index_t> Graph::neighbors(index_t v) const {
  if (v < 0 || v >= n_vertices_) {
    throw HypergraphReorderError("Vertex index out of bounds");
  }
  index_t start = adj_ptr_[v];
  index_t end = adj_ptr_[v + 1];
  return std::span<const index_t>(adj_list_.data() + start, end - start);
}

index_t Graph::degree(index_t v) const {
  if (v < 0 || v >= n_vertices_) {
    throw HypergraphReorderError("Vertex index out of bounds");
  }
  return adj_ptr_[v + 1] - adj_ptr_[v];
}

bool Graph::has_edge(index_t u, index_t v) const {
  if (u < 0 || u >= n_vertices_ || v < 0 || v >= n_vertices_) {
    return false;
  }

  if (degree(u) <= degree(v)) {
    auto neighs = neighbors(u);
    return std::binary_search(neighs.begin(), neighs.end(), v);
  } else {
    auto neighs = neighbors(v);
    return std::binary_search(neighs.begin(), neighs.end(), u);
  }
}

Graph Graph::from_symmetric_matrix(const CSRMatrix &matrix) {
  if (matrix.n_rows() != matrix.n_cols()) {
    throw HypergraphReorderError(
        "Matrix must be square for graph construction");
  }

  index_t n = matrix.n_rows();

  std::vector<index_t> degrees(n, 0);
  index_t n_edges = 0;

  for (index_t i = 0; i < n; ++i) {
    for (auto j : matrix.row_indices(i)) {
      if (i < j) {
        degrees[i]++;
        degrees[j]++;
        n_edges++;
      } else if (i == j) {
        continue;
      } else {
        if (matrix.is_symmetric()) {
          degrees[i]++;
          degrees[j]++;
          n_edges++;
        }
      }
    }
  }

  std::vector<index_t> adj_ptr(n + 1);
  adj_ptr[0] = 0;
  for (index_t i = 0; i < n; ++i) {
    adj_ptr[i + 1] = adj_ptr[i] + degrees[i];
  }

  std::vector<index_t> adj_list(adj_ptr[n]);
  std::vector<index_t> pos(n, 0);

  for (index_t i = 0; i < n; ++i) {
    for (auto j : matrix.row_indices(i)) {
      if (i != j) {
        if (i < j || matrix.is_symmetric()) {
          index_t idx_i = adj_ptr[i] + pos[i]++;
          index_t idx_j = adj_ptr[j] + pos[j]++;
          adj_list[idx_i] = j;
          adj_list[idx_j] = i;
        }
      }
    }
  }

  for (index_t i = 0; i < n; ++i) {
    std::sort(adj_list.begin() + adj_ptr[i], adj_list.begin() + adj_ptr[i + 1]);
  }

  return Graph(n, n_edges, std::move(adj_ptr), std::move(adj_list));
}

bool Graph::validate() const {
  if (adj_ptr_.size() != static_cast<size_t>(n_vertices_ + 1)) return false;
  if (adj_ptr_[n_vertices_] != static_cast<index_t>(adj_list_.size()))
    return false;

  for (index_t i = 0; i < n_vertices_; ++i) {
    if (adj_ptr_[i] > adj_ptr_[i + 1]) return false;
  }

  for (index_t i = 0; i < n_vertices_; ++i) {
    auto neighs = neighbors(i);
    for (size_t k = 0; k < neighs.size(); ++k) {
      if (neighs[k] < 0 || neighs[k] >= n_vertices_) return false;
      if (neighs[k] == i) return false;
      if (k > 0 && neighs[k] <= neighs[k - 1]) return false;
    }
  }

  return true;
}

std::vector<index_t> Graph::compute_degeneracy_ordering() const {
  // Matula-Beck algorithm for degeneracy ordering
  std::vector<index_t> ordering;
  ordering.reserve(n_vertices_);

  std::vector<bool> removed(n_vertices_, false);
  std::vector<index_t> current_degrees(n_vertices_);

  for (index_t i = 0; i < n_vertices_; ++i) {
    current_degrees[i] = degree(i);
  }

  index_t max_degree = 0;
  for (auto d : current_degrees) {
    max_degree = std::max(max_degree, d);
  }

  std::vector<std::set<index_t>> buckets(max_degree + 1);
  for (index_t i = 0; i < n_vertices_; ++i) {
    buckets[current_degrees[i]].insert(i);
  }

  for (index_t iter = 0; iter < n_vertices_; ++iter) {
    index_t min_deg = 0;
    while (min_deg <= max_degree && buckets[min_deg].empty()) {
      min_deg++;
    }

    if (min_deg > max_degree) break;

    index_t v = *buckets[min_deg].begin();
    buckets[min_deg].erase(buckets[min_deg].begin());
    removed[v] = true;
    ordering.push_back(v);

    for (auto u : neighbors(v)) {
      if (!removed[u]) {
        index_t old_deg = current_degrees[u];
        buckets[old_deg].erase(u);
        current_degrees[u]--;
        buckets[current_degrees[u]].insert(u);
      }
    }
  }

  return ordering;
}

Graph Graph::induced_subgraph(const std::vector<index_t> &vertices) const {
  std::unordered_map<index_t, index_t> old_to_new;
  for (size_t i = 0; i < vertices.size(); ++i) {
    old_to_new[vertices[i]] = i;
  }

  index_t new_n = vertices.size();
  std::vector<index_t> new_degrees(new_n, 0);

  index_t new_m = 0;
  for (auto v : vertices) {
    for (auto u : neighbors(v)) {
      if (old_to_new.count(u)) {
        new_degrees[old_to_new[v]]++;
        new_m++;
      }
    }
  }
  new_m /= 2;

  std::vector<index_t> new_adj_ptr(new_n + 1);
  new_adj_ptr[0] = 0;
  for (index_t i = 0; i < new_n; ++i) {
    new_adj_ptr[i + 1] = new_adj_ptr[i] + new_degrees[i];
  }

  std::vector<index_t> new_adj_list(new_adj_ptr[new_n]);
  std::vector<index_t> pos(new_n, 0);

  for (size_t i = 0; i < vertices.size(); ++i) {
    index_t v = vertices[i];
    for (auto u : neighbors(v)) {
      if (old_to_new.count(u)) {
        index_t new_u = old_to_new[u];
        index_t idx = new_adj_ptr[i] + pos[i]++;
        new_adj_list[idx] = new_u;
      }
    }
  }

  for (index_t i = 0; i < new_n; ++i) {
    std::sort(new_adj_list.begin() + new_adj_ptr[i],
              new_adj_list.begin() + new_adj_ptr[i + 1]);
  }

  return Graph(new_n, new_m, std::move(new_adj_ptr), std::move(new_adj_list));
}

}  // namespace hypergraph_reorder
