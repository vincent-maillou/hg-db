#include "hypergraph_reorder/graph.hpp"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace hypergraph_reorder {

Graph::Graph() : n_vertices_(0), n_edges_(0) {}

Graph::Graph(index_t n_vertices) : n_vertices_(n_vertices), n_edges_(0) {
  adj_ptr_.resize(n_vertices + 1, 0);
}

Graph::Graph(
  index_t n_vertices, 
  index_t n_edges, 
  std::vector<index_t> adj_ptr,
  std::vector<index_t> adj_list
) : n_vertices_(n_vertices),
    n_edges_(n_edges),
    adj_ptr_(std::move(adj_ptr)),
    adj_list_(std::move(adj_list))
{
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

  // Compute degrees and number of edges
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

  // Build adjacency structure
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

  // Sort adjacency lists to ensure neighbors are in ascending order
  for (index_t i = 0; i < n; ++i) {
    std::sort(adj_list.begin() + adj_ptr[i], adj_list.begin() + adj_ptr[i + 1]);
  }

  // Canonicalize to a simple graph: deduplicate adjacency entries (a
  // full-storage symmetric CSR inserts every neighbor twice) and drop any
  // diagonal entry that may have slipped in, then recompute the edge count
  // from the surviving adjacency (spec §4.4 / D5).
  std::vector<index_t> canonical_ptr(n + 1, 0);
  std::vector<index_t> canonical_list;
  canonical_list.reserve(adj_list.size());
  for (index_t i = 0; i < n; ++i) {
    auto begin = adj_list.begin() + adj_ptr[i];
    auto end = adj_list.begin() + adj_ptr[i + 1];
    auto last = std::unique(begin, end);
    for (auto it = begin; it != last; ++it) {
      if (*it != i) {
        canonical_list.push_back(*it);
      }
    }
    canonical_ptr[i + 1] = static_cast<index_t>(canonical_list.size());
  }
  index_t total_deg = static_cast<index_t>(canonical_list.size());
  n_edges = total_deg / 2;  // every undirected edge contributes 2 to the sum

  Graph graph(n, n_edges, std::move(canonical_ptr), std::move(canonical_list));

#ifndef NDEBUG
  if (!graph.validate()) {
    throw HypergraphReorderError(
        "Graph::from_symmetric_matrix produced an invalid graph");
  }
#endif

  return graph;
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

index_t Graph::compute_degeneracy_ordering(std::vector<index_t> &ordering) const {
  /* Matula-Beck algorithm for degeneracy ordering (smallest-last ordering).

  Linear-time implementation (spec v5 D4): degree-count buckets over an
  intrusive doubly-linked vertex list (pos[]/prev[]/next[], one head per
  bucket). Each removal pops a vertex of minimum current degree, decrements
  the live degrees of the not-yet-removed neighbors, and relinks each one
  exactly one bucket down. Because the current minimum degree can drop by at
  most 1 per removal, the min-bucket cursor resumes its forward scan one
  bucket below the last removal — total scan work O(n + max_deg), overall
  O(n + m) time and O(n + max_deg) memory.

  The within-bucket extraction order is arbitrary (FIFO here). This is safe:
  the degeneracy value d is tie-break invariant, the *set* of maximal cliques
  found downstream (ELS/Bron-Kerbosch) is invariant to the ordering, and the
  v4.1.0 Phase-2 canonicalization makes the final cover a pure function of
  that candidate set (spec v5 D4).

  References:
  -----------
  Matula-Beck algorithm:
  - Matula, D. W., & Beck, L. L. (1983). Smallest-last ordering and clustering and graph coloring algorithms. Journal of the ACM (JACM), 30(3), 417–427. https://doi.org/10.1145/2402.322385

  Usage of Matula-Beck algorithm for degeneracy ordering:
  - Eppstein, D., Löffler, M., & Strash, D. (2013). Listing All Maximal Cliques in Large Sparse Real-World Graphs. Journal of Experimental Algorithmics (JEA), 18, 3.1:3.1-3.1:3.21. https://doi.org/10.1145/2543629

  */

  const index_t n = n_vertices_;
  ordering.clear();
  ordering.reserve(static_cast<size_t>(n));

  // Live degrees (current degree within the not-yet-removed subgraph).
  std::vector<index_t> current_degrees(n);
  index_t max_degree = 0;
  for (index_t v = 0; v < n; ++v) {
    current_degrees[v] = degree(v);
    max_degree = std::max(max_degree, current_degrees[v]);
  }

  // Intrusive doubly-linked lists, one per degree bucket:
  //   pos[v]        = bucket currently holding v,
  //   prev[v]/next[v] = neighbors of v within that bucket's list,
  //   head[d]       = first vertex of bucket d (-1 = empty).
  std::vector<index_t> pos(n);
  std::vector<index_t> prev(n, -1);
  std::vector<index_t> next(n, -1);
  std::vector<index_t> head(static_cast<size_t>(max_degree) + 1, -1);

  // Append v to the tail of bucket d (FIFO extraction order).
  auto list_append = [&](index_t d, index_t v) {
    pos[v] = d;
    prev[v] = -1;
    next[v] = head[d];
    if (head[d] != -1) {
      prev[head[d]] = v;
    }
    head[d] = v;
  };

  // Detach v from its bucket's list.
  auto list_remove = [&](index_t v) {
    const index_t p = prev[v];
    const index_t nx = next[v];
    if (p != -1) {
      next[p] = nx;
    } else {
      head[pos[v]] = nx;
    }
    if (nx != -1) {
      prev[nx] = p;
    }
    prev[v] = -1;
    next[v] = -1;
  };

  std::vector<bool> removed(n, false);
  for (index_t v = 0; v < n; ++v) {
    list_append(current_degrees[v], v);
  }

  index_t degeneracy = 0;
  // Cursor into the bucket array. Invariant: the current minimum live
  // degree never drops by more than 1 per removal (every surviving vertex
  // had degree >= d before the removal at d; a neighbor of the removed
  // vertex loses exactly 1). So after each removal it suffices to resume
  // scanning at d - 1 — the scan only ever moves forward from there, which
  // keeps the total scan work O(n + max_deg).
  index_t min_bucket = 0;

  for (index_t iter = 0; iter < n; ++iter) {
    if (iter > 0 && min_bucket > 0) {
      --min_bucket;  // the last removal may have opened bucket d - 1
    }
    while (min_bucket <= max_degree && head[min_bucket] == -1) {
      ++min_bucket;
    }
    if (min_bucket > max_degree) break;  // every vertex removed

    // Pop a vertex of minimum current degree (FIFO: the bucket head).
    const index_t v = head[min_bucket];
    list_remove(v);
    removed[v] = true;
    degeneracy = std::max(degeneracy, min_bucket);
    ordering.push_back(v);

    // Decrement the live degree of each surviving neighbor, relinking it
    // exactly one bucket down.
    for (auto u : neighbors(v)) {
      if (removed[u]) continue;
      list_remove(u);
      --current_degrees[u];
      list_append(current_degrees[u], u);
    }
  }

  return degeneracy;
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
