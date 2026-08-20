#include "hypergraph_reorder/clique_cover.hpp"

#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_set>
#include <iterator>

#ifdef USE_OPENMP
#include <omp.h>
#endif

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

  // Build vertex-centric representation (vertex -> cliques)
  std::vector<std::vector<index_t>> vertex_to_cliques(n_vertices);

  for (index_t cid = 0; cid < n_cliques_; ++cid) {
    for (auto v : get_clique(cid)) {
      vertex_to_cliques[v].push_back(cid);
    }
  }

  // Convert to compressed format
  vertex_cliques_ptr_.resize(n_vertices + 1);
  vertex_cliques_ptr_[0] = 0;

  for (index_t v = 0; v < n_vertices; ++v) {
    vertex_cliques_ptr_[v + 1] =
        vertex_cliques_ptr_[v] + vertex_to_cliques[v].size();
  }

  vertex_cliques_.reserve(vertex_cliques_ptr_[n_vertices]);
  for (const auto &vec : vertex_to_cliques) {
    vertex_cliques_.insert(vertex_cliques_.end(), vec.begin(), vec.end());
  }
}

double CliqueCover::avg_clique_size() const {
  if (n_cliques_ == 0) return 0.0;
  return static_cast<double>(clique_members_.size()) / n_cliques_;
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

// ===== CliqueCoverSolver Implementation =====

CliqueCoverSolver::CliqueCoverSolver(const Options &opts) : opts_(opts) {}

CliqueCover CliqueCoverSolver::solve(const Graph &graph) {
  Timer timer;

  std::vector<std::vector<index_t>> all_cliques;

  // Phase 1: Find cliques (strategy depends on graph size and options)
  if (opts_.use_maximal_cliques && !opts_.greedy_only) {
    bool use_bk = graph.n_vertices() <= opts_.max_clique_enum_vertices &&
                  graph.n_edges() <= opts_.max_clique_enum_edges;
    if (use_bk) {
      if (!opts_.suppress_output)
        std::cout << "Finding maximal cliques..." << std::endl;
      all_cliques = find_maximal_cliques(graph);
      if (!opts_.suppress_output)
        std::cout << "Found " << all_cliques.size() << " maximal cliques"
                  << std::endl;

      // Phase 2: Greedily select cliques to cover edges
      if (!opts_.suppress_output)
        std::cout << "Selecting covering cliques..." << std::endl;
      all_cliques = select_covering_cliques(graph, all_cliques);
    } else {
      // For large/dense graphs, use greedy triangle selection (much faster)
      if (!opts_.suppress_output)
        std::cout << "Graph too large for full maximal clique enumeration"
                  << std::endl;
      if (!opts_.suppress_output)
        std::cout << "Using greedy triangle selection instead..." << std::endl;
      auto all_triangles = enumerate_triangles(graph);
      if (!opts_.suppress_output)
        std::cout << "Found " << all_triangles.size() << " triangles"
                  << std::endl;

      // Greedily select triangles to cover edges
      if (!opts_.suppress_output)
        std::cout << "Selecting covering triangles..." << std::endl;
      all_cliques = select_covering_cliques(graph, all_triangles);
    }
  }

  // Phase 3: Cover remaining edges with 2-cliques
  std::vector<bool> covered_edges(graph.n_edges(), false);

  // Mark covered edges
  for (const auto &clique : all_cliques) {
    for (size_t i = 0; i < clique.size(); ++i) {
      for (size_t j = i + 1; j < clique.size(); ++j) {
        index_t u = clique[i];
        index_t v = clique[j];
        if (graph.has_edge(u, v)) {
          // Mark edge as covered (simplified tracking)
        }
      }
    }
  }

  // Add 2-cliques for any remaining edges
  cover_remaining_edges(graph, covered_edges, all_cliques);

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << all_cliques.size() << std::endl;

  stats_.time_clique_cover_ms = timer.elapsed_ms();
  stats_.n_cliques = all_cliques.size();

  return CliqueCover(graph.n_vertices(), std::move(all_cliques));
}

std::vector<std::vector<index_t>> CliqueCoverSolver::enumerate_triangles(
    const Graph &graph) {
  std::vector<std::vector<std::vector<index_t>>> thread_triangles;

#ifdef USE_OPENMP
  int num_threads =
      opts_.num_threads > 0 ? opts_.num_threads : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  thread_triangles.resize(num_threads);
#else
  thread_triangles.resize(1);
#endif

  // Parallel triangle enumeration
#ifdef USE_OPENMP
#pragma omp parallel if (opts_.use_parallel)
#endif
  {
#ifdef USE_OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

#ifdef USE_OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
    for (index_t u = 0; u < graph.n_vertices(); ++u) {
      auto u_neighbors = graph.neighbors(u);

      for (auto v : u_neighbors) {
        if (v <= u) continue;  // Avoid duplicates, ensure u < v

        auto v_neighbors = graph.neighbors(v);

        // Two-pointer intersection of u_neighbors and v_neighbors for w > v.
        // Both lists are sorted, so this is O(deg(u) + deg(v)) per edge.
        auto u_it = std::upper_bound(u_neighbors.begin(), u_neighbors.end(), v);
        auto v_it = std::upper_bound(v_neighbors.begin(), v_neighbors.end(), v);

        while (u_it != u_neighbors.end() && v_it != v_neighbors.end()) {
          if (*u_it == *v_it) {
            thread_triangles[tid].push_back({u, v, *u_it});
            ++u_it;
            ++v_it;
          } else if (*u_it < *v_it) {
            ++u_it;
          } else {
            ++v_it;
          }
        }
      }
    }
  }

  // Merge thread-local results
  std::vector<std::vector<index_t>> all_triangles;
  for (const auto &tt : thread_triangles) {
    all_triangles.insert(all_triangles.end(), tt.begin(), tt.end());
  }

  return all_triangles;
}

std::vector<std::vector<index_t>> CliqueCoverSolver::find_maximal_cliques(
    const Graph &graph) {
  // Compute degeneracy ordering.
  //
  // ordering[i] = vertex at position i in the degeneracy ordering.
  auto ordering = graph.compute_degeneracy_ordering();

  // Build inverse permutation:
  //
  // rank[v] = position of vertex v in the degeneracy ordering.
  //
  // Vertex IDs are not assumed to correspond to the degeneracy ordering.
  std::vector<index_t> rank(graph.n_vertices());

  for (index_t i = 0; i < static_cast<index_t>(ordering.size()); ++i) {
    rank[ordering[i]] = i;
  }

  // Thread-local clique storage
  std::vector<std::vector<std::vector<index_t>>> thread_cliques;

#ifdef USE_OPENMP
  int num_threads =
      opts_.num_threads > 0 ? opts_.num_threads : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  thread_cliques.resize(num_threads);
#else
  thread_cliques.resize(1);
#endif

  // Parallel Bron-Kerbosch over degeneracy ordering
#ifdef USE_OPENMP
#pragma omp parallel if (opts_.use_parallel)
#endif
  {
#ifdef USE_OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

  std::vector<uint64_t> p_marker(graph.n_vertices(), 0);
  uint64_t marker_id = 0;

#ifdef USE_OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (size_t idx = 0; idx < ordering.size(); ++idx) {
      index_t v = ordering[idx];

      // R = {v}
      std::vector<index_t> R = {v};

      // Partition the neighbors of v according to the degeneracy ordering:
      //
      // P = later neighbors
      // X = earlier neighbors
      //
      std::vector<index_t> P;
      std::vector<index_t> X;

      for (auto u : graph.neighbors(v)) {
        if (rank[u] > rank[v]) {
          P.push_back(u);
        } else {
          X.push_back(u);
        }
      }

      // Run Bron-Kerbosch from this starting point.
      bron_kerbosch_pivot(graph, R, P, X, thread_cliques[tid], p_marker, marker_id);
    }
  }

  // Merge thread-local results
  std::vector<std::vector<index_t>> all_cliques;
  for (const auto &tc : thread_cliques) {
    all_cliques.insert(all_cliques.end(), tc.begin(), tc.end());
  }

  return all_cliques;
}

void CliqueCoverSolver::bron_kerbosch_pivot(
    const Graph &graph, 
    std::vector<index_t> &R, 
    std::vector<index_t> &P,
    std::vector<index_t> &X, 
    std::vector<std::vector<index_t>> &cliques,
    std::vector<uint64_t> &p_marker,
    uint64_t &marker_id
  ) {
  // If both sets are empty, R is a maximal clique.
  if (P.empty() && X.empty()) {
    // R is a maximal clique
    if (R.size() >= 2) {  // Only keep cliques with at least 2 vertices
      cliques.push_back(R);
    }
    return;
  }

  if (P.empty()) return;

  // Mark all vertices in P.
  // . this makes membership testing `v in P` O(1)
  ++marker_id;
  const uint64_t current_marker = marker_id;

  for (auto v : P) {
    p_marker[v] = current_marker;
  }

  // Choose pivot from P ∪ X maximizing |N(u) ∩ P|.
  //
  // Since adjacency lists are sorted, no std::set or std::find is needed.
  // We simply scan the adjacency list of each candidate and use the marker
  // array to determine whether each neighbor belongs to P.
  index_t pivot = P.front();
  index_t max_connections = 0;

  auto count_p_neighbors = [&](index_t u) -> index_t {
    index_t connections = 0;

    for (auto v : graph.neighbors(u)) {
      if (p_marker[v] == current_marker) {
        ++connections;
      }
    }

    return connections;
  };

  // Candidates from P.
  for (auto u : P) {
    index_t connections = count_p_neighbors(u);

    if (connections > max_connections) {
      max_connections = connections;
      pivot = u;
    }
  }

  // Candidates from X.
  for (auto u : X) {
    index_t connections = count_p_neighbors(u);

    if (connections > max_connections) {
      max_connections = connections;
      pivot = u;
    }
  }

  // Construct P \ N(pivot).
  //
  // P and graph.neighbors(pivot) are both sorted, so this can be computed
  // with a linear two-pointer traversal.
  std::vector<index_t> candidates;
  const auto pivot_neighbors = graph.neighbors(pivot);

  candidates.reserve(P.size());

  size_t p_idx = 0;
  size_t n_idx = 0;

  while (p_idx < P.size()) {
    const index_t v = P[p_idx];

    while (n_idx < pivot_neighbors.size() &&
           pivot_neighbors[n_idx] < v) {
      ++n_idx;
    }

    if (n_idx == pivot_neighbors.size() ||
        pivot_neighbors[n_idx] != v) {
      candidates.push_back(v);
    }

    ++p_idx;
  }

  // Process candidates.
  for (auto v : candidates) {
    R.push_back(v);

    // Compute P' = P ∩ N(v).
    //
    // Both P and N(v) are sorted.
    std::vector<index_t> P_new;
    const auto v_neighbors = graph.neighbors(v);

    P_new.reserve(std::min(P.size(), v_neighbors.size()));

    size_t i = 0;
    size_t j = 0;

    while (i < P.size() && j < v_neighbors.size()) {
      if (P[i] == v_neighbors[j]) {
        P_new.push_back(P[i]);
        ++i;
        ++j;
      } else if (P[i] < v_neighbors[j]) {
        ++i;
      } else {
        ++j;
      }
    }

    // Compute X' = X ∩ N(v).
    //
    // X is also maintained in sorted order.
    std::vector<index_t> X_new;
    X_new.reserve(std::min(X.size(), v_neighbors.size()));

    i = 0;
    j = 0;

    while (i < X.size() && j < v_neighbors.size()) {
      if (X[i] == v_neighbors[j]) {
        X_new.push_back(X[i]);
        ++i;
        ++j;
      } else if (X[i] < v_neighbors[j]) {
        ++i;
      } else {
        ++j;
      }
    }

    // Recursive call.
    bron_kerbosch_pivot(
        graph, R, P_new, X_new, cliques, p_marker, marker_id);

    R.pop_back();

    // Move v from P to X.
    //
    // We deliberately avoid std::remove/erase here. Since candidates is a
    // snapshot of the original P \ N(pivot), we can locate v using lower_bound.
    auto it = std::lower_bound(P.begin(), P.end(), v);
    if (it != P.end() && *it == v) {
      P.erase(it);
    }

    // X is sorted, so insert v while preserving sorted order.
    auto x_it = std::lower_bound(X.begin(), X.end(), v);
    X.insert(x_it, v);
  }
}

std::vector<std::vector<index_t>> CliqueCoverSolver::select_covering_cliques(
    const Graph &graph,
    const std::vector<std::vector<index_t>> &maximal_cliques) {
  // Sort cliques by size (largest first)
  std::vector<std::pair<index_t, index_t>> clique_sizes;
  for (size_t i = 0; i < maximal_cliques.size(); ++i) {
    clique_sizes.emplace_back(maximal_cliques[i].size(), i);
  }
  std::sort(clique_sizes.rbegin(), clique_sizes.rend());

  // Greedy selection
  std::set<std::pair<index_t, index_t>> covered_edges;
  std::vector<std::vector<index_t>> selected_cliques;

  for (const auto &[size, idx] : clique_sizes) {
    const auto &clique = maximal_cliques[idx];

    // Count uncovered edges in this clique
    index_t uncovered_count = 0;
    for (size_t i = 0; i < clique.size(); ++i) {
      for (size_t j = i + 1; j < clique.size(); ++j) {
        index_t u = std::min(clique[i], clique[j]);
        index_t v = std::max(clique[i], clique[j]);

        if (graph.has_edge(u, v) &&
            covered_edges.find({u, v}) == covered_edges.end()) {
          uncovered_count++;
        }
      }
    }

    // Add clique if it covers new edges
    if (uncovered_count > 0) {
      selected_cliques.push_back(clique);

      // Mark edges as covered
      for (size_t i = 0; i < clique.size(); ++i) {
        for (size_t j = i + 1; j < clique.size(); ++j) {
          index_t u = std::min(clique[i], clique[j]);
          index_t v = std::max(clique[i], clique[j]);
          if (graph.has_edge(u, v)) {
            covered_edges.insert({u, v});
          }
        }
      }
    }
  }

  return selected_cliques;
}

void CliqueCoverSolver::cover_remaining_edges(
    const Graph &graph, const std::vector<bool> &covered_edges,
    std::vector<std::vector<index_t>> &cliques) {
  // Simple approach: iterate all edges and add 2-cliques for uncovered ones
  std::set<std::pair<index_t, index_t>> edge_set;

  // Collect all edges
  for (index_t u = 0; u < graph.n_vertices(); ++u) {
    for (auto v : graph.neighbors(u)) {
      if (u < v) {
        edge_set.insert({u, v});
      }
    }
  }

  // Remove edges already covered by selected cliques
  for (const auto &clique : cliques) {
    for (size_t i = 0; i < clique.size(); ++i) {
      for (size_t j = i + 1; j < clique.size(); ++j) {
        index_t u = std::min(clique[i], clique[j]);
        index_t v = std::max(clique[i], clique[j]);
        edge_set.erase({u, v});
      }
    }
  }

  // Add 2-cliques for remaining edges
  for (const auto &[u, v] : edge_set) {
    cliques.push_back({u, v});
  }

  if (!opts_.suppress_output)
    std::cout << "Added " << edge_set.size() << " 2-cliques for uncovered edges"
              << std::endl;
}

}  // namespace hypergraph_reorder
