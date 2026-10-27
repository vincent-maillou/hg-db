// ecc_bk.cpp - Bron-Kerbosch ECC solver (maximal-clique candidates)
//
// Phase 1 enumerates all maximal cliques (size >= 2) via a degeneracy-
// ordered parallel Bron-Kerbosch with pivoting (Eppstein-Löffler-Strash
// outer loop, Tomita pivot rule); Phase 2/3 are the shared greedy
// selection and 2-clique remainder in ecc_common.cpp.
#include "hypergraph_reorder/ecc.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace hypergraph_reorder {

namespace {

// Bron-Kerbosch algorithm with pivoting
void bron_kerbosch_pivot(const Graph &graph,
                         std::vector<index_t> &R,  // Current clique
                         std::vector<index_t> &P,  // Candidates
                         std::vector<index_t> &X,  // Already processed
                         std::vector<std::vector<index_t>> &cliques,  // Output
                         std::vector<uint64_t> &p_marker,  // Marker array for
                                                           // P membership
                                                           // (timestamping)
                         uint64_t &marker_id  // Current marker ID
) {
  /* Bron-Kerbosch algorithm for maximal cliques enumeration
  - Cariant with pivoting
  
  References:
  -----------
  - Bron, C., & Kerbosch, J. (1973). Algorithm 457: Finding all cliques of an undirected graph. Communications of the ACM, 16(9), 575–577. https://doi.org/10.1145/362342.362367
  - Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time complexity for generating all maximal cliques and computational experiments. Theoretical Computer Science, Computing and Combinatorics, 363(1), 28–42. https://doi.org/10.1016/j.tcs.2006.06.015
  - Cazals, F., & Karande, C. (2008). A note on the problem of reporting maximal cliques. Theoretical Computer Science, 407(1), 564–568. https://doi.org/10.1016/j.tcs.2008.05.010

  
  */


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

    while (n_idx < pivot_neighbors.size() && pivot_neighbors[n_idx] < v) {
      ++n_idx;
    }

    if (n_idx == pivot_neighbors.size() || pivot_neighbors[n_idx] != v) {
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
    bron_kerbosch_pivot(graph, R, P_new, X_new, cliques, p_marker, marker_id);

    // Backtrack
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

// Phase 1b: Find maximal cliques using parallel Bron-Kerbosch over the
// given degeneracy ordering. `ordering` and the inverse permutation `rank`
// are computed by the caller (timed as the degeneracy stage there):
//   ordering[i] = vertex at position i of the degeneracy ordering,
//   rank[v]     = position of vertex v in that ordering.
std::vector<std::vector<index_t>> find_maximal_cliques(
    const Graph &graph, const std::vector<index_t> &ordering,
    const std::vector<index_t> &rank, int num_threads, bool use_parallel) {
  // Thread-local clique storage
  std::vector<std::vector<std::vector<index_t>>> thread_cliques;

#ifdef USE_OPENMP
  int omp_threads =
      num_threads > 0 ? num_threads : omp_get_max_threads();
  omp_set_num_threads(omp_threads);
  thread_cliques.resize(omp_threads);
#else
  (void)num_threads;
  (void)use_parallel;
  thread_cliques.resize(1);
#endif

  // Parallel Bron-Kerbosch over degeneracy ordering
#ifdef USE_OPENMP
#pragma omp parallel if (use_parallel)
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
      bron_kerbosch_pivot(graph, R, P, X, thread_cliques[tid], p_marker,
                          marker_id);
    }
  }

  // Merge thread-local results
  std::vector<std::vector<index_t>> all_cliques;
  for (const auto &tc : thread_cliques) {
    all_cliques.insert(all_cliques.end(), tc.begin(), tc.end());
  }

  return all_cliques;
}

}  // namespace

// Bron-Kerbosch ECC solver: Phase 1 enumerates maximal cliques (size >= 2)
// via degeneracy-ordered parallel Bron-Kerbosch with pivoting; Phase 2/3 are
// the shared greedy selection and 2-clique remainder.
BronKerboschEccSolver::BronKerboschEccSolver(const EccOptions& opts)
    : opts_(opts) {}

CliqueCover BronKerboschEccSolver::solve(const Graph& graph) {
  // Stage 2 of 8 (v5): degeneracy ordering + rank[] construction.
  Timer degeneracy_timer;

  std::vector<index_t> ordering;
  const index_t d = graph.compute_degeneracy_ordering(ordering);

  // rank[v] = position of vertex v in the degeneracy ordering. Vertex IDs
  // are not assumed to correspond to the degeneracy ordering.
  std::vector<index_t> rank(graph.n_vertices());
  for (index_t i = 0; i < static_cast<index_t>(ordering.size()); ++i) {
    rank[ordering[i]] = i;
  }

  stats_.degeneracy = d;
  stats_.time_degeneracy_ms = degeneracy_timer.elapsed_ms();

  // Stage 3 of 8 (v5): parallel Bron-Kerbosch enumeration + thread merge.
  if (!opts_.suppress_output)
    std::cout << "Finding maximal cliques..." << std::endl;
  Timer enumeration_timer;
  std::vector<std::vector<index_t>> all_cliques = find_maximal_cliques(
      graph, ordering, rank, opts_.num_threads, opts_.use_parallel);
  stats_.time_enumeration_ms = enumeration_timer.elapsed_ms();
  if (!opts_.suppress_output)
    std::cout << "Found " << all_cliques.size() << " maximal cliques"
              << std::endl;

  // Capture the enumeration histogram (v5, D2) BEFORE the greedy selection
  // consumes the candidate vector: one O(total members) pass, no retention.
  for (const auto& c : all_cliques) {
    ++stats_.maximal_clique_enumeration_histogram[static_cast<index_t>(c.size())];
  }
  stats_.n_maximal_cliques = static_cast<index_t>(all_cliques.size());

  // Stage 4 of 8 (v5): greedy Phase 2 (selection) + Phase 3 (remainder).
  if (!opts_.suppress_output)
    std::cout << "Selecting covering cliques..." << std::endl;
  Timer greedy_timer;
  std::set<std::pair<index_t, index_t>> covered;
  std::vector<index_t> masses;
  all_cliques =
      detail::select_covering_cliques(graph, all_cliques, covered, masses);

  // Phase 3: Cover remaining edges with 2-cliques
  detail::cover_remaining_edges(graph, covered, all_cliques, masses);
  stats_.time_greedy_ms = greedy_timer.elapsed_ms();

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << all_cliques.size() << std::endl;

  CliqueCover cover(graph.n_vertices(), std::move(all_cliques));
  cover.set_edge_mass(std::move(masses));
  return cover;
}

const Statistics& BronKerboschEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
