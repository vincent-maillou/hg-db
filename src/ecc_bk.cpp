// ecc_bk.cpp - Bron-Kerbosch ECC solver (maximal-clique candidates)
//
// find_maximal_cliques and bron_kerbosch_pivot are moved verbatim from
// clique_cover.cpp (spec §2 / §4.7) — relocated, not rewritten. The
// fix_core_algs branch fixes the algorithm bodies; do not edit them here.
#include "hypergraph_reorder/ecc.hpp"

#include <algorithm>
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
                         std::vector<std::vector<index_t>> &cliques) {
  if (P.empty() && X.empty()) {
    // R is a maximal clique
    if (R.size() >= 2) {  // Only keep cliques with at least 2 vertices
      cliques.push_back(R);
    }
    return;
  }

  if (P.empty()) return;

  // Choose pivot from P ∪ X with maximum degree in P
  index_t pivot = P[0];
  index_t max_connections = 0;

  for (auto u : P) {
    index_t connections = 0;
    for (auto v : graph.neighbors(u)) {
      if (std::find(P.begin(), P.end(), v) != P.end()) {
        connections++;
      }
    }
    if (connections > max_connections) {
      max_connections = connections;
      pivot = u;
    }
  }

  // Try vertices in P \ neighbors(pivot)
  std::set<index_t> pivot_neighbors;
  for (auto u : graph.neighbors(pivot)) {
    pivot_neighbors.insert(u);
  }

  std::vector<index_t> candidates;
  for (auto v : P) {
    if (pivot_neighbors.find(v) == pivot_neighbors.end()) {
      candidates.push_back(v);
    }
  }

  for (auto v : candidates) {
    // R' = R ∪ {v}
    R.push_back(v);

    // P' = P ∩ neighbors(v)
    std::vector<index_t> P_new;
    std::set<index_t> v_neighbors;
    for (auto u : graph.neighbors(v)) {
      v_neighbors.insert(u);
    }
    for (auto u : P) {
      if (v_neighbors.count(u)) {
        P_new.push_back(u);
      }
    }

    // X' = X ∩ neighbors(v)
    std::vector<index_t> X_new;
    for (auto u : X) {
      if (v_neighbors.count(u)) {
        X_new.push_back(u);
      }
    }

    // Recursive call
    bron_kerbosch_pivot(graph, R, P_new, X_new, cliques);

    // Backtrack
    R.pop_back();

    // Move v from P to X
    P.erase(std::remove(P.begin(), P.end(), v), P.end());
    X.push_back(v);
  }
}

// Phase 1b: Find maximal cliques using parallel Bron-Kerbosch
std::vector<std::vector<index_t>> find_maximal_cliques(
    const Graph &graph, int num_threads, bool use_parallel) {
  // Compute degeneracy ordering
  auto ordering = graph.compute_degeneracy_ordering();

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

#ifdef USE_OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (size_t idx = 0; idx < ordering.size(); ++idx) {
      index_t v = ordering[idx];

      // R = {v}
      std::vector<index_t> R = {v};

      // P = neighbors of v that come after v in ordering
      std::vector<index_t> P;
      for (auto u : graph.neighbors(v)) {
        if (u > v) {
          P.push_back(u);
        }
      }

      // X = empty (neighbors before v are already processed)
      std::vector<index_t> X;

      // Run Bron-Kerbosch from this starting point
      bron_kerbosch_pivot(graph, R, P, X, thread_cliques[tid]);
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
  Timer timer;

  if (!opts_.suppress_output)
    std::cout << "Finding maximal cliques..." << std::endl;
  std::vector<std::vector<index_t>> all_cliques =
      find_maximal_cliques(graph, opts_.num_threads, opts_.use_parallel);
  if (!opts_.suppress_output)
    std::cout << "Found " << all_cliques.size() << " maximal cliques"
              << std::endl;

  // Phase 2: Greedily select cliques to cover edges
  if (!opts_.suppress_output)
    std::cout << "Selecting covering cliques..." << std::endl;
  std::set<std::pair<index_t, index_t>> covered;
  all_cliques = detail::select_covering_cliques(graph, all_cliques, covered);

  // Phase 3: Cover remaining edges with 2-cliques
  detail::cover_remaining_edges(graph, covered, all_cliques);

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << all_cliques.size() << std::endl;

  stats_.time_ecc_ms = timer.elapsed_ms();
  stats_.n_cliques = all_cliques.size();

  return CliqueCover(graph.n_vertices(), std::move(all_cliques));
}

const Statistics& BronKerboschEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
