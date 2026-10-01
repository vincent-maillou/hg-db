// ecc_tri.cpp - Triangle-enumeration ECC solver (triangle candidates)
//
// enumerate_triangles is moved verbatim from clique_cover.cpp (spec §2 /
// §4.7) — relocated, not rewritten. The fix_core_algs branch fixes the
// algorithm body; do not edit it here.
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

// Phase 1a: Fast triangle enumeration (for large graphs)
std::vector<std::vector<index_t>> enumerate_triangles(
    const Graph &graph, int num_threads, bool use_parallel) {
  std::vector<std::vector<std::vector<index_t>>> thread_triangles;

#ifdef USE_OPENMP
  int omp_threads =
      num_threads > 0 ? num_threads : omp_get_max_threads();
  omp_set_num_threads(omp_threads);
  thread_triangles.resize(omp_threads);
#else
  (void)num_threads;
  (void)use_parallel;
  thread_triangles.resize(1);
#endif

  // Parallel triangle enumeration
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

}  // namespace

// Triangle-enumeration ECC solver: Phase 1 enumerates all triangles; Phase
// 2/3 are the shared greedy selection and 2-clique remainder.
TriangleEnumEccSolver::TriangleEnumEccSolver(const EccOptions& opts)
    : opts_(opts) {}

CliqueCover TriangleEnumEccSolver::solve(const Graph& graph) {
  Timer timer;

  if (!opts_.suppress_output)
    std::cout << "Enumerating triangles..." << std::endl;
  auto all_triangles =
      enumerate_triangles(graph, opts_.num_threads, opts_.use_parallel);
  if (!opts_.suppress_output)
    std::cout << "Found " << all_triangles.size() << " triangles"
              << std::endl;

  // Greedily select triangles to cover edges
  if (!opts_.suppress_output)
    std::cout << "Selecting covering triangles..." << std::endl;
  std::set<std::pair<index_t, index_t>> covered;
  std::vector<std::vector<index_t>> all_cliques =
      detail::select_covering_cliques(graph, all_triangles, covered);

  // Phase 3: Cover remaining edges with 2-cliques
  detail::cover_remaining_edges(graph, covered, all_cliques);

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << all_cliques.size() << std::endl;

  stats_.time_ecc_ms = timer.elapsed_ms();
  stats_.n_cliques = all_cliques.size();

  return CliqueCover(graph.n_vertices(), std::move(all_cliques));
}

const Statistics& TriangleEnumEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
