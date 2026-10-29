// ecc_tri.cpp - Triangle-enumeration ECC solver (triangle candidates)
//
// Phase 1 enumerates all triangles via parallel two-pointer adjacency
// intersection; Phase 2/3 are the shared greedy selection and 2-clique
// remainder in ecc_common.cpp.
#include "hypergraph_reorder/ecc.hpp"

#include <algorithm>
#include <iostream>
#include <utility>
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
  // Stage 3 of 8 (v5): triangle enumeration + thread merge. (No degeneracy
  // stage for TriEnum: degeneracy stays -1, time_degeneracy_ms stays 0.0,
  // and no enumeration histogram is recorded — owner decision, spec v5 D2.)
  if (!opts_.suppress_output)
    std::cout << "Enumerating triangles..." << std::endl;
  Timer enumeration_timer;
  auto all_triangles =
      enumerate_triangles(graph, opts_.num_threads, opts_.use_parallel);
  stats_.time_enumeration_ms = enumeration_timer.elapsed_ms();
  if (!opts_.suppress_output)
    std::cout << "Found " << all_triangles.size() << " triangles"
              << std::endl;

  // Stage 4 of 8 (v5): greedy Phase 2 (selection) + Phase 3 (remainder),
  // then cover CSR construction — all inside the stage-4 timing window
  // (WP2 D3).
  if (!opts_.suppress_output)
    std::cout << "Selecting covering triangles..." << std::endl;
  Timer greedy_timer;
  std::vector<uint64_t> covered_bits;
  std::vector<index_t> masses;
  std::vector<std::vector<index_t>> all_cliques =
      detail::select_covering_cliques(graph, all_triangles, covered_bits, masses);

  // Phase 3: Cover remaining edges with 2-cliques
  detail::cover_remaining_edges(graph, covered_bits, all_cliques, masses);

  // Cover construction (and the mass move) is timed inside stage 4
  // (WP2 D3): the CSR build is O(n + K) and was the dominant part of the
  // v5 untimed total-minus-stages gap.
  CliqueCover cover(graph.n_vertices(), std::move(all_cliques));
  cover.set_edge_mass(std::move(masses));
  stats_.time_greedy_ms = greedy_timer.elapsed_ms();

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << cover.n_cliques() << std::endl;

  return cover;
}

const Statistics& TriangleEnumEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
