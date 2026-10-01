// ecc_edge.cpp - Edge-cover ECC solver (every edge becomes a 2-clique)
//
// This is the mandatory ablation baseline: it bypasses candidate enumeration
// and Phase 2 entirely, isolating the benefit of higher-order clique
// structure for the partitioner. Equivalent to partitioning the standard
// edge-hypergraph (incidence) formulation through the same CNH machinery.
// Zero enumeration cost.
#include "hypergraph_reorder/ecc.hpp"

#include <iostream>
#include <set>
#include <utility>
#include <vector>

namespace hypergraph_reorder {

// Edge-cover ECC solver: no candidates, Phase 3 over all edges only.
EdgeCoverEccSolver::EdgeCoverEccSolver(const EccOptions& opts)
    : opts_(opts) {}

CliqueCover EdgeCoverEccSolver::solve(const Graph& graph) {
  Timer timer;

  // All edges are uncovered by construction — pass an empty covered set
  // (the pre-fix behavior, deliberately: every edge becomes a 2-clique
  // either way).
  std::set<std::pair<index_t, index_t>> covered;
  std::vector<std::vector<index_t>> all_cliques;
  detail::cover_remaining_edges(graph, covered, all_cliques);

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << all_cliques.size() << std::endl;

  stats_.time_ecc_ms = timer.elapsed_ms();
  stats_.n_cliques = all_cliques.size();

  return CliqueCover(graph.n_vertices(), std::move(all_cliques));
}

const Statistics& EdgeCoverEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
