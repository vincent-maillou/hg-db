// ecc_edge.cpp - Edge-cover ECC solver (every edge becomes a 2-clique)
//
// This is the mandatory ablation baseline: it bypasses candidate enumeration
// and Phase 2 entirely, isolating the benefit of higher-order clique
// structure for the partitioner. Equivalent to partitioning the standard
// edge-hypergraph (incidence) formulation through the same CNH machinery.
// Zero enumeration cost.
#include "hypergraph_reorder/ecc.hpp"

#include <iostream>
#include <utility>
#include <vector>

namespace hypergraph_reorder {

// Edge-cover ECC solver: no candidates, Phase 3 over all edges only.
EdgeCoverEccSolver::EdgeCoverEccSolver(const EccOptions& opts)
    : opts_(opts) {}

CliqueCover EdgeCoverEccSolver::solve(const Graph& graph) {
  // Stage 4 of 8 (v5): the whole ECC work is the Phase-3 remainder walk
  // plus the cover CSR construction (no enumeration, no greedy
  // selection). time_degeneracy_ms and time_enumeration_ms stay 0.0;
  // degeneracy stays -1 (spec v5 D1/D2).
  Timer greedy_timer;

  // All edges are uncovered by construction — pass a zeroed coverage
  // bitmap (the pre-fix behavior, deliberately: every edge becomes a
  // 2-clique either way). Each remainder 2-clique covers exactly one
  // edge, so the recorded mass is all-ones with n_edges entries.
  std::vector<uint64_t> covered_bits(
      (static_cast<size_t>(graph.n_edges()) + 63) / 64, 0);
  std::vector<std::vector<index_t>> all_cliques;
  std::vector<index_t> masses;
  detail::cover_remaining_edges(graph, covered_bits, all_cliques, masses);

  // Cover construction (and the mass move) is timed inside stage 4
  // (WP2 D3).
  CliqueCover cover(graph.n_vertices(), std::move(all_cliques));
  cover.set_edge_mass(std::move(masses));
  stats_.time_greedy_ms = greedy_timer.elapsed_ms();

  if (!opts_.suppress_output)
    std::cout << "Total cliques: " << cover.n_cliques() << std::endl;

  return cover;
}

const Statistics& EdgeCoverEccSolver::get_stats() const { return stats_; }

}  // namespace hypergraph_reorder
