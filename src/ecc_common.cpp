// ecc_common.cpp - Shared Phase-2/Phase-3 ECC machinery + solver factory
//
// select_covering_cliques and cover_remaining_edges are moved verbatim from
// clique_cover.cpp (spec §2 / §4.7); the only sanctioned changes are:
//   1. select_covering_cliques additionally exposes its covered-edge set via
//      an out-parameter so Phase 3 no longer rebuilds coverage from scratch.
//   2. cover_remaining_edges consumes that populated covered-edge set
//      directly (its former covered_edges parameter was unused).
#include "hypergraph_reorder/ecc.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

namespace hypergraph_reorder {

namespace detail {

std::vector<std::vector<index_t>> select_covering_cliques(
    const Graph &graph,
    const std::vector<std::vector<index_t>> &maximal_cliques,
    std::set<std::pair<index_t, index_t>> &covered_edges) {
  // Sort cliques by size (largest first)
  std::vector<std::pair<index_t, index_t>> clique_sizes;
  for (size_t i = 0; i < maximal_cliques.size(); ++i) {
    clique_sizes.emplace_back(maximal_cliques[i].size(), i);
  }
  std::sort(clique_sizes.rbegin(), clique_sizes.rend());

  // Greedy selection
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

void cover_remaining_edges(
    const Graph &graph,
    const std::set<std::pair<index_t, index_t>> &covered_edges,
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
}

}  // namespace detail

std::unique_ptr<EccSolver> make_ecc_solver(const EccOptions &opts) {
  switch (opts.method) {
    case EccMethod::kEdgeCover:
      return std::make_unique<EdgeCoverEccSolver>(opts);
    case EccMethod::kTriEnum:
      return std::make_unique<TriangleEnumEccSolver>(opts);
    case EccMethod::kBk:
    default:
      return std::make_unique<BronKerboschEccSolver>(opts);
  }
}

}  // namespace hypergraph_reorder
