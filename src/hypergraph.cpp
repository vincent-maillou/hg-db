
#include "hypergraph_reorder/hypergraph.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

namespace hypergraph_reorder {

namespace {
// SCALE for vertex_share: cap the total budget W = SCALE * n_nonisolated
// at ~2^30 (int32-safe for MT-KaHyPar) while keeping per-vertex
// resolution up to 2^20 (spec D6).
index_t vertex_share_scale(index_t n_nonisolated) {
  constexpr index_t k1GiB = static_cast<index_t>(1) << 30;
  constexpr index_t k1MiB = static_cast<index_t>(1) << 20;
  return std::max<index_t>(
      1, std::min<index_t>(k1MiB, n_nonisolated > 0 ? k1GiB / n_nonisolated
                                                     : k1GiB));
}
}  // namespace

Hypergraph::Hypergraph()
    : n_nodes_(0), n_nets_(0), total_pins_(0), total_node_weight_(0) {}

std::span<const index_t> Hypergraph::net_pins(index_t net_id) const {
  if (net_id < 0 || net_id >= n_nets_) {
    throw HypergraphReorderError("Net ID out of bounds");
  }
  size_t start = edge_indices_[net_id];
  size_t end = edge_indices_[net_id + 1];
  return std::span<const index_t>(edges_.data() + start, end - start);
}

Hypergraph Hypergraph::from_clique_cover(const CliqueCover &cover,
                                         bool suppress_output,
                                         NodeWeighting weighting) {
  Hypergraph hg;
  hg.n_nodes_ = cover.n_cliques();
  index_t n_vertices = cover.n_vertices();

  if (!suppress_output)
    std::cout << "Constructing CNH: " << hg.n_nodes_ << " nodes, " << n_vertices
              << " vertices" << std::endl;

  // Pre-compute total pins for pre-allocation
  hg.total_pins_ = 0;
  for (index_t v = 0; v < n_vertices; ++v) {
    hg.total_pins_ += cover.get_vertex_cliques(v).size();
  }

  if (!suppress_output)
    std::cout << "Total pins: " << hg.total_pins_ << std::endl;

  // Pre-allocate arrays
  hg.edge_indices_.reserve(n_vertices + 1);
  hg.edges_.reserve(hg.total_pins_);
  hg.edge_weights_.reserve(n_vertices);
  hg.net_to_vertex_.reserve(n_vertices);

  hg.edge_indices_.push_back(0);
  hg.n_nets_ = 0;

  // Single-pass construction using pre-built vertex->cliques mapping
  for (index_t v = 0; v < n_vertices; ++v) {
    auto cliques = cover.get_vertex_cliques(v);

    if (!cliques.empty()) {
      // Add net (hyperedge) for this vertex
      // Deduplicate cliques (cheap safety: a no-op now that the input
      // graph is canonicalized — a vertex can no longer appear twice in
      // the same clique's adjacency, but keep it as a guard)
      std::vector<index_t> unique_cliques(cliques.begin(), cliques.end());
      std::sort(unique_cliques.begin(), unique_cliques.end());
      auto last = std::unique(unique_cliques.begin(), unique_cliques.end());
      unique_cliques.erase(last, unique_cliques.end());

      hg.edges_.insert(hg.edges_.end(), unique_cliques.begin(),
                       unique_cliques.end());
      hg.edge_indices_.push_back(hg.edges_.size());
      hg.edge_weights_.push_back(1);  // Unit weight
      hg.net_to_vertex_.push_back(v);
      hg.n_nets_++;
    }
  }

  // ---- Hypernode (clique) weights (v4) --------------------------------
  // MT-KaHyPar's balance constraint is defined over hypernode weights; the
  // selected scheme decides what a clique "costs" (see ECC.md "Hypernode
  // weighting"). Accumulate in index_t and guard the int32 conversion —
  // MT-KaHyPar's weight type is int (D13).
  const auto guard_int32 = [](index_t w) -> int {
    if (w < 0 || w > static_cast<index_t>(std::numeric_limits<int>::max())) {
      throw HypergraphReorderError(
          "hypernode weight " + std::to_string(w) +
          " exceeds int32 (MT-KaHyPar weight type)");
    }
    return static_cast<int>(w);
  };

  std::vector<index_t> weights64(hg.n_nodes_, 0);

  switch (weighting) {
    case NodeWeighting::kUnit: {
      // Ablation baseline (v3 behavior): balance clique count per part.
      for (index_t c = 0; c < hg.n_nodes_; ++c) {
        weights64[c] = 1;
      }
      break;
    }
    case NodeWeighting::kCliqueSize: {
      // w_i = |C_i|; W = total_pins (vertex incidences, overlaps counted
      // per containing clique).
      const auto &ptr = cover.clique_ptr();
      for (index_t c = 0; c < hg.n_nodes_; ++c) {
        weights64[c] = ptr[c + 1] - ptr[c];
      }
      break;
    }
    case NodeWeighting::kVertexShare: {
      // Each non-isolated vertex holds a budget of SCALE integer units
      // split over its d(v) cliques by largest remainder; the r = SCALE % d
      // remainder units go to the r lowest clique ids (determinism, D12).
      // Per-vertex contribution is exactly SCALE, hence W == SCALE *
      // n_nonisolated exactly (D6). vertex_cliques is ascending and
      // deduplicated by construction.
      const index_t n_nonisolated = hg.n_nets_;  // one net per non-empty vertex
      const index_t scale = vertex_share_scale(n_nonisolated);
      for (index_t v = 0; v < n_vertices; ++v) {
        const auto cliques = cover.get_vertex_cliques(v);
        if (cliques.empty()) continue;  // isolated vertex: no budget
        const index_t d = static_cast<index_t>(cliques.size());
        const index_t base = scale / d;
        const index_t r = scale % d;
        index_t contributed = 0;
        for (index_t i = 0; i < d; ++i) {
          const index_t share = base + (i < r ? 1 : 0);
          weights64[cliques[i]] += share;
          contributed += share;
        }
        assert(contributed == scale && "per-vertex budget must be exact");
      }
      break;
    }
    case NodeWeighting::kEdgeFillMass: {
      // w_i = number of graph edges first-covered by clique i (Phase 2
      // greedy selection), 1 per Phase-3 remainder 2-clique; W == n_edges
      // exactly. The mass array is recorded by every ECC solver; absence
      // is a solver bug, not a user error.
      if (!cover.has_edge_mass() ||
          cover.edge_mass().size() != static_cast<size_t>(hg.n_nodes_)) {
        throw HypergraphReorderError(
            "edge_fill_mass weighting requires an edge-mass array sized "
            "n_cliques (ECC solver did not record edge masses)");
      }
      const auto &mass = cover.edge_mass();
      for (index_t c = 0; c < hg.n_nodes_; ++c) {
        weights64[c] = mass[c];
      }
      break;
    }
    default:
      throw HypergraphReorderError("unknown NodeWeighting");
  }

  hg.node_weights_.resize(hg.n_nodes_);
  hg.total_node_weight_ = 0;
  for (index_t c = 0; c < hg.n_nodes_; ++c) {
    hg.node_weights_[c] = guard_int32(weights64[c]);
    hg.total_node_weight_ += weights64[c];
  }
  // W == SCALE * n_nonisolated exactly for vertex_share (D6).
  assert(weighting != NodeWeighting::kVertexShare ||
         hg.total_node_weight_ ==
             vertex_share_scale(hg.n_nets_) * hg.n_nets_);

  if (!suppress_output)
    std::cout << "CNH constructed: " << hg.n_nets_
              << " nets (efficiency: "
              << (100.0 * static_cast<double>(hg.n_nets_) /
                  static_cast<double>(n_vertices))
              << "%)" << std::endl;

  return hg;
}

}  // namespace hypergraph_reorder
