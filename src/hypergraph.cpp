
#include "hypergraph_reorder/hypergraph.hpp"

#include <algorithm>
#include <iostream>

namespace hypergraph_reorder {

Hypergraph::Hypergraph() : n_nodes_(0), n_nets_(0), total_pins_(0) {}

std::span<const index_t> Hypergraph::net_pins(index_t net_id) const {
  if (net_id < 0 || net_id >= n_nets_) {
    throw HypergraphReorderError("Net ID out of bounds");
  }
  size_t start = edge_indices_[net_id];
  size_t end = edge_indices_[net_id + 1];
  return std::span<const index_t>(edges_.data() + start, end - start);
}

Hypergraph Hypergraph::from_clique_cover(const CliqueCover &cover,
                                         bool suppress_output) {
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
  hg.node_weights_.resize(hg.n_nodes_, 1);  // Unit weights
  hg.net_to_vertex_.reserve(n_vertices);

  hg.edge_indices_.push_back(0);
  hg.n_nets_ = 0;

  // Single-pass construction using pre-built vertex->cliques mapping
  for (index_t v = 0; v < n_vertices; ++v) {
    auto cliques = cover.get_vertex_cliques(v);

    if (!cliques.empty()) {
      // Add net (hyperedge) for this vertex
      // Deduplicate cliques (a vertex might appear in same clique multiple
      // times due to graph structure)
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

  if (!suppress_output)
    std::cout << "CNH constructed: " << hg.n_nets_
              << " nets (efficiency: " << (100.0 * hg.n_nets_ / n_vertices)
              << "%)" << std::endl;

  return hg;
}

}  // namespace hypergraph_reorder
