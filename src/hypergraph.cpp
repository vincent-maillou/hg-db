
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

Hypergraph Hypergraph::from_clique_cover(const CliqueCover &cover) {
  Hypergraph hg;
  hg.n_nodes_ = cover.n_cliques();
  index_t n_vertices = cover.n_vertices();

  std::cout << "Constructing CNH: " << hg.n_nodes_ << " nodes, " << n_vertices
            << " vertices" << std::endl;

  // Pre-compute total pins for pre-allocation
  hg.total_pins_ = 0;
  for (index_t v = 0; v < n_vertices; ++v) {
    hg.total_pins_ += cover.get_vertex_cliques(v).size();
  }

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
  index_t total_duplicates = 0;
  for (index_t v = 0; v < n_vertices; ++v) {
    auto cliques = cover.get_vertex_cliques(v);

    if (!cliques.empty()) {
      // Add net (hyperedge) for this vertex
      // Deduplicate cliques (a vertex might appear in same clique multiple
      // times due to graph structure)
      size_t original_size = cliques.size();
      std::vector<index_t> unique_cliques(cliques.begin(), cliques.end());
      std::sort(unique_cliques.begin(), unique_cliques.end());
      auto last = std::unique(unique_cliques.begin(), unique_cliques.end());
      unique_cliques.erase(last, unique_cliques.end());

      if (unique_cliques.size() < original_size) {
        total_duplicates += (original_size - unique_cliques.size());
        if (v < 10) {  // Debug first few
          std::cout << "Vertex " << v << ": removed "
                    << (original_size - unique_cliques.size())
                    << " duplicates (" << original_size << " -> "
                    << unique_cliques.size() << ")" << std::endl;
        }
      }

      hg.edges_.insert(hg.edges_.end(), unique_cliques.begin(),
                       unique_cliques.end());
      hg.edge_indices_.push_back(hg.edges_.size());
      hg.edge_weights_.push_back(1);  // Unit weight
      hg.net_to_vertex_.push_back(v);
      hg.n_nets_++;
    }
  }

  if (total_duplicates > 0) {
    std::cout << "Removed " << total_duplicates
              << " duplicate clique entries total" << std::endl;
  }

  // Verification: check for duplicates in final hypergraph AND print first few
  // nets
  index_t verification_duplicates = 0;
  for (index_t net = 0; net < std::min(hg.n_nets_, index_t(5)); ++net) {
    size_t start = hg.edge_indices_[net];
    size_t end = hg.edge_indices_[net + 1];
    std::cout << "Net " << net << " (vertex " << hg.net_to_vertex_[net] << "): "
              << "indices[" << start << ":" << end << "] = [";
    for (size_t i = start; i < end && i < start + 10; ++i) {
      std::cout << hg.edges_[i];
      if (i + 1 < end && i + 1 < start + 10) std::cout << ", ";
    }
    if (end - start > 10) std::cout << ", ...";
    std::cout << "]" << std::endl;

    std::vector<index_t> pins(hg.edges_.begin() + start,
                              hg.edges_.begin() + end);
    std::vector<index_t> sorted_pins = pins;
    std::sort(sorted_pins.begin(), sorted_pins.end());
    auto last = std::unique(sorted_pins.begin(), sorted_pins.end());
    size_t unique_count = std::distance(sorted_pins.begin(), last);
    if (unique_count < pins.size()) {
      verification_duplicates++;
      std::cout << "  -> HAS DUPLICATES: " << pins.size() << " -> "
                << unique_count << std::endl;
    }
  }
  if (verification_duplicates > 0) {
    std::cout << "WARNING: Found " << verification_duplicates
              << " nets with duplicate pins in final hypergraph!" << std::endl;
  }

  std::cout << "CNH constructed: " << hg.n_nets_
            << " nets (efficiency: " << (100.0 * hg.n_nets_ / n_vertices)
            << "%)" << std::endl;

  return hg;
}

}  // namespace hypergraph_reorder
