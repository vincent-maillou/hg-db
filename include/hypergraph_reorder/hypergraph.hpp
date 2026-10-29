// hypergraph.hpp - Clique-node hypergraph (CNH) construction
#ifndef HYPERGRAPH_REORDER_HYPERGRAPH_HPP
#define HYPERGRAPH_REORDER_HYPERGRAPH_HPP

#include <span>
#include <vector>

#include "clique_cover.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Clique-Node Hypergraph (CNH) representation
// Notes are basically the cliques from the clique cover and nets are the
// original vertices Net j connects node i iff clique i contains vertex j
class Hypergraph {
 public:
  Hypergraph();

  index_t n_nodes() const { return n_nodes_; }
  index_t n_nets() const { return n_nets_; }
  index_t total_pins() const { return total_pins_; }
  // Total hypernode weight W under the scheme used at construction
  // (stored, not re-summed).
  index_t total_node_weight() const { return total_node_weight_; }

  // Access hypergraph structure (kahyapr compatible format)
  const std::vector<size_t> &edge_indices() const { return edge_indices_; }
  const std::vector<index_t> &edges() const { return edges_; }
  const std::vector<int> &node_weights() const { return node_weights_; }

  // Mapping: net_id -> original vertex_id
  const std::vector<index_t> &net_to_vertex() const { return net_to_vertex_; }

  // Get nodes connected by a net (zero-copy span)
  std::span<const index_t> net_pins(index_t net_id) const;

  // Create CNH from clique cover (optimized with pre-built mapping).
  // The weighting scheme selects what a hypernode (clique) weighs for the
  // partitioner's balance constraint (v4); default kCliqueSize.
  static Hypergraph from_clique_cover(const CliqueCover &cover,
                                      bool suppress_output = false,
                                      NodeWeighting weighting =
                                          NodeWeighting::kCliqueSize);

 private:
  index_t n_nodes_;     // Number of hypernodes (cliques)
  index_t n_nets_;      // Number of hyperedges (original vertices)
  index_t total_pins_;  // Total number of pins
  index_t total_node_weight_;  // W = sum of node weights (v4)

  // Hypergraph in KaHyPar format (zero-copy compatible)
  // Nets are unit-weight by construction: no net-weight array is stored
  // and null is passed at the partition boundary (WP2 D5) — MT-KaHyPar's
  // documented unit-weight path. A future net-weighting feature
  // re-introduces the field only when it has a consumer.
  std::vector<size_t> edge_indices_;  // Size: n_nets + 1
  std::vector<index_t> edges_;        // Size: total_pins
  std::vector<int> node_weights_;     // Size: n_nodes

  // Mapping back to original graph
  std::vector<index_t> net_to_vertex_;  // net_id -> vertex_id
};

// Partition result from hypergraph partitioning
struct HypergraphPartition {
  index_t n_parts;
  std::vector<index_t> node_partition;  // node_id -> part_id

  // Metrics of the achieved CNH partition (unit net weights):
  //   km1_objective: sum over nets of (connectivity - 1)   [connectivity metric]
  //   cut_objective: number of nets spanning >= 2 parts     [cut-net metric]
  // Reported for every run regardless of the optimized objective.
  index_t km1_objective = 0;
  index_t cut_objective = 0;

  // Imbalance actually handed to MT-KaHyPar: the requested value, or the
  // minimum feasible value when a single hypernode exceeded a block's
  // weight budget and epsilon was auto-relaxed (v4, see partitioner).
  double effective_imbalance = 0.0;

  std::vector<index_t> part_sizes;
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_HYPERGRAPH_HPP
