// partitioner.hpp - MT-KaHyPar integration for hypergraph partitioning
#ifndef HYPERGRAPH_REORDER_PARTITIONER_HPP
#define HYPERGRAPH_REORDER_PARTITIONER_HPP

#include <string>
#include <vector>

#include "clique_cover.hpp"
#include "hypergraph.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// MtKahyparPreset and VertexPartition now live in types.hpp.

// Options for the CNH partitioning stage. Declared exactly once; the
// reorderer passes this sub-struct straight through (no field-by-field
// transcription). Runtime concerns (threads, output) are supplied by the
// reorderer from RuntimeOptions.
struct PartitionOptions {
  index_t n_parts = 4;
  double imbalance = 0.03;
  MtKahyparPreset preset = MtKahyparPreset::DEFAULT;
  int seed = 42;  // fixed default for reproducibility; -1 = random
};

// Hypergraph partitioner using MT-KaHyPar
class HypergraphPartitioner {
 public:
  // The partitioning options plus the runtime fields the partitioner needs
  // (taken from RuntimeOptions by the reorderer).
  HypergraphPartitioner(const PartitionOptions& opts, int num_threads,
                        bool suppress_output);
  ~HypergraphPartitioner();

  // Partition hypergraph using MT-KaHyPar
  HypergraphPartition partition(const Hypergraph& hg);

  // Convert CNH partition to vertex separator partition
  VertexPartition create_vertex_partition(
      const HypergraphPartition& cnh_partition, const CliqueCover& cover,
      index_t n_vertices);

 private:
  PartitionOptions opts_;
  int num_threads_;      // 0 = auto-detect
  bool suppress_output_; // collapses the former suppress_partitioner_output
  void* context_;  // Opaque MT-KaHyPar context

  void init_context();
  void cleanup_context();
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_PARTITIONER_HPP
