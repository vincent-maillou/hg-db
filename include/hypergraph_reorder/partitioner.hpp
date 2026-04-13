// partitioner.hpp - MT-KaHyPar integration for hypergraph partitioning
#ifndef HYPERGRAPH_REORDER_PARTITIONER_HPP
#define HYPERGRAPH_REORDER_PARTITIONER_HPP

#include <string>
#include <vector>

#include "clique_cover.hpp"
#include "hypergraph.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// MT-KaHyPar preset configuration
enum class MtKahyparPreset {
  DEFAULT,        // Fast, good quality (default preset)
  QUALITY,        // Higher quality, slower
  DETERMINISTIC,  // Deterministic partitioning
  LARGE_K         // Optimized for large number of parts
};

// Vertex partition (result of converting CNH partition to vertex separator)
struct VertexPartition {
  index_t n_parts;
  std::vector<std::vector<index_t>> parts;  // parts[i] = vertices in part i
  std::vector<index_t> separator;           // Separator vertices

  index_t separator_size() const { return separator.size(); }
  double separator_ratio(index_t n_total) const {
    return static_cast<double>(separator.size()) / n_total;
  }
};

// Hypergraph partitioner using MT-KaHyPar
class HypergraphPartitioner {
 public:
  struct Options {
    index_t n_parts;
    double imbalance;
    MtKahyparPreset preset;  // MT-KaHyPar preset (replaces config file)
    int seed;
    bool suppress_output;
    int num_threads;  // Number of threads for MT-KaHyPar

    Options()
        : n_parts(4),
          imbalance(0.03),
          preset(MtKahyparPreset::DEFAULT),
          seed(-1),
          suppress_output(false),
          num_threads(0) {}  // 0 = auto-detect
  };

  explicit HypergraphPartitioner(const Options& opts = Options());
  ~HypergraphPartitioner();

  // Partition hypergraph using MT-KaHyPar
  HypergraphPartition partition(const Hypergraph& hg);

  // Convert CNH partition to vertex separator partition
  VertexPartition create_vertex_partition(
      const HypergraphPartition& cnh_partition, const CliqueCover& cover,
      index_t n_vertices);

 private:
  Options opts_;
  void* context_;  // Opaque MT-KaHyPar context

  void init_context();
  void cleanup_context();
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_PARTITIONER_HPP
