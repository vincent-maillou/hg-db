// reorderer.hpp - Main API class for symmetric DB reordering
#ifndef HYPERGRAPH_REORDER_REORDERER_HPP
#define HYPERGRAPH_REORDER_REORDERER_HPP

#include <memory>
#include <string>

#include "clique_cover.hpp"
#include "graph.hpp"
#include "hypergraph.hpp"
#include "io.hpp"
#include "ordering.hpp"
#include "partitioner.hpp"
#include "sparse_matrix.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Main reorderer class
class SymmetricDBReorderer {
 public:
  struct Options {
    // Partitioning options
    index_t n_parts;
    double imbalance;
    MtKahyparPreset preset;  // MT-KaHyPar preset (replaces config file)
    int seed;

    // Clique cover options
    bool use_maximal_cliques;
    bool parallel_clique_finding;
    index_t max_clique_enum_vertices;

    // Ordering options
    OrderingMethod ordering_method;
    bool parallel_block_ordering;

    // Performance options
    bool use_openmp;
    int num_threads;
    bool suppress_partitioner_output;
    bool suppress_output;  // Suppress all reorderer progress output

    Options()
        : n_parts(4),
          imbalance(0.03),
          preset(MtKahyparPreset::DEFAULT),
          seed(-1),
          use_maximal_cliques(true),
          parallel_clique_finding(true),
          max_clique_enum_vertices(
              100000),  // Use triangles for graphs > 100k vertices
          ordering_method(OrderingMethod::AMD),
          parallel_block_ordering(true),
          use_openmp(true),
          num_threads(0),  // 0 = auto-detect
          suppress_partitioner_output(false),
          suppress_output(false) {}
  };

  struct Result {
    CSRMatrix reordered_matrix;
    std::vector<index_t> permutation;
    VertexPartition partition;
    Statistics stats;
  };

  explicit SymmetricDBReorderer(const Options& opts = Options());

  // Main pipeline: file input
  Result reorder_from_file(const std::string& matrix_path,
                           MatrixFormat format = MatrixFormat::AUTO);

  // Main pipeline: in-memory matrix
  Result reorder(const CSRMatrix& matrix);

  // Advanced: step-by-step API for fine control
  Graph create_graph(const CSRMatrix& matrix);
  CliqueCover find_clique_cover(const Graph& graph);
  Hypergraph create_hypergraph(const CliqueCover& cover);
  HypergraphPartition partition_hypergraph(const Hypergraph& hg);
  VertexPartition create_vertex_partition(
      const HypergraphPartition& cnh_partition, const CliqueCover& cover,
      index_t n_vertices);
  BlockOrderingResult apply_block_ordering(const CSRMatrix& matrix,
                                           const VertexPartition& partition);
  std::vector<index_t> create_permutation(const VertexPartition& partition,
                                          index_t n_vertices);
  CSRMatrix permute_matrix(const CSRMatrix& matrix,
                           const std::vector<index_t>& perm);

 private:
  Options opts_;

  // Helper components
  std::unique_ptr<HypergraphPartitioner> partitioner_;
  std::unique_ptr<BlockOrderer> block_orderer_;
  std::unique_ptr<CliqueCoverSolver> clique_solver_;
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_REORDERER_HPP
