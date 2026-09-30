// reorderer.hpp - Main API class for symmetric DB reordering
#ifndef HYPERGRAPH_REORDER_REORDERER_HPP
#define HYPERGRAPH_REORDER_REORDERER_HPP

#include <memory>
#include <string>

#include "clique_cover.hpp"
#include "ecc.hpp"
#include "graph.hpp"
#include "hypergraph.hpp"
#include "partitioner.hpp"
#include "sparse_matrix.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Runtime options shared by all pipeline stages.
struct RuntimeOptions {
  int num_threads = 0;          // 0 = auto-detect
  bool suppress_output = true;  // suppresses ALL progress output,
                                // including MT-KaHyPar verbose
};

// Main reorderer class
class SymmetricDBReorderer {
 public:
  // Options are a composition of sub-structs; each sub-struct's defaults
  // are declared exactly once (in its own header) and never re-declared or
  // hand-copied between layers.
  struct Options {
    EccOptions ecc;
    PartitionOptions partition;
    RuntimeOptions runtime;
  };

  // Consumer contract (spec D12): the result is permutation-only — HG-DB
  // never permutes the matrix itself; applying the permutation is the
  // caller's job.
  struct Result {
    std::vector<index_t> permutation;
    VertexPartition partition;
    Statistics stats;
  };

  explicit SymmetricDBReorderer(const Options& opts = Options());

  // Main pipeline: in-memory matrix
  Result reorder(const CSRMatrix& matrix);

 private:
  Options opts_;

  // Helper components, built once in the constructor
  std::unique_ptr<HypergraphPartitioner> partitioner_;
  std::unique_ptr<EccSolver> ecc_solver_;
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_REORDERER_HPP
