// ordering.hpp - AMD/MMD integration for block ordering
#ifndef HYPERGRAPH_REORDER_ORDERING_HPP
#define HYPERGRAPH_REORDER_ORDERING_HPP

#include <vector>

#include "partitioner.hpp"
#include "sparse_matrix.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Block ordering result
struct BlockOrderingResult {
  VertexPartition reordered_partition;  // Partition with reordered blocks
  index_t successful_blocks = 0;
  index_t failed_blocks = 0;
  OrderingMethod method_used = OrderingMethod::NONE;
};

// Block ordering algorithms
class BlockOrderer {
 public:
  struct Options {
    OrderingMethod method;
    bool use_parallel;
    int num_threads;

    Options()
        : method(OrderingMethod::AMD), use_parallel(true), num_threads(-1) {}
  };

  explicit BlockOrderer(const Options& opts = Options());

  // Apply ordering to all blocks in partition
  BlockOrderingResult apply_ordering(const CSRMatrix& matrix,
                                     const VertexPartition& partition);

 private:
  Options opts_;

  // Compute ordering for a single block submatrix
  std::vector<index_t> compute_block_ordering(
      const CSRMatrix& matrix, const std::vector<index_t>& block_vertices);

  // AMD ordering
  std::vector<index_t> compute_amd_ordering(const CSRMatrix& block_matrix);

  // CHOLMOD-based ordering (METIS, NESDIS, etc.)
  std::vector<index_t> compute_cholmod_ordering(const CSRMatrix& block_matrix,
                                                int cholmod_method);
};

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_ORDERING_HPP
