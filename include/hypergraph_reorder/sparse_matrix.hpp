#ifndef HYPERGRAPH_REORDER_SPARSE_MATRIX_HPP
#define HYPERGRAPH_REORDER_SPARSE_MATRIX_HPP

#include <memory>
#include <span>
#include <vector>

#include "types.hpp"

namespace hypergraph_reorder {

// CSR (Compressed Sparse Row) matrix
class CSRMatrix {
 public:
  // Default constructor
  CSRMatrix();

  // Constructor from CSR arrays
  CSRMatrix(index_t n_rows, index_t n_cols, std::vector<index_t> row_ptr,
            std::vector<index_t> col_idx, std::vector<value_t> values = {},
            bool is_symmetric = false);

  // Move constructor
  CSRMatrix(CSRMatrix&& other) noexcept = default;
  CSRMatrix& operator=(CSRMatrix&& other) noexcept = default;

  // Delete copy constructor (use clone() if copy is needed)
  CSRMatrix(const CSRMatrix&) = delete;
  CSRMatrix& operator=(const CSRMatrix&) = delete;

  // Dimensions
  index_t n_rows() const { return n_rows_; }
  index_t n_cols() const { return n_cols_; }
  index_t nnz() const { return nnz_; }
  bool is_symmetric() const { return is_symmetric_; }
  bool pattern_only() const { return pattern_only_; }

  // Access CSR arrays
  const std::vector<index_t>& row_ptr() const { return row_ptr_; }
  const std::vector<index_t>& col_idx() const { return col_idx_; }
  const std::vector<value_t>& values() const { return values_; }

  // Get row data (zero-copy spans)
  std::span<const index_t> row_indices(index_t row) const;
  std::span<const value_t> row_values(index_t row) const;

  // Clone the matrix
  CSRMatrix clone() const;

  // Validate matrix structure
  bool validate() const;

  // Create from coordinate format (COO)
  static CSRMatrix from_coo(index_t n_rows, index_t n_cols,
                            const std::vector<index_t>& rows,
                            const std::vector<index_t>& cols,
                            const std::vector<value_t>& vals = {},
                            bool is_symmetric = false);

 private:
  index_t n_rows_;
  index_t n_cols_;
  index_t nnz_;
  bool is_symmetric_;
  bool pattern_only_;

  std::vector<index_t> row_ptr_;  // Size: n_rows + 1
  std::vector<index_t> col_idx_;  // Size: nnz
  std::vector<value_t> values_;   // Size: nnz (or 0 if pattern_only)
};

// Apply permutation to matrix (create permuted copy)
CSRMatrix permute_matrix(const CSRMatrix& matrix,
                         const std::vector<index_t>& perm);

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_SPARSE_MATRIX_HPP
