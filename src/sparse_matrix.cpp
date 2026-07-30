// sparse_matrix.cpp - Implementation of CSR sparse matrix
#include "hypergraph_reorder/sparse_matrix.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace hypergraph_reorder {

CSRMatrix::CSRMatrix()
    : n_rows_(0),
      n_cols_(0),
      nnz_(0),
      is_symmetric_(false),
      pattern_only_(false) {}

CSRMatrix::CSRMatrix(index_t n_rows, index_t n_cols,
                     std::vector<index_t> row_ptr, std::vector<index_t> col_idx,
                     std::vector<value_t> values, bool is_symmetric)
    : n_rows_(n_rows),
      n_cols_(n_cols),
      is_symmetric_(is_symmetric),
      row_ptr_(std::move(row_ptr)),
      col_idx_(std::move(col_idx)),
      values_(std::move(values)) {
  if (row_ptr_.size() != static_cast<size_t>(n_rows + 1)) {
    throw HypergraphReorderError("Invalid row_ptr size");
  }

  nnz_ = row_ptr_[n_rows];

  if (col_idx_.size() != static_cast<size_t>(nnz_)) {
    throw HypergraphReorderError("col_idx size doesn't match nnz");
  }

  pattern_only_ = values_.empty();
  if (!pattern_only_ && values_.size() != static_cast<size_t>(nnz_)) {
    throw HypergraphReorderError("values size doesn't match nnz");
  }
}

std::span<const index_t> CSRMatrix::row_indices(index_t row) const {
  if (row < 0 || row >= n_rows_) {
    throw HypergraphReorderError("Row index out of bounds");
  }
  index_t start = row_ptr_[row];
  index_t end = row_ptr_[row + 1];
  return std::span<const index_t>(col_idx_.data() + start, end - start);
}

std::span<const value_t> CSRMatrix::row_values(index_t row) const {
  if (pattern_only_) {
    return std::span<const value_t>();
  }
  if (row < 0 || row >= n_rows_) {
    throw HypergraphReorderError("Row index out of bounds");
  }
  index_t start = row_ptr_[row];
  index_t end = row_ptr_[row + 1];
  return std::span<const value_t>(values_.data() + start, end - start);
}

CSRMatrix CSRMatrix::clone() const {
  return CSRMatrix(n_rows_, n_cols_, std::vector<index_t>(row_ptr_),
                   std::vector<index_t>(col_idx_),
                   std::vector<value_t>(values_), is_symmetric_);
}

bool CSRMatrix::validate() const {
  if (row_ptr_.size() != static_cast<size_t>(n_rows_ + 1)) return false;
  if (col_idx_.size() != static_cast<size_t>(nnz_)) return false;
  if (!pattern_only_ && values_.size() != static_cast<size_t>(nnz_))
    return false;

  // Check row_ptr is non-decreasing
  for (index_t i = 0; i < n_rows_; ++i) {
    if (row_ptr_[i] > row_ptr_[i + 1]) return false;
  }

  // Check column indices are in range
  for (index_t col : col_idx_) {
    if (col < 0 || col >= n_cols_) return false;
  }

  return true;
}

CSRMatrix CSRMatrix::from_coo(index_t n_rows, index_t n_cols,
                              const std::vector<index_t> &rows,
                              const std::vector<index_t> &cols,
                              const std::vector<value_t> &vals,
                              bool is_symmetric) {
  if (rows.size() != cols.size()) {
    throw HypergraphReorderError("COO rows and cols size mismatch");
  }

  bool has_values = !vals.empty();
  if (has_values && vals.size() != rows.size()) {
    throw HypergraphReorderError("COO values size mismatch");
  }

  index_t nnz = rows.size();

  // Count nonzeros per row
  std::vector<index_t> row_nnz(n_rows, 0);
  for (auto r : rows) {
    if (r < 0 || r >= n_rows) {
      throw HypergraphReorderError("COO row index out of bounds");
    }
    row_nnz[r]++;
  }

  // Build row_ptr
  std::vector<index_t> row_ptr(n_rows + 1);
  row_ptr[0] = 0;
  for (index_t i = 0; i < n_rows; ++i) {
    row_ptr[i + 1] = row_ptr[i] + row_nnz[i];
  }

  // Fill col_idx and values
  std::vector<index_t> col_idx(nnz);
  std::vector<value_t> values;
  if (has_values) values.resize(nnz);

  std::vector<index_t> row_pos(n_rows, 0);
  for (index_t k = 0; k < nnz; ++k) {
    index_t i = rows[k];
    index_t j = cols[k];

    if (j < 0 || j >= n_cols) {
      throw HypergraphReorderError("COO col index out of bounds");
    }

    index_t pos = row_ptr[i] + row_pos[i]++;
    col_idx[pos] = j;
    if (has_values) values[pos] = vals[k];
  }

  // Sort each row by column index
  for (index_t i = 0; i < n_rows; ++i) {
    index_t start = row_ptr[i];
    index_t end = row_ptr[i + 1];

    if (has_values) {
      std::vector<std::pair<index_t, value_t>> entries;
      for (index_t k = start; k < end; ++k) {
        entries.emplace_back(col_idx[k], values[k]);
      }
      std::sort(entries.begin(), entries.end());
      for (size_t k = 0; k < entries.size(); ++k) {
        col_idx[start + k] = entries[k].first;
        values[start + k] = entries[k].second;
      }
    } else {
      std::sort(col_idx.begin() + start, col_idx.begin() + end);
    }
  }

  return CSRMatrix(n_rows, n_cols, std::move(row_ptr), std::move(col_idx),
                   std::move(values), is_symmetric);
}

CSRMatrix permute_matrix(const CSRMatrix &matrix,
                         const std::vector<index_t> &perm) {
  index_t n = matrix.n_rows();
  if (n != matrix.n_cols()) {
    throw HypergraphReorderError("Can only permute square matrices");
  }
  if (static_cast<index_t>(perm.size()) != n) {
    throw HypergraphReorderError("Permutation size doesn't match matrix size");
  }

  // Create inverse permutation
  std::vector<index_t> inv_perm(n);
  for (index_t i = 0; i < n; ++i) {
    if (perm[i] < 0 || perm[i] >= n) {
      throw HypergraphReorderError("Invalid permutation");
    }
    inv_perm[perm[i]] = i;
  }

  // Build permuted matrix in COO format
  std::vector<index_t> rows, cols;
  std::vector<value_t> vals;

  rows.reserve(matrix.nnz());
  cols.reserve(matrix.nnz());
  if (!matrix.pattern_only()) vals.reserve(matrix.nnz());

  for (index_t old_i = 0; old_i < n; ++old_i) {
    index_t new_i = inv_perm[old_i];
    auto old_cols = matrix.row_indices(old_i);
    auto old_vals = matrix.pattern_only() ? std::span<const value_t>()
                                          : matrix.row_values(old_i);

    for (size_t k = 0; k < old_cols.size(); ++k) {
      index_t old_j = old_cols[k];
      index_t new_j = inv_perm[old_j];

      rows.push_back(new_i);
      cols.push_back(new_j);
      if (!matrix.pattern_only()) vals.push_back(old_vals[k]);
    }
  }

  // Convert to CSR
  return CSRMatrix::from_coo(n, n, rows, cols, vals, matrix.is_symmetric());
}

}  // namespace hypergraph_reorder
