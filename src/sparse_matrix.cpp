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

CSRMatrixView CSRMatrix::get_block(
    const std::vector<index_t> &row_indices) const {
  return CSRMatrixView(this, row_indices);
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

void CSRMatrix::make_symmetric() {
  if (n_rows_ != n_cols_) {
    throw HypergraphReorderError("Can only make square matrices symmetric");
  }
  if (is_symmetric_) return;  // Already symmetric

  // Count nonzeros in symmetric matrix
  std::vector<index_t> row_nnz(n_rows_, 0);

  for (index_t i = 0; i < n_rows_; ++i) {
    for (auto j : row_indices(i)) {
      row_nnz[i]++;
      if (i != j) {
        row_nnz[j]++;  // Transpose entry
      }
    }
  }

  // Build new CSR structure
  std::vector<index_t> new_row_ptr(n_rows_ + 1);
  new_row_ptr[0] = 0;
  for (index_t i = 0; i < n_rows_; ++i) {
    new_row_ptr[i + 1] = new_row_ptr[i] + row_nnz[i];
  }

  index_t new_nnz = new_row_ptr[n_rows_];
  std::vector<index_t> new_col_idx(new_nnz);
  std::vector<value_t> new_values;
  if (!pattern_only_) {
    new_values.resize(new_nnz);
  }

  // Fill new structure
  std::vector<index_t> row_pos(n_rows_, 0);
  for (index_t i = 0; i < n_rows_; ++i) {
    auto cols = row_indices(i);
    auto vals = pattern_only_ ? std::span<const value_t>() : row_values(i);

    for (size_t k = 0; k < cols.size(); ++k) {
      index_t j = cols[k];
      value_t v = pattern_only_ ? 1.0 : vals[k];

      // Add (i,j)
      index_t pos = new_row_ptr[i] + row_pos[i]++;
      new_col_idx[pos] = j;
      if (!pattern_only_) new_values[pos] = v;

      // Add (j,i) if not diagonal
      if (i != j) {
        pos = new_row_ptr[j] + row_pos[j]++;
        new_col_idx[pos] = i;
        if (!pattern_only_) new_values[pos] = v;
      }
    }
  }

  // Sort each row by column index
  for (index_t i = 0; i < n_rows_; ++i) {
    index_t start = new_row_ptr[i];
    index_t end = new_row_ptr[i + 1];

    if (!pattern_only_) {
      // Sort with values
      std::vector<std::pair<index_t, value_t>> entries;
      for (index_t k = start; k < end; ++k) {
        entries.emplace_back(new_col_idx[k], new_values[k]);
      }
      std::sort(entries.begin(), entries.end());
      for (size_t k = 0; k < entries.size(); ++k) {
        new_col_idx[start + k] = entries[k].first;
        new_values[start + k] = entries[k].second;
      }
    } else {
      // Sort column indices only
      std::sort(new_col_idx.begin() + start, new_col_idx.begin() + end);
    }
  }

  // Update matrix
  row_ptr_ = std::move(new_row_ptr);
  col_idx_ = std::move(new_col_idx);
  values_ = std::move(new_values);
  nnz_ = new_nnz;
  is_symmetric_ = true;
}

void CSRMatrix::to_pattern_only() {
  if (pattern_only_) return;
  values_.clear();
  values_.shrink_to_fit();
  pattern_only_ = true;
}

CSRMatrix CSRMatrix::extract_upper_triangle() const {
  if (n_rows_ != n_cols_) {
    throw HypergraphReorderError(
        "Can only extract upper triangle from square matrix");
  }

  // Count upper triangle entries
  index_t upper_nnz = 0;
  for (index_t i = 0; i < n_rows_; ++i) {
    for (auto j : row_indices(i)) {
      if (j >= i) upper_nnz++;
    }
  }

  // Build upper triangle
  std::vector<index_t> upper_row_ptr(n_rows_ + 1);
  std::vector<index_t> upper_col_idx;
  std::vector<value_t> upper_values;

  upper_row_ptr[0] = 0;
  upper_col_idx.reserve(upper_nnz);
  if (!pattern_only_) upper_values.reserve(upper_nnz);

  for (index_t i = 0; i < n_rows_; ++i) {
    auto cols = row_indices(i);
    auto vals = pattern_only_ ? std::span<const value_t>() : row_values(i);

    for (size_t k = 0; k < cols.size(); ++k) {
      if (cols[k] >= i) {
        upper_col_idx.push_back(cols[k]);
        if (!pattern_only_) upper_values.push_back(vals[k]);
      }
    }
    upper_row_ptr[i + 1] = upper_col_idx.size();
  }

  return CSRMatrix(n_rows_, n_cols_, std::move(upper_row_ptr),
                   std::move(upper_col_idx), std::move(upper_values), false);
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

// ===== CSRMatrixView Implementation =====

CSRMatrixView::CSRMatrixView(const CSRMatrix *parent,
                             const std::vector<index_t> &row_indices)
    : parent_(parent), row_map_(row_indices) {
  n_rows_ = row_indices.size();

  // Build column mapping
  std::unordered_set<index_t> col_set;
  for (auto row : row_indices) {
    for (auto col : parent->row_indices(row)) {
      col_set.insert(col);
    }
  }

  col_map_.assign(col_set.begin(), col_set.end());
  std::sort(col_map_.begin(), col_map_.end());
  n_cols_ = col_map_.size();
}

index_t CSRMatrixView::nnz() const {
  index_t count = 0;
  for (auto row : row_map_) {
    count += parent_->row_indices(row).size();
  }
  return count;
}

std::span<const index_t> CSRMatrixView::row_indices(index_t local_row) const {
  if (local_row < 0 || local_row >= n_rows_) {
    throw HypergraphReorderError("View row index out of bounds");
  }
  index_t parent_row = row_map_[local_row];
  return parent_->row_indices(parent_row);
}

std::span<const value_t> CSRMatrixView::row_values(index_t local_row) const {
  if (local_row < 0 || local_row >= n_rows_) {
    throw HypergraphReorderError("View row index out of bounds");
  }
  index_t parent_row = row_map_[local_row];
  return parent_->row_values(parent_row);
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
