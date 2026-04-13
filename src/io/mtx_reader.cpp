#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "hypergraph_reorder/io.hpp"

namespace hypergraph_reorder {
namespace mtx {

CSRMatrix read(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open file: " + filepath);
  }

  std::string line;
  bool is_matrix, is_coordinate, is_pattern, is_symmetric;
  index_t n_rows, n_cols, nnz;

  if (!std::getline(file, line)) {
    throw HypergraphReorderError("Failed to read MTX banner");
  }

  std::transform(line.begin(), line.end(), line.begin(), ::tolower);
  is_matrix = line.find("matrix") != std::string::npos;
  is_coordinate = line.find("coordinate") != std::string::npos;
  is_pattern = line.find("pattern") != std::string::npos;
  is_symmetric = line.find("symmetric") != std::string::npos;

  if (!is_matrix || !is_coordinate) {
    throw HypergraphReorderError("Only coordinate matrix format supported");
  }

  while (std::getline(file, line)) {
    if (line.empty() || line[0] != '%') break;
  }

  auto trim = [](const std::string& str) {
    auto start = std::find_if_not(str.begin(), str.end(), [](unsigned char c) {
      return std::isspace(c);
    });
    auto end = std::find_if_not(str.rbegin(), str.rend(), [](unsigned char c) {
                 return std::isspace(c);
               }).base();
    return (start < end) ? std::string(start, end) : std::string();
  };

  std::istringstream iss(trim(line));
  iss >> n_rows >> n_cols >> nnz;

  if (iss.fail()) {
    throw HypergraphReorderError("Failed to parse MTX dimensions");
  }

  std::vector<index_t> rows, cols;
  std::vector<value_t> vals;

  rows.reserve(nnz);
  cols.reserve(nnz);
  if (!is_pattern) vals.reserve(nnz);

  index_t entries_read = 0;

  while (std::getline(file, line) && entries_read < nnz) {
    line = trim(line);
    if (line.empty()) continue;

    std::istringstream iss(line);
    index_t i, j;
    value_t v = 1.0;

    iss >> i >> j;
    if (iss.fail()) continue;

    if (!is_pattern) {
      iss >> v;
      if (iss.fail()) v = 1.0;
    }

    rows.push_back(i - 1);
    cols.push_back(j - 1);
    if (!is_pattern) vals.push_back(v);

    entries_read++;
  }

  file.close();

  if (entries_read != nnz) {
    throw HypergraphReorderError("MTX file has fewer entries than declared");
  }

  if (is_symmetric) {
    size_t original_size = rows.size();
    for (size_t k = 0; k < original_size; ++k) {
      if (rows[k] != cols[k]) {
        rows.push_back(cols[k]);
        cols.push_back(rows[k]);
        if (!is_pattern) vals.push_back(vals[k]);
      }
    }
  }

  return CSRMatrix::from_coo(n_rows, n_cols, rows, cols, vals, is_symmetric);
}

void write(const std::string& filepath, const CSRMatrix& matrix) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open file for writing: " + filepath);
  }

  file << "%%MatrixMarket matrix coordinate ";
  if (matrix.pattern_only()) {
    file << "pattern ";
  } else {
    file << "real ";
  }
  if (matrix.is_symmetric()) {
    file << "symmetric\n";
  } else {
    file << "general\n";
  }

  index_t entries_to_write = 0;
  if (matrix.is_symmetric()) {
    for (index_t i = 0; i < matrix.n_rows(); ++i) {
      for (auto j : matrix.row_indices(i)) {
        if (j >= i) entries_to_write++;
      }
    }
  } else {
    entries_to_write = matrix.nnz();
  }

  file << matrix.n_rows() << " " << matrix.n_cols() << " " << entries_to_write
       << "\n";

  file << std::scientific << std::setprecision(16);

  for (index_t i = 0; i < matrix.n_rows(); ++i) {
    auto cols = matrix.row_indices(i);
    auto vals = matrix.pattern_only() ? std::span<const value_t>()
                                      : matrix.row_values(i);

    for (size_t k = 0; k < cols.size(); ++k) {
      index_t j = cols[k];

      if (matrix.is_symmetric() && j < i) continue;

      file << (i + 1) << " " << (j + 1);

      if (!matrix.pattern_only()) {
        file << " " << vals[k];
      }

      file << "\n";
    }
  }

  file.close();
}

}  // namespace mtx
}  // namespace hypergraph_reorder
