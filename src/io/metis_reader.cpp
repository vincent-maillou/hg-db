#include <algorithm>
#include <fstream>
#include <sstream>

#include "hypergraph_reorder/io.hpp"

namespace hypergraph_reorder {
namespace metis {

CSRMatrix read(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open METIS file: " + filepath);
  }

  std::string line;

  while (std::getline(file, line)) {
    if (!line.empty() && line[0] != '%') break;
  }

  std::istringstream iss(line);
  index_t n_vertices, n_edges;
  int fmt = 0, ncon = 0;

  iss >> n_vertices >> n_edges;
  if (!iss.eof()) iss >> fmt;
  if (!iss.eof()) iss >> ncon;

  if (iss.fail() || n_vertices <= 0) {
    throw HypergraphReorderError("Invalid METIS header");
  }

  bool has_edge_weights = (fmt % 10) == 1;
  bool has_vertex_weights = (fmt / 10) == 1;

  std::vector<index_t> row_ptr(n_vertices + 1);
  std::vector<index_t> col_idx;
  std::vector<value_t> values;

  row_ptr[0] = 0;
  col_idx.reserve(2 * n_edges);

  for (index_t i = 0; i < n_vertices; ++i) {
    if (!std::getline(file, line)) {
      throw HypergraphReorderError(
          "METIS file has fewer vertices than declared");
    }

    std::istringstream line_stream(line);

    if (has_vertex_weights) {
      for (int c = 0; c < (ncon > 0 ? ncon : 1); ++c) {
        value_t vw;
        line_stream >> vw;
      }
    }

    index_t neighbor;
    while (line_stream >> neighbor) {
      if (neighbor < 1 || neighbor > n_vertices) {
        throw HypergraphReorderError("Invalid METIS neighbor index");
      }

      col_idx.push_back(neighbor - 1);

      if (has_edge_weights) {
        value_t ew;
        line_stream >> ew;
      }
    }

    row_ptr[i + 1] = col_idx.size();
  }

  file.close();

  for (index_t i = 0; i < n_vertices; ++i) {
    std::sort(col_idx.begin() + row_ptr[i], col_idx.begin() + row_ptr[i + 1]);
  }

  return CSRMatrix(n_vertices, n_vertices, std::move(row_ptr),
                   std::move(col_idx), std::vector<value_t>(), true);
}

void write(const std::string& filepath, const CSRMatrix& matrix) {
  if (matrix.n_rows() != matrix.n_cols()) {
    throw HypergraphReorderError("METIS format requires square matrix");
  }

  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open file for writing: " + filepath);
  }

  index_t n_edges = 0;
  for (index_t i = 0; i < matrix.n_rows(); ++i) {
    for (auto j : matrix.row_indices(i)) {
      if (j > i) n_edges++;
    }
  }

  file << matrix.n_rows() << " " << n_edges << " 0\n";

  for (index_t i = 0; i < matrix.n_rows(); ++i) {
    bool first = true;
    for (auto j : matrix.row_indices(i)) {
      if (i != j) {
        if (!first) file << " ";
        file << (j + 1);
        first = false;
      }
    }
    file << "\n";
  }

  file.close();
}

}  // namespace metis
}  // namespace hypergraph_reorder
