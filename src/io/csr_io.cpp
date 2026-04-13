#include <cstring>
#include <fstream>

#include "hypergraph_reorder/io.hpp"

namespace hypergraph_reorder {
namespace csr_binary {

// Binary CSR format:
// Header (fixed size):
//   - magic number (8 bytes): "CSRMAT\0\0"
//   - version (4 bytes): 1
//   - n_rows (8 bytes)
//   - n_cols (8 bytes)
//   - nnz (8 bytes)
//   - is_symmetric (1 byte)
//   - pattern_only (1 byte)
//   - padding (6 bytes)
// Data:
//   - row_ptr array (8 * (n_rows + 1) bytes)
//   - col_idx array (8 * nnz bytes)
//   - values array (8 * nnz bytes, if not pattern_only)

struct BinaryHeader {
  char magic[8];
  int32_t version;
  int64_t n_rows;
  int64_t n_cols;
  int64_t nnz;
  uint8_t is_symmetric;
  uint8_t pattern_only;
  uint8_t padding[6];
};

CSRMatrix read(const std::string& filepath) {
  std::ifstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open binary CSR file: " + filepath);
  }

  BinaryHeader header;
  file.read(reinterpret_cast<char*>(&header), sizeof(header));

  if (file.fail() || std::memcmp(header.magic, "CSRMAT\0\0", 8) != 0) {
    throw HypergraphReorderError("Invalid binary CSR file format");
  }

  if (header.version != 1) {
    throw HypergraphReorderError("Unsupported binary CSR version");
  }

  std::vector<index_t> row_ptr(header.n_rows + 1);
  file.read(reinterpret_cast<char*>(row_ptr.data()),
            (header.n_rows + 1) * sizeof(index_t));

  std::vector<index_t> col_idx(header.nnz);
  file.read(reinterpret_cast<char*>(col_idx.data()),
            header.nnz * sizeof(index_t));

  std::vector<value_t> values;
  if (!header.pattern_only) {
    values.resize(header.nnz);
    file.read(reinterpret_cast<char*>(values.data()),
              header.nnz * sizeof(value_t));
  }

  if (file.fail()) {
    throw HypergraphReorderError("Error reading binary CSR data");
  }

  file.close();

  return CSRMatrix(header.n_rows, header.n_cols, std::move(row_ptr),
                   std::move(col_idx), std::move(values), header.is_symmetric);
}

void write(const std::string& filepath, const CSRMatrix& matrix) {
  std::ofstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open file for writing: " + filepath);
  }

  BinaryHeader header;
  std::memcpy(header.magic, "CSRMAT\0\0", 8);
  header.version = 1;
  header.n_rows = matrix.n_rows();
  header.n_cols = matrix.n_cols();
  header.nnz = matrix.nnz();
  header.is_symmetric = matrix.is_symmetric() ? 1 : 0;
  header.pattern_only = matrix.pattern_only() ? 1 : 0;
  std::memset(header.padding, 0, sizeof(header.padding));

  file.write(reinterpret_cast<const char*>(&header), sizeof(header));

  file.write(reinterpret_cast<const char*>(matrix.row_ptr().data()),
             (matrix.n_rows() + 1) * sizeof(index_t));

  file.write(reinterpret_cast<const char*>(matrix.col_idx().data()),
             matrix.nnz() * sizeof(index_t));

  if (!matrix.pattern_only()) {
    file.write(reinterpret_cast<const char*>(matrix.values().data()),
               matrix.nnz() * sizeof(value_t));
  }

  if (file.fail()) {
    throw HypergraphReorderError("Error writing binary CSR data");
  }

  file.close();
}

}  // namespace csr_binary

MatrixFormat detect_format(const std::string& filepath) {
  if (filepath.size() >= 4 && filepath.substr(filepath.size() - 4) == ".mtx") {
    return MatrixFormat::MTX;
  }
  if (filepath.size() >= 6 &&
      filepath.substr(filepath.size() - 6) == ".graph") {
    return MatrixFormat::METIS;
  }
  if (filepath.size() >= 4 && filepath.substr(filepath.size() - 4) == ".csr") {
    return MatrixFormat::CSR_BINARY;
  }

  std::ifstream file(filepath, std::ios::binary);
  if (!file.is_open()) {
    throw HypergraphReorderError("Cannot open file: " + filepath);
  }

  char magic[8];
  file.read(magic, 8);
  file.close();

  if (std::memcmp(magic, "CSRMAT\0\0", 8) == 0) {
    return MatrixFormat::CSR_BINARY;
  }

  if (magic[0] == '%' && magic[1] == '%') {
    return MatrixFormat::MTX;
  }

  return MatrixFormat::METIS;
}

CSRMatrix read_matrix(const std::string& filepath, MatrixFormat format) {
  if (format == MatrixFormat::AUTO) {
    format = detect_format(filepath);
  }

  switch (format) {
    case MatrixFormat::MTX:
      return mtx::read(filepath);
    case MatrixFormat::METIS:
      return metis::read(filepath);
    case MatrixFormat::CSR_BINARY:
      return csr_binary::read(filepath);
    default:
      throw HypergraphReorderError("Unknown matrix format");
  }
}

void write_matrix(const std::string& filepath, const CSRMatrix& matrix,
                  MatrixFormat format) {
  if (format == MatrixFormat::AUTO) {
    format = detect_format(filepath);
  }

  switch (format) {
    case MatrixFormat::MTX:
      mtx::write(filepath, matrix);
      break;
    case MatrixFormat::METIS:
      metis::write(filepath, matrix);
      break;
    case MatrixFormat::CSR_BINARY:
      csr_binary::write(filepath, matrix);
      break;
    default:
      throw HypergraphReorderError("Unknown matrix format");
  }
}

}  // namespace hypergraph_reorder
