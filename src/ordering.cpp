// ordering.cpp - AMD/MMD ordering implementation
#include "hypergraph_reorder/ordering.hpp"

#include <amd.h>
#include <cholmod.h>

#include <algorithm>
#include <iostream>
#include <unordered_map>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace hypergraph_reorder {

BlockOrderer::BlockOrderer(const Options& opts) : opts_(opts) {}

BlockOrderingResult BlockOrderer::apply_ordering(
    const CSRMatrix& matrix, const VertexPartition& partition) {
  if (opts_.method == OrderingMethod::NONE ||
      opts_.method == OrderingMethod::NATURAL) {
    std::cout << "Block ordering disabled, using natural ordering" << std::endl;
    BlockOrderingResult result;
    result.reordered_partition = partition;  // No reordering
    result.method_used = opts_.method;
    return result;
  }

  std::cout << "Applying " << ordering_method_to_string(opts_.method)
            << " ordering to " << partition.n_parts << " blocks" << std::endl;

  BlockOrderingResult result;
  result.reordered_partition = partition;
  result.method_used = opts_.method;

#ifdef USE_OPENMP
  if (opts_.use_parallel) {
    int num_threads =
        opts_.num_threads > 0 ? opts_.num_threads : omp_get_max_threads();
    omp_set_num_threads(num_threads);
    std::cout << "Using " << num_threads << " threads for parallel ordering"
              << std::endl;
  }
#endif

  // Apply ordering to each block in parallel
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic) if (opts_.use_parallel)
#endif
  for (index_t part_id = 0; part_id < partition.n_parts; ++part_id) {
    const auto& block_vertices = partition.parts[part_id];

    if (block_vertices.size() <= 1) {
      // Block too small for reordering
      continue;
    }

    try {
      // Compute ordering for this block
      auto local_perm = compute_block_ordering(matrix, block_vertices);

      // Apply permutation to block vertices
      std::vector<index_t> reordered_block(block_vertices.size());
      for (size_t i = 0; i < local_perm.size(); ++i) {
        reordered_block[i] = block_vertices[local_perm[i]];
      }

      // Update partition (thread-safe since each thread writes to different
      // block)
      result.reordered_partition.parts[part_id] = std::move(reordered_block);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
      result.successful_blocks++;

    } catch (const std::exception& e) {
      std::cerr << "Warning: Ordering failed for block " << part_id << ": "
                << e.what() << std::endl;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
      result.failed_blocks++;
    }
  }

  std::cout << "Block ordering completed: " << result.successful_blocks
            << " successful, " << result.failed_blocks << " failed"
            << std::endl;

  return result;
}

std::vector<index_t> BlockOrderer::compute_block_ordering(
    const CSRMatrix& matrix, const std::vector<index_t>& block_vertices) {
  // Extract block submatrix
  std::vector<index_t> sorted_vertices = block_vertices;
  std::sort(sorted_vertices.begin(), sorted_vertices.end());

  index_t block_size = sorted_vertices.size();

  // Build index mapping
  std::unordered_map<index_t, index_t> old_to_new;
  for (size_t i = 0; i < sorted_vertices.size(); ++i) {
    old_to_new[sorted_vertices[i]] = i;
  }

  // Extract submatrix in COO format
  std::vector<index_t> rows, cols;
  std::vector<value_t> vals;

  for (index_t i = 0; i < block_size; ++i) {
    index_t global_row = sorted_vertices[i];
    auto global_cols = matrix.row_indices(global_row);
    auto global_vals = matrix.pattern_only() ? std::span<const value_t>()
                                             : matrix.row_values(global_row);

    for (size_t k = 0; k < global_cols.size(); ++k) {
      index_t global_col = global_cols[k];

      if (old_to_new.count(global_col)) {
        rows.push_back(i);
        cols.push_back(old_to_new[global_col]);
        if (!matrix.pattern_only()) vals.push_back(global_vals[k]);
      }
    }
  }

  // Create block matrix
  CSRMatrix block_matrix = CSRMatrix::from_coo(
      block_size, block_size, rows, cols, vals, matrix.is_symmetric());

  // Compute ordering based on method
  switch (opts_.method) {
    case OrderingMethod::AMD:
    case OrderingMethod::CAMD:  // Use AMD for CAMD (simplified)
      return compute_amd_ordering(block_matrix);

    case OrderingMethod::METIS:
      return compute_cholmod_ordering(block_matrix, CHOLMOD_METIS);

    case OrderingMethod::NESDIS:
      return compute_cholmod_ordering(block_matrix, CHOLMOD_NESDIS);

    case OrderingMethod::COLAMD:
      return compute_cholmod_ordering(block_matrix, CHOLMOD_COLAMD);

    default:
      throw HypergraphReorderError("Unsupported ordering method");
  }
}

std::vector<index_t> BlockOrderer::compute_amd_ordering(
    const CSRMatrix& block_matrix) {
  index_t n = block_matrix.n_rows();

  // Convert to CSC format (AMD expects column-oriented)
  std::vector<int64_t> Ap(n + 1);
  std::vector<int64_t> Ai;

  // Count column nonzeros
  std::vector<int64_t> col_counts(n, 0);
  for (index_t i = 0; i < n; ++i) {
    for (auto j : block_matrix.row_indices(i)) {
      col_counts[j]++;
    }
  }

  // Build CSC structure
  Ap[0] = 0;
  for (index_t j = 0; j < n; ++j) {
    Ap[j + 1] = Ap[j] + col_counts[j];
  }

  Ai.resize(Ap[n]);
  std::vector<int64_t> col_pos(n, 0);

  for (index_t i = 0; i < n; ++i) {
    for (auto j : block_matrix.row_indices(i)) {
      index_t pos = Ap[j] + col_pos[j]++;
      Ai[pos] = i;
    }
  }

  // Sort each column
  for (index_t j = 0; j < n; ++j) {
    std::sort(Ai.begin() + Ap[j], Ai.begin() + Ap[j + 1]);
  }

  // Call AMD
  std::vector<int64_t> perm(n);
  double control[AMD_CONTROL], info[AMD_INFO];

  amd_l_defaults(control);
  int status = amd_l_order(n, Ap.data(), Ai.data(), perm.data(), control, info);

  if (status != AMD_OK && status != AMD_OK_BUT_JUMBLED) {
    throw HypergraphReorderError("AMD ordering failed");
  }

  // Convert to index_t
  return std::vector<index_t>(perm.begin(), perm.end());
}

std::vector<index_t> BlockOrderer::compute_cholmod_ordering(
    const CSRMatrix& block_matrix, int cholmod_method) {
  // Initialize CHOLMOD
  cholmod_common c;
  cholmod_l_start(&c);

  index_t n = block_matrix.n_rows();

  // Build CHOLMOD sparse matrix (CSC format)
  cholmod_sparse* A = static_cast<cholmod_sparse*>(
      cholmod_l_allocate_sparse(n, n, block_matrix.nnz(),
                                1,  // sorted
                                1,  // packed
                                0,  // stype (symmetric, lower triangle)
                                CHOLMOD_REAL,  // xtype
                                &c));

  // Convert CSR to CSC
  int64_t* Ap = static_cast<int64_t*>(A->p);
  int64_t* Ai = static_cast<int64_t*>(A->i);
  double* Ax = static_cast<double*>(A->x);

  // Count column nonzeros
  std::vector<int64_t> col_counts(n, 0);
  for (index_t i = 0; i < n; ++i) {
    for (auto j : block_matrix.row_indices(i)) {
      col_counts[j]++;
    }
  }

  Ap[0] = 0;
  for (index_t j = 0; j < n; ++j) {
    Ap[j + 1] = Ap[j] + col_counts[j];
  }

  std::vector<int64_t> col_pos(n, 0);
  for (index_t i = 0; i < n; ++i) {
    auto cols = block_matrix.row_indices(i);
    auto vals = block_matrix.pattern_only() ? std::span<const value_t>()
                                            : block_matrix.row_values(i);

    for (size_t k = 0; k < cols.size(); ++k) {
      index_t j = cols[k];
      int64_t pos = Ap[j] + col_pos[j]++;
      Ai[pos] = i;
      if (!block_matrix.pattern_only()) {
        Ax[pos] = vals[k];
      } else {
        Ax[pos] = 1.0;
      }
    }
  }

  // Set ordering method
  c.nmethods = 1;
  c.method[0].ordering = cholmod_method;
  c.supernodal = CHOLMOD_SIMPLICIAL;

  // Analyze
  cholmod_factor* L = cholmod_l_analyze(A, &c);

  if (!L || !L->Perm) {
    cholmod_l_free_factor(&L, &c);
    cholmod_l_free_sparse(&A, &c);
    cholmod_l_finish(&c);
    throw HypergraphReorderError("CHOLMOD ordering failed");
  }

  // Extract permutation
  std::vector<index_t> perm(n);
  int64_t* perm_ptr = static_cast<int64_t*>(L->Perm);
  for (index_t i = 0; i < n; ++i) {
    perm[i] = perm_ptr[i];
  }

  // Cleanup
  cholmod_l_free_factor(&L, &c);
  cholmod_l_free_sparse(&A, &c);
  cholmod_l_finish(&c);

  return perm;
}

}  // namespace hypergraph_reorder
