#include "hypergraph_reorder/reorderer.hpp"

#include <iostream>

namespace hypergraph_reorder {

SymmetricDBReorderer::SymmetricDBReorderer(const Options& opts) : opts_(opts) {
  partitioner_ = std::make_unique<HypergraphPartitioner>(
      opts_.partition, opts_.runtime.num_threads, opts_.runtime.suppress_output);

  ecc_solver_ = make_ecc_solver(opts_.ecc);
}

SymmetricDBReorderer::Result SymmetricDBReorderer::reorder(
    const CSRMatrix& matrix) {
  Result result;
  Timer timer;
  Timer total_timer;

  const bool quiet = opts_.runtime.suppress_output;

  // [1/6] Create graph
  if (!quiet) std::cout << "\n[1/6] Creating standard graph" << std::endl;
  timer.reset();
  Graph graph = Graph::from_symmetric_matrix(matrix);
  result.stats.time_graph_ms = timer.elapsed_ms();
  result.stats.n_vertices = graph.n_vertices();
  result.stats.n_edges = graph.n_edges();
  if (!quiet) {
    std::cout << "Graph: " << graph.n_vertices() << " vertices, "
              << graph.n_edges() << " edges" << std::endl;
    std::cout << "Time: " << result.stats.time_graph_ms << " ms" << std::endl;
  }

  // [2/6] Find edge-clique cover
  if (!quiet) std::cout << "\n[2/6] Finding edge-clique cover" << std::endl;
  timer.reset();
  CliqueCover cover = ecc_solver_->solve(graph);
  result.stats.time_ecc_ms = timer.elapsed_ms();
  result.stats.ecc_method_used = opts_.ecc.method;
  result.stats.n_cliques = cover.n_cliques();
  result.stats.max_clique_size = cover.max_clique_size();
  result.stats.avg_clique_size = cover.avg_clique_size();
  result.stats.clique_order_counts = cover.clique_order_counts();
  if (!quiet) {
    std::cout << "Clique cover: " << cover.n_cliques() << " cliques"
              << std::endl;
    std::cout << "Time: " << result.stats.time_ecc_ms << " ms" << std::endl;
  }

  // [3/6] Create clique-node hypergraph
  if (!quiet) std::cout << "\n[3/6] Creating clique-node hypergraph" << std::endl;
  timer.reset();
  Hypergraph hg = Hypergraph::from_clique_cover(cover, quiet);
  result.stats.time_cnh_ms = timer.elapsed_ms();
  result.stats.n_hypernodes = hg.n_nodes();
  result.stats.n_hyperedges = hg.n_nets();
  result.stats.total_pins = hg.total_pins();
  if (!quiet)
    std::cout << "Time: " << result.stats.time_cnh_ms << " ms" << std::endl;

  // [4/6] Partition hypergraph
  if (!quiet) std::cout << "\n[4/6] Partitioning hypergraph" << std::endl;
  timer.reset();
  HypergraphPartition hg_partition = partitioner_->partition(hg);
  result.stats.time_partition_ms = timer.elapsed_ms();
  result.stats.km1_objective = hg_partition.objective;
  if (!quiet)
    std::cout << "Time: " << result.stats.time_partition_ms << " ms"
              << std::endl;

  // [5/6] Create vertex separator (map hypergraph partition back to vertices)
  if (!quiet) std::cout << "\n[5/6] Creating vertex separator" << std::endl;
  timer.reset();
  VertexPartition vertex_partition =
      partitioner_->create_vertex_partition(hg_partition, cover,
                                            graph.n_vertices());
  result.stats.time_separator_ms = timer.elapsed_ms();
  result.stats.n_parts = vertex_partition.n_parts;
  result.stats.separator_size = vertex_partition.separator_size();
  result.stats.separator_ratio =
      vertex_partition.separator_ratio(graph.n_vertices());
  for (const auto& part : vertex_partition.parts) {
    result.stats.part_sizes.push_back(part.size());
  }
  if (!quiet)
    std::cout << "Time: " << result.stats.time_separator_ms << " ms"
              << std::endl;

  // [6/6] Build the DB permutation: diagonal blocks first, separator last.
  // (Permutation-only contract, spec D12 — the matrix itself is never
  // permuted inside HG-DB.)
  if (!quiet) std::cout << "\n[6/6] Building permutation" << std::endl;
  timer.reset();
  result.partition = vertex_partition;
  {
    std::vector<index_t> perm;
    perm.reserve(graph.n_vertices());

    for (const auto& part : vertex_partition.parts) {
      perm.insert(perm.end(), part.begin(), part.end());
    }

    perm.insert(perm.end(), vertex_partition.separator.begin(),
                vertex_partition.separator.end());

    if (static_cast<index_t>(perm.size()) != graph.n_vertices()) {
      throw HypergraphReorderError("Permutation size mismatch");
    }

    result.permutation = std::move(perm);
  }
  result.stats.time_permutation_ms = timer.elapsed_ms();
  if (!quiet)
    std::cout << "Time: " << result.stats.time_permutation_ms << " ms"
              << std::endl;

  result.stats.n_rows = matrix.n_rows();
  result.stats.n_cols = matrix.n_cols();
  result.stats.nnz = matrix.nnz();

  result.stats.time_total_ms = total_timer.elapsed_ms();

  return result;
}

}  // namespace hypergraph_reorder
