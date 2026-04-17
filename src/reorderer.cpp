#include "hypergraph_reorder/reorderer.hpp"

#include <iostream>

namespace hypergraph_reorder {

SymmetricDBReorderer::SymmetricDBReorderer(const Options& opts) : opts_(opts) {
  HypergraphPartitioner::Options part_opts;
  part_opts.n_parts = opts.n_parts;
  part_opts.imbalance = opts.imbalance;
  part_opts.preset = opts.preset;
  part_opts.seed = opts.seed;
  part_opts.suppress_output = opts.suppress_partitioner_output;
  part_opts.num_threads = opts.num_threads;
  partitioner_ = std::make_unique<HypergraphPartitioner>(part_opts);

  BlockOrderer::Options ord_opts;
  ord_opts.method = opts.ordering_method;
  ord_opts.use_parallel = opts.parallel_block_ordering && opts.use_openmp;
  ord_opts.num_threads = opts.num_threads;
  block_orderer_ = std::make_unique<BlockOrderer>(ord_opts);

  CliqueCoverSolver::Options clique_opts;
  clique_opts.use_maximal_cliques = opts.use_maximal_cliques;
  clique_opts.use_parallel = opts.parallel_clique_finding && opts.use_openmp;
  clique_opts.max_clique_enum_vertices = opts.max_clique_enum_vertices;
  clique_opts.max_clique_enum_edges = opts.max_clique_enum_edges;
  clique_opts.num_threads = opts.num_threads;
  clique_solver_ = std::make_unique<CliqueCoverSolver>(clique_opts);
}

SymmetricDBReorderer::Result SymmetricDBReorderer::reorder_from_file(
    const std::string& matrix_path, MatrixFormat format) {
  Timer total_timer;

  if (!opts_.suppress_output) {
    std::cout << "========================================" << std::endl;
    std::cout << "Symmetric DB Form Reordering" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\n[1/8] Loading matrix from " << matrix_path << std::endl;
  }
  Timer timer;
  CSRMatrix matrix = read_matrix(matrix_path, format);
  double time_load = timer.elapsed_ms();
  if (!opts_.suppress_output) {
    std::cout << "Loaded: " << matrix.n_rows() << " x " << matrix.n_cols()
              << " with " << matrix.nnz() << " nonzeros" << std::endl;
    std::cout << "Time: " << time_load << " ms" << std::endl;
  }

  Result result = reorder(matrix);

  result.stats.time_load_ms = time_load;
  result.stats.time_total_ms = total_timer.elapsed_ms();

  if (!opts_.suppress_output) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Total time: " << result.stats.time_total_ms << " ms"
              << std::endl;
    std::cout << "========================================" << std::endl;
  }

  return result;
}

SymmetricDBReorderer::Result SymmetricDBReorderer::reorder(
    const CSRMatrix& matrix) {
  Result result;
  Timer timer;

  if (!opts_.suppress_output)
    std::cout << "\n[2/8] Creating standard graph" << std::endl;
  timer.reset();
  Graph graph = create_graph(matrix);
  result.stats.time_graph_construction_ms = timer.elapsed_ms();
  result.stats.n_vertices = graph.n_vertices();
  result.stats.n_edges = graph.n_edges();
  if (!opts_.suppress_output) {
    std::cout << "Graph: " << graph.n_vertices() << " vertices, "
              << graph.n_edges() << " edges" << std::endl;
    std::cout << "Time: " << result.stats.time_graph_construction_ms << " ms"
              << std::endl;
  }

  if (!opts_.suppress_output)
    std::cout << "\n[3/8] Finding edge-clique cover" << std::endl;
  timer.reset();
  CliqueCover cover = find_clique_cover(graph);
  result.stats.time_clique_cover_ms = timer.elapsed_ms();
  result.stats.n_cliques = cover.n_cliques();
  result.stats.max_clique_size = cover.max_clique_size();
  result.stats.avg_clique_size = cover.avg_clique_size();
  if (!opts_.suppress_output) {
    std::cout << "Clique cover: " << cover.n_cliques() << " cliques"
              << std::endl;
    std::cout << "Time: " << result.stats.time_clique_cover_ms << " ms"
              << std::endl;
  }

  if (!opts_.suppress_output)
    std::cout << "\n[4/8] Creating clique-node hypergraph" << std::endl;
  timer.reset();
  Hypergraph hg = create_hypergraph(cover);
  result.stats.time_hypergraph_construction_ms = timer.elapsed_ms();
  result.stats.n_hypernodes = hg.n_nodes();
  result.stats.n_hyperedges = hg.n_nets();
  result.stats.total_pins = hg.total_pins();
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_hypergraph_construction_ms
              << " ms" << std::endl;

  if (!opts_.suppress_output)
    std::cout << "\n[5/8] Partitioning hypergraph" << std::endl;
  timer.reset();
  HypergraphPartition hg_partition = partition_hypergraph(hg);
  result.stats.time_partitioning_ms = timer.elapsed_ms();
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_partitioning_ms << " ms"
              << std::endl;

  if (!opts_.suppress_output)
    std::cout << "\n[6/8] Creating vertex separator" << std::endl;
  timer.reset();
  VertexPartition vertex_partition =
      create_vertex_partition(hg_partition, cover, graph.n_vertices());
  result.stats.time_separator_construction_ms = timer.elapsed_ms();
  result.stats.n_parts = vertex_partition.n_parts;
  result.stats.separator_size = vertex_partition.separator_size();
  result.stats.separator_ratio =
      vertex_partition.separator_ratio(graph.n_vertices());
  for (const auto& part : vertex_partition.parts) {
    result.stats.part_sizes.push_back(part.size());
  }
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_separator_construction_ms
              << " ms" << std::endl;

  if (!opts_.suppress_output)
    std::cout << "\n[7/8] Applying block ordering" << std::endl;
  timer.reset();
  BlockOrderingResult ordering_result =
      apply_block_ordering(matrix, vertex_partition);
  result.stats.time_block_ordering_ms = timer.elapsed_ms();
  result.stats.blocks_ordered = ordering_result.successful_blocks;
  result.stats.blocks_failed = ordering_result.failed_blocks;
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_block_ordering_ms << " ms"
              << std::endl;

  result.partition = ordering_result.reordered_partition;

  if (!opts_.suppress_output)
    std::cout << "\n[8/8] Permuting matrix" << std::endl;
  timer.reset();
  result.permutation = create_permutation(result.partition, graph.n_vertices());
  result.reordered_matrix = permute_matrix(matrix, result.permutation);
  result.stats.time_permutation_ms = timer.elapsed_ms();
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_permutation_ms << " ms"
              << std::endl;

  result.stats.n_rows = matrix.n_rows();
  result.stats.n_cols = matrix.n_cols();
  result.stats.nnz = matrix.nnz();

  return result;
}

Graph SymmetricDBReorderer::create_graph(const CSRMatrix& matrix) {
  return Graph::from_symmetric_matrix(matrix);
}

CliqueCover SymmetricDBReorderer::find_clique_cover(const Graph& graph) {
  return clique_solver_->solve(graph);
}

Hypergraph SymmetricDBReorderer::create_hypergraph(const CliqueCover& cover) {
  return Hypergraph::from_clique_cover(cover);
}

HypergraphPartition SymmetricDBReorderer::partition_hypergraph(
    const Hypergraph& hg) {
  return partitioner_->partition(hg);
}

VertexPartition SymmetricDBReorderer::create_vertex_partition(
    const HypergraphPartition& cnh_partition, const CliqueCover& cover,
    index_t n_vertices) {
  return partitioner_->create_vertex_partition(cnh_partition, cover,
                                               n_vertices);
}

BlockOrderingResult SymmetricDBReorderer::apply_block_ordering(
    const CSRMatrix& matrix, const VertexPartition& partition) {
  return block_orderer_->apply_ordering(matrix, partition);
}

std::vector<index_t> SymmetricDBReorderer::create_permutation(
    const VertexPartition& partition, index_t n_vertices) {
  std::vector<index_t> perm;
  perm.reserve(n_vertices);

  for (const auto& part : partition.parts) {
    perm.insert(perm.end(), part.begin(), part.end());
  }

  perm.insert(perm.end(), partition.separator.begin(),
              partition.separator.end());

  if (static_cast<index_t>(perm.size()) != n_vertices) {
    throw HypergraphReorderError("Permutation size mismatch");
  }

  return perm;
}

CSRMatrix SymmetricDBReorderer::permute_matrix(
    const CSRMatrix& matrix, const std::vector<index_t>& perm) {
  return ::hypergraph_reorder::permute_matrix(matrix, perm);
}

}  // namespace hypergraph_reorder
