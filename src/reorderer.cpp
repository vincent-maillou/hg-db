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

  CliqueCoverSolver::Options clique_opts;
  clique_opts.use_maximal_cliques = opts.use_maximal_cliques;
  clique_opts.use_parallel = opts.parallel_clique_finding && opts.use_openmp;
  clique_opts.max_clique_enum_vertices = opts.max_clique_enum_vertices;
  clique_opts.max_clique_enum_edges = opts.max_clique_enum_edges;
  clique_opts.num_threads = opts.num_threads;
  clique_opts.suppress_output = opts.suppress_output;
  clique_solver_ = std::make_unique<CliqueCoverSolver>(clique_opts);
}

SymmetricDBReorderer::Result SymmetricDBReorderer::reorder(
    const CSRMatrix& matrix) {
  Result result;
  Timer timer;

  // [1/6] Create graph
  if (!opts_.suppress_output)
    std::cout << "\n[1/6] Creating standard graph" << std::endl;
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

  // [2/6] Find edge-clique cover
  if (!opts_.suppress_output)
    std::cout << "\n[2/6] Finding edge-clique cover" << std::endl;
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

  // [3/6] Create clique-node hypergraph
  if (!opts_.suppress_output)
    std::cout << "\n[3/6] Creating clique-node hypergraph" << std::endl;
  timer.reset();
  Hypergraph hg = create_hypergraph(cover, opts_.suppress_output);
  result.stats.time_hypergraph_construction_ms = timer.elapsed_ms();
  result.stats.n_hypernodes = hg.n_nodes();
  result.stats.n_hyperedges = hg.n_nets();
  result.stats.total_pins = hg.total_pins();
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_hypergraph_construction_ms
              << " ms" << std::endl;

  // [4/6] Partition hypergraph
  if (!opts_.suppress_output)
    std::cout << "\n[4/6] Partitioning hypergraph" << std::endl;
  timer.reset();
  HypergraphPartition hg_partition = partition_hypergraph(hg);
  result.stats.time_partitioning_ms = timer.elapsed_ms();
  if (!opts_.suppress_output)
    std::cout << "Time: " << result.stats.time_partitioning_ms << " ms"
              << std::endl;

  // [5/6] Create vertex separator (map hypergraph partition back to vertices)
  if (!opts_.suppress_output)
    std::cout << "\n[5/6] Creating vertex separator" << std::endl;
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

  // [6/6] Build permutation and permute matrix
  // (Block ordering is handled externally by Parallax's ComposedOrdering)
  if (!opts_.suppress_output)
    std::cout << "\n[6/6] Building permutation" << std::endl;
  timer.reset();
  result.partition = vertex_partition;
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

Hypergraph SymmetricDBReorderer::create_hypergraph(const CliqueCover& cover,
                                                    bool suppress_output) {
  return Hypergraph::from_clique_cover(cover, suppress_output);
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
