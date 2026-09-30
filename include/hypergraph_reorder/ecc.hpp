// ecc.hpp - Edge-clique cover (ECC) strategy interface
#ifndef HYPERGRAPH_REORDER_ECC_HPP
#define HYPERGRAPH_REORDER_ECC_HPP

#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "clique_cover.hpp"
#include "graph.hpp"
#include "types.hpp"

namespace hypergraph_reorder {

// Options for the edge-clique cover stage. Defaults are the single source of
// truth mirrored 1:1 by the C API's hgr_default_options.
struct EccOptions {
  EccMethod method = EccMethod::kBk;  // default ECC method (spec D3)
  bool use_parallel = true;           // OpenMP for enumeration + selection
  int num_threads = 0;                // 0 = auto-detect
  bool suppress_output = true;
};

// Strategy interface: graph in, edge-clique cover out.
class EccSolver {
 public:
  virtual ~EccSolver() = default;
  virtual CliqueCover solve(const Graph& graph) = 0;
};

// Concrete strategies (defined in src/ecc_edge.cpp, src/ecc_tri.cpp,
// src/ecc_bk.cpp).

// Edge-cover ECC solver: every edge becomes a 2-clique; no candidate
// enumeration, no greedy selection (ablation baseline).
class EdgeCoverEccSolver : public EccSolver {
 public:
  explicit EdgeCoverEccSolver(const EccOptions& opts);
  CliqueCover solve(const Graph& graph) override;
  const Statistics& get_stats() const;

 private:
  EccOptions opts_;
  Statistics stats_;
};

// Triangle-enumeration ECC solver: Phase 1 enumerates all triangles
// (parallel two-pointer intersection), then shared Phase 2/3.
class TriangleEnumEccSolver : public EccSolver {
 public:
  explicit TriangleEnumEccSolver(const EccOptions& opts);
  CliqueCover solve(const Graph& graph) override;
  const Statistics& get_stats() const;

 private:
  EccOptions opts_;
  Statistics stats_;
};

// Bron-Kerbosch ECC solver: Phase 1 enumerates all maximal cliques (size
// >= 2) via degeneracy-ordered parallel Bron-Kerbosch with pivoting, then
// shared Phase 2/3. Exponential worst case; no size guard (spec D2).
class BronKerboschEccSolver : public EccSolver {
 public:
  explicit BronKerboschEccSolver(const EccOptions& opts);
  CliqueCover solve(const Graph& graph) override;
  const Statistics& get_stats() const;

 private:
  EccOptions opts_;
  Statistics stats_;
};

// Factory (defined in ecc_common.cpp)
std::unique_ptr<EccSolver> make_ecc_solver(const EccOptions& opts);

// Shared Phase-2/Phase-3 machinery used by TriangleEnumEccSolver and
// BronKerboschEccSolver (definitions in ecc_common.cpp).
namespace detail {

// Phase 2: Greedy clique selection to cover edges (largest first).
// Returns the selected cliques and, via *covered, the set of edges they
// cover (consumed by Phase 3).
std::vector<std::vector<index_t>> select_covering_cliques(
    const Graph& graph,
    const std::vector<std::vector<index_t>>& maximal_cliques,
    std::set<std::pair<index_t, index_t>>& covered);

// Phase 3: Cover remaining (uncovered) edges with 2-cliques.
void cover_remaining_edges(
    const Graph& graph,
    const std::set<std::pair<index_t, index_t>>& covered_edges,
    std::vector<std::vector<index_t>>& cliques);

}  // namespace detail

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_ECC_HPP
