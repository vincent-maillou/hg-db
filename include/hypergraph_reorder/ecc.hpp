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
  bool use_parallel = true;           // OpenMP for Phase-1 enumeration only
                                     // (the greedy selection is sequential)
  int num_threads = 0;                // 0 = auto-detect
  bool suppress_output = true;
};

// Strategy interface: graph in, edge-clique cover out.
//
// Timing/statistics ownership (v5): each solver populates ONLY the
// solver-owned Statistics fields — degeneracy, n_maximal_cliques,
// maximal_clique_enumeration_histogram, time_degeneracy_ms,
// time_enumeration_ms, time_greedy_ms — exposed via get_stats(). The
// reorderer merges them right after solve() and keeps owning the
// cover-level fields (n_cliques, max_clique_size, avg_clique_size,
// clique_order_counts). Cover-stage wall time is the sum of the three
// solver-owned timings (no monolithic ECC timing key anymore).
class EccSolver {
 public:
  virtual ~EccSolver() = default;
  virtual CliqueCover solve(const Graph& graph) = 0;
  virtual const Statistics& get_stats() const = 0;
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

// Phase 2: classical greedy set-cover selection (Chvátal 1979; Johnson
// 1974). At every step accepts the candidate clique with the maximum
// number of currently-uncovered edges; ties break by larger clique
// first, then by canonical content key (ascending member list,
// lexicographic) ascending. Candidates are canonicalized up front
// (members ascending, list sorted by size desc / key asc), so the result
// is a pure function of the candidate set — independent of the Phase-1
// OpenMP schedule (deterministic). Complexity: O(E_Sigma log K) where
// E_Sigma = total candidate-edge incidences, K = candidates; coverage is
// tracked with a bitmap plus an ordered active set of current marginals.
// Returns the selected cliques and, via *covered, the set of edges they
// cover (consumed by Phase 3). Via *selected_mass, records each accepted
// clique's marginal at acceptance (the number of edges it first-covered
// — the edge_fill_mass hypernode weight); index-aligned with the
// returned selected cliques.
std::vector<std::vector<index_t>> select_covering_cliques(
    const Graph& graph,
    const std::vector<std::vector<index_t>>& maximal_cliques,
    std::set<std::pair<index_t, index_t>>& covered,
    std::vector<index_t>& selected_mass);

// Phase 3: Cover remaining (uncovered) edges with 2-cliques. Walks all
// graph edges in ascending (u, v) order (= ascending edge-id order) and
// skips those present in *covered_edges (the Phase-2 coverage, consumed
// as-is — never re-derived). Appends 1 per added 2-clique to *masses,
// keeping it index-aligned with cliques. For BK input this is provably a
// no-op; it is the safety net for TriEnum (triangle-free edges) and the
// whole pipeline for EdgeCover (empty covered set).
void cover_remaining_edges(
    const Graph& graph,
    const std::set<std::pair<index_t, index_t>>& covered_edges,
    std::vector<std::vector<index_t>>& cliques,
    std::vector<index_t>& masses);

}  // namespace detail

}  // namespace hypergraph_reorder

#endif  // HYPERGRAPH_REORDER_ECC_HPP
