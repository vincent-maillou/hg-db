// ecc_common.cpp - Shared Phase-2/Phase-3 ECC machinery + solver factory
//
// Phase 2 (select_covering_cliques) implements the classical greedy
// set-cover rule (Chvátal 1979; Johnson 1974): at every step it accepts
// the candidate clique with the maximum number of currently-uncovered
// edges, which earns the H(s') <= ln(s') + 1 approximation guarantee
// over the candidate family. Ties break by larger clique first, then by
// canonical content key ascending (spec
// dev-code/hgdb-greedy-phase2-spec.md, D1-D3). Candidates are
// canonicalized up front, so the selection — and hence the entire ECC
// stage — is a pure function of the input graph, independent of the
// OpenMP schedule of Phase 1 (D3).
//
// Phase 3 (cover_remaining_edges) appends a 2-clique (mass 1) for every
// edge the greedy left uncovered: provably none for BK, the safety net
// for TriEnum (triangle-free edges); EdgeCover runs it alone over all
// edges (D5).
#include "hypergraph_reorder/ecc.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

namespace hypergraph_reorder {

namespace detail {

namespace {

// Ordered active-set entry. `rank` is the candidate's position in the
// canonicalized candidate list, whose order is already (size desc,
// canonical key asc) — so ordering by (marginal desc, rank asc) is
// exactly the spec D2 tie-break: maximum current marginal first, then
// the larger clique, then the lexicographically smallest member list.
struct ActiveEntry {
  index_t marginal;
  index_t rank;
};

struct ActiveEntryLess {
  bool operator()(const ActiveEntry &a, const ActiveEntry &b) const {
    if (a.marginal != b.marginal) return a.marginal > b.marginal;
    return a.rank < b.rank;
  }
};

// CSR slot of the arc u -> v (v != u). The graph is canonical (adjacency
// lists sorted ascending, deduplicated), so this is an O(log deg(u))
// lower_bound with no allocation.
index_t arc_slot(const Graph &graph, index_t u, index_t v) {
  const auto neighbors = graph.neighbors(u);
  const auto it = std::lower_bound(neighbors.begin(), neighbors.end(), v);
  return graph.adj_ptr()[u] + (it - neighbors.begin());
}

}  // namespace

std::vector<std::vector<index_t>> select_covering_cliques(
    const Graph &graph,
    const std::vector<std::vector<index_t>> &maximal_cliques,
    std::set<std::pair<index_t, index_t>> &covered_edges,
    std::vector<index_t> &selected_mass) {
  // --- Canonicalize candidates (D3) ---------------------------------------
  // Member lists ascending; candidate order (size desc, canonical key
  // asc). The position in this order ("rank") is the deterministic
  // tie-break key; from here on the result is a pure function of the
  // candidate *set*, independent of Phase-1 enumeration / merge order.
  std::vector<std::vector<index_t>> cliques = maximal_cliques;
  for (auto &clique : cliques) {
    std::sort(clique.begin(), clique.end());
  }
  std::sort(cliques.begin(), cliques.end(),
            [](const std::vector<index_t> &a, const std::vector<index_t> &b) {
              if (a.size() != b.size()) return a.size() > b.size();
              return a < b;  // lexicographic on ascending member lists
            });

  const index_t n = graph.n_vertices();
  const index_t m = graph.n_edges();
  const index_t k_candidates = static_cast<index_t>(cliques.size());

  // --- Dense undirected edge ids (D4) --------------------------------------
  // Edges numbered 0..m-1 ascending by (u, v); edge_id_at[p] is the id of
  // the undirected edge at CSR slot p. The graph is canonical (simple,
  // symmetric, no self-loops), so pass 1 numbers every u < v arc in
  // ascending (u, v) order and pass 2 mirrors the ids to the transpose
  // arcs.
  const std::vector<index_t> &adj_ptr = graph.adj_ptr();
  const std::vector<index_t> &adj_list = graph.adj_list();
  std::vector<index_t> edge_id_at(adj_list.size(), -1);
  index_t next_edge_id = 0;
  for (index_t u = 0; u < n; ++u) {
    for (index_t p = adj_ptr[u]; p < adj_ptr[u + 1]; ++p) {
      if (u < adj_list[p]) {
        edge_id_at[p] = next_edge_id++;
      }
    }
  }
  for (index_t u = 0; u < n; ++u) {
    for (index_t p = adj_ptr[u]; p < adj_ptr[u + 1]; ++p) {
      const index_t v = adj_list[p];
      if (v < u) {
        edge_id_at[p] = edge_id_at[arc_slot(graph, v, u)];
      }
    }
  }
  assert(next_edge_id == m);  // canonical graph: arcs pair up exactly

  // Edge id of an intra-clique pair (u, v): every pair of clique members
  // is adjacent by construction of the Phase-1 producers (maximal cliques
  // / triangles), so the lookup always hits and no has_edge check is
  // needed.
  auto edge_id = [&](index_t u, index_t v) -> index_t {
    return edge_id_at[arc_slot(graph, std::min(u, v), std::max(u, v))];
  };

  // --- Per-clique edge lists + edge->cliques incidence (D4) -----------------
  // One pass over all candidate pairs: each candidate's edge-id list
  // (computed once) and, CSR-style, for every edge the candidates
  // containing it (drives the marginal decrements in the greedy loop).
  std::vector<std::vector<index_t>> clique_edges(k_candidates);
  std::vector<index_t> incidence_count(m, 0);
  for (index_t c = 0; c < k_candidates; ++c) {
    const auto &clique = cliques[c];
    clique_edges[c].reserve(clique.size() * (clique.size() - 1) / 2);
    for (size_t i = 0; i < clique.size(); ++i) {
      for (size_t j = i + 1; j < clique.size(); ++j) {
        const index_t e = edge_id(clique[i], clique[j]);
        clique_edges[c].push_back(e);
        ++incidence_count[e];
      }
    }
  }
  std::vector<index_t> incidence_ptr(m + 1, 0);
  for (index_t e = 0; e < m; ++e) {
    incidence_ptr[e + 1] = incidence_ptr[e] + incidence_count[e];
  }
  std::vector<index_t> incidence(incidence_ptr[m]);
  {
    std::vector<index_t> cursor(incidence_ptr.begin(),
                                incidence_ptr.end() - 1);
    for (index_t c = 0; c < k_candidates; ++c) {
      for (const index_t e : clique_edges[c]) {
        incidence[cursor[e]++] = c;
      }
    }
  }

  // --- Coverage bitmap + ordered active set (D4) ----------------------------
  std::vector<uint64_t> covered_bits((static_cast<size_t>(m) + 63) / 64, 0);
  auto is_covered = [&](index_t e) -> uint64_t {
    return (covered_bits[static_cast<size_t>(e) >> 6] >> (e & 63)) & 1ULL;
  };
  auto mark_covered = [&](index_t e) {
    covered_bits[static_cast<size_t>(e) >> 6] |= 1ULL << (e & 63);
  };

  std::vector<index_t> marginal(k_candidates, 0);
  std::vector<char> is_active(k_candidates, 0);
  std::set<ActiveEntry, ActiveEntryLess> active;
  for (index_t c = 0; c < k_candidates; ++c) {
    marginal[c] = static_cast<index_t>(clique_edges[c].size());
    if (marginal[c] > 0) {
      is_active[c] = 1;
      active.insert({marginal[c], c});
    }
  }

  // --- Classical greedy selection (D1) ---------------------------------------
  // While edges remain uncovered and an active candidate exists: accept
  // the maximal-marginal candidate, record its edge mass (edges it
  // first-covered), mark its edges covered, and decrement the marginals
  // of every other still-active candidate on those edges. Zero-marginal
  // candidates are dropped (a set's uncovered-edge count never grows).
  std::vector<std::vector<index_t>> selected_cliques;
  index_t uncovered = m;
  while (uncovered > 0 && !active.empty()) {
    const ActiveEntry best = *active.begin();
    active.erase(active.begin());
    const index_t c = best.rank;
    is_active[c] = 0;

    selected_cliques.push_back(std::move(cliques[c]));

    index_t mass = 0;
    for (const index_t e : clique_edges[c]) {
      if (is_covered(e)) continue;
      mark_covered(e);
      --uncovered;
      ++mass;
      for (index_t p = incidence_ptr[e]; p < incidence_ptr[e + 1]; ++p) {
        const index_t d = incidence[p];
        if (!is_active[d]) continue;
        active.erase({marginal[d], d});
        --marginal[d];
        if (marginal[d] == 0) {
          is_active[d] = 0;
        } else {
          active.insert({marginal[d], d});
        }
      }
    }
    selected_mass.push_back(mass);
  }

  // --- Export the coverage for Phase 3 ---------------------------------------
  // Ascending (u, v) walk — the same order the edge ids were assigned in
  // — so the exported set is sorted, as Phase 3's merge assumes.
  for (index_t u = 0; u < n; ++u) {
    for (index_t p = adj_ptr[u]; p < adj_ptr[u + 1]; ++p) {
      if (u < adj_list[p] && is_covered(edge_id_at[p])) {
        covered_edges.insert({u, adj_list[p]});
      }
    }
  }

#ifndef NDEBUG
  // Debug cross-check (spec §9 risk mitigation): the incremental
  // bookkeeping must agree with a naive recount — coverage rebuilt from
  // the selected cliques' pairs, and the mass sum equal to the number of
  // covered edges.
  {
    std::vector<uint64_t> check_bits(covered_bits.size(), 0);
    auto mark_check = [&](index_t e) {
      check_bits[static_cast<size_t>(e) >> 6] |= 1ULL << (e & 63);
    };
    index_t mass_sum = 0;
    for (size_t s = 0; s < selected_cliques.size(); ++s) {
      const auto &clique = selected_cliques[s];
      for (size_t i = 0; i < clique.size(); ++i) {
        for (size_t j = i + 1; j < clique.size(); ++j) {
          mark_check(edge_id(clique[i], clique[j]));
        }
      }
      mass_sum += selected_mass[s];
    }
    index_t covered_count = 0;
    for (index_t e = 0; e < m; ++e) {
      covered_count += is_covered(e) ? 1 : 0;
    }
    assert(check_bits == covered_bits);
    assert(mass_sum == covered_count);
  }
#endif

  return selected_cliques;
}

void cover_remaining_edges(
    const Graph &graph,
    const std::set<std::pair<index_t, index_t>> &covered_edges,
    std::vector<std::vector<index_t>> &cliques,
    std::vector<index_t> &masses) {
  // Walk all graph edges in ascending (u, v) order — the order in which
  // Phase 2 assigned edge ids — and skip the covered ones with a
  // two-pointer merge over the (sorted) covered set. Coverage is consumed
  // as-is; it is never re-derived from the selected cliques.
  auto covered_it = covered_edges.begin();
  for (index_t u = 0; u < graph.n_vertices(); ++u) {
    for (const index_t v : graph.neighbors(u)) {
      if (v <= u) continue;
      if (covered_it != covered_edges.end() && covered_it->first == u &&
          covered_it->second == v) {
        ++covered_it;
        continue;
      }
      // Each remainder 2-clique covers exactly one edge, so its edge
      // mass is 1 (keeps masses index-aligned with cliques).
      cliques.push_back({u, v});
      masses.push_back(1);
    }
  }
}

}  // namespace detail

std::unique_ptr<EccSolver> make_ecc_solver(const EccOptions &opts) {
  switch (opts.method) {
    case EccMethod::kEdgeCover:
      return std::make_unique<EdgeCoverEccSolver>(opts);
    case EccMethod::kTriEnum:
      return std::make_unique<TriangleEnumEccSolver>(opts);
    case EccMethod::kBk:
    default:
      return std::make_unique<BronKerboschEccSolver>(opts);
  }
}

}  // namespace hypergraph_reorder
