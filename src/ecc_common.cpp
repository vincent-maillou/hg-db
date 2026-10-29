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

// Indexed 4-ary max-heap over the active candidates, ordered by
// ActiveEntryLess (marginal descending, rank ascending), with a
// candidate -> slot map (`heap_slot`, -1 = inactive).
//
// Bit-identity argument (WP3 spec D1/D4): `rank` is unique, so the
// comparator is a strict total order, so the maximum element is unique
// and the extraction sequence is independent of the heap internals —
// the same guarantee the std::set active set gave, which is what makes
// this a drop-in replacement. Every operation is O(log K) with ZERO
// allocation per greedy event: the std::set paid a red-black-tree node
// allocation/free pair for every erase+reinsert, i.e. for E_Sigma
// incidence decrements, which dominated the stage at benchmark scale.
//
// 4-ary (not binary): shallower tree, better cache behavior on large K;
// the choice is invisible to output under the total-order argument.
class ActiveHeap {
 public:
  explicit ActiveHeap(index_t k_candidates)
      : heap_slot_(static_cast<size_t>(k_candidates), -1) {}

  bool empty() const { return heap_.empty(); }

  bool contains(index_t rank) const { return slot(rank) != -1; }

  // Insert a candidate with the given marginal.
  void push(index_t marginal, index_t rank) {
    heap_.push_back({marginal, rank});
    slot(rank) = static_cast<index_t>(heap_.size()) - 1;
    sift_up(static_cast<index_t>(heap_.size()) - 1);
  }

  // Remove and return the maximum element (unique under the strict
  // total order — exactly what `*active.begin()` returned).
  ActiveEntry extract_max() {
    const ActiveEntry best = at(0);
    slot(best.rank) = -1;
    const index_t last = static_cast<index_t>(heap_.size()) - 1;
    if (last > 0) {
      at(0) = at(last);
      slot(at(0).rank) = 0;
    }
    heap_.pop_back();
    if (last > 0) sift_down(0);
    return best;
  }

  // Decrement a member's marginal by one and restore the heap order.
  // The key only ever decreases (a set's uncovered-edge count never
  // grows), so the parent relation upward remains valid and the element
  // can only move DOWN — an in-place sift-down, no allocation. A
  // candidate whose marginal reaches zero leaves the active set.
  void decrement(index_t rank) {
    const index_t i = slot(rank);
    if (--at(i).marginal == 0) {
      erase(rank);
    } else {
      sift_down(i);
    }
  }

  // Remove an arbitrary member (position known via the slot map): swap
  // in the last element, pop, then re-order the moved element in its
  // new position (it may need to go up or down).
  void erase(index_t rank) {
    const index_t i = slot(rank);
    slot(rank) = -1;
    const index_t last = static_cast<index_t>(heap_.size()) - 1;
    if (i < last) {
      at(i) = at(last);
      slot(at(i).rank) = i;
      heap_.pop_back();
      const index_t parent = (i - 1) / 4;
      if (i > 0 && precedes(at(i), at(parent))) {
        sift_up(i);
      } else {
        sift_down(i);
      }
    } else {
      heap_.pop_back();
    }
  }

#ifndef NDEBUG
  // Full O(K) invariant check (debug builds only; fixture-scale runs):
  // heap order under ActiveEntryLess, slot-map consistency, and the
  // active count equal to the heap size.
  void assert_invariants() const {
    for (index_t i = 0; i < size(); ++i) {
      assert(slot(at(i).rank) == i);
      for (index_t child = 4 * i + 1;
           child < 4 * i + 5 && child < size();
           ++child) {
        assert(precedes(at(i), at(child)));
      }
    }
    index_t n_active = 0;
    for (const index_t s : heap_slot_) {
      if (s != -1) ++n_active;
    }
    assert(n_active == size());
  }
#endif

 private:
  static bool precedes(const ActiveEntry &a, const ActiveEntry &b) {
    return ActiveEntryLess{}(a, b);
  }

  const ActiveEntry &at(index_t i) const {
    return heap_[static_cast<size_t>(i)];
  }
  ActiveEntry &at(index_t i) { return heap_[static_cast<size_t>(i)]; }
  index_t slot(index_t rank) const {
    return heap_slot_[static_cast<size_t>(rank)];
  }
  index_t &slot(index_t rank) {
    return heap_slot_[static_cast<size_t>(rank)];
  }
  index_t size() const { return static_cast<index_t>(heap_.size()); }

  void sift_up(index_t i) {
    while (i > 0) {
      const index_t parent = (i - 1) / 4;
      if (!precedes(at(i), at(parent))) break;
      std::swap(at(i), at(parent));
      slot(at(i).rank) = i;
      slot(at(parent).rank) = parent;
      i = parent;
    }
  }

  void sift_down(index_t i) {
    const index_t n = size();
    while (true) {
      const index_t first = 4 * i + 1;
      if (first >= n) break;
      index_t best = i;
      for (index_t child = first; child < first + 4 && child < n; ++child) {
        if (precedes(at(child), at(best))) best = child;
      }
      if (best == i) break;
      std::swap(at(i), at(best));
      slot(at(i).rank) = i;
      slot(at(best).rank) = best;
      i = best;
    }
  }

  std::vector<ActiveEntry> heap_;                   // the heap
  std::vector<index_t> heap_slot_;                  // candidate -> position
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
    std::vector<uint64_t> &covered_bits,
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

  // --- Per-clique edge lists + edge->cliques incidence (D4, WP3 D3) ---------
  // Each candidate's edge-id list lives in ONE flat CSR-style buffer
  // (clique_edges_ptr/clique_edges_flat) instead of a vector-of-vectors:
  // pass 1 sizes each candidate combinatorially — s·(s−1)/2 edges for a
  // clique of size s (0 for s < 2), no pair enumeration needed — and
  // pass 2 enumerates the pairs in the same (i < j, ascending members)
  // order as before, so the per-candidate edge-id sequences — and hence
  // the incidence CSR fill and the greedy marginals — are identical to
  // the vector-of-vectors version. One allocation for the member array,
  // one for the prefix; zero per-candidate allocations.
  std::vector<index_t> clique_edges_ptr(k_candidates + 1, 0);
  for (index_t c = 0; c < k_candidates; ++c) {
    const size_t s = cliques[c].size();
    clique_edges_ptr[c + 1] =
        clique_edges_ptr[c] + static_cast<index_t>(s >= 2 ? s * (s - 1) / 2 : 0);
  }
  std::vector<index_t> clique_edges_flat(clique_edges_ptr[k_candidates]);
  std::vector<index_t> incidence_count(m, 0);
  for (index_t c = 0; c < k_candidates; ++c) {
    const auto &clique = cliques[c];
    index_t cursor = clique_edges_ptr[c];
    for (size_t i = 0; i < clique.size(); ++i) {
      for (size_t j = i + 1; j < clique.size(); ++j) {
        const index_t e = edge_id(clique[i], clique[j]);
        clique_edges_flat[cursor++] = e;
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
      for (index_t p = clique_edges_ptr[c]; p < clique_edges_ptr[c + 1]; ++p) {
        incidence[cursor[clique_edges_flat[p]]++] = c;
      }
    }
  }

  // --- Coverage bitmap + indexed active heap (D4, WP3 D1/D2) -----------------
  // The bitmap is the exported coverage (WP3 D2): Phase 3 consumes it
  // directly, reconstructing edge ids by walk order — no intermediate
  // set of covered edges is ever built (the v5.2 export walk inserted
  // every covered edge into a node-based std::set, ≈m allocations on a
  // BK cover, plus a node-destroying teardown at solve() scope exit
  // that landed outside every timed stage).
  covered_bits.assign((static_cast<size_t>(m) + 63) / 64, 0);
  auto is_covered = [&](index_t e) -> uint64_t {
    return (covered_bits[static_cast<size_t>(e) >> 6] >> (e & 63)) & 1ULL;
  };
  auto mark_covered = [&](index_t e) {
    covered_bits[static_cast<size_t>(e) >> 6] |= 1ULL << (e & 63);
  };

  // Active set = indexed 4-ary max-heap under ActiveEntryLess. The slot
  // map is the single activity flag; a candidate's current marginal
  // lives in its heap entry. Initial marginals come from the Phase-1
  // CSR spans. Zero-marginal candidates never enter (a candidate with
  // no edges can never cover anything).
  ActiveHeap active_heap(k_candidates);
  for (index_t c = 0; c < k_candidates; ++c) {
    const index_t marginal = clique_edges_ptr[c + 1] - clique_edges_ptr[c];
    if (marginal > 0) {
      active_heap.push(marginal, c);
    }
  }
#ifndef NDEBUG
  active_heap.assert_invariants();
#endif

  // --- Classical greedy selection (D1) ---------------------------------------
  // While edges remain uncovered and an active candidate exists: accept
  // the maximal-marginal candidate, record its edge mass (edges it
  // first-covered), mark its edges covered, and decrement the marginals
  // of every other still-active candidate on those edges. Zero-marginal
  // candidates are dropped (a set's uncovered-edge count never grows).
  std::vector<std::vector<index_t>> selected_cliques;
  index_t uncovered = m;
  while (uncovered > 0 && !active_heap.empty()) {
    const ActiveEntry best = active_heap.extract_max();
    const index_t c = best.rank;

    selected_cliques.push_back(std::move(cliques[c]));

    index_t mass = 0;
    for (index_t p = clique_edges_ptr[c]; p < clique_edges_ptr[c + 1]; ++p) {
      const index_t e = clique_edges_flat[p];
      if (is_covered(e)) continue;
      mark_covered(e);
      --uncovered;
      ++mass;
      for (index_t q = incidence_ptr[e]; q < incidence_ptr[e + 1]; ++q) {
        const index_t d = incidence[q];
        if (active_heap.contains(d)) {
          active_heap.decrement(d);
        }
      }
    }
    selected_mass.push_back(mass);
#ifndef NDEBUG
    active_heap.assert_invariants();
#endif
  }

  // The coverage export is the bitmap itself (WP3 D2) — nothing more to
  // do here.

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
    const std::vector<uint64_t> &covered_bits,
    std::vector<std::vector<index_t>> &cliques,
    std::vector<index_t> &masses) {
  // Walk all graph edges in ascending (u, v) order — exactly the order
  // in which Phase 2 assigned edge ids (its pass 1 numbers every u < v
  // arc over ascending CSR slots, which is precisely this walk) — so a
  // running counter reconstructs each edge's id with no lookup table
  // (WP3 D2; the id domain here is Phase 2's id domain by construction).
  // Per arc: one bitmap test; uncovered edges become remainder 2-cliques.
  // Coverage is consumed as-is; it is never re-derived from the selected
  // cliques.
  index_t e = 0;
  for (index_t u = 0; u < graph.n_vertices(); ++u) {
    for (const index_t v : graph.neighbors(u)) {
      if (v <= u) continue;
      const bool covered =
          (covered_bits[static_cast<size_t>(e) >> 6] >> (e & 63)) & 1ULL;
      ++e;  // every u < v arc consumes an id, covered or not
      if (covered) continue;
      // Each remainder 2-clique covers exactly one edge, so its edge
      // mass is 1 (keeps masses index-aligned with cliques).
      cliques.push_back({u, v});
      masses.push_back(1);
    }
  }
  assert(e == graph.n_edges());  // the walk enumerated every edge exactly once
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
