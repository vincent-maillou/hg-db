# Edge-Clique Cover (ECC) Methods

This document describes the edge-clique cover methods implemented in HG-DB
v2, which are the core algorithmic contribution of the library.

## Overview

The ECC computes a set of cliques (complete subgraphs) that cover every edge
of the input graph. These cliques are then used to construct a
**Clique-Node Hypergraph (CNH)**, where:

- **Hypernodes** = cliques from the cover
- **Hyperedges (nets)** = original vertices
- Net *j* connects node *i* iff clique *i* contains vertex *j*

The CNH is partitioned by MT-KaHyPar, and the partition is mapped back to
original vertices to produce a doubly-bordered (DB) permutation.

## The three ECC methods

ECC method selection is an **explicit, first-class choice** (`ecc_method`);
there is no hidden size-based heuristic and **no guard** against
combinatorial blow-up — the caller is responsible for choosing an
appropriate method for their graph (D2).

### `EdgeCover` (`HGR_ECC_EDGE_COVER`)

Every edge becomes a 2-clique. No candidate enumeration, no greedy
selection — Phase 3 alone over all edges. Equivalent to partitioning the
standard edge-hypergraph (incidence) formulation through the same CNH
machinery. **This is the mandatory ablation baseline**: it isolates the
benefit of higher-order clique structure for the partitioner. Zero
enumeration cost; the resulting histogram is exactly `{2: n_edges}`.

### `TriEnum` (`HGR_ECC_TRI_ENUM`)

1. For each vertex *u*, iterate over edges *(u, v)* where *v > u*.
2. Use **two-pointer intersection** of N(u) and N(v) restricted to *w > v*
   to find common neighbors.
3. Each common neighbor *w* produces triangle {u, v, w}.
4. Parallelized via OpenMP over vertices (`schedule(dynamic, 64)`).
5. Thread-local storage, merged at end.

A scalable baseline; candidates are limited to triangles (order-3 cliques).

### `BK` (`HGR_ECC_BK`) — the default

1. Compute a **degeneracy ordering** using the Matula-Beck algorithm
   (bucket-sorted degrees, O(n + m)).
2. For each vertex *v* in the ordering, spawn an independent
   **Bron-Kerbosch with pivoting** (Tomita et al. optimization):
   - R = {v} (current clique)
   - P = neighbors of *v* that appear after *v* in the ordering (candidates)
   - X = ∅ (already processed)
3. **Pivot selection**: choose the vertex in P ∪ X with the maximum number
   of neighbors in P.
4. Only explore vertices in P \ N(pivot) — the Tomita optimization.
5. Only report cliques of size ≥ 2.
6. Parallelized via OpenMP over the degeneracy ordering
   (`schedule(dynamic, 1)`).

Best CNH quality; **exponential worst case** — accepted explicitly (no
guard).

## Shared Phase 2: Greedy Covering Selection

Used by `TriEnum` and `BK` (not `EdgeCover`):

1. Sort all candidate cliques (triangles or maximal cliques) by **size
   descending**.
2. For each clique, count how many **uncovered edges** it contains.
3. If `uncovered_count > 0`, add the clique to the cover and mark its edges
   as covered.
4. Edge coverage is tracked using `std::set<std::pair<index_t, index_t>>`.
   This is correct but O(log E) per insertion — a potential bottleneck for
   very large graphs.

The populated covered-edge set is passed forward to Phase 3.

## Shared Phase 3: Cover Remaining Edges

1. Iterate all graph edges.
2. Remove edges already covered by selected cliques (the Phase-2 covered
   set; coverage is no longer rebuilt from scratch).
3. Add a **2-clique** {u, v} for every remaining uncovered edge.

This guarantees that every edge is covered by at least one clique. For
`EdgeCover` this phase runs alone (with an empty covered set — every edge
is uncovered by construction).

## Old → new mapping

| Old configuration | Effective old behavior | New equivalent |
|---|---|---|
| `use_maximal_cliques=true`, small graph (≤ 5000 vertices, ≤ 100000 edges) | BK | `ecc_method="BK"` |
| `use_maximal_cliques=true`, large graph | TriEnum (hidden heuristic) | `ecc_method="TriEnum"` |
| `use_maximal_cliques=false` *or* `greedy_only=true` | all-2-clique cover (undocumented) | `ecc_method="EdgeCover"` |

Note that the old flags `use_maximal_cliques=false` and `greedy_only=true`
never selected "triangles only" as documented — both bypassed candidate
enumeration entirely and produced the trivial all-2-clique cover.

## `EccOptions`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `method` | `kBk` | ECC method: `kEdgeCover` / `kTriEnum` / `kBk` |
| `use_parallel` | `true` | OpenMP parallelism for enumeration + selection |
| `num_threads` | `0` | Thread count (0 = auto-detect) |
| `suppress_output` | `true` | Suppress ECC progress output |

The deleted flags `use_maximal_cliques`, `greedy_only`,
`max_clique_enum_vertices`, and `max_clique_enum_edges` had no equivalent —
their behavior is replaced by the explicit `method` choice.

## Clique-Node Hypergraph Construction

After the clique cover is computed:

1. For each vertex *v*, collect all cliques containing *v* → this forms a
   hyperedge (net).
2. Deduplicate cliques per vertex via `sort + unique` (a cheap safety no-op
   now that the input graph is canonicalized).
3. Skip vertices with zero cliques (isolated vertices).
4. All hypernodes (cliques) have **unit weight** (1).
5. All hyperedges (nets) have **unit weight** (1).
6. The `net_to_vertex_` mapping enables reverse lookup from hyperedge to
   original vertex.

## Partitioning and the km1 objective

The CNH is partitioned by MT-KaHyPar minimizing the **connectivity (km1)**
metric: the sum over nets of (connectivity λ − 1). Every net (original
vertex) whose cliques span λ ≥ 2 parts contributes λ − 1. This is **not**
the separator size: a net split across parts puts its vertex in the
separator, and vertices in no net (isolated ones) land in the separator
without affecting km1.

Configuration levers:

- `n_parts` — more parts → more cut nets → larger km1 and typically a
  larger separator.
- `imbalance` — tighter balance → larger km1 and separator.
- `preset` — `quality`/`deterministic` trade km1 quality vs time and
  reproducibility.
- `seed` — run-to-run variance of km1 at fixed configuration (default 42
  for reproducibility; −1 for random).

## Clique-order histogram

`Statistics::clique_order_counts` (C++ `std::map<index_t, index_t>`,
exposed via the C `hgr_get_clique_order_counts` accessor) maps clique
**order** (size: 2 for an edge, 3 for a triangle, …) to the number of
cliques of that order in the **final selected cover** (after greedy
selection and the 2-clique remainder — not the raw candidate
enumeration). Populated for all three methods; `EdgeCover` yields the
degenerate histogram `{2: n_edges}`.

## Data Structures

### CliqueCover

Dual CSR representation:

- **Clique-centric**: `clique_ptr_[c]` / `clique_members_` — "what vertices
  are in this clique?"
- **Vertex-centric**: `vertex_cliques_ptr_[v]` / `vertex_cliques_` — "which
  cliques contain this vertex?"

Both provide zero-copy `std::span` access.

## Performance Considerations

- **BK** is exponential in the worst case; bounded in practice by the
  degeneracy ordering. No size guard — caller's responsibility (D2).
- **TriEnum** is O(m^(3/2)) in the worst case but practical for sparse
  graphs.
- **EdgeCover** is O(m) — pure bookkeeping.
- **Greedy covering** with `std::set` edge tracking is O(k log E) where k
  is the number of candidate cliques. For very large graphs with many
  candidates, this can dominate runtime.
- **CNH construction** is O(total_pins) with a sort+unique per vertex.
