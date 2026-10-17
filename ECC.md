# Edge-Clique Cover (ECC) Methods

This document describes the edge-clique cover methods implemented in HG-DB
v4.1.0, which are the core algorithmic contribution of the library.

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
   **Bron-Kerbosch with pivoting** (Tomita et al. optimization) over the
   later neighborhood (Eppstein–Löffler–Strash outer loop):
   - R = {v} (current clique)
   - P = neighbors of *v* that appear after *v* in the ordering (candidates)
   - X = neighbors of *v* that appear before *v* in the ordering
     (already processed)
3. **Pivot selection**: choose the vertex in P ∪ X with the maximum number
   of neighbors in P.
4. Only explore vertices in P \ N(pivot) — the Tomita optimization.
5. Only report cliques of size ≥ 2.
6. Parallelized via OpenMP over the degeneracy ordering
   (`schedule(dynamic, 1)`).

Best CNH quality; **exponential worst case** — accepted explicitly (no
guard).

## Shared Phase 2: Classical Greedy Set-Cover Selection

Used by `TriEnum` and `BK` (not `EdgeCover`). Implements the **classical
greedy set-cover algorithm** (Chvátal 1979; Johnson 1974): the universe is
the set of all *m* graph edges, the sets are the candidates' edge families
(every pair of clique members is a graph edge by construction of both
Phase-1 producers, so no membership checks are needed), and at every step
the candidate with the **maximum number of currently-uncovered edges** is
accepted. This per-step re-ranking is what earns the
H(s') ≤ ln s' + 1 approximation guarantee with respect to the best cover
over the same candidate family. (It replaced the v4.0.0 static
largest-first accept-if-positive pass, which forfeited the redundancy
elimination — it accepted 86–91% of all enumerated maximal cliques.)

1. **Canonicalize** the Phase-1 candidates: member lists ascending; the
   candidate list sorted by (size descending, canonical content key —
   the ascending member list compared lexicographically — ascending).
   The position in this order ("rank") is the deterministic tie-break
   key, so the selection is a pure function of the candidate *set*,
   independent of the OpenMP schedule of Phase 1.
2. **Tie-break** (equal current marginal): larger clique first (serves
   the blocked-computation interest in large cliques), then the
   canonical content key ascending. Any consistent tie-break preserves
   the approximation guarantee.
3. **Loop**: while uncovered edges remain and an active candidate
   exists, accept the max-marginal candidate; record its edge mass (the
   number of edges it first-covered — the `edge_fill_mass` hypernode
   weight); mark its edges covered; decrement the current marginals of
   every other still-active candidate containing one of those edges
   (zero-marginal candidates are dropped — a set's uncovered-edge count
   never grows). Early exit when coverage completes; a trailing
   zero-marginal pass cannot exist.

Data structures (v4.1.0): dense undirected edge ids numbered 0..m-1
ascending by (u, v); per-candidate edge-id lists (one pass over pairs);
a CSR-style edge→candidates incidence list; a coverage bitmap over
edges; and an ordered active set keyed by (marginal desc, size desc,
canonical rank asc) with erase+reinsert on each decrement. Selection
cost is O(E_Σ log K), where E_Σ = Σ_C (|C| choose 2) candidate-edge
incidences and K = candidates — no `has_edge` binary searches, no
per-edge `std::set` node allocations.

The final coverage is exported from the bitmap in ascending (u, v)
order into the covered-edge set passed forward to Phase 3.

## Shared Phase 3: Cover Remaining Edges

1. Walk all graph edges in ascending (u, v) order — the same order in
   which Phase 2 assigned edge ids — skipping the covered ones via a
   linear merge against the (sorted) Phase-2 covered-edge set. Coverage
   is consumed as-is; it is never re-derived from the selected cliques.
2. Add a **2-clique** {u, v} (mass 1) for every remaining uncovered
   edge.

This guarantees that every edge is covered by at least one clique. For
`BK` the phase is provably a no-op (every graph edge is contained in
some maximal clique); it is the safety net for `TriEnum` (triangle-free
edges); and for `EdgeCover` it runs alone with an empty covered set —
every edge is uncovered by construction.

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
| `use_parallel` | `true` | OpenMP for Phase-1 **enumeration only** (the greedy selection is inherently sequential) |
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
4. Hypernode (clique) weights follow the selected **node weighting scheme**
   (v4) — see [Hypernode weighting](#hypernode-weighting); `unit` is the
   v3 behavior.
5. All hyperedges (nets) have **unit weight** (1).
6. The `net_to_vertex_` mapping enables reverse lookup from hyperedge to
   original vertex.

## Hypernode weighting (v4)

MT-KaHyPar's balance constraint is defined over **hypernode weights**:

```
sum_{i in block} w_i  <=  (1 + eps) * W / k,      W = sum_i w_i
```

The `node_weighting` option selects what a hypernode (clique) weighs —
i.e. what the partitioner *balances*. With the v3 default (`unit`),
`imbalance` bounds **clique count per part**, which is only a weak proxy
for the quantity the doubly-bordered output cares about: **vertices per
diagonal block**. A part's vertex set is derived after partitioning (a
vertex is interior to part *p* iff *all* its cliques landed in *p*, else
it goes to the separator), so block sizes are a non-linear function of
the partition. On clique-heterogeneous covers the divergence is large
(measured on the spec fixture: 2.35x vertex imbalance at eps = 0.03).

### The four schemes

Let `C_i` be clique *i* of the final cover, `d(v)` the number of cliques
containing vertex *v*, `n_nonisolated` the number of vertices with
`d(v) >= 1` (= `n_hyperedges` of the CNH).

| Mode | Weight formula | `W = sum w_i` | Part weight ≈ |
|------|----------------|----------------|----------------|
| `unit` | `w_i = 1` | `n_cliques` | clique count (v3 behavior) |
| `clique_size` (default) | `w_i = \|C_i\|` | `total_pins` | vertex incidences (overlap over-counts) |
| `vertex_share` | `w_i = sum_{v in C_i} s(v,i)` — each vertex splits SCALE integer units over its d(v) cliques (largest remainder, ties to lowest clique ids) | `SCALE * n_nonisolated` (exact) | interior vertices + fractional separator share |
| `edge_fill_mass` | `w_i` = #edges first-covered by `C_i` (Phase 2), 1 per Phase-3 2-clique | `n_edges` (exact) | exclusive edge mass (off-diagonal nnz proxy) |

Per-scheme notes:

- **`unit`** is the ablation baseline and reproduces v3 results.
- **`clique_size`** is the v4 default; fixes the demonstrated vertex-
  imbalance defect. `EdgeCover` covers are provably unaffected by the
  flip (constant order 2 → identical feasible set as `unit` up to a
  factor 2).
- **`vertex_share`**: SCALE = `max(1, min(2^20, 2^30 / n_nonisolated))`,
  so `W <= ~2^30` (int32-safe; MT-KaHyPar's weight type is `int`).
  Per-vertex contribution is *exactly* SCALE — zero rounding drift.
  Degradations (documented): SCALE = 1 (only when `n_nonisolated > 2^30`)
  degenerates to "each vertex gives its single unit to its lowest-id
  clique" — not equivalent to `clique_size`; when `d(v) > SCALE` only
  SCALE of v's cliques receive a unit. The invariants below still hold.
- **`edge_fill_mass`**: attributes each edge to exactly one clique (the
  Phase-2 `uncovered_count` already computed during greedy selection;
  Phase-3 remainder 2-cliques get mass 1). The mass array is recorded by
  all three ECC solvers; a missing/mis-sized array is a solver bug and
  throws, not a user error.

### Interior-vertex bound (hard invariant for weighted modes)

An interior vertex of part *p* has *all* its cliques in *p*, hence it
contributes >= 1 to part *p*'s weight under every weighted scheme:

- `unit`: 1 (per containing clique — the bound applies to clique counts
  instead, which is the defect);
- `clique_size`: `d(v) >= 1` pins;
- `vertex_share`: exactly SCALE >= 1;
- `edge_fill_mass`: every incident edge's owning clique is in `p`, so
  >= deg(v) >= 1.

With part weights <= `(1 + eps_eff) * W / k` it follows that

```
max_i part_size_i  <=  (1 + effective_imbalance) * W / k
```

for every weighted mode — a guaranteed bound, not a heuristic (pinned by
`test_interior_bound_weighted_modes`).

### Epsilon auto-relax

A single hypernode heavier than `(1 + eps) * W / k` makes the instance
infeasible. Rather than failing, the partitioner relaxes eps to the
minimum feasible value `k * w_max / W - 1` with an unconditional stderr
warning and re-issues the partitioning parameters (a single code path —
epsilon is always re-set from `partition()`). The value actually used is
reported as `effective_imbalance` (>= the requested imbalance). Caveat:
the `LARGE_K` preset may apply different block-weight semantics
internally; the relaxation remains a necessary-condition guard there.

### Determinism

Since v4.1.0 the ECC stage is deterministic: the Phase-2 selection is a
pure function of the canonicalized candidate set, so the cover, the
edge masses, and the CNH are identical across thread counts and OpenMP
schedules (`EccOptions::use_parallel` gates only Phase-1 enumeration,
whose arrival order is canonicalized away). `vertex_share` remainders go
to the lowest clique ids — the greedy acceptance order, a *defined*
order rather than an incidental one; `edge_fill_mass` attribution
likewise follows the greedy acceptance order and is fully deterministic
(no process/thread-count caveats). End-to-end result determinism
additionally requires MT-KaHyPar's `deterministic` preset; a fixed seed
alone does not pin the partition search.

### Mode interactions

| Combination | Behavior |
|---|---|
| `EdgeCover` + `clique_size` | ≡ `unit` up to constant factor 2 — identical feasible set, identical partition expected under a deterministic preset |
| `EdgeCover` + `edge_fill_mass` | ≡ `unit` exactly (every mass is 1) |
| `EdgeCover` + `vertex_share` | each vertex's 2-cliques receive SCALE/2 each (+1 to the lower-id one for odd SCALE) |
| Triangle-free graph + any mode | BK/TriEnum covers degenerate to all-2-cliques → `clique_size ≡ unit×2` (why the v3→v4 default flip is invisible on triangle-free fixtures) |
| `vertex_share` at SCALE = 1 | lowest-id-clique attribution — **not** equivalent to `clique_size` |

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
- **Greedy selection** (v4.1.0) is O(E_Σ log K) with small constants
  (E_Σ = total candidate-edge incidences, K = candidates): bitmap
  coverage, dense edge ids, incidence-list marginal decrements, ordered
  active set. It replaced v4.0.0's scan carrying two binary searches
  per pair plus a `std::set` insert/find per edge.
- **Phase 3 remainder** is a linear O(n + m) bitmap/merge walk.
- **CNH construction** is O(total_pins) with a sort+unique per vertex.
