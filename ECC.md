# Edge-Clique Cover (ECC) Algorithm

This document describes the ECC algorithm implemented in HG-DB, which is the
core algorithmic contribution of the library.

## Overview

The ECC algorithm finds a set of cliques (complete subgraphs) that cover every
edge of the input graph. These cliques are then used to construct a
**Clique-Node Hypergraph (CNH)**, where:

- **Hypernodes** = cliques from the cover
- **Hyperedges (nets)** = original vertices
- Net *j* connects node *i* iff clique *i* contains vertex *j*

The CNH is partitioned by MT-KaHyPar, and the partition is mapped back to
original vertices to produce a doubly-bordered (DB) permutation.

## Algorithm Phases

### Phase 1: Candidate Clique Enumeration

Two strategies are available, selected based on graph size:

#### Strategy A: Parallel Bron-Kerbosch (small graphs)

Used when `n_vertices ≤ max_clique_enum_vertices` (default 5000) **AND**
`n_edges ≤ max_clique_enum_edges` (default 100000).

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

#### Strategy B: Triangle Enumeration (large graphs)

Used when the graph exceeds the Bron-Kerbosch thresholds.

1. For each vertex *u*, iterate over edges *(u, v)* where *v > u*.
2. Use **two-pointer intersection** of N(u) and N(v) restricted to *w > v*
   to find common neighbors.
3. Each common neighbor *w* produces triangle {u, v, w}.
4. Parallelized via OpenMP over vertices (`schedule(dynamic, 64)`).
5. Thread-local storage, merged at end.

### Phase 2: Greedy Covering Selection

1. Sort all candidate cliques (maximal cliques or triangles) by **size
   descending**.
2. For each clique, count how many **uncovered edges** it contains.
3. If `uncovered_count > 0`, add the clique to the cover and mark its edges
   as covered.
4. Edge coverage is tracked using `std::set<std::pair<index_t, index_t>>`.
   This is correct but O(log E) per insertion — a potential bottleneck for
   very large graphs.

### Phase 3: Cover Remaining Edges

1. Iterate all graph edges.
2. Remove edges already covered by selected cliques.
3. Add a **2-clique** {u, v} for every remaining uncovered edge.

This guarantees that every edge is covered by at least one clique.

## Clique-Node Hypergraph Construction

After the clique cover is computed:

1. For each vertex *v*, collect all cliques containing *v* → this forms a
   hyperedge (net).
2. Deduplicate cliques per vertex via `sort + unique`.
3. Skip vertices with zero cliques (isolated vertices).
4. All hypernodes (cliques) have **unit weight** (1).
5. All hyperedges (nets) have **unit weight** (1).
6. The `net_to_vertex_` mapping enables reverse lookup from hyperedge to
   original vertex.

## Data Structures

### CliqueCover

Dual CSR representation:

- **Clique-centric**: `clique_ptr_[c]` / `clique_members_` — "what vertices
  are in this clique?"
- **Vertex-centric**: `vertex_cliques_ptr_[v]` / `vertex_cliques_` — "which
  cliques contain this vertex?"

Both provide zero-copy `std::span` access.

### CliqueCoverSolver::Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_maximal_cliques` | `true` | Use BK (true) or triangles only (false) |
| `use_parallel` | `true` | Enable OpenMP parallelism |
| `max_clique_enum_vertices` | 5000 | Max vertices for BK |
| `max_clique_enum_edges` | 100000 | Max edges for BK |
| `greedy_only` | `false` | Skip BK entirely, use triangles only |
| `num_threads` | -1 | Thread count (-1 = auto) |

## Performance Considerations

- **Bron-Kerbosch** is exponential in worst case but bounded by the
  degeneracy ordering and the 5000/100000 thresholds.
- **Triangle enumeration** is O(m^(3/2)) in worst case but practical for
  sparse graphs.
- **Greedy covering** with `std::set` edge tracking is O(k log E) where k
  is the number of candidate cliques. For very large graphs with many
  candidates, this can dominate runtime.
- **CNH construction** is O(total_pins) with a sort+unique per vertex.
