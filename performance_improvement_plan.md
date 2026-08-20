# Clique Cover / Graph / Hypergraph Performance Improvement Plan

## 1. Purpose

This document updates the original performance-improvement list in light of the Bron–Kerbosch changes already implemented.

The current implementation has already addressed the four most important Bron–Kerbosch issues:

1. **Degeneracy ordering is now represented explicitly by `rank[v]`.**
2. **The outer Bron–Kerbosch initialization now uses**
   - `P = later neighbors`
   - `X = earlier neighbors`.
3. **Pivot membership tests no longer use repeated `std::find`;** a marker array gives O(1) membership.
4. **Pivot selection considers `P ∪ X`.**
5. **`std::set` has been removed from the recursive Bron–Kerbosch path**, with sorted adjacency lists and two-pointer intersections used for `P ∩ N(v)`, `X ∩ N(v)`, and `P \ N(pivot)`.

These changes should be treated as **completed** and should not be reimplemented.

The remaining work should now focus on the edge-cover phase, graph construction/ordering, allocation behavior, and only then further Bron–Kerbosch micro-optimizations.

---

# 2. Current implementation assessment

## 2.1 Bron–Kerbosch

The current implementation is substantially better than the original one.

The important structural properties are now correct:

```cpp
auto ordering = graph.compute_degeneracy_ordering();

std::vector<index_t> rank(graph.n_vertices());

for (index_t i = 0; i < ordering.size(); ++i)
    rank[ordering[i]] = i;
```

and:

```cpp
if (rank[u] > rank[v])
    P.push_back(u);
else
    X.push_back(u);
```

The recursive implementation also uses a marker array:

```cpp
std::vector<uint64_t> p_marker(graph.n_vertices(), 0);
```

so membership testing during pivot selection is O(1).

The pivot is now selected from both `P` and `X`, and the intersections use sorted adjacency lists.

### Remaining BK issue

The principal remaining BK issue is **allocation and copying inside recursion**.

Every recursive call still creates:

```cpp
std::vector<index_t> P_new;
std::vector<index_t> X_new;
std::vector<index_t> candidates;
```

and copies the relevant subsets into them.

This is not an algorithmic correctness problem, and it should be addressed only after the edge-cover and graph-level issues below have been measured.

---

# 3. Priority overview

| Priority | Improvement | Area | Expected impact | Recommendation |
|---|---|---|---|---|
| ★★★★★ | Replace edge `std::set` with edge IDs + coverage array | Clique cover | Very high | Implement next |
| ★★★★★ | Reuse one coverage state across all cover phases | Clique cover | High | Implement together with edge IDs |
| ★★★★☆ | Remove unnecessary `has_edge()` calls inside cliques | Clique cover | Medium, very easy | Implement with edge IDs |
| ★★★★☆ | Improve degeneracy-order computation | Graph | Potentially high | Profile, then implement |
| ★★★★☆ | Reduce BK recursive allocations | Clique enumeration | High on difficult BK instances | Next BK optimization |
| ★★★☆☆ | Store candidate clique edge IDs | Clique cover | High if greedy selection dominates | Implement if profiling justifies memory |
| ★★★☆☆ | Dynamic greedy gains | Clique cover | Better cover quality, more work | Optional algorithmic improvement |
| ★★★☆☆ | Degree-oriented triangle enumeration | Triangle enumeration | Medium/high on skewed graphs | Good next optimization |
| ★★★☆☆ | Compact triangle representation | Triangle enumeration | Medium | Easy follow-up |
| ★★★☆☆ | Remove redundant sorting/deduplication in hypergraph construction | Hypergraph | Low/medium | Implement after invariants are explicit |
| ★★☆☆☆ | Reduce `CliqueCover` construction allocations | Clique cover representation | Medium on large covers | Later |
| ★★☆☆☆ | Improve OpenMP scheduling | Parallel enumeration | Workload-dependent | Benchmark, do not assume |
| ★★☆☆☆ | Dense remapping instead of `unordered_map` in induced subgraphs | Graph | Medium | Easy if this path is hot |
| ★★☆☆☆ | Galloping intersection/search for triangles | Triangle enumeration | Workload-dependent | Later |

The most important change from the original list is that **edge coverage is now the main remaining structural bottleneck outside BK**.

---

# 4. Improvement 5 — Remove repeated `graph.has_edge()` inside clique loops

## Current situation

`select_covering_cliques()` currently performs:

```cpp
for (size_t i = 0; i < clique.size(); ++i) {
    for (size_t j = i + 1; j < clique.size(); ++j) {
        ...
        if (graph.has_edge(u, v) &&
            covered_edges.find({u, v}) == covered_edges.end()) {
```

and later performs the same test again while marking edges.

For a clique, every pair of vertices is an edge by definition.

Therefore, once the input invariant is:

> every candidate in `maximal_cliques` is actually a clique,

the `graph.has_edge()` calls are redundant.

The same is true for the triangle candidates.

### Recommended change

Do not perform:

```cpp
graph.has_edge(u, v)
```

inside clique pair loops.

Instead, obtain the edge ID directly for `(u,v)`.

This should be implemented together with the edge-ID redesign described below.

### Priority

**★★★★☆ — easy and essentially free once edge IDs exist.**

---

# 5. Improvement 6 — Replace `std::set<pair<index_t,index_t>>` with edge IDs

## Current situation

The greedy cover currently uses:

```cpp
std::set<std::pair<index_t, index_t>> covered_edges;
```

Every coverage operation therefore performs tree-based lookup/insertion.

This is unnecessary because the graph already has a fixed set of `m = graph.n_edges()` edges.

## Recommended design

Give every undirected graph edge a stable integer ID:

```text
(u, v)  --->  edge_id
```

and maintain:

```cpp
std::vector<bool> covered;
```

or, preferably if the graph is not extremely large:

```cpp
std::vector<uint8_t> covered;
```

depending on memory/performance measurements.

Then:

```cpp
if (!covered[eid]) {
    ...
}
```

is O(1).

## Graph API

A natural extension is:

```cpp
index_t edge_id(index_t u, index_t v) const;
```

and potentially:

```cpp
std::span<const index_t> edge_ids(index_t v) const;
```

where the edge IDs are aligned with the existing sorted adjacency list.

For example:

```text
adj_list_[p]       = neighbor
adj_edge_ids_[p]   = corresponding undirected edge ID
```

The same edge ID appears in both endpoint adjacency lists.

This preserves the current CSR structure while giving every adjacency entry an O(1) associated edge ID once its adjacency position is known.

An implementation can initially use binary search in the sorted adjacency list to find the position of `(u,v)`. If this becomes hot, the graph representation can be extended so that the clique-generation path can obtain edge IDs more directly.

## Important design point

Do **not** introduce a general-purpose:

```cpp
std::unordered_map<pair<index_t,index_t>, index_t>
```

just to replace the `std::set`.

That would remove logarithmic tree operations but introduce hashing, random memory access, and potentially large memory overhead.

The graph is already a compressed sparse representation, so the edge-ID representation should preferably remain CSR-oriented.

### Priority

**★★★★★ — one of the highest-value remaining changes.**

---

# 6. Improvement 7/8 — One coverage state for the entire edge-cover phase

## Current situation

`solve()` creates:

```cpp
std::vector<bool> covered_edges(graph.n_edges(), false);
```

but this array is never modified.

Instead, `select_covering_cliques()` creates a completely separate:

```cpp
std::set<std::pair<index_t, index_t>> covered_edges;
```

This creates duplicated coverage logic.

## Recommended design

Coverage should have a single owner.

Conceptually:

```cpp
std::vector<uint8_t> covered(graph.n_edges(), false);
```

Then:

```text
find candidates
       |
       v
select_covering_cliques()
       |
       | updates covered[]
       v
cover_remaining_edges()
       |
       | reads covered[]
       v
final clique cover
```

A possible interface is:

```cpp
select_covering_cliques(
    graph,
    candidates,
    covered_edges);
```

followed by:

```cpp
cover_remaining_edges(
    graph,
    covered_edges,
    cliques);
```

Alternatively, the coverage state can be encapsulated in a small internal helper class.

The important point is that **there must be exactly one coverage representation**.

### Priority

**★★★★★ — implement together with edge IDs.**

---

# 7. Improvement 9 — Replace static clique-size ordering with a better greedy criterion

## Current situation

Candidates are currently sorted once:

```cpp
std::sort(clique_sizes.rbegin(), clique_sizes.rend());
```

and then a clique is selected if it covers at least one new edge.

This is cheap, but it does not implement the usual greedy set-cover heuristic.

A large clique can become almost useless after earlier cliques have covered most of its edges.

## Better criterion

For each candidate clique `C`, define:

```text
gain(C) = number of currently uncovered edges in C
```

Then a natural greedy rule is:

```text
select clique with maximum gain
```

or, if clique cost is not uniform:

```text
select maximum gain / cost
```

Since the current objective assigns equal cost to each selected clique, `gain` is the natural first criterion.

## Important tradeoff

A dynamically updated priority queue can improve the resulting cover, but it introduces bookkeeping.

Therefore this is **not automatically a performance optimization**.

It trades:

```text
less computation + potentially worse cover
```

for:

```text
more computation + potentially better cover
```

This should be treated as an **algorithmic-quality option**, not simply a low-level optimization.

### Priority

**★★★☆☆ — optional.**

Implement only after measuring whether the number/size of resulting cliques matters to the downstream reordering.

---

# 8. Improvement 10 — Dynamic greedy gains

A stronger version of the previous improvement maintains candidate gains dynamically.

A useful implementation strategy is a **lazy priority queue**:

1. Initially compute each candidate's gain.
2. Put candidates in a max-heap.
3. Pop the candidate with the largest stored gain.
4. Recompute its actual gain against the current `covered[]`.
5. If the stored gain is stale, update it and reinsert it.
6. Otherwise select it.

This avoids explicitly updating every candidate whenever an edge becomes covered.

## Why this is attractive

The edge coverage state is monotonic:

```text
uncovered -> covered
```

so stale gains can only decrease.

This makes lazy heap evaluation natural.

## Caveat

If there are very many maximal cliques, the heap itself may become expensive.

Therefore this should be benchmarked against the current static size ordering.

### Priority

**★★★☆☆.**

---

# 9. Improvement 11 — Triangle enumeration

The current triangle implementation is already structurally sound:

- adjacency lists are sorted,
- `upper_bound()` skips irrelevant prefixes,
- intersections use two-pointer traversal,
- duplicates are avoided.

It should therefore not be rewritten prematurely.

## 9.1 Degree-based orientation

The current orientation is:

```cpp
u < v
```

This is primarily a duplicate-elimination rule.

A stronger orientation is based on:

```text
(degree(v), vertex_id)
```

rather than vertex ID alone.

For every edge `(u,v)`, orient it from the lower-ranked endpoint to the higher-ranked endpoint under:

```text
(degree, id)
```

Then triangle enumeration intersects forward-neighbor lists.

This is the standard degeneracy/degree-oriented approach used by efficient triangle enumeration algorithms.

The main benefit is that forward adjacency lists are substantially smaller on high-degree vertices, particularly for graphs with skewed degree distributions.

### Priority

**★★★☆☆ — potentially significant for sparse, power-law-like graphs.**

This should be implemented after edge coverage and graph-ordering improvements unless triangle enumeration is already a measured bottleneck.

---

# 10. Improvement 11b — Galloping search for triangle intersections

Two-pointer intersection is excellent when the lists have comparable sizes.

When:

```text
|A| << |B|
```

a galloping/exponential search can reduce the amount of work required to search `B`.

A hybrid implementation can choose between:

```text
two-pointer merge
```

and:

```text
binary/exponential search
```

based on the relative list sizes.

This is a classic constant-factor optimization, but it adds branching and implementation complexity.

### Recommendation

Do not implement initially.

Benchmark degree-oriented intersections first.

### Priority

**★★☆☆☆.**

---

# 11. Improvement 12 — Reduce Bron–Kerbosch recursive allocations

## Current situation

The current recursive function creates:

```cpp
std::vector<index_t> candidates;
std::vector<index_t> P_new;
std::vector<index_t> X_new;
```

at every recursion level.

This is now the most important remaining optimization inside BK itself.

The algorithm can spend substantial time allocating, copying, and freeing these vectors even though the graph operations themselves have become much cheaper.

## Possible next implementation

Move toward an in-place Tomita-style representation using reusable work arrays.

The general idea is:

```text
one large working buffer
+
index/range boundaries
+
in-place partitioning
```

rather than:

```text
allocate P_new
allocate X_new
copy
recurse
destroy
```

This is substantially more complicated than the current implementation.

## Recommended staged approach

Do not immediately rewrite BK.

First profile:

- number of recursive calls,
- total sizes of `P_new`,
- total sizes of `X_new`,
- total candidates generated,
- time spent in allocation/copying.

If these dominate, implement reusable work buffers.

### Priority

**★★★★☆ for difficult BK instances, but after edge coverage.**

---

# 12. Improvement 13 — OpenMP scheduling

The current outer loop uses:

```cpp
#pragma omp for schedule(dynamic, 1)
```

This is a reasonable choice because Bron–Kerbosch work is highly irregular.

However, `dynamic,1` has relatively high scheduling overhead.

Potential alternatives include:

```cpp
schedule(dynamic, 4)
schedule(dynamic, 8)
schedule(guided)
```

or explicit OpenMP tasks.

## Important qualification

There is no universally superior scheduling policy here.

A larger chunk size improves scheduling overhead but can worsen load balance because one degeneracy-order vertex can generate vastly more work than another.

Therefore this should be benchmarked rather than changed on principle.

## Additional issue

The call:

```cpp
omp_set_num_threads(num_threads);
```

changes OpenMP's global/default thread setting.

If the solver is used as a library, it may be preferable to avoid changing global OpenMP state and instead control parallelism through the parallel region or a documented application-level policy.

### Priority

**★★☆☆☆.**

---

# 13. Improvement 14 — Compact triangle representation

The current triangle representation is:

```cpp
std::vector<std::vector<index_t>>
```

so every triangle is a separate dynamically allocated vector.

A triangle has fixed cardinality, so a much better internal representation is:

```cpp
using Triangle = std::array<index_t, 3>;
```

This removes per-triangle heap allocations and improves locality.

The public API can still convert the final result to:

```cpp
std::vector<std::vector<index_t>>
```

if required.

This is particularly attractive because triangle enumeration can generate a very large number of objects.

### Priority

**★★★☆☆.**

---

# 14. Improvement 15 — Store clique edge IDs

This is a potentially important extension of the edge-ID redesign.

Currently, when evaluating a clique, the code repeatedly generates:

```text
all vertex pairs in the clique
```

and then determines whether each edge is covered.

If edge IDs are available, each candidate clique can instead carry:

```cpp
std::vector<index_t> edge_ids;
```

where:

```text
edge_ids = all edges contained in the clique
```

Then the greedy phase becomes:

```cpp
for (auto eid : clique.edge_ids) {
    if (!covered[eid])
        ++gain;
}
```

and selection becomes:

```cpp
for (auto eid : clique.edge_ids)
    covered[eid] = true;
```

## Benefit

The edge IDs are computed once instead of regenerating the same `(u,v)` pairs every time a candidate is examined.

This becomes especially useful for dynamic greedy gains, because the same candidate can be evaluated many times.

## Memory tradeoff

The total storage is:

```text
sum_C |E(C)|
```

where:

```text
|E(C)| = |C| choose 2
```

For large maximal cliques this can be substantial.

Therefore this optimization should be enabled only if profiling shows that repeated clique-edge traversal dominates.

### Priority

**★★★☆☆ — highly useful if dynamic greedy selection is adopted.**

---

# 15. Improvement 16 — Degeneracy ordering itself can be improved

## Current situation

`Graph::compute_degeneracy_ordering()` uses:

```cpp
std::vector<std::set<index_t>> buckets(max_degree + 1);
```

and moves vertices between sets as their current degree changes.

This is correct in spirit, but it introduces tree allocations and logarithmic operations for every degree update.

Since degeneracy ordering is required before every full BK enumeration, it can become a measurable preprocessing cost.

## Better options

### Option A — priority queue

A relatively simple improvement is a min-heap containing:

```cpp
(current_degree, vertex)
```

with lazy deletion of stale entries.

This is still O((n+m) log n), but often has better implementation simplicity and memory behavior than a vector of many `std::set`s.

### Option B — bucket queue

A true bucket-based implementation can achieve near-linear behavior:

```text
O(n + m)
```

for the degeneracy ordering.

It requires more careful data structures, typically intrusive linked lists or arrays storing bucket positions.

### Recommendation

Do not jump directly to the most complicated bucket implementation.

First benchmark the current ordering time against:

```text
BK time
triangle enumeration time
cover-selection time
```

If degeneracy ordering is a significant fraction of runtime, replace the `std::set` buckets.

### Priority

**★★★★☆ if preprocessing is significant; otherwise ★★★☆☆.**

---

# 16. Improvement 17 — `Graph::induced_subgraph()`

The current implementation uses:

```cpp
std::unordered_map<index_t, index_t> old_to_new;
```

for vertex remapping.

Since graph vertices are already dense integer IDs:

```text
0 ... n-1
```

a dense array is usually a better fit:

```cpp
std::vector<index_t> old_to_new(n_vertices_, invalid);
```

or a marker/remapping array.

This avoids hash-table allocation and random access.

The transformation then becomes essentially:

```cpp
old_to_new[v] = new_id;
```

followed by direct array lookup.

### Priority

**★★☆☆☆**, unless `induced_subgraph()` is used heavily.

---

# 17. Improvement 18 — `CliqueCover` construction allocates many small vectors

The constructor currently creates:

```cpp
std::vector<std::vector<index_t>> vertex_to_cliques(n_vertices);
```

and fills it before converting it into CSR-like storage.

This creates one vector object per graph vertex and potentially many individual allocations.

The final representation is already compressed:

```cpp
vertex_cliques_ptr_
vertex_cliques_
```

so the temporary vector-of-vectors is not strictly necessary.

## Better approach

Use two passes:

### Pass 1

Count the number of clique memberships for every vertex.

### Pass 2

Build:

```cpp
vertex_cliques_ptr_
vertex_cliques_
```

directly.

This gives the same final representation without the intermediate vector-of-vectors.

### Priority

**★★☆☆☆ to ★★★☆☆**, depending on cover size.

---

# 18. Improvement 19 — Remove redundant sorting in `Hypergraph::from_clique_cover`

The current code does:

```cpp
std::vector<index_t> unique_cliques(cliques.begin(), cliques.end());
std::sort(unique_cliques.begin(), unique_cliques.end());
auto last = std::unique(unique_cliques.begin(), unique_cliques.end());
```

However, `CliqueCover` constructs `vertex_to_cliques` by traversing:

```cpp
for (index_t cid = 0; cid < n_cliques_; ++cid)
```

and appending `cid`.

Therefore, under the current representation, each vertex's clique list is already sorted by clique ID.

The sort is consequently redundant.

The remaining question is whether duplicate clique IDs can occur.

If the project invariant is:

> each clique contains a vertex at most once,

then duplicates cannot arise from a valid `CliqueCover`, and the entire temporary vector/sort/unique operation can be removed.

## Recommended design

Make the invariant explicit:

```text
Every clique contains unique vertex IDs.
Every vertex-to-clique list is sorted by clique ID.
```

Then `from_clique_cover()` can directly append the span.

If robustness against externally constructed invalid clique lists is important, validate the invariant once when constructing `CliqueCover` rather than repeatedly sorting during hypergraph construction.

### Priority

**★★★☆☆ — simple and directly relevant to `hypergraph.cpp`.**

---

# 19. Potential correctness issue to verify in `Graph::from_symmetric_matrix`

This is not merely a performance optimization and should be checked before benchmarking.

The current code treats:

```cpp
matrix.is_symmetric()
```

specially while iterating over every stored matrix entry.

If `is_symmetric()` means the usual situation where both `(i,j)` and `(j,i)` are physically present in CSR, then the current logic appears to process the same undirected edge twice.

In particular, the degree construction contains:

```cpp
if (i < j) {
    ...
} else if (i == j) {
    continue;
} else {
    if (matrix.is_symmetric()) {
        degrees[i]++;
        degrees[j]++;
        n_edges++;
    }
}
```

and the adjacency construction contains:

```cpp
if (i != j) {
    if (i < j || matrix.is_symmetric()) {
        ...
    }
}
```

For a conventionally stored symmetric matrix, both the upper and lower entries would satisfy the condition.

### Action

Before performance work, verify the exact semantics of:

```cpp
CSRMatrix::is_symmetric()
```

and whether the input matrix stores:

1. both triangular parts, or
2. only one triangular part with symmetric semantics.

If both parts are stored, this construction needs correction.

### Priority

**★★★★★ if confirmed.**

Correctness must precede optimization.

---

# 20. Suggested edge-ID architecture

The edge-ID redesign is the most useful structural change remaining.

A good target architecture is:

```text
Graph
 ├── adj_ptr_
 ├── adj_list_
 └── adj_edge_ids_
          |
          +---- undirected edge ID
```

For every undirected edge:

```text
u ---- v
 \      /
  \    /
   edge_id = e
```

both adjacency entries refer to the same `e`.

Then the cover algorithm becomes:

```text
candidate clique
      |
      v
edge IDs of clique
      |
      v
covered[e]
      |
      +---- gain
      |
      +---- mark covered
```

This eliminates the current repeated dependence on:

```cpp
std::set<std::pair<index_t,index_t>>
graph.has_edge(u,v)
```

and makes the edge-centric nature of the algorithm explicit.

---

# 21. Recommended implementation order

## Phase 0 — Validate invariants and correctness

Before optimization:

- verify `Graph::from_symmetric_matrix()`;
- verify every clique contains unique vertices;
- verify every candidate clique really is a clique;
- verify adjacency lists are sorted;
- verify `n_edges()` equals the number of undirected edges;
- verify edge IDs are unique and consistent at both endpoints.

Add debug-only assertions where useful.

---

## Phase 1 — Edge coverage redesign

Implement:

1. stable undirected edge IDs;
2. graph edge-ID lookup;
3. `covered[eid]`;
4. one shared coverage array;
5. remove `std::set<pair<...>>`;
6. remove `graph.has_edge()` from clique pair loops;
7. make `cover_remaining_edges()` a single graph-edge scan.

Target structure:

```text
select_covering_cliques()
        |
        | updates covered[]
        v
cover_remaining_edges()
        |
        | scans graph edges once
        v
final cover
```

This should be the **next major implementation step**.

---

## Phase 2 — Measure the new baseline

Record separately:

```text
degeneracy ordering
triangle enumeration
BK enumeration
greedy clique selection
remaining-edge covering
CliqueCover construction
Hypergraph construction
total solve time
```

Also record:

```text
number of maximal cliques
number of selected cliques
maximum clique size
average clique size
number of uncovered edges
number of 2-cliques added
```

Without this breakdown, later optimizations risk optimizing the wrong component.

---

## Phase 3 — Graph preprocessing

If degeneracy ordering is significant:

1. replace `vector<set>` buckets with a heap-based implementation;
2. benchmark;
3. if still important, consider a true bucket queue.

Also optimize `induced_subgraph()` if it is on a hot path.

---

## Phase 4 — Triangle enumeration

If triangle enumeration is significant:

1. degree-oriented forward adjacency;
2. compact triangle storage;
3. benchmark;
4. only then consider galloping search.

---

## Phase 5 — Greedy cover quality

If the quality of the final clique cover matters:

1. compare static clique-size ordering;
2. compare maximum-current-gain;
3. test lazy priority-queue gains;
4. measure both runtime and resulting:
   - number of cliques,
   - number of 2-cliques,
   - total clique memberships,
   - downstream reordering quality.

Do not assume a more sophisticated greedy rule is better overall until the downstream objective confirms it.

---

## Phase 6 — BK allocation optimization

If BK remains dominant after the structural changes:

1. measure recursive allocation/copy volume;
2. introduce reusable work buffers;
3. move toward in-place partitioning;
4. benchmark against the current vector-based recursion.

This is the point where a more substantial Tomita-style implementation becomes justified.

---

## Phase 7 — Representation cleanup

Finally:

- remove temporary vector-of-vectors in `CliqueCover`;
- remove redundant sorting/deduplication in `Hypergraph`;
- use compact triangle storage;
- review OpenMP scheduling;
- remove unnecessary includes.

---

# 22. Benchmarking plan

Every optimization should be evaluated on at least three graph regimes:

### A. Low-degeneracy sparse graphs

Expected behavior:

- BK should be relatively manageable;
- degeneracy ordering should be cheap;
- triangle orientation may help significantly.

### B. High-degree / skewed graphs

Expected behavior:

- triangle enumeration can become dominant;
- degree orientation should have larger benefits;
- edge coverage may contain large cliques.

### C. Dense graphs near the BK threshold

Expected behavior:

- BK recursion can dominate;
- recursive allocation reduction becomes important;
- pivot quality and degeneracy ordering matter strongly.

For every benchmark, record:

```text
N
M
degeneracy
number of triangles
number of maximal cliques
maximum clique size
average maximal-clique size
selected cliques
2-cliques added
runtime by phase
peak memory
```

The most important comparison is not only total runtime but **where the runtime moves after each optimization**.

---

# 23. Updated priority list

| Priority | Item | Status |
|---|---|---|
| ★★★★★ | Correct degeneracy ordering / `rank[]` | **DONE** |
| ★★★★★ | Correct `P`/`X` initialization | **DONE** |
| ★★★★★ | O(1) BK membership marker | **DONE** |
| ★★★★★ | Pivot from `P ∪ X` | **DONE** |
| ★★★★★ | Remove recursive `std::set`/`std::find` | **DONE** |
| ★★★★★ | Verify/fix symmetric-matrix graph construction | **CHECK NOW** |
| ★★★★★ | Edge IDs + coverage array | **NEXT** |
| ★★★★★ | Single shared coverage state | **NEXT** |
| ★★★★☆ | Remove `has_edge()` inside clique loops | **NEXT** |
| ★★★★☆ | Improve degeneracy-order data structure | **PROFILE → IMPLEMENT** |
| ★★★★☆ | Reduce BK recursive allocations | **PROFILE → IMPLEMENT** |
| ★★★☆☆ | Degree-oriented triangle enumeration | **LATER** |
| ★★★☆☆ | Compact triangle storage | **LATER** |
| ★★★☆☆ | Store candidate clique edge IDs | **IF COVER SELECTION IS HOT** |
| ★★★☆☆ | Dynamic/lazy greedy gains | **OPTIONAL QUALITY IMPROVEMENT** |
| ★★★☆☆ | Remove redundant hypergraph sorting/deduplication | **LATER** |
| ★★☆☆☆ | `CliqueCover` CSR construction without temporary vectors | **LATER** |
| ★★☆☆☆ | Dense remapping in `induced_subgraph()` | **LATER** |
| ★★☆☆☆ | OpenMP scheduling experiments | **BENCHMARK** |
| ★★☆☆☆ | Galloping triangle intersection | **LATE OPTIMIZATION** |

---

# 24. Final recommended roadmap

The project should now proceed in this order:

```text
                 ┌─────────────────────────────┐
                 │ Validate graph/clique       │
                 │ invariants and symmetry      │
                 └──────────────┬──────────────┘
                                │
                                v
                 ┌─────────────────────────────┐
                 │ Edge IDs + shared coverage  │
                 │ state                       │
                 └──────────────┬──────────────┘
                                │
                                v
                 ┌─────────────────────────────┐
                 │ Remove has_edge/set/edge    │
                 │ reconstruction              │
                 └──────────────┬──────────────┘
                                │
                                v
                 ┌─────────────────────────────┐
                 │ Establish phase-level       │
                 │ performance baseline        │
                 └──────────────┬──────────────┘
                                │
                ┌───────────────┴────────────────┐
                v                                v
       ┌──────────────────┐             ┌──────────────────┐
       │ Graph ordering   │             │ Triangle path    │
       │ if hot           │             │ if hot           │
       └────────┬─────────┘             └────────┬─────────┘
                │                                │
                └───────────────┬────────────────┘
                                v
                 ┌─────────────────────────────┐
                 │ Profile BK allocation      │
                 └──────────────┬──────────────┘
                                │
                                v
                 ┌─────────────────────────────┐
                 │ In-place / reusable-buffer  │
                 │ BK if justified             │
                 └──────────────┬──────────────┘
                                │
                                v
                 ┌─────────────────────────────┐
                 │ Optional cover-quality and  │
                 │ representation improvements │
                 └─────────────────────────────┘
```

The central principle is now:

> **Do not optimize individual loops in isolation. First make the representation edge-centric and eliminate repeated work; then profile the resulting implementation before introducing more sophisticated algorithms.**

The Bron–Kerbosch implementation has already crossed that first structural threshold. The next major structural threshold is the **edge-ID/coverage redesign**.