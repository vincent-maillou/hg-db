# HG-DB — Hypergraph-based Doubly-Bordered Reordering

This library implements a doubly-bordered (DB) matrix reordering via
edge-clique cover (ECC), clique-node hypergraph (CNH) construction, and
MT-KaHyPar partitioning.  It is used as a C/C++ dependency of the
[parallax](https://github.com/vincent-maillou/parallax) Python package.

## Dependencies

- GCC 14+ (C++20)
- CMake 3.16+
- [MT-KaHyPar](https://github.com/kahypar/mt-kahypar) (and its dependencies)

## Building

This library is **not** meant to be built standalone.  It is built
automatically by the parallax CMake build system:

```bash
cd parallax/
cmake -B build -DPARALLAX_PROFILE=default
cmake --build build -j$(nproc)
```

On an HPC cluster without system TBB, the parallax profile should set:

```cmake
set(PARALLAX_MTKAHYPAR_EXTRA_ARGS
    "-DKAHYPAR_DOWNLOAD_TBB=ON"
    "-DKAHYPAR_DISABLE_HWLOC=ON"
    "-DCMAKE_EXE_LINKER_FLAGS='-lpthread'")
```

The parallax build system handles TBB installation and RPATH setup
automatically.

## Python Interface

This library is exposed via `parallax.ordering.general.HGDB`:

```python
from parallax.ordering.general import HGDB, HGDBConfig

config = HGDBConfig(n_parts=8)
hgdb = HGDB(config=config)
result = hgdb.order(sparse_matrix)
```

The result contains a `DoublyBordered` structural descriptor with `General`
leaves for each diagonal block and the global separator.  Per-block refinement
is delegated to parallax's `ComposedOrdering`:

```python
from parallax.ordering.composed import ComposedOrdering
from parallax.ordering.general import AMD

composed = ComposedOrdering(HGDB(config))
composed.add_refinement(AMD())
result = composed.order(sparse_matrix)
```

## Algorithm

1. Convert CSR matrix to an undirected graph.
2. Compute an edge-clique cover (ECC) — parallel Bron-Kerbosch with pivoting
   for small graphs, parallel triangle enumeration for larger graphs.
3. Build a clique-node hypergraph (CNH) from the ECC.
4. Partition the CNH with MT-KaHyPar.
5. Detect vertex separators from the partition.
6. Construct the DB permutation: diagonal blocks first, separator last.

See [`ECC.md`](ECC.md) for details on the clique-cover algorithm.
