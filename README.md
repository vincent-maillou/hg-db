# hypergraph-reordering

## Dependencies

- GCC 14+ (C++20)
- CMake 3.16+
- SuiteSparse
- BLAS and LAPACK
- [MT-KaHyPar](https://github.com/kahypar/mt-kahypar) (and it's dependencies; included as a submodule, builds automatically)

## Building

```bash
git clone --recurse-submodules https://github.com/ohcomely/hypergraph-reordering.git
cd hypergraph-reordering/
mkdir build && cd build
cmake ..
make -j$(nproc)
```

On an HPC cluster without system TBB:

```bash
cmake .. -DKAHYPAR_DOWNLOAD_TBB=ON -DCMAKE_EXE_LINKER_FLAGS="-lpthread"
make -j$(nproc)
```

If SuiteSparse is not in a standard location, set `-DSUITESPARSE_ROOT=/path/to/suitesparse`. If using conda, `$CONDA_PREFIX` is picked up automatically.

## Usage

```
hypergraph_reorder <input_matrix> <k> [options]
```

`k` is the number of diagonal blocks. The input format is detected from the file extension (`.mtx` for Matrix Market, `.graph` for METIS, otherwise binary CSR).

**Options:**

| Flag | Description |
|---|---|
| `--preset <name>` | MT-KaHyPar preset: `default`, `quality`, `deterministic`, `large_k` |
| `--threads <n>` | Number of threads (default: auto) |
| `--quiet` | Suppress partitioner output |
| `--output-perm` | Write permutation to `<input>.perm` |
| `--output-graph` | Write METIS graph to `<input>.graph` |

The reordered matrix is written to `<input>_reordered.mtx`.
