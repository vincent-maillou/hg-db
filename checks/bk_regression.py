#!/usr/bin/env python3
"""BK correctness regression check for the HG-DB edge-clique cover.

Standalone ctypes script (no pytest, no parallax import). Locates the library
built by the parallax superproject at external/builds/hgdb/lib/ relative to
the parallax root (three levels up from this file's directory).

ABI note (v5): this script used to carry a v2-era ctypes mirror. Against any
build >= v3, hgr_get_statistics wrote the full newer struct PAST the end of
the smaller v2 buffer — an out-of-bounds write / heap corruption (the script
has been unsafe since v3). The mirrors below are now the v5 layout
(SOVERSION 5) and must move in lockstep with include/hypergraph_reorder.h
and src/parallax/ordering/general/hgdb.py::_HGRStatistics.

What it checks
--------------
1. Structural invariants of the ECC result for all three methods
   (BK, TriEnum, EdgeCover) on small fixtures:
     - the returned permutation is a valid permutation of range(n),
     - sum(clique-order histogram) == n_cliques,
     - max(clique orders) == max_clique_size,
     - all clique orders are >= 2,
     - avg_clique_size is consistent with the histogram,
     - EdgeCover produces exactly one 2-clique per graph edge.
2. BK ground-truth pinning: on fixtures whose degeneracy ordering differs
   from the vertex numbering (this is exactly where the historical
   `u > v` P-construction with X = {} was wrong), the exact BK candidate
   set (all maximal cliques, nothing else) feeds the Phase-2/3 greedy cover.
   The expected aggregate statistics below were derived from an independent
   ground-truth simulation (plain maximal-clique enumeration + the exact C++
   greedy semantics: size-DESC / index-DESC sort order, ascending remainder)
   and confirmed against the fixed library.
3. v5 observables (BK only): degeneracy d, n_maximal_cliques and the
   maximal-clique-enumeration histogram, pinned against the same
   independent ground truth (bitmask brute-force maximal cliques +
   tie-break-free smallest-last degeneracy). Cross-checks:
   sum(enum_hist) == n_maximal_cliques and, for every order >= 3, the
   enumeration count dominates the cover count (Phase 3 only ever adds
   2-cliques).

Why the BK pins matter (red/green)
----------------------------------
The pre-fix BK (root P = higher-ID neighbors, X empty, as found at the
v2.0.0 merge point) does NOT miss any maximal clique, but it ADDS spurious
non-maximal cliques (sub-maximal R terminating with P = X = {} because X is
never seeded at the root). Those extra candidates shift the greedy
tie-breaks, which observably changes the clique-order histogram, so the
enumeration pins (n_maximal_cliques / enum_hist) and the cover pins fail on
a pre-fix build. (Historical v2-era measurement at the merge commit: F1
n_cliques 9, F2 10, F3 8 with histogram {2: 1, 3: 2, 4: 5} - the numbers
were relative to the then-current pin values.)

Note: a plain edge-coverage assertion ("every edge in some final clique")
cannot discriminate the two versions - Phase 3 (cover_remaining_edges)
patches every uncovered edge with a 2-clique, and an exhaustive audit over
all graphs with n <= 6 (plus randomized n = 7..12) found no graph on which
the pre-fix code produces an invalid or non-covering final cover.

The n_cliques(BK) <= n_cliques(TriEnum) relation does NOT hold robustly on
these fixtures (F1: BK 10 > TriEnum 9; F2: BK 11 > TriEnum 10) and is
therefore not asserted; it is reported for information only.

The per-fixture `degeneracy_ordering` entries below are INFORMATIONAL only:
the v5.0.0 Matula-Beck rewrite is linear-time and free in its within-bucket
extraction order, so the ordering it produces may differ; the degeneracy
value, the maximal-clique set, and (via the v4.1.0 Phase-2
canonicalization) the final cover are invariant — those are what is pinned.

Usage (from the parallax repo root):
    LD_PRELOAD=/usr/pack/gcc-15.3.0-af/lib64/libstdc++.so.6 \
        python external/HG-DB/checks/bk_regression.py
"""

import ctypes
import os
import sys

HGR_ECC_EDGE_COVER = 0
HGR_ECC_TRI_ENUM = 1
HGR_ECC_BK = 2

METHOD_NAMES = {
    HGR_ECC_EDGE_COVER: "EdgeCover",
    HGR_ECC_TRI_ENUM: "TriEnum",
    HGR_ECC_BK: "BK",
}

# ---------------------------------------------------------------------------
# Fixtures
#
# Small graphs whose degeneracy ordering clearly differs from the vertex
# numbering (computed with the exact Matula-Beck / std::set-bucket semantics
# of Graph::compute_degeneracy_ordering). This is precisely the situation in
# which the pre-fix `u > v` P-construction with X = {} was wrong.
# ---------------------------------------------------------------------------
FIXTURES = {
    "F1_n9": {
        "n": 9,
        "edges": [
            (0, 3),
            (0, 6),
            (0, 7),
            (0, 8),
            (1, 3),
            (1, 4),
            (1, 5),
            (1, 6),
            (1, 7),
            (1, 8),
            (2, 6),
            (2, 8),
            (3, 6),
            (3, 7),
            (4, 8),
            (5, 6),
            (5, 7),
            (6, 8),
        ],
        "degeneracy_ordering": [2, 4, 5, 7, 0, 3, 1, 6, 8],
    },
    "F2_n10": {
        "n": 10,
        "edges": [
            (0, 4),
            (1, 2),
            (1, 4),
            (1, 6),
            (1, 7),
            (1, 8),
            (2, 4),
            (2, 6),
            (2, 8),
            (2, 9),
            (3, 4),
            (3, 7),
            (3, 9),
            (4, 5),
            (4, 7),
            (4, 9),
            (5, 6),
            (6, 9),
        ],
        "degeneracy_ordering": [0, 5, 8, 3, 7, 1, 4, 2, 6, 9],
    },
    "F3_n9": {
        "n": 9,
        "edges": [
            (0, 1),
            (0, 3),
            (0, 7),
            (0, 8),
            (1, 2),
            (1, 4),
            (1, 5),
            (1, 7),
            (1, 8),
            (2, 4),
            (2, 5),
            (2, 6),
            (2, 7),
            (4, 5),
            (4, 6),
            (4, 7),
            (4, 8),
            (5, 6),
            (5, 8),
            (6, 7),
            (6, 8),
        ],
        "degeneracy_ordering": [3, 0, 7, 1, 2, 4, 5, 6, 8],
    },
}

# Expected aggregate statistics (single-threaded, deterministic).
# BK values are the ground-truth pins; TriEnum/EdgeCover are context pins.
#
# PIN PROVENANCE (all values): independent ground truth computed offline —
# brute-force bitmask maximal-clique enumeration (each subset maximal iff
# complete and no single-vertex superset is complete), brute-force triangle
# enumeration for TriEnum, a tie-break-free pure-Python smallest-last
# degeneracy (the degeneracy value is tie-break invariant), and a pure-
# Python simulation of the CURRENT greedy semantics (v4.1.0 classical rule:
# argmax current marginal, ties by marginal desc / rank asc with rank =
# size desc / canonical key asc; Phase 3 appends 2-cliques in ascending
# (u, v) order). The library must reproduce this exactly.
#
# HISTORY: the cover pins below (n_cliques / max_clique_size / hist) were
# originally v2-era values confirmed against a v2-era library. The v4.1.0
# greedy change (static largest-first -> classical max-current-marginal)
# changed the covers on F2/F3, but the script was unrunnable against any
# build >= v3 (its stale v2 mirror suffered the out-of-bounds write
# described above), so the stale cover pins went unnoticed until this v5
# refresh re-derived them from the current-greedy ground truth. The v5
# observables (degeneracy / n_maximal_cliques / enum_hist) are new pins.
EXPECTED = {
    "F1_n9": {
        "n_edges": 18,
        HGR_ECC_BK: {
            "n_cliques": 8,
            "max_clique_size": 3,
            "hist": {3: 8},
            "degeneracy": 3,
            "n_maximal_cliques": 10,
            "enum_hist": {3: 10},
        },
        HGR_ECC_TRI_ENUM: {"n_cliques": 8, "max_clique_size": 3, "hist": {3: 8}},
        HGR_ECC_EDGE_COVER: {"n_cliques": 18, "max_clique_size": 2, "hist": {2: 18}},
    },
    "F2_n10": {
        "n_edges": 18,
        HGR_ECC_BK: {
            "n_cliques": 10,
            "max_clique_size": 3,
            "hist": {2: 3, 3: 7},
            "degeneracy": 3,
            "n_maximal_cliques": 11,
            "enum_hist": {2: 3, 3: 8},
        },
        HGR_ECC_TRI_ENUM: {"n_cliques": 10, "max_clique_size": 3, "hist": {2: 3, 3: 7}},
        HGR_ECC_EDGE_COVER: {"n_cliques": 18, "max_clique_size": 2, "hist": {2: 18}},
    },
    "F3_n9": {
        "n_edges": 21,
        HGR_ECC_BK: {
            "n_cliques": 6,
            "max_clique_size": 4,
            "hist": {2: 1, 3: 2, 4: 3},
            "degeneracy": 4,
            "n_maximal_cliques": 9,
            "enum_hist": {2: 1, 3: 2, 4: 6},
        },
        HGR_ECC_TRI_ENUM: {"n_cliques": 10, "max_clique_size": 3, "hist": {2: 1, 3: 9}},
        HGR_ECC_EDGE_COVER: {"n_cliques": 21, "max_clique_size": 2, "hist": {2: 21}},
    },
}


# ---------------------------------------------------------------------------
# ctypes mirrors of the v5 C API (hypergraph_reorder.h, SOVERSION 5).
# MIRROR DISCIPLINE: copy the field order of hgr_statistics_t verbatim; a
# stale mirror against a newer .so lets hgr_get_statistics write past this
# buffer (this script's own v2-era mirror did exactly that — see the header
# comment). Keep in lockstep with include/hypergraph_reorder.h and
# src/parallax/ordering/general/hgdb.py::_HGRStatistics.
# ---------------------------------------------------------------------------
class EccOptions(ctypes.Structure):
    _fields_ = [("method", ctypes.c_int), ("use_parallel", ctypes.c_int)]


class PartitionOptions(ctypes.Structure):
    _fields_ = [
        ("n_parts", ctypes.c_int64),
        ("imbalance", ctypes.c_double),
        ("preset", ctypes.c_int),
        ("seed", ctypes.c_int),
        ("objective", ctypes.c_int),  # v3
        ("node_weighting", ctypes.c_int),  # v4
    ]


class RuntimeOptions(ctypes.Structure):
    _fields_ = [("num_threads", ctypes.c_int), ("suppress_output", ctypes.c_int)]


class Options(ctypes.Structure):
    _fields_ = [
        ("ecc", EccOptions),
        ("partition", PartitionOptions),
        ("runtime", RuntimeOptions),
    ]


class Statistics(ctypes.Structure):
    _fields_ = [
        # Matrix
        ("n_rows", ctypes.c_int64),
        ("n_cols", ctypes.c_int64),
        ("nnz", ctypes.c_int64),
        # Graph
        ("n_vertices", ctypes.c_int64),
        ("n_edges", ctypes.c_int64),
        ("degeneracy", ctypes.c_int64),  # v5: -1 non-BK; 0 is a real value
        # ECC
        ("ecc_method_used", ctypes.c_int),
        ("n_cliques", ctypes.c_int64),
        ("max_clique_size", ctypes.c_int64),
        ("avg_clique_size", ctypes.c_double),
        ("n_maximal_cliques", ctypes.c_int64),  # v5: 0 for non-BK
        # CNH
        ("n_hypernodes", ctypes.c_int64),
        ("n_hyperedges", ctypes.c_int64),
        ("total_pins", ctypes.c_int64),
        # Partition
        ("n_parts", ctypes.c_int64),
        ("separator_size", ctypes.c_int64),
        ("separator_ratio", ctypes.c_double),
        ("km1_objective", ctypes.c_int64),
        ("cut_objective", ctypes.c_int64),  # v3
        # Node weighting (v4)
        ("node_weighting_used", ctypes.c_int),
        ("total_node_weight", ctypes.c_int64),
        ("effective_imbalance", ctypes.c_double),
        ("achieved_block_imbalance", ctypes.c_double),
        # Timings (v5): 8 stages + total, ms
        ("time_ingestion_ms", ctypes.c_double),
        ("time_degeneracy_ms", ctypes.c_double),
        ("time_enumeration_ms", ctypes.c_double),
        ("time_greedy_ms", ctypes.c_double),
        ("time_cnh_ms", ctypes.c_double),
        ("time_partition_ms", ctypes.c_double),
        ("time_separator_ms", ctypes.c_double),
        ("time_permutation_ms", ctypes.c_double),
        ("time_total_ms", ctypes.c_double),
    ]


def find_library():
    here = os.path.dirname(os.path.abspath(__file__))
    # .../external/HG-DB/checks/bk_regression.py -> parallax root is 3 up
    parallax_root = os.path.abspath(os.path.join(here, "..", "..", ".."))
    cand = os.path.join(
        parallax_root, "external", "builds", "hgdb", "lib", "libhypergraph_reorder.so"
    )
    if not os.path.exists(cand):
        raise SystemExit(
            "library not found at %s (build it first: "
            "cmake --build build --target parallax_hgdb)" % cand
        )
    return cand


def load_library(path):
    lib = ctypes.CDLL(path, mode=1)  # RTLD_LAZY

    lib.hgr_version.argtypes = [ctypes.POINTER(ctypes.c_int)]

    lib.hgr_create.restype = ctypes.c_void_p
    lib.hgr_create.argtypes = [ctypes.POINTER(Options)]

    lib.hgr_free.argtypes = [ctypes.c_void_p]

    lib.hgr_reorder_csr.restype = ctypes.c_int
    lib.hgr_reorder_csr.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int64,
        ctypes.c_int64,
        ctypes.POINTER(ctypes.c_int64),
        ctypes.POINTER(ctypes.c_int64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_void_p),
    ]

    lib.hgr_get_permutation.restype = ctypes.c_int
    lib.hgr_get_permutation.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
        ctypes.POINTER(ctypes.c_int64),
    ]

    lib.hgr_get_statistics.restype = ctypes.c_int
    lib.hgr_get_statistics.argtypes = [ctypes.c_void_p, ctypes.POINTER(Statistics)]

    lib.hgr_get_clique_order_counts.restype = ctypes.c_int
    lib.hgr_get_clique_order_counts.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int64),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
    ]

    lib.hgr_get_maximal_clique_enumeration_histogram.restype = ctypes.c_int
    lib.hgr_get_maximal_clique_enumeration_histogram.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int64),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
    ]

    lib.hgr_result_free.argtypes = [ctypes.c_void_p]

    lib.hgr_get_last_error.restype = ctypes.c_char_p
    return lib


def edges_to_csr(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    for a in adj:
        a.sort()
    indptr = [0]
    indices = []
    for a in adj:
        indices.extend(a)
        indptr.append(len(indices))
    return indptr, indices


def run_case(lib, n, edges, method):
    """Run the full pipeline single-threaded (deterministic ECC ordering)."""
    indptr, indices = edges_to_csr(n, edges)

    opts = Options()
    lib.hgr_default_options(ctypes.byref(opts))
    opts.ecc.method = method
    opts.ecc.use_parallel = 0  # serial enumeration: deterministic
    opts.runtime.num_threads = 1
    opts.runtime.suppress_output = 1
    opts.partition.n_parts = 2

    reorderer = lib.hgr_create(ctypes.byref(opts))
    if not reorderer:
        raise RuntimeError("hgr_create failed: %s" % lib.hgr_get_last_error().decode())

    result = ctypes.c_void_p()
    rc = lib.hgr_reorder_csr(
        reorderer,
        ctypes.c_int64(n),
        ctypes.c_int64(len(indices)),
        (ctypes.c_int64 * len(indptr))(*indptr),
        (ctypes.c_int64 * len(indices))(*indices),
        None,
        ctypes.byref(result),
    )
    if rc != 0:
        lib.hgr_free(reorderer)
        raise RuntimeError(
            "hgr_reorder_csr failed (rc=%d): %s"
            % (rc, lib.hgr_get_last_error().decode())
        )

    stats = Statistics()
    lib.hgr_get_statistics(result, ctypes.byref(stats))

    hn = ctypes.c_int64()
    orders = ctypes.POINTER(ctypes.c_int64)()
    counts = ctypes.POINTER(ctypes.c_int64)()
    lib.hgr_get_clique_order_counts(
        result, ctypes.byref(hn), ctypes.byref(orders), ctypes.byref(counts)
    )
    hist = {int(orders[i]): int(counts[i]) for i in range(hn.value)} if hn.value else {}

    # v5: maximal-clique-enumeration histogram (BK only; *n == 0 otherwise)
    en = ctypes.c_int64()
    eorders = ctypes.POINTER(ctypes.c_int64)()
    ecounts = ctypes.POINTER(ctypes.c_int64)()
    lib.hgr_get_maximal_clique_enumeration_histogram(
        result, ctypes.byref(en), ctypes.byref(eorders), ctypes.byref(ecounts)
    )
    enum_hist = (
        {int(eorders[i]): int(ecounts[i]) for i in range(en.value)} if en.value else {}
    )

    perm_n = ctypes.c_int64()
    perm = ctypes.POINTER(ctypes.c_int64)()
    lib.hgr_get_permutation(result, ctypes.byref(perm), ctypes.byref(perm_n))
    perm_list = [int(perm[i]) for i in range(perm_n.value)]

    lib.hgr_result_free(result)
    lib.hgr_free(reorderer)

    return stats, hist, enum_hist, perm_list


def check(cond, msg, failures):
    if cond:
        print("  ok   %s" % msg)
    else:
        print("  FAIL %s" % msg)
        failures.append(msg)


def main():
    lib = load_library(find_library())

    version = (ctypes.c_int * 3)()
    lib.hgr_version(version)
    print("library: hgr_version %d.%d.%d" % tuple(version))
    if tuple(version) < (5, 0, 0):
        raise SystemExit(
            "this script's ctypes mirrors are the v5 layout; running it "
            "against an older .so would corrupt memory (hgr_statistics_t "
            "is caller-allocated). Rebuild first: "
            "cmake --build build --target parallax_hgdb"
        )
    print("")

    failures = []

    for fname, fx in FIXTURES.items():
        n = fx["n"]
        edges = fx["edges"]
        print(
            "fixture %s: n=%d, %d edges, degeneracy ordering %s (differs from vertex numbering: %s)"
            % (
                fname,
                n,
                len(edges),
                fx["degeneracy_ordering"],
                "yes" if fx["degeneracy_ordering"] != list(range(n)) else "NO",
            )
        )

        for method in (HGR_ECC_BK, HGR_ECC_TRI_ENUM, HGR_ECC_EDGE_COVER):
            mname = METHOD_NAMES[method]
            stats, hist, enum_hist, perm = run_case(lib, n, edges, method)
            exp = EXPECTED[fname][method]

            print(
                " %s: n_cliques=%d max=%d hist=%s"
                % (mname, stats.n_cliques, stats.max_clique_size, hist)
            )

            # --- structural invariants (any build must satisfy these) ---
            check(
                sorted(perm) == list(range(n)),
                "%s: permutation valid" % mname,
                failures,
            )
            check(
                stats.n_edges == EXPECTED[fname]["n_edges"],
                "%s: n_edges == %d" % (mname, EXPECTED[fname]["n_edges"]),
                failures,
            )
            check(
                sum(hist.values()) == stats.n_cliques,
                "%s: sum(hist) == n_cliques" % mname,
                failures,
            )
            check(
                max(hist) == stats.max_clique_size,
                "%s: max(hist orders) == max_clique_size" % mname,
                failures,
            )
            check(min(hist) >= 2, "%s: all clique orders >= 2" % mname, failures)
            avg_from_hist = sum(k * c for k, c in hist.items()) / stats.n_cliques
            check(
                abs(avg_from_hist - stats.avg_clique_size) < 1e-9,
                "%s: avg_clique_size consistent with histogram" % mname,
                failures,
            )
            if method == HGR_ECC_EDGE_COVER:
                check(
                    hist == {2: len(edges)},
                    "EdgeCover: exactly one 2-clique per edge",
                    failures,
                )

            # --- ground-truth pins (discriminate fixed vs pre-fix BK) ---
            check(
                stats.n_cliques == exp["n_cliques"],
                "%s: n_cliques == %d (got %d)"
                % (mname, exp["n_cliques"], stats.n_cliques),
                failures,
            )
            check(
                stats.max_clique_size == exp["max_clique_size"],
                "%s: max_clique_size == %d (got %d)"
                % (mname, exp["max_clique_size"], stats.max_clique_size),
                failures,
            )
            check(
                hist == exp["hist"],
                "%s: histogram == %s (got %s)" % (mname, exp["hist"], hist),
                failures,
            )

            if method == HGR_ECC_BK:
                # --- v5 observables (independent ground-truth pins) ---
                check(
                    stats.degeneracy == exp["degeneracy"],
                    "BK: degeneracy == %d (got %d)"
                    % (exp["degeneracy"], stats.degeneracy),
                    failures,
                )
                check(
                    stats.n_maximal_cliques == exp["n_maximal_cliques"],
                    "BK: n_maximal_cliques == %d (got %d)"
                    % (exp["n_maximal_cliques"], stats.n_maximal_cliques),
                    failures,
                )
                check(
                    enum_hist == exp["enum_hist"],
                    "BK: enum_hist == %s (got %s)" % (exp["enum_hist"], enum_hist),
                    failures,
                )
                # cross-checks (any build must satisfy these)
                check(
                    sum(enum_hist.values()) == stats.n_maximal_cliques,
                    "BK: sum(enum_hist) == n_maximal_cliques",
                    failures,
                )
                check(
                    all(k >= 2 for k in enum_hist),
                    "BK: all enumerated clique orders >= 2",
                    failures,
                )
                check(
                    max(enum_hist) == stats.max_clique_size,
                    "BK: max enumerated order == max_clique_size (greedy's "
                    "first pick is a maximum-order candidate)",
                    failures,
                )
                dominance = all(
                    hist.get(k, 0) <= enum_hist.get(k, 0) for k in hist if k >= 3
                )
                check(
                    dominance,
                    "BK: cover count <= enumeration count for every order >= 3",
                    failures,
                )
            else:
                # non-BK methods never compute the degeneracy: -1 sentinel
                # (never 0 — 0 is a real degeneracy), no enumeration.
                check(
                    stats.degeneracy == -1,
                    "%s: degeneracy == -1 (not-computed sentinel)" % mname,
                    failures,
                )
                check(
                    stats.n_maximal_cliques == 0 and enum_hist == {},
                    "%s: no enumeration observables" % mname,
                    failures,
                )
        print("")

    # Informational only - does not hold robustly on these fixtures.
    print(
        "informational: n_cliques(BK) vs n_cliques(TriEnum) relation varies by "
        "fixture (F1: 8 vs 8, F2: 10 vs 10, F3: 6 vs 10); not asserted."
    )
    print("")

    if failures:
        print("RESULT: FAIL (%d failed check(s))" % len(failures))
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
