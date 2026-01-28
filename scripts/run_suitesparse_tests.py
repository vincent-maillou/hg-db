#!/usr/bin/env python3
"""
Test script for running hypergraph reordering on SuiteSparse matrices.

This script:
1. Loads matrices from a local directory ($HPCVAULT/matrices or --matrix-dir)
2. Runs the hypergraph_reorder CLI on each matrix
3. Collects and reports timing statistics
"""

import sys
import subprocess
import time
import csv
import argparse
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


# Test matrices organized by group
MATRICES = {
    "CFD": [
        "copter2",
        "ex10",
        "parabolic_fem",
        "3Dspectralwave",
    ],
    "Chemistry": [
        "CO",
        "crystk03",
        "Ga19As19H42",
        "Si87H76",
        "SiH4",
    ],
    "Circuit": [
        "circuit_3",
        "G2_circuit",
        "memplus",
        "rajat06",
        "rajat10",
        "rajat15",
    ],
    "Graphs": [
        "144",
        "598a",
        "auto",
        "ca-AstroPh",
        "citationCiteseer",
        "delaunay_n20",
        "dictionary28",
    ],
    "Optimization": [
        "boyd2",
        "c-55",
        "c-70",
        "gupta2",
        "gupta3",
        "lpl1",
        "ncvxbqp1",
        "apache1",
        "apache2",
        "nasasrb",
    ],
}


@dataclass
class TestResult:
    """Result of a single matrix test."""
    matrix_name: str
    group: str
    n_rows: int = 0
    n_cols: int = 0
    nnz: int = 0
    n_parts: int = 0
    preset: str = ""
    time_total_ms: float = 0.0
    separator_size: int = 0
    separator_ratio: float = 0.0
    # Fill-in metrics (computed by cmpfillin)
    fillin_original: int = 0
    fillin_reordered: int = 0
    fillin_reduction: float = 0.0
    nnz_l_original: int = 0
    nnz_l_reordered: int = 0
    ops_original: int = 0
    ops_reordered: int = 0
    success: bool = False
    error_message: str = ""
    output_path: str = ""


def write_identity_permutation(perm_path: Path, n: int) -> None:
    """Write an identity permutation file (0, 1, 2, ..., n-1)."""
    with open(perm_path, 'w') as f:
        for i in range(n):
            f.write(f"{i}\n")


def run_cmpfillin(cmpfillin_path: Path, graph_path: Path, perm_path: Path, timeout: int = 300) -> dict:
    """
    Run the cmpfillin program with METIS graph and permutation file.

    cmpfillin expects: <GraphFile> <PermFile>
    - GraphFile: METIS format graph
    - PermFile: One integer per line (inverse permutation, 0-indexed)

    Returns dict with keys: nnz_l, ops (or empty dict on failure)
    """
    try:
        result = subprocess.run(
            [str(cmpfillin_path), str(graph_path), str(perm_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        if result.returncode != 0:
            return {}

        stats = {}
        for line in result.stdout.split('\n'):
            # Parse output format: "Nonzeros: X.XXXe+YY Operation Count: Z.ZZZe+WW"
            if "Nonzeros:" in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "Nonzeros:":
                        # Parse scientific notation (e.g., 1.234e+05)
                        try:
                            stats['nnz_l'] = int(float(parts[i + 1]))
                        except (IndexError, ValueError):
                            pass
                    elif p == "Count:":
                        try:
                            stats['ops'] = int(float(parts[i + 1]))
                        except (IndexError, ValueError):
                            pass

        return stats
    except Exception:
        return {}


def find_matrix(matrix_name: str, matrix_dir: Path) -> Optional[Path]:
    """Find a matrix .mtx file in the local matrix directory."""
    mtx_path = matrix_dir / f"{matrix_name}.mtx"
    if mtx_path.exists():
        return mtx_path
    print(f"  ERROR: {mtx_path} not found")
    return None


def run_reorder(
    cli_path: Path,
    matrix_path: Path,
    n_parts: int,
    preset: str = "default",
    threads: int = 0,
    quiet: bool = True,
    timeout: int = 600,
    output_graph: bool = False,
) -> tuple[bool, str, float]:
    """
    Run the hypergraph_reorder CLI on a matrix.

    Returns: (success, output, elapsed_time_ms)
    """
    cmd = [str(cli_path), str(matrix_path), str(n_parts), "--preset", preset]
    if threads > 0:
        cmd.extend(["--threads", str(threads)])
    if quiet:
        cmd.append("--quiet")
    if output_graph:
        cmd.append("--output-graph")

    env = dict()
    # Set LD_LIBRARY_PATH to include the build directory for shared libraries
    build_dir = cli_path.parent
    env["LD_LIBRARY_PATH"] = f"{build_dir}:{build_dir / 'external' / 'mt-kahypar' / 'lib'}"

    start_time = time.perf_counter()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env={**dict(__import__('os').environ), **env},
        )
        elapsed_ms = (time.perf_counter() - start_time) * 1000

        output = result.stdout + result.stderr
        success = result.returncode == 0
        return success, output, elapsed_ms

    except subprocess.TimeoutExpired:
        elapsed_ms = (time.perf_counter() - start_time) * 1000
        return False, f"Timeout after {timeout}s", elapsed_ms
    except Exception as e:
        elapsed_ms = (time.perf_counter() - start_time) * 1000
        return False, str(e), elapsed_ms


def parse_statistics(output: str) -> dict:
    """Parse statistics from CLI output."""
    stats = {}

    lines = output.split('\n')
    for line in lines:
        # Parse matrix dimensions
        if line.startswith("Matrix:"):
            # Matrix: 144649 x 144649, 1074393 nonzeros
            parts = line.split()
            if len(parts) >= 4:
                stats['n_rows'] = int(parts[1])
                stats['n_cols'] = int(parts[3].rstrip(','))
                stats['nnz'] = int(parts[4])

        # Parse partition info
        elif "Partition:" in line:
            # Partition: 8 parts, separator size: 1234 (5.67%)
            parts = line.split()
            for i, p in enumerate(parts):
                if p == "parts,":
                    stats['n_parts'] = int(parts[i-1])
                elif p == "size:":
                    stats['separator_size'] = int(parts[i+1])
                elif p.startswith('(') and p.endswith('%)'):
                    stats['separator_ratio'] = float(p[1:-2]) / 100

        # Parse total time
        elif "Total time:" in line:
            # Total time: 123.45 ms
            parts = line.split()
            for i, p in enumerate(parts):
                if p == "ms":
                    stats['time_total_ms'] = float(parts[i-1])

    return stats


def run_tests(
    matrices: dict[str, list[str]],
    matrix_dir: Path,
    cli_path: Path,
    cmpfillin_path: Path,
    output_dir: Path,
    n_parts: int = 8,
    preset: str = "default",
    threads: int = 0,
    timeout: int = 600,
) -> list[TestResult]:
    """Run tests on all specified matrices."""
    results = []

    total_matrices = sum(len(m) for m in matrices.values())
    current = 0

    # Cache for original matrix fill-in (computed once per matrix)
    # Keys: matrix_name -> {'nnz_l': int, 'ops': int, 'graph_path': Path, 'n_vertices': int}
    original_fillin_cache: dict[str, dict] = {}

    for group, matrix_list in matrices.items():
        print(f"\n{'='*60}")
        print(f"Group: {group}")
        print(f"{'='*60}")

        for matrix_name in matrix_list:
            current += 1
            print(f"\n[{current}/{total_matrices}] Testing {matrix_name}...")

            result = TestResult(
                matrix_name=matrix_name,
                group=group,
                n_parts=n_parts,
                preset=preset,
            )

            # Find matrix
            matrix_path = find_matrix(matrix_name, matrix_dir)
            if matrix_path is None:
                result.error_message = "Matrix file not found"
                results.append(result)
                continue

            # Run reordering with --output-graph to generate .graph and .perm files
            print(f"  Running reordering (k={n_parts}, preset={preset})...", end=" ", flush=True)
            success, output, elapsed_ms = run_reorder(
                cli_path=cli_path,
                matrix_path=matrix_path,
                n_parts=n_parts,
                preset=preset,
                threads=threads,
                timeout=timeout,
                output_graph=True,  # Generate .graph and .perm files
            )

            if not success:
                print(f"FAILED")
                result.error_message = output[:200] if output else "Unknown error"
                results.append(result)
                continue

            print(f"OK ({elapsed_ms:.1f} ms)")
            stats = parse_statistics(output)
            result.success = True
            result.n_rows = stats.get('n_rows', 0)
            result.n_cols = stats.get('n_cols', 0)
            result.nnz = stats.get('nnz', 0)
            result.time_total_ms = stats.get('time_total_ms', elapsed_ms)
            result.separator_size = stats.get('separator_size', 0)
            result.separator_ratio = stats.get('separator_ratio', 0.0)
            result.output_path = str(matrix_path.with_suffix('')) + "_reordered.mtx"

            # Paths to generated files
            base_path = matrix_path.with_suffix('')
            graph_path = base_path.with_suffix('.graph')
            perm_path = base_path.with_suffix('.perm')

            # Check that graph and perm files were generated
            if not graph_path.exists() or not perm_path.exists():
                print(f"  WARNING: Graph or perm file not generated")
                results.append(result)
                continue

            # Get number of vertices from stats
            n_vertices = result.n_rows

            # Compute fill-in for original matrix (with identity permutation)
            if matrix_name not in original_fillin_cache:
                print(f"  Computing fill-in (original)...", end=" ", flush=True)
                # Create identity permutation file
                identity_perm_path = base_path.with_suffix('.identity_perm')
                write_identity_permutation(identity_perm_path, n_vertices)

                orig_stats = run_cmpfillin(cmpfillin_path, graph_path, identity_perm_path)
                if orig_stats:
                    original_fillin_cache[matrix_name] = orig_stats
                    print(f"OK (nnz_l={orig_stats.get('nnz_l', 0)})")
                else:
                    original_fillin_cache[matrix_name] = {}
                    print("FAILED")

            orig_stats = original_fillin_cache.get(matrix_name, {})
            result.nnz_l_original = orig_stats.get('nnz_l', 0)
            result.ops_original = orig_stats.get('ops', 0)

            # Compute fill-in for reordered matrix (with generated permutation)
            print(f"  Computing fill-in (reordered)...", end=" ", flush=True)
            reord_stats = run_cmpfillin(cmpfillin_path, graph_path, perm_path)
            if reord_stats:
                result.nnz_l_reordered = reord_stats.get('nnz_l', 0)
                result.ops_reordered = reord_stats.get('ops', 0)

                # Calculate fill-in reduction based on nnz_l
                if result.nnz_l_original > 0:
                    result.fillin_reduction = (
                        (result.nnz_l_original - result.nnz_l_reordered)
                        / result.nnz_l_original
                    )
                print(f"OK (nnz_l={result.nnz_l_reordered}, reduction={result.fillin_reduction*100:.1f}%)")
            else:
                print("FAILED")

            results.append(result)

    return results


def print_summary(results: list[TestResult]):
    """Print a summary of test results."""
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]

    print(f"\nTotal: {len(results)} matrices")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if successful:
        print(f"\n{'Matrix':<20} {'Group':<10} {'k':<4} {'Time(ms)':<10} {'NNZ_L Orig':<12} {'NNZ_L Reord':<12} {'Reduction':<10}")
        print("-" * 90)
        for r in sorted(successful, key=lambda x: (x.matrix_name, x.n_parts)):
            reduction_str = f"{r.fillin_reduction*100:.1f}%" if r.nnz_l_original > 0 else "N/A"
            print(f"{r.matrix_name:<20} {r.group:<10} {r.n_parts:<4} {r.time_total_ms:<10.1f} {r.nnz_l_original:<12} {r.nnz_l_reordered:<12} {reduction_str:<10}")

    if failed:
        print(f"\nFailed matrices:")
        for r in failed:
            print(f"  - {r.matrix_name}: {r.error_message[:60]}")


def save_results_csv(results: list[TestResult], output_path: Path):
    """Save results to a CSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'matrix_name', 'group', 'n_rows', 'n_cols', 'nnz',
            'n_parts', 'preset', 'time_total_ms',
            'separator_size', 'separator_ratio',
            'fillin_original', 'fillin_reordered', 'fillin_reduction',
            'nnz_l_original', 'nnz_l_reordered',
            'ops_original', 'ops_reordered',
            'success', 'error_message'
        ])
        for r in results:
            writer.writerow([
                r.matrix_name, r.group, r.n_rows, r.n_cols, r.nnz,
                r.n_parts, r.preset, r.time_total_ms,
                r.separator_size, r.separator_ratio,
                r.fillin_original, r.fillin_reordered, r.fillin_reduction,
                r.nnz_l_original, r.nnz_l_reordered,
                r.ops_original, r.ops_reordered,
                r.success, r.error_message
            ])
    print(f"\nResults saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Test hypergraph reordering on SuiteSparse matrices"
    )
    parser.add_argument(
        "--matrices", "-m",
        nargs="+",
        help="Specific matrices to test (default: all)"
    )
    parser.add_argument(
        "--groups", "-g",
        nargs="+",
        choices=list(MATRICES.keys()),
        help="Groups of matrices to test (default: all)"
    )
    parser.add_argument(
        "--n-parts", "-k",
        type=int,
        nargs="+",
        default=[2, 4, 8, 16],
        help="Number of diagonal blocks/parts (default: 2 4 8 16)"
    )
    parser.add_argument(
        "--preset", "-p",
        choices=["default", "quality", "deterministic", "large_k"],
        default="default",
        help="MT-KaHyPar preset (default: default)"
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=0,
        help="Number of threads (default: auto-detect)"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Timeout per matrix in seconds (default: 600)"
    )
    parser.add_argument(
        "--output-csv", "-o",
        type=Path,
        help="Output CSV file for results"
    )
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        default=Path(os.environ.get("HPCVAULT", ""), "matrices"),
        help="Directory containing .mtx files (default: $HPCVAULT/matrices)"
    )

    args = parser.parse_args()

    # Setup paths
    project_root = Path(__file__).resolve().parent.parent
    cli_path = project_root / "build" / "hypergraph_reorder"
    cmpfillin_path = project_root.parent / "cmpfillin" / "cmpfillin"
    matrix_dir = args.matrix_dir
    output_dir = project_root / "results"

    # Verify matrix directory exists
    if not matrix_dir.is_dir():
        print(f"Error: Matrix directory not found at {matrix_dir}")
        print("Set $HPCVAULT or use --matrix-dir")
        return 1

    # Determine which matrices to test
    if args.matrices:
        # Specific matrices requested
        matrices_to_test = {"Custom": args.matrices}
    elif args.groups:
        # Specific groups requested
        matrices_to_test = {g: MATRICES[g] for g in args.groups}
    else:
        # All matrices
        matrices_to_test = MATRICES

    # Verify CLI exists
    if not cli_path.exists():
        print(f"Error: CLI not found at {cli_path}")
        print("Please build the project first: cd build && make")
        return 1

    # Verify cmpfillin exists
    if not cmpfillin_path.exists():
        print(f"Error: cmpfillin not found at {cmpfillin_path}")
        print("Please build cmpfillin first: cd ../cmpfillin && make")
        return 1

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Hypergraph Reordering Test Suite")
    print("=" * 60)
    print(f"CLI: {cli_path}")
    print(f"cmpfillin: {cmpfillin_path}")
    print(f"Matrix dir: {matrix_dir}")
    print(f"Parameters: k={args.n_parts}, preset={args.preset}, threads={args.threads or 'auto'}")

    # Run tests for each k value
    all_results = []
    for n_parts in args.n_parts:
        print(f"\n{'#'*60}")
        print(f"# Running with k = {n_parts}")
        print(f"{'#'*60}")

        results = run_tests(
            matrices=matrices_to_test,
            matrix_dir=matrix_dir,
            cli_path=cli_path,
            cmpfillin_path=cmpfillin_path,
            output_dir=output_dir,
            n_parts=n_parts,
            preset=args.preset,
            threads=args.threads,
            timeout=args.timeout,
        )
        all_results.extend(results)

        # Print summary for this k
        print_summary(results)

        # Save results for this k
        if not args.output_csv:
            k_csv = output_dir / f"results_k{n_parts}_{args.preset}.csv"
            save_results_csv(results, k_csv)

    # Save combined results if output file specified
    if args.output_csv:
        save_results_csv(all_results, args.output_csv)
    else:
        # Also save combined results
        combined_csv = output_dir / f"results_combined_{args.preset}.csv"
        save_results_csv(all_results, combined_csv)

    # Print overall summary
    print(f"\n{'='*80}")
    print("OVERALL SUMMARY")
    print(f"{'='*80}")
    for n_parts in args.n_parts:
        k_results = [r for r in all_results if r.n_parts == n_parts]
        successful = sum(1 for r in k_results if r.success)
        print(f"k={n_parts}: {successful}/{len(k_results)} successful")

    # Return exit code based on success rate
    successful = sum(1 for r in all_results if r.success)
    return 0 if successful == len(all_results) else 1


if __name__ == "__main__":
    sys.exit(main())
