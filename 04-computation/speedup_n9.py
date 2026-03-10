"""
speedup_n9.py — Profile and optimize n=9 path homology computation

The n=8 speedup (S48) went from 95ms to 43.6ms (2.2x) via numpy vectorized Gauss.
At n=9, the computation is ~250ms per tournament. Target: <100ms.

Profiling breakdown needed:
1. Path enumeration (enumerate_all_allowed)
2. Constraint matrix construction (_build_constraint_matrix)
3. Gaussian elimination (_gauss_nullbasis_modp)
4. Boundary matrix construction + rank computation

Speedup ideas:
A. Bitwise path enumeration (replace list-of-lists with int bitmasks)
B. Sparse constraint matrices (most entries are 0)
C. Early termination: if dim(Omega_3) = 0, skip boundary computation
D. Cache adjacency lists
E. Parallel path enumeration for different starting vertices

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp,
    enumerate_all_allowed, enumerate_allowed_paths, build_adj_lists,
    _build_constraint_matrix, _build_boundary_matrix,
    _gauss_nullbasis_modp, _gauss_rank_np,
    boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def profile_n9_detailed(A, n):
    """Profile each step of chain complex computation at n=9."""
    times = {}

    # Step 1: Path enumeration
    t0 = time.perf_counter()
    adj = build_adj_lists(A, n)
    t1 = time.perf_counter()
    times['adj_lists'] = t1 - t0

    ap = {}
    for p in range(min(8, n)):  # max_p=7 for full Betti
        t2 = time.perf_counter()
        ap[p] = enumerate_allowed_paths(A, n, p, adj)
        t3 = time.perf_counter()
        times[f'enum_p{p}'] = t3 - t2

    t4 = time.perf_counter()
    times['total_enum'] = t4 - t0

    # Step 2: For each degree p, constraint matrix + Gauss
    for p in [2, 3, 4, 5]:
        paths = ap.get(p, [])
        if not paths:
            times[f'constraint_p{p}'] = 0
            times[f'gauss_p{p}'] = 0
            continue

        t5 = time.perf_counter()
        P, na_rows, na_cols = _build_constraint_matrix(ap, p, PRIME)
        t6 = time.perf_counter()
        times[f'constraint_p{p}'] = t6 - t5

        if P is not None:
            t7 = time.perf_counter()
            rk, nbasis = _gauss_nullbasis_modp(P, na_rows, na_cols, PRIME)
            t8 = time.perf_counter()
            times[f'gauss_p{p}'] = t8 - t7
            times[f'omega_dim_p{p}'] = na_cols - rk
        else:
            times[f'gauss_p{p}'] = 0
            times[f'omega_dim_p{p}'] = len(paths)

    # Step 3: Boundary matrices + rank
    for p in [3, 4, 5]:
        t9 = time.perf_counter()
        bd, bd_nrows, bd_ncols = _build_boundary_matrix(ap, p, PRIME)
        t10 = time.perf_counter()
        times[f'boundary_p{p}'] = t10 - t9

        if bd is not None:
            # Need omega basis for composition
            P, na_rows, na_cols = _build_constraint_matrix(ap, p, PRIME)
            if P is not None:
                rk, nbasis = _gauss_nullbasis_modp(P, na_rows, na_cols, PRIME)
                if nbasis:
                    basis_np = np.array(nbasis, dtype=np.int64)
                    t11 = time.perf_counter()
                    composed = bd @ basis_np.T % PRIME
                    t12 = time.perf_counter()
                    times[f'compose_p{p}'] = t12 - t11
                    t13 = time.perf_counter()
                    rank_bd = _gauss_rank_np(composed, PRIME)
                    t14 = time.perf_counter()
                    times[f'rank_bd_p{p}'] = t14 - t13

    # Path count summary
    for p in range(8):
        if p in ap:
            times[f'|A_{p}|'] = len(ap[p])

    return times


def fast_enumerate_paths_bitwise(A, n, p, adj=None):
    """Faster path enumeration using integer bitmask for visited set.
    Preallocate arrays instead of list append."""
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    if adj is None:
        adj = build_adj_lists(A, n)

    # Convert adj to tuple-of-tuples for speed
    adj_t = [tuple(adj[i]) for i in range(n)]

    paths = []
    # Iterative DFS with tuple-based stack
    for start in range(n):
        stack = [(start, 1 << start, 0)]  # (current_vertex, visited_mask, depth)
        path_stack = [start]  # Current path being built

        # Use explicit stack with path tracking
        stack2 = [((start,), 1 << start)]
        while stack2:
            path, visited = stack2.pop()
            if len(path) == p + 1:
                paths.append(path)
                continue
            v = path[-1]
            for u in adj_t[v]:
                if not (visited & (1 << u)):
                    stack2.append((path + (u,), visited | (1 << u)))

    return paths


def benchmark_path_enumeration(n, n_trials=20):
    """Benchmark path enumeration variants."""
    print(f"\n--- Path Enumeration Benchmark at n={n} ---")
    rng = np.random.RandomState(12345)

    total_orig = 0
    total_fast = 0

    for trial in range(n_trials):
        A = random_tournament(n, rng)
        adj = build_adj_lists(A, n)

        # Original
        t0 = time.perf_counter()
        for p in range(min(6, n)):
            enumerate_allowed_paths(A, n, p, adj)
        t1 = time.perf_counter()
        total_orig += t1 - t0

        # Fast bitwise
        t2 = time.perf_counter()
        for p in range(min(6, n)):
            fast_enumerate_paths_bitwise(A, n, p, adj)
        t3 = time.perf_counter()
        total_fast += t3 - t2

    avg_orig = total_orig / n_trials * 1000
    avg_fast = total_fast / n_trials * 1000
    print(f"  Original: {avg_orig:.1f}ms avg")
    print(f"  Bitwise:  {avg_fast:.1f}ms avg")
    print(f"  Speedup:  {avg_orig/avg_fast:.2f}x")


def main():
    print("=" * 70)
    print("N=9 SPEEDUP ANALYSIS")
    print("=" * 70)

    # Profile n=8 for baseline
    n = 8
    print(f"\n--- Detailed profiling at n={n} ---")
    rng = np.random.RandomState(42)
    A = random_tournament(n, rng)

    times = profile_n9_detailed(A, n)
    print(f"\n  n={n} profile (single tournament):")
    for key in sorted(times.keys()):
        val = times[key]
        if isinstance(val, float):
            if val > 0.001:
                print(f"    {key}: {val*1000:.1f}ms")
            elif val > 0:
                print(f"    {key}: {val*1000:.3f}ms")
        else:
            print(f"    {key}: {val}")

    # Profile n=9
    n = 9
    print(f"\n--- Detailed profiling at n={n} ---")
    rng = np.random.RandomState(42)
    A = random_tournament(n, rng)

    times = profile_n9_detailed(A, n)
    print(f"\n  n={n} profile (single tournament):")
    for key in sorted(times.keys()):
        val = times[key]
        if isinstance(val, float):
            if val > 0.001:
                print(f"    {key}: {val*1000:.1f}ms")
            elif val > 0:
                print(f"    {key}: {val*1000:.3f}ms")
        else:
            print(f"    {key}: {val}")

    # Full chain complex timing comparison
    print(f"\n--- Full chain complex timing ---")
    for n in [8, 9]:
        rng = np.random.RandomState(42)
        A = random_tournament(n, rng)
        n_trials = 10 if n <= 8 else 5

        total = 0
        for _ in range(n_trials):
            A = random_tournament(n, rng)
            t0 = time.perf_counter()
            full_chain_complex_modp(A, n, max_p=min(n-1, 7))
            t1 = time.perf_counter()
            total += t1 - t0

        print(f"  n={n}: avg {total/n_trials*1000:.1f}ms per tournament")

    # Benchmark path enumeration variants
    benchmark_path_enumeration(8, 20)
    benchmark_path_enumeration(9, 10)

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
