"""
profile_bottleneck.py — Identify performance bottlenecks in tournament_utils

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')

from tournament_utils import (
    random_tournament, enumerate_all_allowed, full_chain_complex_modp,
    _build_constraint_matrix, _gauss_nullbasis_modp, _gauss_rank_np,
    boundary_faces, RANK_PRIME
)

def profile_n8():
    n = 8
    rng = np.random.RandomState(42)
    A = random_tournament(n, rng)
    prime = RANK_PRIME

    # Time each phase
    t0 = time.time()
    for _ in range(10):
        ap = enumerate_all_allowed(A, n, max_p=7)
    t1 = time.time()
    path_time = (t1 - t0) / 10

    # Count paths
    for p in range(8):
        paths = ap.get(p, [])
        print(f"  A_{p}: {len(paths)} paths")

    print(f"\n  Path enumeration: {path_time*1000:.1f}ms")

    # Time constraint building
    t2 = time.time()
    for _ in range(10):
        for p in range(1, 8):
            _build_constraint_matrix(ap, p, prime)
    t3 = time.time()
    constraint_time = (t3 - t2) / 10
    print(f"  Constraint building: {constraint_time*1000:.1f}ms")

    # Time null basis computation
    t4 = time.time()
    for _ in range(10):
        for p in range(1, 8):
            P, nrows, ncols = _build_constraint_matrix(ap, p, prime)
            if P is not None:
                _gauss_nullbasis_modp(P, nrows, ncols, prime)
    t5 = time.time()
    null_time = (t5 - t4) / 10
    print(f"  Null basis: {null_time*1000:.1f}ms")

    # Time boundary rank computation
    t6 = time.time()
    for _ in range(10):
        full_chain_complex_modp(A, n, max_p=7)
    t7 = time.time()
    full_time = (t7 - t6) / 10
    print(f"  Full chain complex: {full_time*1000:.1f}ms")

    # Breakdown by p
    print(f"\n  Per-degree breakdown:")
    for p in range(1, 7):
        paths_p = ap.get(p, [])
        paths_pm1 = ap.get(p-1, [])
        if not paths_p:
            continue

        t_a = time.time()
        for _ in range(10):
            P, nrows, ncols = _build_constraint_matrix(ap, p, prime)
        t_b = time.time()

        omega_dim = 0
        if P is not None:
            rk, nb = _gauss_nullbasis_modp(P, nrows, ncols, prime)
            omega_dim = ncols - rk
        else:
            omega_dim = len(paths_p)

        print(f"    p={p}: |A_p|={len(paths_p)}, dim(Omega_p)={omega_dim}, "
              f"constraint={((t_b-t_a)/10)*1000:.1f}ms")


def profile_path_enum():
    """Profile just path enumeration at different n."""
    print("\n  Path enumeration scaling:")
    for n in [5, 6, 7, 8, 9]:
        rng = np.random.RandomState(42)
        A = random_tournament(n, rng)
        t0 = time.time()
        reps = max(1, min(100, int(1000 / (n**2))))
        for _ in range(reps):
            ap = enumerate_all_allowed(A, n, max_p=n-1)
        t1 = time.time()
        total_paths = sum(len(ap.get(p, [])) for p in range(n))
        print(f"    n={n}: {((t1-t0)/reps)*1000:.1f}ms, {total_paths} total paths")


if __name__ == '__main__':
    print("=" * 70)
    print("PERFORMANCE PROFILING")
    print("=" * 70)
    print("\n--- n=8 breakdown ---")
    profile_n8()
    print("\n--- Path enumeration scaling ---")
    profile_path_enum()
    print("\nDONE.")
