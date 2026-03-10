"""
fast_beta3.py — Optimized beta_3 computation for n=9+

Key optimizations over full_chain_complex_modp:
1. Only compute degrees 2-5 (minimal for beta_3)
2. Use rank-only computation (no null basis extraction when possible)
3. Sparse-aware Gaussian elimination
4. Early termination when beta_3 must be 0

The formula: beta_3 = ker(d_3|Omega_3) - im(d_4|Omega_4)
           = (dim Omega_3 - rank(d_3|Omega_3)) - rank(d_4|Omega_4)

We need:
  - A_2, A_3, A_4, A_5 (allowed paths at degrees 2-5)
  - Omega_3 dimension and basis (for d_3 restriction)
  - rank(d_3|Omega_3) = rank of boundary d_3 restricted to Omega_3
  - Omega_4 dimension and basis (for d_4 restriction)
  - rank(d_4|Omega_4) = rank of boundary d_4 restricted to Omega_4

Alternative approach using block matrix rank:
  rank(d_p|Omega_p) = rank([P_p; d_p]) - rank(P_p)
  where P_p is the constraint matrix for Omega_p.
  This avoids computing the null basis entirely!

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp,
    enumerate_allowed_paths, build_adj_lists,
    boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME


def _build_constraint_rows(paths_p, paths_pm1_set, prime):
    """Build constraint matrix P for Omega_p.
    Returns (P_np, na_count, num_paths) or (None, 0, num_paths) if no constraints."""
    if not paths_p:
        return None, 0, 0

    non_allowed = {}
    na_count = 0
    entries = []
    for j, path in enumerate(paths_p):
        for sign, face in boundary_faces(path):
            if len(set(face)) == len(face) and face not in paths_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
                entries.append((non_allowed[face], j, sign))

    if na_count == 0:
        return None, 0, len(paths_p)

    P = np.zeros((na_count, len(paths_p)), dtype=np.int64)
    for row, col, sign in entries:
        P[row, col] = (P[row, col] + sign) % prime
    return P, na_count, len(paths_p)


def _build_boundary_np(paths_p, paths_pm1, prime):
    """Build boundary matrix d_p: A_p -> A_{p-1}."""
    if not paths_p or not paths_pm1:
        return None
    idx_prev = {path: i for i, path in enumerate(paths_pm1)}
    bd = np.zeros((len(paths_pm1), len(paths_p)), dtype=np.int64)
    for j, path in enumerate(paths_p):
        for sign, face in boundary_faces(path):
            if face in idx_prev:
                bd[idx_prev[face], j] = (bd[idx_prev[face], j] + sign) % prime
    return bd


def fast_beta3_block(A, n, prime=PRIME):
    """Compute beta_3 using block-matrix rank method.

    Key insight: rank(d_p|Omega_p) = rank([P_p; bd_p]) - rank(P_p)
    where P_p is the constraint matrix for Omega_p.

    This avoids computing the null basis entirely!

    beta_3 = (dim Omega_3 - rank(d_3|Omega_3)) - rank(d_4|Omega_4)
    """
    adj = build_adj_lists(A, n)

    # Enumerate only the paths we need: p=1,2,3,4,5
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p, adj)

    paths_sets = {p: set(ap.get(p, [])) for p in range(6)}

    # --- Omega_3 computation ---
    P3, na3, nc3 = _build_constraint_rows(ap[3], paths_sets[2], prime)

    if P3 is None:
        # No constraints: Omega_3 = all A_3 paths
        dim_omega3 = nc3 if nc3 else len(ap.get(3, []))
        # rank(d_3) unrestricted
        bd3 = _build_boundary_np(ap[3], ap[2], prime)
        if bd3 is not None:
            rank_d3 = _gauss_rank_np(bd3.copy(), prime)
        else:
            rank_d3 = 0
    else:
        # Omega_3 dim = nc3 - rank(P3)
        rank_P3 = _gauss_rank_np(P3.copy(), prime)
        dim_omega3 = nc3 - rank_P3

        if dim_omega3 == 0:
            return 0  # No 3-cycles => beta_3 = 0

        # rank(d_3|Omega_3) = rank([P3; bd3]) - rank(P3)
        bd3 = _build_boundary_np(ap[3], ap[2], prime)
        if bd3 is not None:
            # Stack P3 on top of bd3
            block = np.vstack([P3, bd3]) % prime
            rank_block = _gauss_rank_np(block, prime)
            rank_d3 = rank_block - rank_P3
        else:
            rank_d3 = 0

    ker_d3 = dim_omega3 - rank_d3
    if ker_d3 == 0:
        return 0  # No 3-cycles in kernel => beta_3 = 0

    # --- Omega_4 computation ---
    P4, na4, nc4 = _build_constraint_rows(ap.get(4, []), paths_sets[3], prime)

    if not ap.get(4, []):
        return ker_d3  # No 4-paths => rank(d_4) = 0

    if P4 is None:
        dim_omega4 = nc4 if nc4 else len(ap.get(4, []))
        bd4 = _build_boundary_np(ap[4], ap[3], prime)
        if bd4 is not None:
            rank_d4 = _gauss_rank_np(bd4.copy(), prime)
        else:
            rank_d4 = 0
    else:
        rank_P4 = _gauss_rank_np(P4.copy(), prime)
        dim_omega4 = nc4 - rank_P4

        if dim_omega4 == 0:
            return ker_d3  # rank(d_4|Omega_4) = 0

        bd4 = _build_boundary_np(ap[4], ap[3], prime)
        if bd4 is not None:
            block4 = np.vstack([P4, bd4]) % prime
            rank_block4 = _gauss_rank_np(block4, prime)
            rank_d4 = rank_block4 - rank_P4
        else:
            rank_d4 = 0

    return ker_d3 - rank_d4


def fast_beta3_nullbasis(A, n, prime=PRIME):
    """Compute beta_3 using null basis approach (current method, for comparison)."""
    adj = build_adj_lists(A, n)
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p, adj)
    paths_sets = {p: set(ap.get(p, [])) for p in range(6)}

    # Omega_3
    P3, na3, nc3 = _build_constraint_rows(ap[3], paths_sets[2], prime)
    if P3 is None:
        dim_omega3 = len(ap.get(3, []))
        omega3_basis = None
    else:
        rank_P3, nbasis3 = _gauss_nullbasis_modp(P3, na3, nc3, prime)
        dim_omega3 = nc3 - rank_P3
        omega3_basis = np.array(nbasis3, dtype=np.int64) if nbasis3 else None

    if dim_omega3 == 0:
        return 0

    bd3 = _build_boundary_np(ap[3], ap[2], prime)
    if bd3 is not None and omega3_basis is not None:
        composed3 = bd3 @ omega3_basis.T % prime
        rank_d3 = _gauss_rank_np(composed3, prime)
    elif bd3 is not None:
        rank_d3 = _gauss_rank_np(bd3.copy(), prime)
    else:
        rank_d3 = 0

    ker_d3 = dim_omega3 - rank_d3
    if ker_d3 == 0:
        return 0

    # Omega_4
    P4, na4, nc4 = _build_constraint_rows(ap.get(4, []), paths_sets[3], prime)
    if not ap.get(4, []):
        return ker_d3

    if P4 is None:
        dim_omega4 = len(ap.get(4, []))
        omega4_basis = None
    else:
        rank_P4, nbasis4 = _gauss_nullbasis_modp(P4, na4, nc4, prime)
        dim_omega4 = nc4 - rank_P4
        omega4_basis = np.array(nbasis4, dtype=np.int64) if nbasis4 else None

    if dim_omega4 == 0:
        return ker_d3

    bd4 = _build_boundary_np(ap[4], ap[3], prime)
    if bd4 is not None and omega4_basis is not None:
        composed4 = bd4 @ omega4_basis.T % prime
        rank_d4 = _gauss_rank_np(composed4, prime)
    elif bd4 is not None:
        rank_d4 = _gauss_rank_np(bd4.copy(), prime)
    else:
        rank_d4 = 0

    return ker_d3 - rank_d4


def benchmark():
    """Compare block method vs null basis method."""
    print("=" * 70)
    print("FAST BETA_3 BENCHMARK")
    print("=" * 70)

    for n in [7, 8, 9]:
        print(f"\n--- n={n} ---")
        rng = np.random.RandomState(42)
        n_trials = 50 if n <= 8 else 20

        total_block = 0
        total_null = 0
        mismatches = 0

        for trial in range(n_trials):
            A = random_tournament(n, rng)

            t0 = time.perf_counter()
            b3_block = fast_beta3_block(A, n)
            t1 = time.perf_counter()
            total_block += t1 - t0

            t2 = time.perf_counter()
            b3_null = fast_beta3_nullbasis(A, n)
            t3 = time.perf_counter()
            total_null += t3 - t2

            if b3_block != b3_null:
                mismatches += 1
                print(f"  MISMATCH at trial {trial}: block={b3_block}, null={b3_null}")

        avg_block = total_block / n_trials * 1000
        avg_null = total_null / n_trials * 1000
        speedup = avg_null / avg_block if avg_block > 0 else 0

        print(f"  Block method:  {avg_block:.1f}ms avg")
        print(f"  Null method:   {avg_null:.1f}ms avg")
        print(f"  Speedup:       {speedup:.2f}x")
        print(f"  Mismatches:    {mismatches}/{n_trials}")

    # Also compare with full_chain_complex_modp
    print(f"\n--- Comparison with full_chain_complex_modp ---")
    for n in [8, 9]:
        rng = np.random.RandomState(42)
        n_trials = 20 if n <= 8 else 10

        total_fast = 0
        total_full = 0
        mismatches = 0

        for trial in range(n_trials):
            A = random_tournament(n, rng)

            t0 = time.perf_counter()
            b3_fast = fast_beta3_block(A, n)
            t1 = time.perf_counter()
            total_fast += t1 - t0

            t2 = time.perf_counter()
            res = full_chain_complex_modp(A, n, max_p=5)
            b3_full = res['bettis'].get(3, 0)
            t3 = time.perf_counter()
            total_full += t3 - t2

            if b3_fast != b3_full:
                mismatches += 1
                print(f"  MISMATCH at trial {trial}: fast={b3_fast}, full={b3_full}")

        avg_fast = total_fast / n_trials * 1000
        avg_full = total_full / n_trials * 1000

        print(f"  n={n}: fast_beta3={avg_fast:.1f}ms, full_chain={avg_full:.1f}ms, "
              f"speedup={avg_full/avg_fast:.2f}x, mismatches={mismatches}")

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    benchmark()
