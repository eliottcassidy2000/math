"""
gauss_speedup.py — Explore faster Gaussian elimination for tournament homology

Current bottleneck: _gauss_rank_np with PRIME=2^31-1 uses int64 numpy which is slow.

Ideas:
1. Use smaller prime (e.g., 2^23-1 or 2^16-1) — products fit in float64
2. Use float64 matmul with multi-precision correction
3. Batch row operations more aggressively
4. Use scipy sparse operations

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_allowed_paths, build_adj_lists,
    boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    RANK_PRIME
)


def _gauss_rank_np_small_prime(M_np, prime):
    """Gaussian elimination using a smaller prime.
    If prime < 2^26, products fit in float64 mantissa (52 bits).
    We can use float64 matmul for the elimination step."""
    nrows, ncols = M_np.shape
    mat = M_np.astype(np.float64) % prime
    rank = 0
    for col in range(ncols):
        nonzero = np.where(mat[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            mat[[rank, pivot]] = mat[[pivot, rank]]
        inv = pow(int(mat[rank, col]), prime - 2, prime)
        mat[rank] = mat[rank] * inv % prime
        factors = mat[:, col].copy()
        factors[rank] = 0
        nonzero_rows = np.where(factors != 0)[0]
        if len(nonzero_rows) > 0:
            mat[nonzero_rows] = (mat[nonzero_rows] - np.outer(factors[nonzero_rows], mat[rank])) % prime
        rank += 1
    return rank


def _gauss_nullbasis_small_prime(M, nrows, ncols, prime):
    """Gaussian elimination with null basis using smaller prime."""
    if isinstance(M, list):
        mat = np.array(M, dtype=np.float64) % prime
    else:
        mat = np.asarray(M, dtype=np.float64) % prime

    rank = 0
    pivot_cols = []
    for col in range(ncols):
        nonzero = np.where(mat[rank:nrows, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        pivot_cols.append(col)
        if pivot != rank:
            mat[[rank, pivot]] = mat[[pivot, rank]]
        inv = pow(int(mat[rank, col]), prime - 2, prime)
        mat[rank] = mat[rank] * inv % prime
        factors = mat[:nrows, col].copy()
        factors[rank] = 0
        nonzero_rows = np.where(factors != 0)[0]
        if len(nonzero_rows) > 0:
            mat[nonzero_rows] = (mat[nonzero_rows] - np.outer(factors[nonzero_rows], mat[rank])) % prime
        rank += 1

    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    null_basis = []
    for fc in free_cols:
        vec = [0] * ncols
        vec[fc] = 1
        for r, pc in enumerate(pivot_cols):
            vec[pc] = int((prime - mat[r, fc]) % prime)
        null_basis.append(vec)
    return rank, null_basis


SMALL_PRIME = 8388593  # Largest prime < 2^23 = 8388608. Products < 2^46 < 2^52.
# Actually let's use a prime where p^2 < 2^52 => p < 2^26 = 67108864
# Use 67108859 (largest prime < 2^26)
MEDIUM_PRIME = 67108859


def fast_beta3_small_prime(A, n, prime=MEDIUM_PRIME):
    """Compute beta_3 using smaller prime + float64 operations."""
    adj = build_adj_lists(A, n)
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p, adj)
    paths_sets = {p: set(ap.get(p, [])) for p in range(6)}

    def build_constraint(paths_p, paths_pm1_set):
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

    def build_boundary(paths_p, paths_pm1):
        if not paths_p or not paths_pm1:
            return None
        idx_prev = {path: i for i, path in enumerate(paths_pm1)}
        bd = np.zeros((len(paths_pm1), len(paths_p)), dtype=np.int64)
        for j, path in enumerate(paths_p):
            for sign, face in boundary_faces(path):
                if face in idx_prev:
                    bd[idx_prev[face], j] = (bd[idx_prev[face], j] + sign) % prime
        return bd

    # Omega_3
    P3, na3, nc3 = build_constraint(ap[3], paths_sets[2])
    if P3 is None:
        dim_omega3 = len(ap.get(3, []))
        omega3_basis = None
    else:
        rank_P3, nbasis3 = _gauss_nullbasis_small_prime(P3, na3, nc3, prime)
        dim_omega3 = nc3 - rank_P3
        omega3_basis = np.array(nbasis3, dtype=np.float64) if nbasis3 else None

    if dim_omega3 == 0:
        return 0

    bd3 = build_boundary(ap[3], ap[2])
    if bd3 is not None and omega3_basis is not None:
        composed3 = (bd3.astype(np.float64) @ omega3_basis.T) % prime
        rank_d3 = _gauss_rank_np_small_prime(composed3.astype(np.int64), prime)
    elif bd3 is not None:
        rank_d3 = _gauss_rank_np_small_prime(bd3.copy(), prime)
    else:
        rank_d3 = 0

    ker_d3 = dim_omega3 - rank_d3
    if ker_d3 == 0:
        return 0

    # Omega_4
    if not ap.get(4, []):
        return ker_d3

    P4, na4, nc4 = build_constraint(ap.get(4, []), paths_sets[3])
    if P4 is None:
        dim_omega4 = len(ap.get(4, []))
        omega4_basis = None
    else:
        rank_P4, nbasis4 = _gauss_nullbasis_small_prime(P4, na4, nc4, prime)
        dim_omega4 = nc4 - rank_P4
        omega4_basis = np.array(nbasis4, dtype=np.float64) if nbasis4 else None

    if dim_omega4 == 0:
        return ker_d3

    bd4 = build_boundary(ap[4], ap[3])
    if bd4 is not None and omega4_basis is not None:
        composed4 = (bd4.astype(np.float64) @ omega4_basis.T) % prime
        rank_d4 = _gauss_rank_np_small_prime(composed4.astype(np.int64), prime)
    elif bd4 is not None:
        rank_d4 = _gauss_rank_np_small_prime(bd4.copy(), prime)
    else:
        rank_d4 = 0

    return ker_d3 - rank_d4


def benchmark():
    """Compare different prime sizes."""
    print("=" * 70)
    print("GAUSS SPEEDUP BENCHMARK")
    print("=" * 70)

    from fast_beta3 import fast_beta3_nullbasis

    for n in [8, 9]:
        print(f"\n--- n={n} ---")
        rng = np.random.RandomState(42)
        n_trials = 30 if n <= 8 else 15

        total_big = 0
        total_small = 0
        mismatches = 0

        for trial in range(n_trials):
            A = random_tournament(n, rng)

            t0 = time.perf_counter()
            b3_big = fast_beta3_nullbasis(A, n, RANK_PRIME)
            t1 = time.perf_counter()
            total_big += t1 - t0

            t2 = time.perf_counter()
            b3_small = fast_beta3_small_prime(A, n, MEDIUM_PRIME)
            t3 = time.perf_counter()
            total_small += t3 - t2

            if b3_big != b3_small:
                mismatches += 1
                print(f"  MISMATCH at trial {trial}: big={b3_big}, small={b3_small}")

        avg_big = total_big / n_trials * 1000
        avg_small = total_small / n_trials * 1000
        speedup = avg_big / avg_small if avg_small > 0 else 0

        print(f"  Big prime (2^31-1):    {avg_big:.1f}ms avg")
        print(f"  Medium prime (2^26):   {avg_small:.1f}ms avg")
        print(f"  Speedup:               {speedup:.2f}x")
        print(f"  Mismatches:            {mismatches}/{n_trials}")

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    benchmark()
