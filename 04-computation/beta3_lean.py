"""
beta3_lean.py — Memory-efficient beta_3 computation + search

Uses int32 arrays with small prime (< 2^15) to save 50% memory.
Also uses sparse operations where possible.

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, enumerate_allowed_paths, build_adj_lists,
    boundary_faces, RANK_PRIME
)

# Small prime: products fit in int32 (p < 2^15 = 32768, p^2 < 2^30 < 2^31)
SMALL_PRIME = 32749  # Largest prime < 2^15


def _gauss_rank_i32(M_np, prime):
    """Gaussian elimination mod prime using int32 arrays."""
    nrows, ncols = M_np.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M_np[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M_np[[rank, pivot]] = M_np[[pivot, rank]]
        inv = pow(int(M_np[rank, col]), prime - 2, prime)
        M_np[rank] = M_np[rank] * inv % prime
        factors = M_np[:, col].copy()
        factors[rank] = 0
        nonzero_rows = np.where(factors != 0)[0]
        if len(nonzero_rows) > 0:
            # Cast to int64 for the multiply to avoid overflow, then back
            M_np[nonzero_rows] = ((M_np[nonzero_rows].astype(np.int64) -
                                   np.outer(factors[nonzero_rows].astype(np.int64),
                                           M_np[rank].astype(np.int64))) % prime).astype(np.int32)
        rank += 1
    return rank


def _gauss_nullbasis_i32(P_np, nrows, ncols, prime):
    """Gaussian elimination with null basis using int32."""
    mat = P_np.astype(np.int32) % prime
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
            mat[nonzero_rows] = ((mat[nonzero_rows].astype(np.int64) -
                                  np.outer(factors[nonzero_rows].astype(np.int64),
                                          mat[rank].astype(np.int64))) % prime).astype(np.int32)
        rank += 1

    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    null_basis = []
    for fc in free_cols:
        vec = np.zeros(ncols, dtype=np.int32)
        vec[fc] = 1
        for r, pc in enumerate(pivot_cols):
            vec[pc] = int((prime - mat[r, fc]) % prime)
        null_basis.append(vec)
    return rank, null_basis


def fast_beta3_lean(A, n, prime=SMALL_PRIME):
    """Memory-efficient beta_3 computation using int32."""
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
        P = np.zeros((na_count, len(paths_p)), dtype=np.int32)
        for row, col, sign in entries:
            P[row, col] = (P[row, col] + sign) % prime
        return P, na_count, len(paths_p)

    def build_boundary(paths_p, paths_pm1):
        if not paths_p or not paths_pm1:
            return None
        idx_prev = {path: i for i, path in enumerate(paths_pm1)}
        bd = np.zeros((len(paths_pm1), len(paths_p)), dtype=np.int32)
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
        rank_P3, nbasis3 = _gauss_nullbasis_i32(P3, na3, nc3, prime)
        dim_omega3 = nc3 - rank_P3
        omega3_basis = np.array(nbasis3, dtype=np.int32) if nbasis3 else None
        del P3  # Free memory

    if dim_omega3 == 0:
        return 0

    bd3 = build_boundary(ap[3], ap[2])
    if bd3 is not None and omega3_basis is not None:
        composed3 = (bd3.astype(np.int64) @ omega3_basis.astype(np.int64).T % prime).astype(np.int32)
        rank_d3 = _gauss_rank_i32(composed3, prime)
        del composed3
    elif bd3 is not None:
        rank_d3 = _gauss_rank_i32(bd3.copy(), prime)
    else:
        rank_d3 = 0
    del bd3, omega3_basis  # Free memory

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
        rank_P4, nbasis4 = _gauss_nullbasis_i32(P4, na4, nc4, prime)
        dim_omega4 = nc4 - rank_P4
        omega4_basis = np.array(nbasis4, dtype=np.int32) if nbasis4 else None
        del P4  # Free memory

    if dim_omega4 == 0:
        return ker_d3

    bd4 = build_boundary(ap[4], ap[3])
    if bd4 is not None and omega4_basis is not None:
        composed4 = (bd4.astype(np.int64) @ omega4_basis.astype(np.int64).T % prime).astype(np.int32)
        rank_d4 = _gauss_rank_i32(composed4, prime)
        del composed4
    elif bd4 is not None:
        rank_d4 = _gauss_rank_i32(bd4.copy(), prime)
    else:
        rank_d4 = 0
    del bd4, omega4_basis

    return ker_d3 - rank_d4


def main():
    """Search for beta_3 values and study structural properties."""
    import gc

    print("=" * 70)
    print("LEAN BETA_3 SEARCH + LES ANALYSIS")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(12345)

    b3_2_tours = []
    b3_1_tours = []
    b3_dist = {}

    t0 = time.time()
    for trial in range(5000):
        A = random_tournament(n, rng)
        gc.collect()  # Explicit GC to manage memory
        b3 = fast_beta3_lean(A, n)

        b3_dist[b3] = b3_dist.get(b3, 0) + 1

        if b3 == 2 and len(b3_2_tours) < 4:
            b3_2_tours.append((trial, A.copy()))
            scores = sorted([int(sum(A[i])) for i in range(n)])
            print(f"  beta_3=2 at trial {trial}, scores={scores}")
        elif b3 == 1 and len(b3_1_tours) < 2:
            b3_1_tours.append((trial, A.copy()))

        if (trial + 1) % 500 == 0:
            elapsed = time.time() - t0
            gc.collect()
            print(f"  {trial+1}/5000, {elapsed:.1f}s, b3 dist: {dict(sorted(b3_dist.items()))}")

    print(f"\nb3 distribution: {dict(sorted(b3_dist.items()))}")

    # Now study the LES for b3=2 tournaments
    if b3_2_tours:
        print(f"\n{'='*60}")
        print(f"LES ANALYSIS OF {len(b3_2_tours)} BETA_3=2 TOURNAMENTS")
        print(f"{'='*60}")

        for trial, A in b3_2_tours:
            scores = sorted([int(sum(A[i])) for i in range(n)])
            print(f"\n  Trial {trial}, scores={scores}")

            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
                gc.collect()
                b3_v = fast_beta3_lean(np.array(A_sub, dtype=np.int8), n-1)
                out_deg = int(sum(A[v]))

                tag = "GOOD" if b3_v == 0 else "BAD"
                print(f"    v={v} (out={out_deg}): b3(T\\v)={b3_v} [{tag}]")

    # Also check beta_3=1 patterns
    if b3_1_tours:
        print(f"\n{'='*60}")
        print(f"LES ANALYSIS OF {len(b3_1_tours)} BETA_3=1 TOURNAMENTS (comparison)")
        print(f"{'='*60}")

        for trial, A in b3_1_tours[:1]:
            scores = sorted([int(sum(A[i])) for i in range(n)])
            print(f"\n  Trial {trial}, scores={scores}")

            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
                gc.collect()
                b3_v = fast_beta3_lean(np.array(A_sub, dtype=np.int8), n-1)
                out_deg = int(sum(A[v]))

                tag = "GOOD" if b3_v == 0 else "BAD"
                print(f"    v={v} (out={out_deg}): b3(T\\v)={b3_v} [{tag}]")

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
