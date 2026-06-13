#!/usr/bin/env python3
"""beta2_unfiltered_large.py - Test unfiltered cone at n=9,10

The UNFILTERED single-vertex cone (using ALL T' 2-paths, not just
those where both cone 3-paths are allowed) ALWAYS works at n=7,8.

Test at n=9 and beyond. Also check: does the resulting 3-chain
actually live in Omega_3 (at n=7,8)?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def build_cone_B_unfiltered(A, n, v, paths2, path2_idx):
    """Build unfiltered cone B matrix."""
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    n_prime = len(others)
    paths2_Tp = enumerate_allowed_paths(B_sub, n_prime, 2)
    Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

    if not Tp_orig:
        return np.zeros((len(paths2), 0)), Tp_orig

    B = np.zeros((len(paths2), len(Tp_orig)))
    for j, (a, b, c) in enumerate(Tp_orig):
        terms = [
            ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
            ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
        ]
        for path, coeff in terms:
            if path in path2_idx:
                B[path2_idx[path], j] += coeff
    return B, Tp_orig


def compute_swap_cycles(A, n, v, paths2, path2_idx):
    P = [a for a in range(n) if a != v and A[v][a] == 1]
    Q = [b for b in range(n) if b != v and A[b][v] == 1]
    arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
    if len(arcs_PQ) < 2:
        return []

    m = len(arcs_PQ)
    rows = []
    for a in P:
        row = [0] * m
        for j, (a2, b2) in enumerate(arcs_PQ):
            if a2 == a:
                row[j] = 1
        if any(r != 0 for r in row):
            rows.append(row)
    for b in Q:
        row = [0] * m
        for j, (a2, b2) in enumerate(arcs_PQ):
            if b2 == b:
                row[j] = 1
        if any(r != 0 for r in row):
            rows.append(row)
    if not rows:
        return []

    C_mat = np.array(rows, dtype=float)
    Sc = np.linalg.svd(C_mat, compute_uv=False)
    rank_C = sum(s > 1e-8 for s in Sc)
    ker_dim = m - rank_C
    if ker_dim == 0:
        return []

    _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
    ker_basis = Vt[rank_C:]

    swap_vecs = []
    for ki in range(ker_basis.shape[0]):
        M_vec = ker_basis[ki]
        z = np.zeros(len(paths2))
        for j, (a, b) in enumerate(arcs_PQ):
            coeff = M_vec[j]
            if abs(coeff) < 1e-12:
                continue
            if (a, b, v) in path2_idx and (v, a, b) in path2_idx:
                z[path2_idx[(a, b, v)]] += coeff
                z[path2_idx[(v, a, b)]] -= coeff
        if np.max(np.abs(z)) > 1e-12:
            swap_vecs.append(z)
    return swap_vecs


# ============================================================
# PART 1: Unfiltered cone at n=9
# ============================================================
print("=" * 70)
print("PART 1: UNFILTERED CONE at n=9 (200 random)")
print("=" * 70)

n = 9
t0 = time.time()
total_ok = 0
total_fail = 0
total_swaps = 0

for trial in range(200):
    A = random_tournament(n)
    paths2 = enumerate_allowed_paths(A, n, 2)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

    any_fail = False
    for v in range(n):
        swaps = compute_swap_cycles(A, n, v, paths2, path2_idx)
        if not swaps:
            continue

        B, _ = build_cone_B_unfiltered(A, n, v, paths2, path2_idx)
        if B.shape[1] == 0:
            any_fail = True
            break

        for z in swaps:
            total_swaps += 1
            alpha, _, _, _ = np.linalg.lstsq(B, z, rcond=None)
            err = np.max(np.abs(B @ alpha - z))
            if err > 1e-6:
                any_fail = True
                break
        if any_fail:
            break

    if any_fail:
        total_fail += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  FAIL trial {trial}: scores={scores}")
    else:
        total_ok += 1

    if (trial + 1) % 40 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/200 ({elapsed:.0f}s) ok={total_ok} "
              f"fail={total_fail} swaps={total_swaps}")

elapsed = time.time() - t0
print(f"\nn=9 ({elapsed:.0f}s): unfiltered cone "
      f"OK={total_ok}/200, FAIL={total_fail}")
print(f"  Total swap cycles tested: {total_swaps}")


# ============================================================
# PART 2: Unfiltered cone rank surplus at n=7,8,9
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: RANK SURPLUS ANALYSIS")
print("=" * 70)

for n in [7, 8, 9]:
    random.seed(42)
    num_trials = 200 if n <= 8 else 100
    surpluses = []

    for trial in range(num_trials):
        A = random_tournament(n)
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        for v in range(n):
            swaps = compute_swap_cycles(A, n, v, paths2, path2_idx)
            if not swaps:
                continue

            B, _ = build_cone_B_unfiltered(A, n, v, paths2, path2_idx)
            rank_B = np.linalg.matrix_rank(B, tol=1e-8)

            swap_mat = np.column_stack(swaps)
            swap_dim = np.linalg.matrix_rank(swap_mat, tol=1e-8)

            surpluses.append(rank_B - swap_dim)

    if surpluses:
        print(f"\nn={n}: {len(surpluses)} cases")
        print(f"  Surplus (rank_B - swap_dim): "
              f"min={min(surpluses)}, max={max(surpluses)}, "
              f"mean={np.mean(surpluses):.1f}")
        print(f"  rank(B) >= swap_dim always? "
              f"{'YES' if min(surpluses) >= 0 else 'NO'}")


# ============================================================
# PART 3: Direct beta_2 at n=10 (small sample)
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: DIRECT BETA_2 at n=10 (50 random)")
print("=" * 70)

n = 10
t0 = time.time()
b2_zero = 0
b2_nonzero = 0

for trial in range(50):
    A = random_tournament(n)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        b2_zero += 1
        continue

    D2 = build_full_boundary_matrix(
        [tuple(p) for p in paths2], [tuple(p) for p in paths1]
    )
    D2_om = D2 @ omega2
    r2 = np.linalg.matrix_rank(D2_om, tol=1e-8)
    ker_d2 = dim_O2 - r2

    if ker_d2 == 0:
        b2_zero += 1
        continue

    if dim_O3 == 0:
        b2_nonzero += 1
        continue

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    D3_om = D3 @ omega3

    # Get ker(d2) basis
    _, _, Vt2 = np.linalg.svd(D2_om, full_matrices=True)
    ker_basis = omega2 @ Vt2[r2:].T

    combined = np.hstack([D3_om, ker_basis])
    r_comb = np.linalg.matrix_rank(combined, tol=1e-8)
    r_d3 = np.linalg.matrix_rank(D3_om, tol=1e-8)

    if r_comb == r_d3:
        b2_zero += 1
    else:
        b2_nonzero += 1
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  n={n} NONZERO beta_2 at trial {trial}! scores={scores}")

    if (trial + 1) % 10 == 0:
        elapsed = time.time() - t0
        print(f"  n={n}: {trial + 1}/50 ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=10 ({elapsed:.0f}s): beta_2=0: {b2_zero}/50, "
      f"beta_2>0: {b2_nonzero}")


# ============================================================
# PART 4: Dimension growth analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: DIMENSION ANALYSIS")
print("=" * 70)

for n in range(5, 10):
    random.seed(42)
    A = random_tournament(n)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    D2 = build_full_boundary_matrix(
        [tuple(p) for p in paths2], [tuple(p) for p in paths1]
    )
    D2_om = D2 @ omega2 if dim_O2 > 0 else np.zeros((len(paths1), 0))
    r2 = np.linalg.matrix_rank(D2_om, tol=1e-8) if dim_O2 > 0 else 0

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    D3_om = D3 @ omega3 if dim_O3 > 0 else np.zeros((len(paths2), 0))
    r3 = np.linalg.matrix_rank(D3_om, tol=1e-8) if dim_O3 > 0 else 0

    ker_d2 = dim_O2 - r2

    # T' path counts
    v = 0
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    n_prime = len(others)
    paths2_Tp = enumerate_allowed_paths(B_sub, n_prime, 2)

    print(f"n={n}: |A_2|={len(paths2)}, dim(Om2)={dim_O2}, "
          f"rank(d2)={r2}, ker(d2)={ker_d2}, "
          f"dim(Om3)={dim_O3}, rank(d3)={r3}, "
          f"|T'_paths|={len(paths2_Tp)}")


print("\n\nDone.")
