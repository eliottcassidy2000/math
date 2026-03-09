#!/usr/bin/env python3
"""beta2_omega3_check.py - Check if cone-from-T' filling is in Omega_3

The B*alpha=z system always has a solution, but we need the filling
w = sum alpha_{abc} * [(v,a,b,c) + (a,b,c,v)]
to be in Omega_3, not just A_3.

Two approaches:
1. Project alpha onto Omega_3 and solve there
2. Check when the naive filling is already in Omega_3

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


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


random.seed(42)

# ============================================================
# PART 1: Detailed check at n=5
# ============================================================
print("=" * 70)
print("PART 1: OMEGA_3 CHECK AT n=5 (exhaustive)")
print("=" * 70)

n = 5
total_tested = 0
naive_in_omega = 0
projected_works = 0
total_fail = 0

for A in all_tournaments_gen(n):
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

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
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C
        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        ker_basis = Vt[rank_C:]

        # Get T' paths
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        num_Tp = len(Tp_paths_orig)

        # Build B matrix
        B_fill = np.zeros((len(paths2), num_Tp))
        for j, (a, b, c) in enumerate(Tp_paths_orig):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

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

            if np.max(np.abs(z)) < 1e-12:
                continue

            total_tested += 1

            # Method 1: Naive solve B*alpha = z, build w, check if in Omega_3
            alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
            w_full = np.zeros(len(paths3))
            for j, (a, b, c) in enumerate(Tp_paths_orig):
                if abs(alpha[j]) < 1e-12:
                    continue
                front = (v, a, b, c)
                back = (a, b, c, v)
                if front in path3_idx:
                    w_full[path3_idx[front]] += alpha[j]
                if back in path3_idx:
                    w_full[path3_idx[back]] += alpha[j]

            # Check d_3(w) = z
            check_d3 = D3 @ w_full
            d3_err = np.max(np.abs(check_d3 - z))

            # Check w in Omega_3
            if dim_O3 > 0:
                proj, _, _, _ = np.linalg.lstsq(omega3, w_full, rcond=None)
                w_projected = omega3 @ proj
                omega_err = np.max(np.abs(w_full - w_projected))
            else:
                omega_err = np.max(np.abs(w_full))

            if omega_err < 1e-8:
                naive_in_omega += 1

            # Method 2: Use Omega_3 basis directly
            if dim_O3 > 0:
                D3_omega = D3 @ omega3
                w_om, _, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
                err_om = np.max(np.abs(D3_omega @ w_om - z))
                if err_om < 1e-6:
                    projected_works += 1
                else:
                    total_fail += 1
            else:
                total_fail += 1

print(f"n=5: {total_tested} swap cycles tested")
print(f"  Naive filling in Omega_3: {naive_in_omega}/{total_tested}")
print(f"  Omega_3-projected filling works: {projected_works}/{total_tested}")
print(f"  Failures: {total_fail}")


# ============================================================
# PART 2: The CORRECT approach - d_3 restricted to Omega_3
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: CORRECT APPROACH - solve in Omega_3 directly")
print("=" * 70)
print("Instead of cone-from-T', solve d_3|_{Omega_3} * c = z")
print("This is what beta_2=0 actually requires.")
print("The cone-from-T' insight: im(d_3|_{Omega_3}) >= swap cycle space")

n = 5
print(f"\nn={n}: Checking im(d_3|_Omega_3) vs ker(d_2|_Omega_2)")

total_tours = 0
beta2_zero_count = 0

for A in all_tournaments_gen(n):
    total_tours += 1
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        beta2_zero_count += 1
        continue

    D2 = build_full_boundary_matrix(
        [tuple(p) for p in paths2], [tuple(p) for p in paths1]
    )
    D2_om = D2 @ omega2

    U2, S2, _ = np.linalg.svd(D2_om, full_matrices=False)
    r2 = sum(s > 1e-8 for s in S2)
    ker_d2_dim = dim_O2 - r2

    if ker_d2_dim == 0:
        beta2_zero_count += 1
        continue

    if dim_O3 == 0:
        # ker(d_2) != 0 but Omega_3 = 0 => beta_2 > 0
        continue

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    D3_om = D3 @ omega3

    # im(d_3|_Omega_3) lives in A_2 space
    # ker(d_2|_Omega_2) also in A_2 space (via omega2 embedding)
    # beta_2 = 0 iff ker(d_2) subset im(d_3)

    # Get ker(d_2) basis in A_2
    _, _, Vt2 = np.linalg.svd(D2_om, full_matrices=True)
    ker_d2_basis = omega2 @ Vt2[r2:].T  # columns are ker(d_2) in A_2

    # Check each ker(d_2) vector is in im(d_3)
    all_exact = True
    for k in range(ker_d2_dim):
        z = ker_d2_basis[:, k]
        c, _, _, _ = np.linalg.lstsq(D3_om, z, rcond=None)
        err = np.max(np.abs(D3_om @ c - z))
        if err > 1e-6:
            all_exact = False
            break

    if all_exact:
        beta2_zero_count += 1

print(f"  beta_2 = 0: {beta2_zero_count}/{total_tours}")


# ============================================================
# PART 3: Why does cone-from-T' work? Column space analysis
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: COLUMN SPACE RELATIONSHIP")
print("=" * 70)
print("For each v, the B_fill columns live in a specific subspace of A_2.")
print("Key question: does col(B_fill) intersected with Omega_2 cover")
print("all swap cycles?")

n = 5
for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx > 20:
        break

    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
        if len(arcs_PQ) < 2:
            continue

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
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C
        if ker_dim == 0:
            continue

        # Check: are B columns in image of d_3|_{Omega_3}?
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_paths_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        B_fill = np.zeros((len(paths2), len(Tp_paths_orig)))
        for j, (a, b, c) in enumerate(Tp_paths_orig):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        if dim_O3 > 0:
            D3 = build_full_boundary_matrix(
                [tuple(p) for p in paths3], [tuple(p) for p in paths2]
            )
            D3_om = D3 @ omega3
            rank_D3 = np.linalg.matrix_rank(D3_om, tol=1e-8)

            # Check each B column is in im(D3_om)
            cols_in_im = 0
            for j in range(B_fill.shape[1]):
                col = B_fill[:, j]
                if np.max(np.abs(col)) < 1e-12:
                    cols_in_im += 1
                    continue
                c, _, _, _ = np.linalg.lstsq(D3_om, col, rcond=None)
                err = np.max(np.abs(D3_om @ c - col))
                if err < 1e-6:
                    cols_in_im += 1

            print(f"  T#{tidx} v={v}: B cols in im(d3|_Om3): "
                  f"{cols_in_im}/{B_fill.shape[1]}, "
                  f"rank(d3|_Om3)={rank_D3}, dim_Om3={dim_O3}")

print("\n\nDone.")
