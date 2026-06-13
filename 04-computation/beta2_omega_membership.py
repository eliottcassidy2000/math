#!/usr/bin/env python3
"""beta2_omega_membership.py - Verify cone filling is automatically in Omega_3

Key question: when we solve B*alpha = z and construct
w = sum alpha_{abc} * [(v,a,b,c) + (a,b,c,v)],
is w AUTOMATICALLY in Omega_3, or do we need to project?

Also: check if the B columns are themselves in Omega_2 (as 2-chains).

Part 1: Omega_3 membership of cone filling (n=5 exhaustive, n=6 sampled)
Part 2: Omega_2 membership of B columns
Part 3: Direct beta_2 computation for sanity check

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


# ============================================================
# PART 1: Omega_3 membership
# ============================================================
print("=" * 70)
print("PART 1: IS CONE FILLING AUTOMATICALLY IN OMEGA_3?")
print("=" * 70)

for n in [5, 6]:
    total_tested = 0
    auto_in_omega = 0
    not_in_omega = 0
    max_omega_err = 0

    gen = all_tournaments_gen(n) if n <= 6 else None
    items = list(gen) if n <= 6 else [random_tournament(n) for _ in range(500)]

    for A in items:
        paths1 = enumerate_allowed_paths(A, n, 1)
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)
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

            others = [x for x in range(n) if x != v]
            B_sub, vlist = get_induced(A, n, others)
            paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
            Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
            valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                        if A[v][a] == 1 and A[c][v] == 1]

            if not valid_Tp:
                continue

            B_fill = np.zeros((len(paths2), len(valid_Tp)))
            for j, (a, b, c) in enumerate(valid_Tp):
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

                alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
                err = np.max(np.abs(B_fill @ alpha - z))
                if err > 1e-6:
                    continue

                # Build 3-chain
                w_full = np.zeros(len(paths3))
                for j, (a, b, c) in enumerate(valid_Tp):
                    if abs(alpha[j]) < 1e-12:
                        continue
                    front = (v, a, b, c)
                    back = (a, b, c, v)
                    if front in path3_idx:
                        w_full[path3_idx[front]] += alpha[j]
                    if back in path3_idx:
                        w_full[path3_idx[back]] += alpha[j]

                # Check Omega_3 membership
                if dim_O3 > 0:
                    proj, _, _, _ = np.linalg.lstsq(omega3, w_full, rcond=None)
                    om_err = np.max(np.abs(w_full - omega3 @ proj))
                else:
                    om_err = np.max(np.abs(w_full))

                max_omega_err = max(max_omega_err, om_err)

                if om_err < 1e-8:
                    auto_in_omega += 1
                else:
                    not_in_omega += 1

    print(f"\nn={n}: {total_tested} swap cycles")
    print(f"  Auto in Omega_3: {auto_in_omega}/{total_tested}")
    print(f"  NOT in Omega_3: {not_in_omega}/{total_tested}")
    print(f"  Max Omega_3 error: {max_omega_err:.2e}")


# ============================================================
# PART 2: B columns in Omega_2?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: ARE B COLUMNS IN OMEGA_2?")
print("=" * 70)

n = 5
total_cols = 0
cols_in_om2 = 0

for A in all_tournaments_gen(n):
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
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

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        if not valid_Tp:
            continue

        B_fill = np.zeros((len(paths2), len(valid_Tp)))
        for j, (a, b, c) in enumerate(valid_Tp):
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    B_fill[path2_idx[path], j] += coeff

        if dim_O2 > 0:
            for j in range(B_fill.shape[1]):
                col = B_fill[:, j]
                total_cols += 1
                if np.max(np.abs(col)) < 1e-12:
                    cols_in_om2 += 1
                    continue
                proj, _, _, _ = np.linalg.lstsq(omega2, col, rcond=None)
                err = np.max(np.abs(col - omega2 @ proj))
                if err < 1e-8:
                    cols_in_om2 += 1

print(f"\nn=5: B columns in Omega_2: {cols_in_om2}/{total_cols}")
print(f"  (Each B column = d_3[(v,a,b,c)+(a,b,c,v)] restricted to v-paths)")


# ============================================================
# PART 3: Direct beta_2 at n=7,8
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: DIRECT BETA_2 COMPUTATION (n=7,8 sampled)")
print("=" * 70)

for n in [7, 8]:
    beta2_zero = 0
    beta2_nonzero = 0
    num_trials = 500 if n <= 7 else 200

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        paths1 = enumerate_allowed_paths(A, n, 1)
        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)

        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

        if dim_O2 == 0:
            beta2_zero += 1
            continue

        D2 = build_full_boundary_matrix(
            [tuple(p) for p in paths2], [tuple(p) for p in paths1]
        )
        D2_om = D2 @ omega2
        r2 = np.linalg.matrix_rank(D2_om, tol=1e-8)
        ker_d2 = dim_O2 - r2

        if ker_d2 == 0:
            beta2_zero += 1
            continue

        if dim_O3 == 0:
            beta2_nonzero += 1
            continue

        D3 = build_full_boundary_matrix(
            [tuple(p) for p in paths3], [tuple(p) for p in paths2]
        )
        D3_om = D3 @ omega3
        im_d3 = np.linalg.matrix_rank(D3_om, tol=1e-8)

        # Get ker(d2) basis
        _, _, Vt2 = np.linalg.svd(D2_om, full_matrices=True)
        ker_basis = omega2 @ Vt2[r2:].T

        # Check ker(d2) subset im(d3)
        combined = np.hstack([D3_om, ker_basis])
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        rank_d3 = np.linalg.matrix_rank(D3_om, tol=1e-8)

        if rank_combined == rank_d3:
            beta2_zero += 1
        else:
            beta2_nonzero += 1
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  n={n} NONZERO beta_2! trial={trial}, scores={scores}")
            print(f"    ker(d2)={ker_d2}, im(d3)={im_d3}, "
                  f"combined_rank={rank_combined}")

        if (trial + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): beta_2=0: {beta2_zero}/{num_trials}, "
          f"beta_2>0: {beta2_nonzero}")


print("\n\nDone.")
