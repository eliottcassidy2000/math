#!/usr/bin/env python3
"""beta2_cone_failure.py - Investigate cone failure at n=8

At n=8, 1 out of 1000 random tournaments has filtered cone failure.
Score (3,3,3,3,4,4,4,4). Diagnose:
- Which vertex v fails?
- Is the failure because valid_Tp is too small (rank deficient)?
- Does the UNfiltered (full B) still work?
- Does the direct Omega_3 approach work (it should, since beta_2=0)?

Also: systematically check how often the filtered cone fails at n=7,8

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


# First, reproduce the failure
print("=" * 70)
print("REPRODUCING n=8 CONE FAILURE")
print("=" * 70)

random.seed(42)
n = 8
failure_A = None
failure_trial = -1

for trial in range(1000):
    A = random_tournament(n)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}
    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )

    failed = False
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
            # Check for nonzero swap cycles
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
                    failed = True
                    break
            if failed:
                break
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

        rank_B = np.linalg.matrix_rank(B_fill, tol=1e-8)

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

            alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
            B_err = np.max(np.abs(B_fill @ alpha - z))

            if B_err > 1e-6:
                failed = True
                failure_A = [row[:] for row in A]
                failure_trial = trial

                scores = tuple(sorted([sum(row) for row in A]))
                print(f"\nFAILURE at trial {trial}, v={v}")
                print(f"  Scores: {scores}")
                print(f"  P={P}, Q={Q}")
                print(f"  #arcs_PQ={len(arcs_PQ)}, ker_dim={ker_dim}")
                print(f"  #valid_Tp={len(valid_Tp)}, #all_Tp={len(Tp_all)}")
                print(f"  rank_B={rank_B}")
                print(f"  B_err={B_err:.2e}")

                # Try with ALL T' paths (unfiltered)
                B_full = np.zeros((len(paths2), len(Tp_all)))
                for j, (a, b, c) in enumerate(Tp_all):
                    terms = [
                        ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                        ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
                    ]
                    for path, coeff in terms:
                        if path in path2_idx:
                            B_full[path2_idx[path], j] += coeff

                alpha_full, _, _, _ = np.linalg.lstsq(B_full, z, rcond=None)
                B_full_err = np.max(np.abs(B_full @ alpha_full - z))
                rank_B_full = np.linalg.matrix_rank(B_full, tol=1e-8)
                print(f"\n  UNFILTERED: rank_B={rank_B_full}, err={B_full_err:.2e}")

                # Try Omega_3 direct
                omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
                dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
                if dim_O3 > 0:
                    D3_om = D3 @ omega3
                    c_om, _, _, _ = np.linalg.lstsq(D3_om, z, rcond=None)
                    om_err = np.max(np.abs(D3_om @ c_om - z))
                    print(f"  OMEGA_3 direct: dim_O3={dim_O3}, err={om_err:.2e}")

                # Try each vertex as cone vertex
                print(f"\n  Per-vertex cone status:")
                for vv in range(n):
                    PP = [a for a in range(n) if a != vv and A[vv][a] == 1]
                    QQ = [b for b in range(n) if b != vv and A[b][vv] == 1]

                    oo = [x for x in range(n) if x != vv]
                    Bs, vl = get_induced(A, n, oo)
                    p2t = enumerate_allowed_paths(Bs, len(oo), 2)
                    ta = [tuple(vl[x] for x in p) for p in p2t]
                    vt = [(a, b, c) for (a, b, c) in ta
                          if A[vv][a] == 1 and A[c][vv] == 1]

                    print(f"    v={vv}: d_out={len(PP)}, d_in={len(QQ)}, "
                          f"#all_Tp={len(ta)}, #valid_Tp={len(vt)}")

                break

        if failed:
            break

    if failed:
        break

if not failed:
    print("No failure found (different random seed path?)")


# ============================================================
# PART 2: Multi-vertex cone - combine fillings from ALL vertices
# ============================================================
if failure_A is not None:
    print(f"\n{'=' * 70}")
    print("PART 2: MULTI-VERTEX CONE FOR FAILURE CASE")
    print("=" * 70)
    print("Can we combine cone fillings from multiple vertices to fill z?")

    A = failure_A
    n = 8
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}
    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    # Collect ALL valid cone columns from ALL vertices
    all_B_cols = []
    col_labels = []

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        for a, b, c in valid_Tp:
            col = np.zeros(len(paths2))
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    col[path2_idx[path]] += coeff
            all_B_cols.append(col)
            col_labels.append((v, a, b, c))

    B_multi = np.column_stack(all_B_cols) if all_B_cols else np.zeros((len(paths2), 0))
    rank_multi = np.linalg.matrix_rank(B_multi, tol=1e-8)
    print(f"\n  Multi-vertex B: {B_multi.shape[1]} columns, rank={rank_multi}")

    # Find the swap cycle that failed
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

            # Try multi-vertex cone
            alpha_m, _, _, _ = np.linalg.lstsq(B_multi, z, rcond=None)
            err_m = np.max(np.abs(B_multi @ alpha_m - z))
            print(f"\n  v={v}, swap {ki}: multi-cone err={err_m:.2e}")

            if err_m < 1e-6:
                # Which vertices contribute?
                used_vertices = set()
                for j, (vv, a, b, c) in enumerate(col_labels):
                    if abs(alpha_m[j]) > 1e-10:
                        used_vertices.add(vv)
                print(f"  Uses cone vertices: {sorted(used_vertices)}")

            # Also try Omega_3 direct (should always work)
            if dim_O3 > 0:
                D3_om = D3 @ omega3
                c_om, _, _, _ = np.linalg.lstsq(D3_om, z, rcond=None)
                om_err = np.max(np.abs(D3_om @ c_om - z))
                print(f"  Omega_3 direct: err={om_err:.2e}")


# ============================================================
# PART 3: Failure rate vs n
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: SINGLE-VERTEX CONE FAILURE RATE")
print("=" * 70)

for n in [7, 8, 9]:
    random.seed(42)
    num_trials = 500 if n <= 8 else 100
    single_ok = 0
    single_fail = 0
    multi_ok = 0
    multi_fail = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        any_fail = False
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
                        any_fail = True
                        break
                if any_fail:
                    break
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

                alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
                if np.max(np.abs(B_fill @ alpha - z)) > 1e-6:
                    any_fail = True
                    break

            if any_fail:
                break

        if any_fail:
            single_fail += 1
        else:
            single_ok += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): single-vertex cone "
          f"OK={single_ok}/{num_trials}, FAIL={single_fail}/{num_trials} "
          f"({100*single_fail/num_trials:.1f}%)")


print("\n\nDone.")
