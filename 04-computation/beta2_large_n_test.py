#!/usr/bin/env python3
"""beta2_large_n_test.py - Test beta_2=0 and cone construction at n=8,9

Part 1: Direct beta_2 computation (n=8: 1000, n=9: 200)
Part 2: Filtered cone at n=8 (1000 random)
Part 3: Filtered cone at n=9 (100 random)
Part 4: Omega_3 membership verification at n=7,8

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import Counter
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


def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def test_beta2_direct(A, n):
    """Direct beta_2 computation. Returns beta_2."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        return 0

    D2 = build_full_boundary_matrix(
        [tuple(p) for p in paths2], [tuple(p) for p in paths1]
    )
    D2_om = D2 @ omega2
    r2 = np.linalg.matrix_rank(D2_om, tol=1e-8)
    ker_d2 = dim_O2 - r2

    if ker_d2 == 0:
        return 0

    if dim_O3 == 0:
        return ker_d2

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    D3_om = D3 @ omega3

    # im(d_3) in A_2 space
    combined_rank = np.linalg.matrix_rank(
        np.hstack([D3_om,
                   omega2 @ np.linalg.svd(D2_om, full_matrices=True)[2][r2:].T]),
        tol=1e-8
    )
    im_d3_rank = np.linalg.matrix_rank(D3_om, tol=1e-8)

    return combined_rank - im_d3_rank


def test_filtered_cone_single(A, n):
    """Test filtered cone for one tournament. Returns (all_ok, num_tested, rank_data)."""
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )

    all_ok = True
    num_tested = 0
    rank_data = []

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
            # Check if there are actual nonzero swap cycles
            has_nonzero = False
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
                    has_nonzero = True
                    break
            if has_nonzero:
                all_ok = False
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

            num_tested += 1

            alpha, _, _, _ = np.linalg.lstsq(B_fill, z, rcond=None)
            B_err = np.max(np.abs(B_fill @ alpha - z))

            if B_err > 1e-6:
                all_ok = False
            else:
                # Verify d3(w) = z
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

                d3w = D3 @ w_full
                d3_err = np.max(np.abs(d3w - z))
                if d3_err > 1e-6:
                    all_ok = False

        rank_data.append({
            'swap_dim': ker_dim, 'num_valid': len(valid_Tp),
            'rank_B': rank_B
        })

    return all_ok, num_tested, rank_data


# ============================================================
# PART 1: Direct beta_2 at n=8, n=9
# ============================================================
print("=" * 70)
print("PART 1: DIRECT BETA_2 COMPUTATION")
print("=" * 70)

for n, num_trials in [(8, 1000), (9, 200)]:
    t0 = time.time()
    b2_zero = 0
    b2_nonzero = 0

    for trial in range(num_trials):
        b2 = test_beta2_direct(random_tournament(n), n)
        if b2 == 0:
            b2_zero += 1
        else:
            b2_nonzero += 1
            print(f"  n={n} NONZERO beta_2={b2} at trial {trial}!")

        if (trial + 1) % (num_trials // 5) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"zero={b2_zero} nonzero={b2_nonzero}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): beta_2=0: {b2_zero}/{num_trials}, "
          f"beta_2>0: {b2_nonzero}")

# Paley tournaments
for p in [7, 11]:
    A = paley_tournament(p)
    b2 = test_beta2_direct(A, p)
    print(f"  Paley T_{p}: beta_2 = {b2}")


# ============================================================
# PART 2: Filtered cone at n=8
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: FILTERED CONE at n=8 (1000 random)")
print("=" * 70)

n = 8
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0
surplus_data = []

for trial in range(1000):
    A = random_tournament(n)
    ok, tested, rdata = test_filtered_cone_single(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL trial {trial}: scores={scores}")

    for rd in rdata:
        surplus_data.append(rd['rank_B'] - rd['swap_dim'])

    if (trial + 1) % 200 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/1000 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

elapsed = time.time() - t0
print(f"\nn=8 ({elapsed:.0f}s): {total_tested} swap cycles")
print(f"  OK: {total_ok}/1000, Failures: {total_fail}")
if surplus_data:
    print(f"  Surplus min={min(surplus_data)}, max={max(surplus_data)}, "
          f"mean={np.mean(surplus_data):.1f}")


# ============================================================
# PART 3: Filtered cone at n=9
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: FILTERED CONE at n=9 (100 random)")
print("=" * 70)

n = 9
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for trial in range(100):
    A = random_tournament(n)
    ok, tested, rdata = test_filtered_cone_single(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL trial {trial}: scores={scores}")

    if (trial + 1) % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {trial + 1}/100 ({elapsed:.0f}s) tested={total_tested} "
              f"ok={total_ok} fail={total_fail}")

elapsed = time.time() - t0
print(f"\nn=9 ({elapsed:.0f}s): {total_tested} swap cycles")
print(f"  OK: {total_ok}/100, Failures: {total_fail}")


# ============================================================
# PART 4: Omega_3 auto-membership at n=7,8
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: OMEGA_3 AUTO-MEMBERSHIP at n=7,8")
print("=" * 70)

for n, num_trials in [(7, 200), (8, 100)]:
    t0 = time.time()
    total_tested = 0
    auto_ok = 0
    auto_fail = 0

    for trial in range(num_trials):
        A = random_tournament(n)
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

            for ki in range(min(ker_basis.shape[0], 3)):
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
                if np.max(np.abs(B_fill @ alpha - z)) > 1e-6:
                    continue

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

                if dim_O3 > 0:
                    proj, _, _, _ = np.linalg.lstsq(omega3, w_full, rcond=None)
                    om_err = np.max(np.abs(w_full - omega3 @ proj))
                else:
                    om_err = np.max(np.abs(w_full))

                if om_err < 1e-8:
                    auto_ok += 1
                else:
                    auto_fail += 1

        if (trial + 1) % (num_trials // 5) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"tested={total_tested}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {total_tested} fillings checked")
    print(f"  Auto in Omega_3: {auto_ok}/{total_tested}")
    print(f"  NOT in Omega_3: {auto_fail}/{total_tested}")


print("\n\nDone.")
