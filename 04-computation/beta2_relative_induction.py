#!/usr/bin/env python3
"""beta2_relative_induction.py - Clean inductive proof via relative homology

STRATEGY: beta_2(T) = 0 by induction on n.
  Base: n <= 4 trivially (no 2-cycles).
  Step: For tournament T on n vertices, pick vertex v.
    LES: 0 = H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)
    By induction, H_2(T\v) = 0. So H_2(T) injects into H_2(T,T\v).
    If H_2(T,T\v) = 0 for SOME v, then H_2(T) = 0.

This script computes H_2(T,T\v) via rank-nullity on relative complex.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
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


def compute_h2_rel(A, n, v):
    """Compute dim H_2(T, T\v) directly.

    Uses: beta_2(T) = dim(Z_2(T)) - rk(d_3(T))
    And the LES relation. But more directly:

    H_2(T,T\v) = ker(d_2^rel: R_2 -> R_1) / im(d_3^rel: R_3 -> R_2)

    R_p = Omega_p(T) / Omega_p(T\v)

    SHORTCUT: Use the LES itself.
    0 -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)

    If we know beta_2(T), beta_1(T), beta_1(T\v), and the connecting map delta,
    we can compute H_2(T,T\v).

    Actually: dim H_2(T,T\v) = beta_2(T) + dim(ker(i_*: H_1(T\v) -> H_1(T)))
    from the LES (where i_* = inclusion-induced map).

    Since beta_2(T) should be 0:
    dim H_2(T,T\v) = dim(ker(i_*)).

    So: H_2(T,T\v) = 0 iff i_*: H_1(T\v) -> H_1(T) is injective.
    """
    # Compute beta_2(T)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)
    paths0 = [(i,) for i in range(n)]

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        return {'beta2': 0, 'h2_rel': 0, 'b1_T': 0, 'b1_Tv': 0, 'ker_istar': 0}

    D2 = build_full_boundary_matrix([tuple(p) for p in paths2], [tuple(p) for p in paths1])
    D2_om = D2 @ omega2
    sv2 = np.linalg.svd(D2_om, compute_uv=False)
    rk_d2 = int(sum(s > 1e-8 for s in sv2))
    z2_dim = dim_O2 - rk_d2

    if paths3:
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
        if dim_O3 > 0:
            D3 = build_full_boundary_matrix([tuple(p) for p in paths3], [tuple(p) for p in paths2])
            D3_om = D3 @ omega3
            sv3 = np.linalg.svd(D3_om, compute_uv=False)
            rk_d3 = int(sum(s > 1e-8 for s in sv3))
        else:
            rk_d3 = 0
    else:
        rk_d3 = 0

    beta2 = z2_dim - rk_d3

    # Compute beta_1(T)
    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 > 0:
        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        sv1 = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv1))
        b1_T = dim_O1 - rk_d1 - rk_d2
    else:
        b1_T = 0

    # Compute beta_1(T\v)
    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    n_prime = len(others)

    paths1_Tv = enumerate_allowed_paths(B_sub, n_prime, 1)
    paths0_Tv = [(i,) for i in range(n_prime)]
    paths2_Tv = enumerate_allowed_paths(B_sub, n_prime, 2)

    if paths1_Tv:
        omega1_Tv = compute_omega_basis(B_sub, n_prime, 1, paths1_Tv, paths0_Tv)
        dim_O1_Tv = omega1_Tv.shape[1] if omega1_Tv.ndim == 2 else 0
    else:
        dim_O1_Tv = 0

    if dim_O1_Tv > 0:
        D1_Tv = build_full_boundary_matrix([tuple(p) for p in paths1_Tv], paths0_Tv)
        D1_Tv_om = D1_Tv @ omega1_Tv
        sv1_Tv = np.linalg.svd(D1_Tv_om, compute_uv=False)
        rk_d1_Tv = int(sum(s > 1e-8 for s in sv1_Tv))

        if paths2_Tv:
            omega2_Tv = compute_omega_basis(B_sub, n_prime, 2, paths2_Tv, paths1_Tv)
            dim_O2_Tv = omega2_Tv.shape[1] if omega2_Tv.ndim == 2 else 0
            if dim_O2_Tv > 0:
                D2_Tv = build_full_boundary_matrix([tuple(p) for p in paths2_Tv],
                                                   [tuple(p) for p in paths1_Tv])
                D2_Tv_om = D2_Tv @ omega2_Tv
                sv2_Tv = np.linalg.svd(D2_Tv_om, compute_uv=False)
                rk_d2_Tv = int(sum(s > 1e-8 for s in sv2_Tv))
            else:
                rk_d2_Tv = 0
        else:
            rk_d2_Tv = 0

        b1_Tv = dim_O1_Tv - rk_d1_Tv - rk_d2_Tv
    else:
        b1_Tv = 0

    # Compute ker(i_*: H_1(T\v) -> H_1(T))
    if b1_Tv == 0:
        ker_istar = 0
    else:
        # Get H_1(T\v) cycle basis
        _, _, Vt_Tv = np.linalg.svd(D1_Tv_om, full_matrices=True)
        # Cycles in Omega_1(T\v) coords
        z1_Tv_om = Vt_Tv[rk_d1_Tv:].T  # Omega_1_Tv coords

        # Map to A_1(T\v) coords
        z1_Tv_A = omega1_Tv @ z1_Tv_om  # columns in A_1(T\v)

        # Remove boundaries (im d_2 in T\v)
        if rk_d2_Tv > 0:
            # Get im(d_2_Tv) in A_1(T\v)
            im_d2_Tv = D2_Tv_om  # already in A_1(T\v) coords
            # Project out im(d_2_Tv)
            U_d2, S_d2, _ = np.linalg.svd(im_d2_Tv, full_matrices=True)
            proj = np.eye(im_d2_Tv.shape[0]) - U_d2[:, :rk_d2_Tv] @ U_d2[:, :rk_d2_Tv].T
            z1_Tv_H = proj @ z1_Tv_A
            # Rank of H_1(T\v) representatives
            rk_H1Tv = np.linalg.matrix_rank(z1_Tv_H, tol=1e-8)
        else:
            z1_Tv_H = z1_Tv_A
            rk_H1Tv = b1_Tv

        # Embed into A_1(T) coords
        path1_T_idx = {tuple(p): i for i, p in enumerate(paths1)}
        embed = np.zeros((len(paths1), len(paths1_Tv)))
        for j, p_Tv in enumerate(paths1_Tv):
            p_T = tuple(vlist[x] for x in p_Tv)
            if p_T in path1_T_idx:
                embed[path1_T_idx[p_T], j] = 1

        z1_in_T = embed @ z1_Tv_A  # columns are Z_1(T\v) cycles in A_1(T)

        # Check if these are boundaries in T (i.e., in im(d_2|Omega_2(T)))
        if dim_O2 > 0:
            im_d2_T = D2_om  # in A_1(T) coords
        else:
            im_d2_T = np.zeros((len(paths1), 0))

        ker_istar = 0
        for k in range(z1_in_T.shape[1]):
            z = z1_in_T[:, k]
            if np.max(np.abs(z)) < 1e-12:
                ker_istar += 1
                continue
            if im_d2_T.shape[1] > 0:
                c, _, _, _ = np.linalg.lstsq(im_d2_T, z, rcond=None)
                err = np.max(np.abs(im_d2_T @ c - z))
                if err < 1e-6:
                    ker_istar += 1

    # LES: H_2(T,T\v) = beta_2(T) + ker_istar
    # (by exactness: 0 -> H_2(T) -> H_2(T,T\v) -> ker(i_*) -> 0)
    h2_rel = beta2 + ker_istar

    return {
        'beta2': beta2, 'h2_rel': h2_rel,
        'b1_T': b1_T, 'b1_Tv': b1_Tv,
        'ker_istar': ker_istar,
        'z2_dim': z2_dim, 'rk_d3': rk_d3,
    }


# ============================================================
# Part 1: n=5 exhaustive — find best vertex for each tournament
# ============================================================
print("=" * 70)
print("RELATIVE HOMOLOGY ANALYSIS")
print("=" * 70)

for n in [5, 6]:
    t0 = time.time()

    if n <= 6:
        gen = list(all_tournaments_gen(n))

    total_tours = 0
    has_good_v = 0  # has SOME v with H_2^rel = 0 (i_* injective)
    all_v_bad = 0   # ALL v have H_2^rel > 0
    total_pairs = 0
    inj_ok = 0
    inj_fail = 0
    ker_istar_vals = []

    for tidx, A in enumerate(gen):
        total_tours += 1
        found_good = False

        for v in range(n):
            result = compute_h2_rel(A, n, v)
            total_pairs += 1

            if result['ker_istar'] == 0:
                inj_ok += 1
                found_good = True
            else:
                inj_fail += 1
                ker_istar_vals.append(result['ker_istar'])

        if found_good:
            has_good_v += 1
        else:
            all_v_bad += 1
            if all_v_bad <= 3:
                scores = sorted([sum(row) for row in A])
                print(f"  ALL BAD: T#{tidx}, scores={scores}")

        if (tidx + 1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {tidx+1}/{len(gen)} ({elapsed:.0f}s) "
                  f"good={has_good_v} bad={all_v_bad}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {total_tours} tournaments")
    print(f"  (T,v) with i_* injective: {inj_ok}/{total_pairs}")
    print(f"  (T,v) with i_* NOT injective: {inj_fail}/{total_pairs}")
    print(f"  Tournaments with SOME good v: {has_good_v}/{total_tours}")
    print(f"  Tournaments with NO good v: {all_v_bad}/{total_tours}")
    if ker_istar_vals:
        from collections import Counter
        dist = Counter(ker_istar_vals)
        print(f"  ker(i_*) distribution: {dict(sorted(dist.items()))}")


# ============================================================
# Part 2: Can we always find a good vertex? n=7,8
# ============================================================
print(f"\n{'=' * 70}")
print("n=7,8 SAMPLED")
print("=" * 70)

for n in [7, 8]:
    random.seed(42)
    num_trials = 200 if n == 7 else 100
    has_good = 0
    no_good = 0
    t0 = time.time()

    for trial in range(num_trials):
        A = random_tournament(n)
        found = False
        for v in range(n):
            d_out = sum(A[v])
            # Skip sources and sinks (they always give H_2^rel = 0 trivially)
            result = compute_h2_rel(A, n, v)
            if result['ker_istar'] == 0:
                found = True
                break

        if found:
            has_good += 1
        else:
            no_good += 1

        if (trial + 1) % (num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{num_trials} ({elapsed:.0f}s) "
                  f"good={has_good} bad={no_good}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {num_trials} tournaments")
    print(f"  Has SOME good v: {has_good}/{num_trials}")
    print(f"  ALL v bad: {no_good}/{num_trials}")


print("\n\nDone.")
