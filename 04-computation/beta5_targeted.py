#!/usr/bin/env python3
"""
beta5_targeted.py -- Targeted beta_5 computation for n=9 tournaments.

Only computes beta_5 (not beta_0 through beta_4), which avoids
needing to enumerate high-dimensional paths.

To compute beta_5 we need:
  Omega_4, Omega_5, Omega_6
  => allowed paths at p=3,4,5,6

For a tournament on 9 vertices, these have ~414, 1072, 2252, 3582 paths.
This is large but feasible.

kind-pasteur-2026-03-08-S34
"""

import numpy as np
from itertools import permutations
from collections import defaultdict
import random
import time
import sys

def enumerate_allowed_paths(A, n, p):
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    paths = []
    for perm in permutations(range(n), p + 1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space


def compute_beta5_only(A, n):
    """Compute only beta_5 for a tournament. Returns beta_5 value."""
    # We need allowed paths at p = 3, 4, 5, 6
    t0 = time.time()
    allowed = {}
    for p in [3, 4, 5, 6]:
        allowed[p] = enumerate_allowed_paths(A, n, p)
    t_enum = time.time() - t0

    # Compute Omega_4, Omega_5, Omega_6
    omega = {}
    omega[4] = compute_omega_basis(A, n, 4, allowed[4], allowed[3])
    omega[5] = compute_omega_basis(A, n, 5, allowed[5], allowed[4])
    omega[6] = compute_omega_basis(A, n, 6, allowed[6], allowed[5])

    # beta_5 = dim ker(del_5: Omega_5 -> Omega_4) - dim im(del_6: Omega_6 -> Omega_5)
    dim_omega_5 = omega[5].shape[1] if omega[5].ndim == 2 else 0
    if dim_omega_5 == 0:
        return 0, t_enum

    # del_5: A_5 -> A_4, restricted to Omega_5
    bd_5 = build_full_boundary_matrix(allowed[5], allowed[4])
    bd_5_omega = bd_5 @ omega[5]
    if bd_5_omega.shape[0] > 0 and bd_5_omega.shape[1] > 0:
        S5 = np.linalg.svd(bd_5_omega, compute_uv=False)
        rank_5 = sum(s > 1e-8 for s in S5)
    else:
        rank_5 = 0

    ker_dim = dim_omega_5 - rank_5

    # im(del_6: Omega_6 -> A_5)
    dim_omega_6 = omega[6].shape[1] if omega[6].ndim == 2 else 0
    if dim_omega_6 > 0:
        bd_6 = build_full_boundary_matrix(allowed[6], allowed[5])
        bd_6_omega = bd_6 @ omega[6]
        S6 = np.linalg.svd(bd_6_omega, compute_uv=False)
        im_dim = sum(s > 1e-8 for s in S6)
    else:
        im_dim = 0

    beta_5 = max(0, ker_dim - im_dim)
    return beta_5, t_enum


def compute_beta_range(A, n, dim_low, dim_high):
    """Compute beta_{dim_low} through beta_{dim_high}."""
    # Need allowed paths from (dim_low - 1) to (dim_high + 1)
    allowed = {}
    for p in range(max(0, dim_low - 1), dim_high + 2):
        allowed[p] = enumerate_allowed_paths(A, n, p)

    # Need Omega from dim_low to dim_high + 1
    omega = {}
    for p in range(max(0, dim_low - 1), dim_high + 2):
        if p == 0:
            omega[p] = np.eye(len(allowed[0]))
        elif p in allowed and (p-1) in allowed:
            omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        else:
            omega[p] = np.zeros((0, 0))

    betti = []
    for p in range(dim_low, dim_high + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue

        bd_p = build_full_boundary_matrix(allowed[p], allowed.get(p-1, []))
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p

        dim_omega_p1 = omega.get(p+1, np.zeros((0,0)))
        if hasattr(dim_omega_p1, 'ndim'):
            dim_omega_p1_val = dim_omega_p1.shape[1] if dim_omega_p1.ndim == 2 else 0
        else:
            dim_omega_p1_val = 0
        if dim_omega_p1_val > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0

        beta_p = max(0, ker_dim - im_dim)
        betti.append(beta_p)

    return betti


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def score_sequence(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))


def main():
    n = 9

    print("=" * 70)
    print(f"TARGETED BETA_5 CHECK FOR n={n} TOURNAMENTS")
    print("=" * 70)

    # Timing test: beta_5 only vs full
    print("\nTiming comparison:")
    random.seed(42)
    A_test = random_tournament(n)

    t0 = time.time()
    b5, t_enum = compute_beta5_only(A_test, n)
    dt_targeted = time.time() - t0
    print(f"  beta_5 only: {dt_targeted:.2f}s (enum: {t_enum:.2f}s), beta_5 = {b5}")

    t0 = time.time()
    betti_45 = compute_beta_range(A_test, n, 4, 5)
    dt_range = time.time() - t0
    print(f"  beta_4,5: {dt_range:.2f}s, beta = {betti_45}")

    # Determine sample size
    dt = dt_targeted
    if dt < 5:
        NUM_SAMPLES = 500
    elif dt < 15:
        NUM_SAMPLES = 200
    elif dt < 30:
        NUM_SAMPLES = 100
    else:
        NUM_SAMPLES = 50

    print(f"\n  Per-tournament time: {dt:.1f}s")
    print(f"  Planning {NUM_SAMPLES} random + 20 top-H = {NUM_SAMPLES + 20} total")
    print(f"  Estimated: {(NUM_SAMPLES + 20) * dt / 60:.1f} min")
    sys.stdout.flush()

    # Phase 1: H survey
    print(f"\n--- Phase 1: H-value survey (20000 random) ---")
    random.seed(12345)
    h_survey = []
    t_start = time.time()
    for _ in range(20000):
        A = random_tournament(n)
        H = ham_path_count(A, n)
        h_survey.append((H, A))
    print(f"  Done in {time.time()-t_start:.1f}s")

    h_values = [x[0] for x in h_survey]
    print(f"  H range: [{min(h_values)}, {max(h_values)}]")
    print(f"  H mean: {np.mean(h_values):.1f}")

    h_survey.sort(key=lambda x: -x[0])
    top_unique = []
    seen = set()
    for H, A in h_survey:
        key = str(A)
        if key not in seen:
            seen.add(key)
            top_unique.append((H, A))
        if len(top_unique) >= 20:
            break

    # Also include some LOW-H and MEDIUM-H for variety
    random.shuffle(h_survey)
    for H, A in h_survey:
        key = str(A)
        if key not in seen and H < 500:
            seen.add(key)
            top_unique.append((H, A))
            if len(top_unique) >= 25:
                break

    print(f"\n  Top 20 by H:")
    for i in range(min(20, len(top_unique))):
        H, A = top_unique[i]
        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        print(f"    #{i+1}: H={H}, t3={t3}, score={sc}")
    sys.stdout.flush()

    # Phase 2: Compute beta_4 and beta_5 for top-H tournaments
    print(f"\n--- Phase 2: beta_4,5 for top-H + selected tournaments ---")
    all_results = []
    beta5_positive = []

    for i, (H, A) in enumerate(top_unique):
        t0 = time.time()
        betti_45 = compute_beta_range(A, n, 4, 5)
        elapsed = time.time() - t0
        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        all_results.append((H, betti_45, t3, sc))

        flag = " <-- BETA_5 > 0!" if betti_45[1] > 0 else ""
        print(f"    #{i+1}: H={H}, t3={t3}, [b4,b5]={betti_45} ({elapsed:.1f}s){flag}")

        if betti_45[1] > 0:
            beta5_positive.append((H, betti_45, t3, sc))
        sys.stdout.flush()

    # Phase 3: Random sample for beta_5
    print(f"\n--- Phase 3: beta_5 for {NUM_SAMPLES} random tournaments ---")
    random.seed(99999)
    t_start = time.time()
    b5_count = 0
    b4_count = 0

    for trial in range(NUM_SAMPLES):
        A = random_tournament(n)
        H = ham_path_count(A, n)

        t0 = time.time()
        b5, _ = compute_beta5_only(A, n)
        elapsed = time.time() - t0

        sc = score_sequence(A, n)
        t3 = count_3cycles(A, n)

        if b5 > 0:
            b5_count += 1
            beta5_positive.append((H, [None, b5], t3, sc))
            print(f"    trial {trial+1}: H={H}, t3={t3}, beta_5={b5} <-- FOUND!")

        if (trial + 1) % 25 == 0:
            elapsed_total = time.time() - t_start
            rate = (trial + 1) / elapsed_total
            eta = (NUM_SAMPLES - trial - 1) / rate
            print(f"    {trial+1}/{NUM_SAMPLES} ({rate:.2f}/s, ETA {eta:.0f}s), "
                  f"b5>0: {b5_count}/{trial+1}")
            sys.stdout.flush()

    total_time = time.time() - t_start
    total_computed = len(top_unique) + NUM_SAMPLES

    # ============================================================
    # RESULTS
    # ============================================================
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\nTotal tournaments computed: {total_computed}")

    print(f"\n--- beta_4, beta_5 for top-H tournaments ---")
    print(f"  {'H':>6s}  {'t3':>3s}  {'b4':>3s}  {'b5':>3s}  {'score'}")
    for H, betti_45, t3, sc in all_results:
        print(f"  {H:>6d}  {t3:>3d}  {betti_45[0]:>3d}  {betti_45[1]:>3d}  {sc}")

    print(f"\n--- beta_5 > 0 instances: {len(beta5_positive)} ---")
    if beta5_positive:
        for H, betti, t3, sc in beta5_positive[:30]:
            print(f"  H={H}, t3={t3}, betti_45={betti}, score={sc}")
        print(f"\n  => beta_5 > 0 CONFIRMED at n=9")
        print(f"  => Pattern beta_{{n-4}} > 0 holds for n=6,7,8,9")
    else:
        print("  NONE FOUND in {total_computed} samples")

    # beta_4 distribution for top-H
    print(f"\n--- beta_4 distribution (top-H) ---")
    b4_dist = defaultdict(int)
    for H, betti_45, t3, sc in all_results:
        b4_dist[betti_45[0]] += 1
    for b4 in sorted(b4_dist.keys()):
        print(f"  beta_4 = {b4}: {b4_dist[b4]} tournaments")

    print("\nDone.")


if __name__ == "__main__":
    main()
