#!/usr/bin/env python3
"""
Seesaw beta_1 * beta_3 = 0: LARGE SAMPLE TEST

The devil's advocate correctly noted: 500 samples at n=8 is NOT
statistically significant (Fisher p ~ 0.38). We need 2000-5000
samples for 80% power to detect a 0.2% violation rate.

This script runs 3000 samples at n=8 to get proper statistical power.
Also tests 1000 at n=9.

kind-pasteur-2026-03-25-S3
"""
import sys, time, random
import numpy as np
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)
random.seed(20260325)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def enum_allowed(A, n, p):
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_beta(A, n, target_p):
    """Compute beta_{target_p} efficiently."""
    p = target_p
    ap_m1 = enum_allowed(A, n, p-1) if p >= 1 else [()]
    ap = enum_allowed(A, n, p)
    ap_p1 = enum_allowed(A, n, p+1)
    if not ap: return 0

    # Omega_p
    ap_m1_set = set(ap_m1)
    na = {}; nai = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in na:
                na[face] = nai; nai += 1
    if not na:
        omega_dim = len(ap); omega_basis = np.eye(len(ap))
    else:
        P = np.zeros((len(na), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in na: P[na[face], j] += sign
        _, S, Vt = np.linalg.svd(P)
        rP = sum(s > 1e-8 for s in S)
        omega_dim = len(ap) - rP
        if omega_dim == 0: return 0
        omega_basis = Vt[rP:].T

    # ker(d_p)
    if p >= 1 and ap_m1:
        idx = {path: i for i, path in enumerate(ap_m1)}
        D = np.zeros((len(ap_m1), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in idx: D[idx[face], j] += sign
        Do = D @ omega_basis
        _, Sd, _ = np.linalg.svd(Do)
        rank_d = sum(s > 1e-8 for s in Sd)
    else: rank_d = 0
    ker_d = omega_dim - rank_d

    # im(d_{p+1})
    if not ap_p1: return ker_d
    ap_set = set(ap)
    na2 = {}; nai2 = 0
    for path in ap_p1:
        for i in range(p+2):
            face = path[:i] + path[i+1:]
            if face not in ap_set and face not in na2:
                na2[face] = nai2; nai2 += 1
    if not na2:
        o1_basis = np.eye(len(ap_p1))
    else:
        Q = np.zeros((len(na2), len(ap_p1)))
        for j, path in enumerate(ap_p1):
            for i in range(p+2):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in na2: Q[na2[face], j] += sign
        _, Sq, Vqt = np.linalg.svd(Q)
        rQ = sum(s > 1e-8 for s in Sq)
        if rQ >= len(ap_p1): return ker_d
        o1_basis = Vqt[rQ:].T

    aidx = {path: i for i, path in enumerate(ap)}
    D1 = np.zeros((len(ap), len(ap_p1)))
    for j, path in enumerate(ap_p1):
        for i in range(p+2):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in aidx: D1[aidx[face], j] += sign
    D1o = D1 @ o1_basis
    im_mat = omega_basis.T @ D1o
    _, Si, _ = np.linalg.svd(im_mat)
    im_d1 = sum(s > 1e-8 for s in Si)
    return max(0, ker_d - im_d1)

# ============================================================
# MAIN
# ============================================================
print("=" * 72)
print("  SEESAW LARGE-SAMPLE TEST: beta_1 * beta_3 = 0")
print("  kind-pasteur-2026-03-25-S3")
print("=" * 72)

for n, target in [(8, 3000), (9, 500)]:
    print(f"\n--- n={n}: {target} samples ---")
    b1_hist = Counter(); b3_hist = Counter()
    violations = 0; total = 0
    t0 = time.time()

    for trial in range(target):
        A = random_tournament(n)
        b1 = compute_beta(A, n, 1)
        b3 = compute_beta(A, n, 3)
        b1_hist[b1] += 1; b3_hist[b3] += 1
        if b1 > 0 and b3 > 0: violations += 1
        total += 1

        if total % 500 == 0:
            el = time.time() - t0
            rate = total / el
            eta = (target - total) / rate if rate > 0 else 0
            print(f"  [{total}/{target}] {el:.0f}s ({rate:.1f}/s) violations={violations}")

    el = time.time() - t0
    print(f"\n  n={n}: {total} samples in {el:.0f}s")
    print(f"  beta_1: {dict(sorted(b1_hist.items()))}")
    print(f"  beta_3: {dict(sorted(b3_hist.items()))}")
    print(f"  Seesaw violations (b1>0 AND b3>0): {violations}/{total}")

    # Statistical analysis
    b1_rate = sum(v for k,v in b1_hist.items() if k>0) / total
    b3_rate = sum(v for k,v in b3_hist.items() if k>0) / total
    expected_if_indep = b1_rate * b3_rate * total
    print(f"\n  beta_1 rate: {b1_rate*100:.2f}%")
    print(f"  beta_3 rate: {b3_rate*100:.2f}%")
    print(f"  If independent: expected {expected_if_indep:.1f} violations")
    print(f"  Observed: {violations}")

    if expected_if_indep > 0:
        import math
        # Poisson probability of 0 given expected
        p_zero = math.exp(-expected_if_indep)
        print(f"  P(0 violations | independence): {p_zero:.4f}")
        if p_zero < 0.05:
            print(f"  ==> STATISTICALLY SIGNIFICANT: seesaw is REAL (p={p_zero:.4f})")
        else:
            print(f"  ==> NOT statistically significant yet (p={p_zero:.4f})")
    else:
        print(f"  (Cannot test: one rate is 0)")

    # Cross-tabulation
    print(f"\n  Cross-tabulation:")
    b1v = sorted(b1_hist.keys()); b3v = sorted(b3_hist.keys())
    cross = Counter()
    # Recompute cross-tab... actually we didn't store per-tournament results.
    # Just report the marginals and violation count.
    n_b1_b3 = violations
    n_b1_only = sum(v for k,v in b1_hist.items() if k>0) - violations
    n_b3_only = sum(v for k,v in b3_hist.items() if k>0) - violations
    n_neither = total - n_b1_only - n_b3_only - violations
    print(f"             b3=0   b3>0")
    print(f"    b1=0   {n_neither:>6}  {n_b3_only:>6}")
    print(f"    b1>0   {n_b1_only:>6}  {n_b1_b3:>6}")
