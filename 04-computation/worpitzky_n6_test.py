#!/usr/bin/env python3
"""
TEST KEY PROPERTIES AT n=6 (sampling) AND n=7

Properties to test:
1. Unimodality of F(T,x)
2. Log-concavity of F_k
3. Real-rootedness of F(T,x)
4. F_k(T) = F_{n-1-k}(T^op) complement duality
"""
import random
import numpy as np
from math import comb, factorial
from itertools import permutations

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def complement_tournament(A, n):
    return [[1-A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

def forward_edge_poly_dp(A, n):
    """DP computation of F(T,x) = sum F_k x^k via bitmask DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = [0] * n  # polynomial coefficients indexed by ascent count
        dp[(1 << v, v)][0] = 1  # initially 0 ascents

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            counts = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    new_mask = mask | (1 << u)
                    key = (new_mask, u)
                    if key not in dp:
                        dp[key] = [0] * n

                    is_ascent = 1 if v < u else 0
                    for a in range(n):
                        if counts[a] == 0: continue
                        new_a = a + is_ascent
                        if new_a < n:
                            dp[key][new_a] += counts[a]

    full = (1 << n) - 1
    F = [0] * n
    for v in range(n):
        if (full, v) in dp:
            for a in range(n):
                F[a] += dp[(full, v)][a]
    return F

NUM_SAMPLES = 500

for n in [6, 7, 8]:
    print(f"\n{'='*50}")
    print(f"n={n}: Testing {NUM_SAMPLES} random tournaments")
    print(f"{'='*50}")

    unimodal_count = 0
    lc_count = 0
    real_root_count = 0
    complement_ok = 0
    total = 0
    d = n - 1

    lc_failures = []
    non_real_Fs = []

    for trial in range(NUM_SAMPLES):
        A = random_tournament(n)
        F = forward_edge_poly_dp(A, n)
        H = sum(F)

        # 1. Unimodality
        peak = max(range(n), key=lambda k: F[k])
        is_uni = all(F[k] <= F[k+1] for k in range(peak)) and \
                 all(F[k] >= F[k+1] for k in range(peak, d))
        if is_uni: unimodal_count += 1

        # 2. Log-concavity
        is_lc = True
        for k in range(1, d):
            if F[k]**2 < F[k-1] * F[k+1]:
                is_lc = False
                if len(lc_failures) < 3:
                    lc_failures.append(F)
                break
        if is_lc: lc_count += 1

        # 3. Real-rootedness
        coeffs = list(reversed(F))
        while coeffs and coeffs[0] == 0:
            coeffs.pop(0)
        if len(coeffs) > 1:
            roots = np.roots(coeffs)
            all_real = all(abs(r.imag) < 1e-6 for r in roots)
        else:
            all_real = True
        if all_real: real_root_count += 1
        elif len(non_real_Fs) < 3:
            non_real_Fs.append(F)

        # 4. Complement duality
        A_op = complement_tournament(A, n)
        F_op = forward_edge_poly_dp(A_op, n)
        ok = all(F[k] == F_op[d-k] for k in range(n))
        if ok: complement_ok += 1

        total += 1

        if trial < 3:
            print(f"  Trial {trial}: F={F}, H={H}")

    print(f"\n  Results ({total} samples):")
    print(f"    Unimodal: {unimodal_count}/{total} ({100*unimodal_count/total:.1f}%)")
    print(f"    Log-concave: {lc_count}/{total} ({100*lc_count/total:.1f}%)")
    print(f"    Real roots: {real_root_count}/{total} ({100*real_root_count/total:.1f}%)")
    print(f"    Complement duality: {complement_ok}/{total}")

    if lc_failures:
        print(f"\n  LC failures: {lc_failures}")
    if non_real_Fs:
        print(f"  Non-real root examples: {non_real_Fs}")

print("\n\nDone.")
