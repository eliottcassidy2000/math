#!/usr/bin/env python3
"""
ocr_exact_rational.py — kind-pasteur-2026-03-21-S16

BREAKTHROUGH: OCR(5) = 18/19 and OCR(6) = 12/13 EXACTLY.

These are exact rational numbers. This means OCR is a COMBINATORIAL quantity,
not just a statistical measure. There must be a formula.

18/19: 19 is prime. 18 = 2*9 = 2*3^2.
12/13: 13 is prime. 12 = 4*3 = 2^2*3.

Denominators: 19, 13. Difference = 6 = C(4,2). Sum = 32.
1/19 = 0.0526..., 1/13 = 0.0769...

Let's check: does 1-OCR = 1/(2m-1) where m = C(n,2)?
n=5: m=10, 2m-1=19. 1/19. YES!
n=6: m=15, 2m-1=29. 1/29 = 0.0345. NO (actual is 1/13).

What about related to the number of score classes?
n=5: 9 classes. 18/19 = (2*9)/(2*9+1). Hmm, 18 = 2*9. 19 = 2*9+1.
n=6: 22 classes. 12/13... 2*22 = 44. Doesn't match.

Actually: 1-OCR(5) = 1/19. 1-OCR(6) = 1/13.
Let's think about what these MEAN.

R^2 = [Cov(S2, H)]^2 / [Var(S2) * Var(H)]
    = [sum (S2_i - mean(S2))(H_i - mean(H))]^2 / [sum (S2_i-mean(S2))^2 * sum (H_i-mean(H))^2]

Since c3 = C(n,3) - sum C(s_i,2), and sum C(s_i,2) = (sum s_i^2 - sum s_i)/2
And S2 = sum (s_i - (n-1)/2)^2 = sum s_i^2 - n*(n-1)^2/4

So S2 = 2*sum C(s_i,2) + sum s_i - n*(n-1)^2/4
     = 2*[C(n,3) - c3] + n(n-1)/2 - n(n-1)^2/4

This means S2 = -2*c3 + const. And c3 determines S2 and vice versa (Rao's formula).

So R^2(S2, H) = R^2(c3, H).

And at n<=4: H = 1 + 2*c3 exactly. So R^2 = 1.
At n>=5: H = 1 + 2*c3 + RESIDUAL(5-cycles, 7-cycles, ...).

The key: Var(residual) / Var(H) = 1 - R^2(c3, H).

At n=5: Var(residual)/Var(H) = 1/19.
At n=6: Var(residual)/Var(H) = 1/13.

What is the residual? At n=5: residual = H - (1+2c3) comes from 5-cycle directed counts.
Actually from OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
And 1+2*c3 only counts 3-cycle contributions.

Let's compute the exact variance decomposition.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from math import comb, gcd
from fractions import Fraction
from math import factorial

def build_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_5cycles(A, n):
    """Count all directed 5-cycles (fixing start vertex, counting orientations)."""
    count = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[a]][verts[b]] for b in range(5)] for a in range(5)]
        for perm in permutations(range(5)):
            if all(sub[perm[k]][perm[(k+1)%5]] for k in range(5)):
                count += 1
    return count // 5  # Cyclic rotations

print("=" * 72)
print("  OCR EXACT RATIONAL INVESTIGATION")
print("  kind-pasteur-2026-03-21-S16")
print("=" * 72)

# ============================================================
# Exact computation at n=5 using Fraction arithmetic
# ============================================================

for n in [5, 6]:
    print(f"\n{'='*72}")
    print(f"  EXACT ANALYSIS AT n={n}")
    print(f"{'='*72}")

    m = n*(n-1)//2
    total = 1 << m
    print(f"  m = {m}, total tournaments = {total}")

    all_H = []
    all_S2 = []
    all_c3 = []
    all_c5d = []  # directed 5-cycle count

    for bits in range(total):
        A = build_tournament(n, bits)
        scores = [sum(A[i]) for i in range(n)]
        S2 = sum((s - (n-1)/2)**2 for s in scores)
        c3 = comb(n,3) - sum(s*(s-1)//2 for s in scores)
        H = ham_paths_dp(A, n)

        # Directed 5-cycle count
        c5d = 0
        if n >= 5:
            for verts in combinations(range(n), 5):
                sub = [[A[verts[a]][verts[b]] for b in range(5)] for a in range(5)]
                for perm in permutations(range(5)):
                    if all(sub[perm[k]][perm[(k+1)%5]] for k in range(5)):
                        c5d += 1
            c5d //= 5

        all_H.append(H)
        all_S2.append(S2)
        all_c3.append(c3)
        all_c5d.append(c5d)

    H = np.array(all_H, dtype=np.int64)
    S2 = np.array(all_S2, dtype=np.int64)
    c3 = np.array(all_c3, dtype=np.int64)
    c5d = np.array(all_c5d, dtype=np.int64)

    # Exact sums (use Python integers to avoid overflow)
    N = total
    sum_H = int(np.sum(H))
    sum_S2 = int(np.sum(S2))
    sum_H2 = int(np.sum(H**2))
    sum_S2_2 = int(np.sum(S2**2))
    sum_HS2 = int(np.sum(H * S2))

    print(f"  sum(H) = {sum_H}")
    print(f"  sum(S2) = {sum_S2}")
    print(f"  sum(H^2) = {sum_H2}")
    print(f"  sum(S2^2) = {sum_S2_2}")
    print(f"  sum(H*S2) = {sum_HS2}")

    # R^2 = [N*sum(HS2) - sum(H)*sum(S2)]^2 / [(N*sum(S2^2) - sum(S2)^2)(N*sum(H^2) - sum(H)^2)]
    numer_cov = N * sum_HS2 - sum_H * sum_S2
    var_S2 = N * sum_S2_2 - sum_S2**2
    var_H = N * sum_H2 - sum_H**2

    r2_numer = numer_cov**2
    r2_denom = var_S2 * var_H

    g = gcd(r2_numer, r2_denom)
    p, q = r2_numer // g, r2_denom // g

    print(f"\n  R^2(S2, H) = {p} / {q}")
    print(f"  1 - R^2 = {q - p} / {q}")
    print(f"  OCR = {p/q:.10f}")

    # Mean H
    mean_H = Fraction(sum_H, N)
    print(f"\n  Mean H = {mean_H} = {float(mean_H):.6f}")
    print(f"  Expected random = {n}! / 2^{n-1} = {factorial(n) // (1 << (n-1))}")

    # Var decomposition using c3
    sum_c3 = int(np.sum(c3))
    sum_c3_2 = int(np.sum(c3**2))
    sum_Hc3 = int(np.sum(H * c3))

    # R^2(c3, H)
    numer_c3 = N * sum_Hc3 - sum_H * sum_c3
    var_c3 = N * sum_c3_2 - sum_c3**2

    r2_c3_n = numer_c3**2
    r2_c3_d = var_c3 * var_H
    g2 = gcd(r2_c3_n, r2_c3_d)

    print(f"\n  R^2(c3, H) = {r2_c3_n//g2} / {r2_c3_d//g2}")
    print(f"  (Should equal R^2(S2, H) since S2 = -2c3 + const)")

    # Now: what does the residual correlate with?
    residual = H - (1 + 2*c3)
    sum_res = int(np.sum(residual))
    sum_res2 = int(np.sum(residual**2))
    var_res = N * sum_res2 - sum_res**2

    print(f"\n  Residual = H - (1+2c3):")
    print(f"  sum(res) = {sum_res}")
    print(f"  Var(res)*N^2 = {var_res}")
    print(f"  Var(res)/Var(H) = {Fraction(var_res, var_H)}")

    # Does the residual correlate with c5d (directed 5-cycles)?
    sum_c5d = int(np.sum(c5d))
    sum_c5d2 = int(np.sum(c5d**2))
    sum_res_c5d = int(np.sum(residual * c5d))

    if n == 5:
        print(f"\n  Directed 5-cycle count c5d:")
        print(f"    sum(c5d) = {sum_c5d}")
        print(f"    c5d values: {sorted(set(c5d.tolist()))}")
        print(f"    Correlation(res, c5d):")
        numer_rc5 = N * sum_res_c5d - sum_res * sum_c5d
        var_c5d = N * sum_c5d2 - sum_c5d**2
        if var_c5d > 0 and var_res > 0:
            r2_rc5 = Fraction(numer_rc5**2, var_c5d * var_res)
            print(f"    R^2(residual, c5d) = {r2_rc5} = {float(r2_rc5):.6f}")
        else:
            print(f"    var_c5d = {var_c5d}, var_res = {var_res}")

        # Is residual = 2*c5d exactly?
        check = residual - 2*c5d
        if np.all(check == 0):
            print(f"    EXACT: residual = 2*c5d!")
        else:
            # Check residual = f(c5d, alpha2)
            print(f"    residual - 2*c5d values: {sorted(set(check.tolist()))}")

    # Compute alpha_2 (independent pairs of directed odd cycles)
    # At n=5: odd cycles = 3-cycles + 5-cycles
    # alpha_2 = # independent pairs (disjoint vertex sets)
    print(f"\n  Computing full alpha_1, alpha_2 for {n}...")

    all_alpha1 = []
    all_alpha2 = []
    for idx in range(min(total, total)):
        bits_val = idx
        A_t = build_tournament(n, bits_val)

        # Enumerate all directed odd cycles
        odd_cycles = []

        # 3-cycles
        for i, j, k in combinations(range(n), 3):
            if (A_t[i][j] and A_t[j][k] and A_t[k][i]):
                odd_cycles.append(frozenset({i, j, k}))
            elif (A_t[i][k] and A_t[k][j] and A_t[j][i]):
                odd_cycles.append(frozenset({i, j, k}))

        # 5-cycles (for n>=5): each vertex set can have multiple directed cycles
        if n >= 5:
            for verts in combinations(range(n), 5):
                sub = [[A_t[verts[a]][verts[b]] for b in range(5)] for a in range(5)]
                # Count Hamiltonian directed cycles (fixing vertex 0 as start)
                for perm in permutations(range(1, 5)):
                    cycle = (0,) + perm
                    if all(sub[cycle[k]][cycle[(k+1)%5]] for k in range(5)):
                        odd_cycles.append(frozenset(verts))
                        # Note: same frozenset may appear multiple times
                        # That's correct: each directed cycle is a separate vertex of Omega

        alpha1 = len(odd_cycles)

        # Count independent pairs (vertex-disjoint)
        alpha2 = 0
        for a in range(len(odd_cycles)):
            for b in range(a+1, len(odd_cycles)):
                if len(odd_cycles[a] & odd_cycles[b]) == 0:
                    alpha2 += 1

        all_alpha1.append(alpha1)
        all_alpha2.append(alpha2)

    alpha1 = np.array(all_alpha1, dtype=np.int64)
    alpha2 = np.array(all_alpha2, dtype=np.int64)

    # Check OCF: H = 1 + 2*alpha1 + 4*alpha2
    H_ocf = 1 + 2*alpha1 + 4*alpha2
    mismatches = np.sum(H_ocf != H)
    if mismatches == 0:
        print(f"  OCF H = 1 + 2*alpha1 + 4*alpha2 EXACT at n={n}! (0 mismatches)")
    else:
        # Need alpha3 too
        max_diff = np.max(np.abs(H - H_ocf))
        print(f"  OCF needs higher alpha: {mismatches}/{total} mismatches, max diff={max_diff}")

    # Now: residual = H - (1+2c3) = 2*(alpha1 - c3) + 4*alpha2 + ...
    # alpha1 - c3 = number of directed 5-cycles (at n=5)
    delta_alpha = alpha1 - c3
    print(f"  alpha1 - c3 (= directed 5-cycles at n={n}):")
    print(f"    values: {sorted(set(delta_alpha.tolist()))}")

    res_exact = residual
    # Does residual = 2*(alpha1-c3) + 4*alpha2?
    res_check = 2*delta_alpha + 4*alpha2
    if np.all(res_exact == res_check):
        print(f"  CONFIRMED: residual = 2*(alpha1-c3) + 4*alpha2 EXACTLY")
    else:
        diff = res_exact - res_check
        print(f"  residual - [2*(alpha1-c3) + 4*alpha2] = {sorted(set(diff.tolist()))}")

    # Final: Var(H)/Var(c3) decomposition
    print(f"\n  FINAL DECOMPOSITION:")
    print(f"  H = 1 + 2*c3 + 2*(alpha1-c3) + 4*alpha2 + ...")
    print(f"      = [score-determined] + [5-cycle excess] + [disjoint pairs] + ...")
    print(f"  R^2(c3, H) = {p}/{q} = {p/q:.10f}")
    print(f"  1 - R^2 = {q-p}/{q} = {(q-p)/q:.10f}")

from math import factorial

print("\n" + "=" * 72)
print("  SUMMARY: THE EXACT OCR SEQUENCE")
print("=" * 72)
print()
print("  n=3: OCR = 1 = 1/1")
print("  n=4: OCR = 1 = 1/1")
print("  n=5: OCR = 18/19")
print("  n=6: OCR = 12/13")
print("  n=7: OCR = ? (need exact computation)")
print("  n=8: OCR = ? (need exact or high-precision sampling)")
print()
print("  Denominators: 1, 1, 19, 13, ?, ?")
print("  Numerators:   1, 1, 18, 12, ?, ?")
print()
print("  19 = 2*C(5,2) - 1 = 2*10 - 1")
print("  13 = C(6,2) - 2 = 15 - 2")
print("  Or: 19 = 4*5-1, 13 = 2*7-1")
print("  Or: 1-OCR = 1/19, 1/13")
print("  Differences of denominators: 13-19 = -6 = -C(4,2)")
