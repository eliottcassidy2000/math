#!/usr/bin/env python3
"""
Derive sigma((3,2)) at general n.

sigma((3,2)) = (n-7)! * sum_{ordered (G1 of 4, G2 of 3) disjoint from [n]} H(G1)*dp2(G2)

where H(4-set) = 1+2*c3(4-set), dp2(3-set) = 1+2*cyc(3-set).

Expand: (1+2*c3(G1))*(1+2*cyc(G2)) = 1 + 2*c3(G1) + 2*cyc(G2) + 4*c3(G1)*cyc(G2)

The first three terms are easy. The cross-term c3(G1)*cyc(G2) requires understanding
how 3-cycles in G1 relate to cyclic triples G2 when G1 and G2 are disjoint.

Analysis of cross-term:
  sum_{disjoint (G1,G2)} c3(G1)*cyc(G2)
  = sum_{G2 cyclic} sum_{G1 of 4 from [n]\G2} c3(G1)
  = sum_{G2 cyclic} [(# 3-cycles in [n]\G2) * (n-3-3)]
    ... no, need to expand: sum_{G1 of 4 from S} c3(G1) where S=[n]\G2, |S|=n-3.
  = sum over 3-cycles C1 in S of (|S|-3) = (n-6) * #{3-cycles disjoint from G2}
  So: sum cross = (n-6) * sum_{G2 cyclic} #{3-cycles disjoint from G2}
                = (n-6) * [t3*(t3-1) - #{overlapping pairs}]... no
                = (n-6) * 2*bc  (since each unordered disjoint pair contributes once per G2)

Wait: for each ordered pair (C1, C2) of 3-cycles with C1 != C2, they are either
disjoint or overlapping. If disjoint, they contribute to bc (unordered: bc counts
unordered pairs). For each disjoint pair {C1,C2}, with C2 as G2: C1 is a 3-cycle
in [n]\G2. The extra vertex for G1 can be any of n-6 vertices not in C1 or C2.

So: sum cross = sum_{unordered disjoint pairs (C1,C2)} (n-6) * 2
               (factor 2: C2 can be either role, and G1=C1 union {v}, G2=C2)
  = 2*(n-6)*bc

Actually let me be more careful. The cross-term is:
sum_{G2 cyclic, G1 of 4 disjoint} c3(G1)

For fixed cyclic G2: sum_{G1 of 4 from [n]\G2} c3(G1).
Each 3-cycle in [n]\G2 appears in exactly (n-3-3) = n-6 choices of G1
(the 4th vertex in G1 can be any of the n-6 vertices in [n]\G2 not in the 3-cycle).

So sum_{G1} c3(G1) = (n-6) * t3([n]\G2)
where t3([n]\G2) = #{3-cycles using only vertices in [n]\G2}.

Now: t3([n]\G2) = t3 - #{3-cycles with at least one vertex in G2}
Let G2 = {a,b,c}.
#{3-cycles touching G2} = #{3-cycles with vertex in {a,b,c}}
= sum_v #{3-cycles through v} - sum_{v<w} #{3-cycles through both v and w} + #{3-cycles through a,b,c}
= (d_3(a)+d_3(b)+d_3(c)) - (shared pairs) + cyc(G2)

where d_3(v) = #{3-cycles through vertex v}.

This is getting complicated. Let me just compute numerically.

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
import random

def ham_count_sub(A, verts):
    k = len(verts)
    if k <= 1: return 1
    sub = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    dp = [[0]*k for _ in range(1 << k)]
    for v in range(k):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(k):
                if mask & (1 << u): continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << k) - 1])

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_5_cycles(A, n):
    count = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
                count += 1
    return count // 5

def count_bc(A, n):
    total = 0
    cyc_triples = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            cyc_triples.append(set(triple))
    for i in range(len(cyc_triples)):
        for j in range(i+1, len(cyc_triples)):
            if cyc_triples[i].isdisjoint(cyc_triples[j]):
                total += 1
    return total

# Compute sigma((3,2)) by brute force for various n
for n in [7, 8, 9]:
    print(f"n={n}:")
    free = n - 7  # 4+3 = 7 vertices used

    data = []
    for trial in range(20):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        t3 = count_3_cycles(A, n)
        t5 = count_5_cycles(A, n)
        bc = count_bc(A, n)

        # Compute sigma((3,2)) = (free)! * sum_{ordered (G1 of 4, G2 of 3) disjoint} H(G1)*H(G2)
        total = 0
        for G1 in combinations(range(n), 4):
            remaining = [v for v in range(n) if v not in G1]
            h1 = ham_count_sub(A, list(G1))
            for G2 in combinations(remaining, 3):
                h2 = ham_count_sub(A, list(G2))
                total += h1 * h2
        sigma = factorial(free) * total

        data.append((t3, t5, bc, sigma))
        if trial < 5:
            print(f"  T{trial}: t3={t3}, t5={t5}, bc={bc}, sigma={sigma}")

    # Regression: sigma = a + b*t3 + c*t5 + d*bc
    import numpy as np
    X = np.array([[1, d[0], d[1], d[2]] for d in data])
    y = np.array([d[3] for d in data])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    residuals = y - X @ coeffs
    max_err = np.max(np.abs(residuals))
    print(f"  Regression: sigma = {coeffs[0]:.1f} + {coeffs[1]:.1f}*t3 + {coeffs[2]:.1f}*t5 + {coeffs[3]:.1f}*bc")
    print(f"  Max error: {max_err:.4f}")

    # Check if integer
    c_int = [round(c) for c in coeffs]
    residuals_int = y - X @ np.array(c_int, dtype=float)
    max_err_int = np.max(np.abs(residuals_int))
    print(f"  Integer: sigma = {c_int[0]} + {c_int[1]}*t3 + {c_int[2]}*t5 + {c_int[3]}*bc")
    print(f"  Max error (int): {max_err_int:.4f}")

    # Normalize by (n-7)!
    if free >= 0:
        print(f"  Normalized by {free}!={factorial(free)}:")
        for i, name in enumerate(["const", "t3", "t5", "bc"]):
            print(f"    {name}: {c_int[i]}/{factorial(free)} = {c_int[i]/factorial(free):.4f}")

print("\n" + "="*70)
print("Cross-term analysis for sigma((3,2))")
print("="*70)
print("""
Expected from n=7:
  sigma((3,2)) = 35 + 10*t3 + 8*bc at n=7 (free=0, so no factorial normalization)

The bc term arises from the cross-term c3(G1)*cyc(G2) in the H(G1)*dp2(G2) expansion.

For general n:
  sigma((3,2)) = (n-7)! * [C(n,4)*C(n-4,3) + 2*(n-3)*t3*C(n-4,3)
                            + 2*t3*C(n-3,4) + 4*X]

where X = sum_{disjoint (G1 of 4, G2 of 3)} c3(G1)*cyc(G2).

The question is: how does X depend on t3 and bc?
""")
