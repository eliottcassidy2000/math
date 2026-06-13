#!/usr/bin/env python3
"""
PERPENDICULARITY MECHANISM

Key findings:
1. n=6 (even): c_0 ≡ 0 identically. Perfect "perpendicularity" (trivially).
2. n=5 (odd): Corr(c_0, H) ≈ 0.24 — NOT perpendicular.
3. n=7 (odd): Corr(c_0, H) ≈ -0.03 — NEARLY perpendicular.

The trend is clear: as n grows, the alternating sum cancellation improves.

WHY? The covariance terms (-1)^k Cov(D_k, D_0) have palindromic symmetry
(term k = -(term n-1-k) for even n-1, but not for odd n-1).

Actually for even n (like n=6): n-1=5 is odd.
  (-1)^k * Cov(D_k, D_0) paired with (-1)^{n-1-k} * Cov(D_{n-1-k}, D_0)
  = (-1)^{n-1-k} * Cov(D_k, D_0) (by palindromicity of D)
  = (-1)^{n-1} * (-1)^{-k} * Cov = -(-1)^{-k} * Cov (when n-1 odd)
  = -(-1)^k * Cov

So the paired terms EXACTLY CANCEL when n is even! This proves c_0 ≡ 0.

For odd n (like n=7): n-1=6 is even.
  (-1)^{n-1-k} = (-1)^{-k} * (-1)^{n-1} = (-1)^{-k} * 1 = (-1)^k
  So the paired terms are EQUAL, not opposite. No cancellation from palindromicity.

But the NEAR-cancellation at n=7 must come from something else.

HYPOTHESIS: The near-cancellation comes from the BINOMIAL STRUCTURE of D_k.
For a random tournament, fwd(P) ≈ Binomial(n-1, 1/2), so
D_k ≈ n! * C(n-1,k)/2^{n-1} + fluctuations.

The fluctuations come from edge correlations between permutations.
As n grows, the ratio (fluctuation)/(mean) shrinks,
and the alternating sum of the fluctuations cancels more and more.

Let me quantify this.
"""
from itertools import permutations
import numpy as np
import random
from math import comb, factorial

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_F(A, n):
    D = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        D[fwd] += 1
    return D

random.seed(42)

# Compute Corr(c_0, H) vs n
print("=== Perpendicularity vs n ===")
print(f"{'n':>3} | {'Corr(c0,H)':>12} | {'Var(c0)':>10} | {'Var(H)':>10} | {'Var(c0)/Var(H)':>14}")
print("-" * 60)

for n in range(3, 8):
    N_samples = 500 if n <= 5 else (300 if n <= 6 else 100)

    c0_list = []
    H_list = []
    for _ in range(N_samples):
        A = random_tournament(n)
        D = compute_F(A, n)
        H = D[n-1]
        c0 = sum((-1)**k * D[k] for k in range(n)) / 2**(n-1)
        c0_list.append(c0)
        H_list.append(H)

    c0_arr = np.array(c0_list)
    H_arr = np.array(H_list)

    if np.var(c0_arr) < 1e-10:
        corr = 0.0
    else:
        corr = np.corrcoef(c0_arr, H_arr)[0,1]

    print(f"{n:>3} | {corr:>12.4f} | {np.var(c0_arr):>10.4f} | {np.var(H_arr):>10.2f} | {np.var(c0_arr)/max(np.var(H_arr),1e-10):>14.6f}")

# ========== EXACT ANALYSIS for small n ==========
print("\n\n=== EXACT analysis: all tournaments ===")
for n in [3, 4, 5]:
    all_T = []
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)

    c0_list = []
    H_list = []

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1

        D = compute_F(A, n)
        H = D[n-1]
        c0 = sum((-1)**k * D[k] for k in range(n)) / 2**(n-1)
        c0_list.append(c0)
        H_list.append(H)

    c0_arr = np.array(c0_list)
    H_arr = np.array(H_list)

    # Exact covariance
    cov = np.mean(c0_arr * H_arr) - np.mean(c0_arr) * np.mean(H_arr)
    var_c0 = np.var(c0_arr)
    var_H = np.var(H_arr)

    if var_c0 < 1e-10 or var_H < 1e-10:
        corr_str = "N/A (zero variance)"
    else:
        corr_str = f"{cov / np.sqrt(var_c0 * var_H):.6f}"

    print(f"  n={n}: Cov(c0,H)={cov:.6f}, Var(c0)={var_c0:.6f}, Corr={corr_str}")
    print(f"    c0 values: {sorted(set(c0_arr))}")
    print(f"    H values: {sorted(set(H_arr))}")

# ========== HYPOTHESIS J: c_0 depends ONLY on odd cycle counts ==========
print("\n\n=== HYPOTHESIS J: c_0 = f(t3, t5, t7, ...) ===")
# c_0 = S/2^{n-1} = signed permanent / 2^{n-1}
# We know S = sum_P prod(2A-1). This is a MULTILINEAR function of edges.
# Each edge appears in many permutations.
# S depends on the tournament's structure in a specific way.
#
# From earlier work: at n=5, c_0 takes values {-1, 0, 1} ∪ {±1/2, ±3/2}
# The signed permanent is related to cycle structure via Redei theory.

# At n=5 (exact), what does c_0 depend on?
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

def count_3cycles(A):
    n = len(A)
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

data = []
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    D = compute_F(A, n)
    H = D[n-1]
    c0 = sum((-1)**k * D[k] for k in range(n)) / 2**(n-1)
    t3 = count_3cycles(A)
    data.append({'H': H, 'c0': c0, 't3': t3})

# Group by (H, t3, c0)
from collections import Counter
triples = Counter((d['H'], d['t3'], d['c0']) for d in data)
print(f"  n=5: (H, t3, c0) distribution:")
for (H, t3, c0), count in sorted(triples.items()):
    print(f"    H={H:>2}, t3={t3:>2}, c0={c0:>5.2f}: count={count:>3}")
