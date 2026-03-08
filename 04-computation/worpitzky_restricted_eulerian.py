#!/usr/bin/env python3
"""
F(T,x) AS RESTRICTED EULERIAN POLYNOMIAL — Final analysis

PROVED: Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)}
Reason: Each permutation sigma is a HP of T iff A[sigma_i][sigma_{i+1}]=1 for all i.
The (n-1) path edges are fixed, the remaining C(n,2)-(n-1) edges are free.
So each permutation appears in exactly 2^{C(n,2)-(n-1)} tournaments.
The number of permutations with k ascents = A(n,k), so
Sum_T F_k(T) = A(n,k) * 2^{C(n,2)-(n-1)}.

CONSEQUENCE: Average H = n! / 2^{n-1}

Now investigate:
1. VARIANCE structure: Var(F_k) across tournaments
2. COVARIANCE: Cov(F_k, F_j)
3. Deviation F(T,x) - (H/n!) * A_n(x): what does the residual look like?
4. Connection to palindromic residuals and OCF
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np
from collections import defaultdict

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[fwd] += 1
    return F

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

random.seed(42)

print("=" * 70)
print("F(T,x) - (H/n!) * A_n(x): RESIDUAL ANALYSIS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")

    E = [eulerian_number(n, k) for k in range(n)]
    nfact = factorial(n)

    all_residuals = []
    all_F = []
    all_H = []

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        all_F.append(F)
        all_H.append(H)

        # Residual: F_k - (H/n!) * A(n,k)
        residual = [F[k] - H * E[k] / nfact for k in range(n)]
        all_residuals.append(residual)

    total = len(all_F)
    all_residuals = np.array(all_residuals)
    all_F_arr = np.array(all_F)

    print(f"  Eulerian numbers A({n},k) = {E}")
    print(f"  n! = {nfact}")
    print(f"  Mean residual: {np.mean(all_residuals, axis=0)}")  # should be 0

    # Variance of each residual component
    var_residual = np.var(all_residuals, axis=0)
    var_F = np.var(all_F_arr, axis=0)
    print(f"  Var(F_k):       {var_F}")
    print(f"  Var(residual_k): {var_residual}")
    print(f"  Fraction of variance in residual: {var_residual / (var_F + 1e-10)}")

    # Is the residual palindromic?
    pal_count = 0
    for r in all_residuals:
        if all(abs(r[k] - r[n-1-k]) < 1e-8 for k in range(n)):
            pal_count += 1
    print(f"  Residual palindromic: {pal_count}/{total}")

    # Sum of residual should be 0 (since sum E_k / n! = 1)
    sum_residual = [sum(r) for r in all_residuals]
    print(f"  All residual sums zero: {all(abs(s) < 1e-8 for s in sum_residual)}")

    # Covariance matrix of F_k
    cov_F = np.cov(all_F_arr.T)
    print(f"\n  Covariance matrix of F:")
    print(cov_F.round(4))

    # Correlation with H
    H_arr = np.array(all_H)
    for k in range(n):
        corr = np.corrcoef(all_F_arr[:, k], H_arr)[0, 1]
        print(f"  Corr(F_{k}, H) = {corr:.4f}")

# ===== THE BIG PICTURE =====
print("\n\n" + "=" * 70)
print("THE BIG PICTURE: F(T,x) = (H/n!) * A_n(x) + residual")
print("=" * 70)
print("""
THEOREM: For any tournament T on n vertices:
  F(T,x) = (H(T)/n!) * A_n(x) + R(T,x)

where:
  - A_n(x) = Eulerian polynomial (descent polynomial of S_n)
  - H(T) = total HP count = F(T,1)
  - R(T,x) is a PALINDROMIC polynomial with R(T,1) = 0

Proof: Both F(T,x) and A_n(x) are palindromic polynomials.
  At x=1: F(T,1) = H, A_n(1) = n!, so R(T,1) = H - H = 0.
  Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)}, so Sum_T R(T,x) = 0.

QUESTION: What determines R(T,x)?
Since F is palindromic, so is H * A_n(x) / n! (since A_n is palindromic).
Therefore R(T,x) is also palindromic.
R has degree n-1 and R(1) = 0, so (x-1) divides R.
Since R is palindromic of even degree (n-1 is even at odd n):
  R(x) = (x-1) * Q(x) where Q has degree n-2.
  Palindromicity: R(x) = x^{n-1} R(1/x)
  So (x-1) Q(x) = x^{n-1} (1/x - 1) Q(1/x) = -x^{n-2}(x-1) Q(1/x)
  Hence Q(x) = -x^{n-2} Q(1/x), i.e., Q is ANTI-palindromic.
""")

# Verify that R is palindromic
for n in [5]:
    E = [eulerian_number(n, k) for k in range(n)]
    nfact = factorial(n)

    pal_count = 0
    total = 0
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        R = [F[k] - H * E[k] / nfact for k in range(n)]

        is_pal = all(abs(R[k] - R[n-1-k]) < 1e-8 for k in range(n))
        if is_pal: pal_count += 1
        total += 1

        # Check if (x-1) divides R
        # R(1) = 0 (we know this). What about R(-1)?
        R_at_m1 = sum(R[k] * (-1)**k for k in range(n))

        if total <= 5:
            print(f"  F={[int(x) for x in F]}, H={H}, R={[round(x,4) for x in R]}")
            print(f"    palindromic={is_pal}, R(1)={sum(R):.6f}, R(-1)={R_at_m1:.6f}")

    print(f"\n  Palindromic residuals: {pal_count}/{total}")

# What are the residuals in terms of cycle counts?
print("\n\n" + "=" * 70)
print("RESIDUAL R(T,x) vs CYCLE COUNTS")
print("=" * 70)

for n in [5]:
    E = [eulerian_number(n, k) for k in range(n)]
    nfact = factorial(n)

    data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        R = [F[k] - H * E[k] / nfact for k in range(n)]

        # Count 3-cycles
        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        data.append((H, t3, scores, R))

    # Check if R is determined by (H, t3)
    by_ht3 = defaultdict(list)
    for H, t3, scores, R in data:
        by_ht3[(H, t3)].append(tuple(round(x, 8) for x in R))

    print(f"\nn={n}: Is R determined by (H, t3)?")
    ambiguous = 0
    for key, Rs in by_ht3.items():
        unique = set(Rs)
        if len(unique) > 1:
            ambiguous += 1
            if ambiguous <= 3:
                print(f"  {key}: {len(unique)} distinct R values")

    print(f"  Ambiguous: {ambiguous}/{len(by_ht3)}")

    # Check if R is determined by (H, t3, scores)
    by_hts = defaultdict(list)
    for H, t3, scores, R in data:
        by_hts[(H, t3, scores)].append(tuple(round(x, 8) for x in R))

    ambiguous2 = 0
    for key, Rs in by_hts.items():
        unique = set(Rs)
        if len(unique) > 1:
            ambiguous2 += 1
    print(f"  With scores: Ambiguous: {ambiguous2}/{len(by_hts)}")

    # Is R proportional to some simple polynomial?
    # R is palindromic with R(1)=0. At n=5, R has degree 4 and 2 free coefficients.
    # R = a*(x^4 - 1) + b*(x^3 - x)? Let's check.
    print(f"\n  Testing R = a*(x^4-1) + b*(x^3-x):")
    for H, t3, scores, R in data[:10]:
        # R_0 = -a, R_1 = -b, R_2 = 0, R_3 = b, R_4 = a
        a_est = R[4]
        b_est = R[3]
        R_pred = [-a_est, -b_est, 0, b_est, a_est]
        match = all(abs(R[k] - R_pred[k]) < 1e-6 for k in range(n))
        if not match:
            print(f"    H={H}, t3={t3}: R={[round(x,4) for x in R]}, pred={[round(x,4) for x in R_pred]}, MISMATCH")

    # So R_2 might not be 0. Let me check.
    R2_values = [R[2] for _, _, _, R in data]
    print(f"\n  R_2 values: min={min(R2_values):.4f}, max={max(R2_values):.4f}")
    print(f"  R_2 = 0 for all? {all(abs(x) < 1e-8 for x in R2_values)}")
