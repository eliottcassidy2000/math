#!/usr/bin/env python3
"""
DEEP DIVE: det(I - uA) as predictor of H(T)

Corr(H, |det(I-A/2)|) ≈ 0.98 (!!!)

This suggests: H(T) ≈ f(det(I - A/2)) for some function f.

Let's find f and understand why this works.

det(I - uA) = sum_{k=0}^n (-u)^k * e_k(A)
where e_k(A) is the k-th elementary symmetric function of eigenvalues of A.

But also: det(I - uA) = sum_{k=0}^n (-u)^k * trace of k-th compound matrix of A
= sum_k (-u)^k * sum over k-subsets S of det(A_S)

For a tournament, det(A_S) for subset S counts... directed subgraph structure.

Actually, expanding det(I - uA):
det(I - uA) = 1 - u*tr(A) + u^2 * (sum_{i<j} det(A_{ij})) - ...

tr(A) = 0 for all tournaments (diagonal is 0).
The u^2 term: sum_{i<j} (A[i][i]*A[j][j] - A[i][j]*A[j][i])
= sum_{i<j} (0 - A[i][j]*A[j][i]) = 0 (since exactly one of A[i][j], A[j][i] is 1)
Actually: A[i][j]*A[j][i] = 0 always (tournament property). So u^2 coeff = 0.

Wait: det(minor of I-uA at rows/cols {i,j}) = (1)(1) - (-uA[i][j])(-uA[j][i])
= 1 - u^2 A[i][j]A[j][i] = 1 (since A[i][j]A[j][i]=0 for tournaments).

So the u^2 coefficient is sum of C(n,2) terms, each equal to 1, minus stuff...

Actually let me just compute det(I-uA) symbolically.

For tournament, A + A^T = J - I (all-ones minus identity).
So A^T = J - I - A.

The eigenvalues of A: since A + A^T = J - I, if Av = λv and v ⊥ 1,
then (A + A^T)v = (J-I)v = -v (since Jv = 0 for v ⊥ 1).
So (λ + λ̄)v = -v, meaning Re(λ) = -1/2.

And A*1 = s (score vector), 1^T A = (n-1)*1^T - s^T.

So for REGULAR tournaments (all scores = (n-1)/2):
A*1 = ((n-1)/2)*1, so (n-1)/2 is an eigenvalue.
All other eigenvalues have Re = -1/2.

For general tournaments, the eigenvalues cluster around Re = -1/2
except for the Perron eigenvalue ≈ (n-1)/2.

det(I - uA) = prod(1 - u*λ_i)
For u = 2/(n-1): the Perron factor ≈ 1 - 2/(n-1) * (n-1)/2 = 0.
So det(I - 2A/(n-1)) ≈ 0 for all tournaments!
The DEVIATION from 0 measures how irregular the tournament is.

But we found Corr(H, |det(I-A/2)|) ≈ 0.98, using u=1/2, not u=2/(n-1).
For n=5, 2/(n-1) = 0.5! So they're the SAME!
For n=7, 2/(n-1) = 1/3 ≠ 1/2, but we still had high correlation.

Let me check: is the optimal u = 2/(n-1)?
"""
import numpy as np
import random
from itertools import permutations

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

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

random.seed(42)

# Find optimal u for each n
print("=== Optimal u for det(I - uA) correlation with H ===")
for n in [5, 6, 7]:
    N = 300 if n <= 6 else 100
    data = []
    for _ in range(N):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)
        A_np = np.array(A, dtype=float)
        data.append((H, A_np))

    H_arr = np.array([d[0] for d in data])

    best_corr = 0
    best_u = 0

    for u_100 in range(5, 200):
        u = u_100 / 100.0
        det_arr = np.array([np.linalg.det(np.eye(n) - u * d[1]) for d in data])
        # Use signed det (not abs)
        if np.std(det_arr) > 1e-10:
            corr = abs(np.corrcoef(H_arr, det_arr)[0,1])
            if corr > best_corr:
                best_corr = corr
                best_u = u

    print(f"  n={n}: best_u={best_u:.2f}, best_corr={best_corr:.4f}, "
          f"2/(n-1)={2/(n-1):.4f}")

    # Check specific values
    for u in [2.0/(n-1), 0.5, 1.0/(n-1), 1.0]:
        det_arr = np.array([np.linalg.det(np.eye(n) - u * d[1]) for d in data])
        if np.std(det_arr) > 1e-10:
            corr = abs(np.corrcoef(H_arr, det_arr)[0,1])
        else:
            corr = 0
        print(f"    u={u:.4f}: |Corr|={corr:.4f}")

# ========== Direct relationship H vs det ==========
print("\n\n=== Direct relationship: H vs det(I - 2A/(n-1)) ===")
for n in [5, 7]:
    N = 300 if n <= 5 else 100
    u = 2.0 / (n - 1)
    data = []
    for _ in range(N):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)
        det_val = np.linalg.det(np.eye(n) - u * np.array(A, dtype=float))
        data.append((H, det_val))

    H_arr = np.array([d[0] for d in data])
    det_arr = np.array([d[1] for d in data])

    print(f"\n  n={n}, u=2/(n-1)={u:.4f}:")
    print(f"  Corr(H, det) = {np.corrcoef(H_arr, det_arr)[0,1]:.4f}")

    # Scatter
    from collections import Counter
    pairs = Counter(zip(H_arr.astype(int), np.round(det_arr, 4)))
    for (h, d), count in sorted(pairs.items())[:20]:
        print(f"    H={h:>3}, det={d:>8.4f}: count={count}")

# ========== DEEP: eigenvalue structure ==========
print("\n\n=== Eigenvalue structure of tournament adjacency matrices ===")
n = 5
all_eigs = []
all_H = []

edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = ham_path_count_dp(A, n)
    eigs = np.linalg.eigvals(np.array(A, dtype=float))
    # Sort by real part descending
    eigs_sorted = sorted(eigs, key=lambda x: -x.real)
    all_eigs.append(eigs_sorted)
    all_H.append(H)

# The Perron eigenvalue (largest real part)
perron = np.array([e[0].real for e in all_eigs])
second_eig = np.array([abs(e[1]) for e in all_eigs])
all_H = np.array(all_H)

print(f"  n={n}: {1<<m} tournaments")
print(f"  Corr(H, Perron_eig) = {np.corrcoef(all_H, perron)[0,1]:.4f}")
print(f"  Corr(H, |2nd_eig|) = {np.corrcoef(all_H, second_eig)[0,1]:.4f}")

# The det(I - A/2) at n=5:
# det = prod(1 - λ_i/2)
# The Perron eigenvalue is 2 (for regular tournament at n=5), so
# 1 - 2/2 = 0. For non-regular tournaments, the Perron eigenvalue
# deviates from 2, and det deviates from 0.

# Group by H and show eigenvalues
from collections import defaultdict
by_H = defaultdict(list)
for H, eigs in zip(all_H, all_eigs):
    by_H[H].append(eigs)

print(f"\n  Eigenvalues by H value:")
for H in sorted(by_H.keys()):
    eigs_list = by_H[H]
    perrons = [e[0].real for e in eigs_list]
    print(f"  H={H:>2}: {len(eigs_list):>3} tournaments, "
          f"Perron in [{min(perrons):.2f}, {max(perrons):.2f}]")
