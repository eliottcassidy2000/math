#!/usr/bin/env python3
"""
Compute H for representative Z_19 circulant tournaments.
2^19 * 19 ~ 10M states per tournament — feasible with optimized DP.
"""

import numpy as np
from itertools import combinations
from math import comb
import time

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

def adjacency_matrix(S, p):
    A = np.zeros((p,p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in S:
                A[i][j] = 1
    return A

def count_hp_fast(A):
    """Optimized Held-Karp for n up to ~20."""
    n = len(A)
    full = (1 << n) - 1
    # Use numpy array for DP
    # dp[mask] is array of length n: dp[mask][v] = #paths ending at v visiting mask
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            # Find neighbors
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c

    return int(np.sum(dp[full]))

def eigenvalues_circulant(S, p):
    omega = np.exp(2j*np.pi/p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

p = 19
m = (p-1)//2  # 9

# Get Paley and a few representative connection sets
paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))
print(f"Paley QR_19 = {sorted(paley_S)}")

# Generate connection sets of varying "distance" from Paley
# Strategy: swap one element at a time
all_circulants = []
for S in combinations(range(1,p), m):
    S_set = set(S)
    if all((p-j)%p not in S_set for j in S_set):
        all_circulants.append(frozenset(S_set))

print(f"Total: {len(all_circulants)} circulant tournaments")

# Compute eigenvalue spectra for all
spec_data = []
for S in all_circulants:
    eigs = eigenvalues_circulant(S, p)
    y_half = [eigs[k].imag for k in range(1, m+1)]
    x = [y**2 for y in y_half]
    e2 = sum(x[i]*x[j] for i in range(m) for j in range(i+1, m))
    spec_data.append({'S': S, 'x': x, 'e2': e2, 'paley': S == paley_S})

# Sort by e2 to get diverse representatives
spec_data.sort(key=lambda d: d['e2'])
paley_idx = next(i for i, d in enumerate(spec_data) if d['paley'])
print(f"Paley position in e2 ranking: {paley_idx+1}/{len(spec_data)} (highest)")

# Pick representatives: Paley, lowest e2, and a few in between
indices = [0, len(spec_data)//4, len(spec_data)//2, 3*len(spec_data)//4, paley_idx]
# Remove duplicates and sort
indices = sorted(set(indices))

print(f"\nComputing H for {len(indices)} representative tournaments...")
print(f"(Each takes ~10-30 seconds for p=19)")

results = []
for idx in indices:
    d = spec_data[idx]
    S = d['S']
    A = adjacency_matrix(S, p)

    t0 = time.time()
    H = count_hp_fast(A)
    t1 = time.time()

    # Compute more e_k
    x = d['x']
    ek = {0: 1.0}
    for k in range(1, min(m+1, 6)):
        ek[k] = sum(np.prod(list(combo)) for combo in combinations(x, k))

    results.append({'S': S, 'H': H, 'e': ek, 'x': x,
                    'paley': d['paley'], 'time': t1-t0})

    print(f"  S={sorted(S)}: H={H}, e2={ek.get(2,0):.2f}, "
          f"e3={ek.get(3,0):.2f}, time={t1-t0:.1f}s"
          f" {'(PALEY)' if d['paley'] else ''}")

# Sort by H descending
results.sort(key=lambda r: r['H'], reverse=True)

print(f"\n{'='*60}")
print(f"RESULTS FOR p = 19")
print(f"{'='*60}")

for r in results:
    print(f"  H={r['H']:>12}, e2={r['e'].get(2,0):>12.2f}, "
          f"e3={r['e'].get(3,0):>12.2f}, e4={r['e'].get(4,0):>12.2f}"
          f" {'(PALEY)' if r['paley'] else ''}")

# Check: does Paley maximize H?
paley_H = next(r['H'] for r in results if r['paley'])
max_H = max(r['H'] for r in results)
print(f"\nPaley H = {paley_H}")
print(f"Max H = {max_H}")
print(f"Paley maximizes H among tested? {paley_H == max_H}")

# If we have enough distinct H values, try fitting
distinct_H = {}
for r in results:
    if r['H'] not in distinct_H:
        distinct_H[r['H']] = r

n_distinct = len(distinct_H)
print(f"\n{n_distinct} distinct H values among {len(results)} tested")

if n_distinct >= 3:
    reps = list(distinct_H.values())
    reps.sort(key=lambda r: r['H'], reverse=True)

    # Try fitting with as many e_k as we have distinct values
    n_e = min(n_distinct - 1, 5)
    A_mat = np.zeros((n_distinct, n_e + 1))
    b_vec = np.array([r['H'] for r in reps], dtype=float)

    for i, r in enumerate(reps):
        A_mat[i, 0] = 1
        for j in range(1, n_e + 1):
            A_mat[i, j] = r['e'].get(j + 1, 0)

    rank = np.linalg.matrix_rank(A_mat)
    print(f"Fitting with {n_e} e_k terms, rank={rank}")

    if rank == n_e + 1 and n_distinct == n_e + 1:
        coeffs = np.linalg.solve(A_mat, b_vec)
        residual = np.max(np.abs(A_mat @ coeffs - b_vec))
        print(f"H = {coeffs[0]:.2f}", end="")
        for j in range(1, n_e+1):
            print(f" + ({coeffs[j]:.4f})*e_{j+1}", end="")
        print(f"  [residual: {residual:.2e}]")

        # Hessian
        v = p/4
        lambda_H = sum(-coeffs[j] * comb(m-2, j) * v**(j-1) for j in range(1, n_e+1) if j <= m-2+1)
        # More carefully:
        lambda_H = 0
        for j in range(1, n_e+1):
            k = j + 1
            if k - 2 <= m - 2:
                curv = -comb(m-2, k-2) * v**(k-2)
                lambda_H += coeffs[j] * curv
        print(f"Hessian lambda_H = {lambda_H:.4f}")
        print(f"Paley is {'LOCAL MAX' if lambda_H < 0 else 'NOT local max'}")
