#!/usr/bin/env python3
"""
p17_walsh_decomposition.py — Compute Walsh degree decomposition at p=17

This completes the degree dominance ratio data:
  p=7:  0.00  (Paley)
  p=11: 0.19  (Paley)
  p=13: 4.42  (Interval)
  p=17: ???   (Interval expected)

Uses optimized Held-Karp DP for all 256 circulant orientations.

Author: opus-2026-03-12-S67
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import time

def count_H_fast(sigma, p):
    """Optimized Held-Karp DP for circulant tournament."""
    m = (p-1)//2
    n = p
    # Build adjacency as a list of neighbor sets for each vertex
    adj_out = [[] for _ in range(n)]
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                adj_out[i].append(j)
            else:
                adj_out[j].append(i)

    # Convert to adjacency matrix for DP
    adj = [0] * n  # bitmask adjacency
    for i in range(n):
        for j in adj_out[i]:
            adj[i] |= (1 << j)

    # DP: dp[mask][v] = # Ham paths ending at v using vertex set mask
    full = (1 << n) - 1
    # Use array instead of dict for speed
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        bc = bin(mask).count('1')
        if bc >= n:
            continue
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            # Extend to each neighbor not in mask
            neighbors = adj[v] & ~mask
            while neighbors:
                u = (neighbors & -neighbors).bit_length() - 1
                dp[mask | (1 << u)][u] += dp[mask][v]
                neighbors &= neighbors - 1

    return sum(dp[full][v] for v in range(n))

def walsh_char(S, sigma):
    prod = 1
    for k in S:
        prod *= sigma[k-1]
    return prod

p = 17
m = (p-1)//2
N = 1 << m
print(f"p={p}, m={m}, {N} orientations")
print(f"Computing H for all {N} orientations...")

t0 = time.time()
pairs = []
for bits in range(N):
    sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
    H = count_H_fast(sigma, p)
    pairs.append((sigma, H))
    if (bits + 1) % 32 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{N} done ({elapsed:.1f}s)")

elapsed = time.time() - t0
print(f"All {N} done in {elapsed:.1f}s")

# Walsh decomposition
print(f"\nComputing Walsh decomposition...")
chords = list(range(1, m+1))

walsh_coeffs = {}
for deg in range(0, m+1, 2):
    for S in combinations(chords, deg):
        val = sum(H * walsh_char(S, sigma) for sigma, H in pairs) / N
        if abs(val) > 1e-6:
            walsh_coeffs[S] = val

var_by_deg = defaultdict(float)
for S, val in walsh_coeffs.items():
    var_by_deg[len(S)] += val**2

total_var = sum(v for d, v in var_by_deg.items() if d > 0)

print(f"\nWalsh variance decomposition at p={p}:")
for deg in sorted(var_by_deg):
    if deg == 0:
        continue
    frac = var_by_deg[deg] / total_var * 100
    print(f"  Degree {deg}: variance = {var_by_deg[deg]:.2f} ({frac:.2f}%)")

ratio_4_2 = var_by_deg.get(4, 0) / var_by_deg.get(2, 1)
print(f"\n  deg-4/deg-2 ratio: {ratio_4_2:.4f}")

if 6 in var_by_deg:
    ratio_6_2 = var_by_deg[6] / var_by_deg[2]
    print(f"  deg-6/deg-2 ratio: {ratio_6_2:.4f}")

if 8 in var_by_deg:
    ratio_8_2 = var_by_deg[8] / var_by_deg[2]
    print(f"  deg-8/deg-2 ratio: {ratio_8_2:.4f}")

# Verify Interval is maximum
sigma_int = tuple([1]*m)
H_int = dict(pairs)[sigma_int]
H_max = max(H for _, H in pairs)
print(f"\nH(Interval) = {H_int}")
print(f"H(max)      = {H_max}")
print(f"Interval is max? {H_int == H_max}")

# Summary table
print(f"\n{'='*70}")
print("COMPLETE DEGREE DOMINANCE DATA")
print(f"{'='*70}")
all_data = {
    7: {'ratio': 0.0000, 'winner': 'Paley', 'deg2': 100.0, 'deg4': 0.0},
    11: {'ratio': 0.1885, 'winner': 'Paley', 'deg2': 84.1, 'deg4': 15.9},
    13: {'ratio': 4.4248, 'winner': 'Interval', 'deg2': 18.4, 'deg4': 81.6},
    17: {'ratio': ratio_4_2, 'winner': 'Interval' if H_int == H_max else 'Other',
         'deg2': var_by_deg.get(2, 0)/total_var*100,
         'deg4': var_by_deg.get(4, 0)/total_var*100}
}

print(f"\n  {'p':>3} | {'deg-2%':>7} | {'deg-4%':>7} | {'ratio 4/2':>10} | {'Winner':>10}")
print(f"  " + "-" * 55)
for p_val in sorted(all_data):
    d = all_data[p_val]
    print(f"  {p_val:>3} | {d['deg2']:>6.1f}% | {d['deg4']:>6.1f}% | "
          f"{d['ratio']:>10.4f} | {d['winner']:>10}")

# Trend analysis
ratios = [all_data[p_val]['ratio'] for p_val in [11, 13, 17]]
print(f"\n  Trend: deg-4/deg-2 ratio at p=11,13,17: {ratios}")
if ratios[2] > ratios[1]:
    print(f"  Ratio is INCREASING → degree-4 dominance strengthens with p")
    print(f"  This strongly supports HYP-480: Interval wins for ALL p >= 13")
elif ratios[2] < ratios[1]:
    print(f"  Ratio DECREASED from p=13 to p=17 — unexpected!")

# Number of distinct H values
n_distinct = len(set(H for _, H in pairs))
print(f"\n  Distinct H values at p={p}: {n_distinct} out of {N}")

print("\nDONE.")
