#!/usr/bin/env python3
"""
Explore H vs χ at n=5 (exhaustive) and n=7 (sample).

opus-2026-03-13-S71b
"""

import sys, time, random
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers
from collections import Counter, defaultdict

def count_hp(adj, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and adj[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# n=5 exhaustive
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)
print(f"n={n}: exhaustive ({2**m} tournaments)")

H_chi = defaultdict(list)
chi_counts = Counter()
t0 = time.time()

for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    betti = path_betti_numbers(adj, n)
    chi = sum((-1)**d * b for d, b in enumerate(betti))
    H = count_hp(adj, n)
    H_chi[H].append(chi)
    chi_counts[chi] += 1

t1 = time.time()
print(f"  Done in {t1-t0:.1f}s")

print(f"\n  χ distribution: {dict(sorted(chi_counts.items()))}")
print(f"\n  H → χ:")
for H in sorted(H_chi.keys()):
    chi_vals = sorted(set(H_chi[H]))
    n_tot = len(H_chi[H])
    marker = " ***" if len(chi_vals) > 1 else ""
    print(f"    H={H:3d}: χ ∈ {chi_vals}, count={n_tot}{marker}")

# n=7 sample
n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)
print(f"\nn={n}: random sample (500)")

H_chi = defaultdict(list)
chi_counts = Counter()
betti_counts = Counter()
t0 = time.time()

for _ in range(500):
    adj = [[0]*n for _ in range(n)]
    for i, j in edges:
        if random.random() < 0.5:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    betti = path_betti_numbers(adj, n)
    chi = sum((-1)**d * b for d, b in enumerate(betti))
    H = count_hp(adj, n)
    H_chi[H].append(chi)
    chi_counts[chi] += 1
    betti_counts[tuple(betti)] += 1

t1 = time.time()
print(f"  Done in {t1-t0:.1f}s")

print(f"\n  χ distribution: {dict(sorted(chi_counts.items()))}")
print(f"\n  β distribution (top 10):")
for betti, cnt in betti_counts.most_common(10):
    chi = sum((-1)**d * b for d, b in enumerate(betti))
    print(f"    β={betti}: count={cnt}, χ={chi}")

print(f"\n  H → χ ambiguities:")
ambig = 0
for H in sorted(H_chi.keys()):
    chi_vals = sorted(set(H_chi[H]))
    if len(chi_vals) > 1:
        ambig += 1
        print(f"    H={H}: χ ∈ {chi_vals}")
if ambig == 0:
    print(f"    None! H determines χ in sample")
else:
    print(f"    {ambig} H values with multiple χ")
