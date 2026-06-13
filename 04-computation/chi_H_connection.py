#!/usr/bin/env python3
"""
Explore the relationship between H (Hamiltonian path count) and χ (Euler characteristic
of path homology) across ALL tournaments.

Key finding at n=7 regular:
  H=189 (Paley): χ=7=p
  H=171: χ=1
  H=175: χ=0

Question: Is there a formula relating H and χ?
- Is χ monotone in H?
- Does χ=p characterize Paley-like tournaments?
- What values can χ take?

opus-2026-03-13-S71b
"""

import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers
import itertools
from collections import Counter, defaultdict
import time

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

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

# Analyze n=5 exhaustively
for n in [5, 6]:
    print(f"\n{'='*60}")
    print(f"n = {n}: All tournaments")
    print(f"{'='*60}")

    H_chi = defaultdict(list)
    chi_counts = Counter()
    t0 = time.time()

    count = 0
    for adj in all_tournaments(n):
        betti = path_betti_numbers(adj, n)
        chi = sum((-1)**d * b for d, b in enumerate(betti))
        H = count_hp(adj, n)
        H_chi[H].append(chi)
        chi_counts[chi] += 1
        count += 1

    t1 = time.time()
    print(f"  {count} tournaments in {t1-t0:.1f}s")

    print(f"\n  χ value distribution:")
    for chi_val in sorted(chi_counts.keys()):
        print(f"    χ={chi_val}: {chi_counts[chi_val]} tournaments")

    print(f"\n  H vs χ (is H→χ a function?):")
    non_functional = 0
    for H in sorted(H_chi.keys()):
        chi_vals = sorted(set(H_chi[H]))
        marker = "" if len(chi_vals) == 1 else " ← MULTIPLE χ!"
        if len(chi_vals) > 1:
            non_functional += 1
        print(f"    H={H:3d}: χ ∈ {chi_vals}, count={len(H_chi[H])}{marker}")

    if non_functional == 0:
        print(f"\n  *** H determines χ at n={n}! ***")
    else:
        print(f"\n  H does NOT determine χ at n={n} ({non_functional} ambiguous H values)")

    # Check: is χ always 0 or 1?
    all_01 = all(c in {0, 1} for c in chi_counts.keys())
    print(f"\n  χ ∈ {{0,1}} only: {all_01}")
    if not all_01:
        print(f"  Other χ values: {[c for c in sorted(chi_counts.keys()) if c not in {0,1}]}")

    # Check: which tournament has max χ?
    max_chi = max(chi_counts.keys())
    print(f"  Max χ = {max_chi}")

# n=7 would take too long exhaustively, sample
print(f"\n{'='*60}")
print(f"n = 7: Random sample (1000 tournaments)")
print(f"{'='*60}")

import random
n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

H_chi = defaultdict(list)
chi_counts = Counter()
betti_counts = Counter()
t0 = time.time()

for _ in range(1000):
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
print(f"  1000 tournaments in {t1-t0:.1f}s")

print(f"\n  χ distribution:")
for chi_val in sorted(chi_counts.keys()):
    print(f"    χ={chi_val}: {chi_counts[chi_val]} tournaments")

print(f"\n  β distribution (top 10):")
for betti, cnt in betti_counts.most_common(10):
    chi = sum((-1)**d * b for d, b in enumerate(betti))
    print(f"    β={betti}: {cnt} ({chi})")

print(f"\n  H vs χ (is H→χ a function?):")
non_functional = 0
for H in sorted(H_chi.keys()):
    chi_vals = sorted(set(H_chi[H]))
    if len(chi_vals) > 1:
        non_functional += 1
        print(f"    H={H:3d}: χ ∈ {chi_vals} ← MULTIPLE!")

if non_functional == 0:
    print(f"    All H values map to unique χ in sample")
else:
    print(f"    {non_functional} H values with multiple χ")
