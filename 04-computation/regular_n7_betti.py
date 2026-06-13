#!/usr/bin/env python3
"""
Compute Betti numbers for each of the 3 H-classes of regular n=7 tournaments.

Question: Does path homology distinguish the three classes?
H=189 (Paley, disj=7), H=171 (disj=10), H=175 (disj=14)

opus-2026-03-13-S71b
"""

import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers, enumerate_allowed_paths
import itertools
import time

def all_tournaments_n7():
    n = 7
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

def is_regular(adj, n):
    return all(sum(adj[i]) == (n-1)//2 for i in range(n))

n = 7
print(f"Computing Betti numbers for regular n={n} tournaments")
print(f"(3 H-classes: 189, 171, 175)\n")

# Find one representative of each H-class
reps = {}
count = 0
for adj in all_tournaments_n7():
    if not is_regular(adj, n):
        continue
    count += 1
    H = count_hp(adj, n)
    if H not in reps:
        reps[H] = adj
        print(f"Found H={H} representative (tournament #{count})")
    if len(reps) == 3:
        break

# Compute Betti numbers for each
for H in sorted(reps.keys()):
    adj = reps[H]
    t0 = time.time()
    betti = path_betti_numbers(adj, n)
    t1 = time.time()

    chi = sum((-1)**d * b for d, b in enumerate(betti))
    omega_dims = []
    for p in range(n):
        allowed = enumerate_allowed_paths(adj, n, p)
        # This gives |A_p|, not |Omega_p| directly. But Betti gives us the info.
        omega_dims.append(len(allowed))

    print(f"\nH={H}: β = {betti}, χ = {chi}, time = {t1-t0:.1f}s")
    print(f"  |A_d| = {omega_dims}")
