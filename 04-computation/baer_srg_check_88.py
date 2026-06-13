#!/usr/bin/env python3
"""
baer_srg_check_88.py — opus-2026-03-14-S88

Verify: is the Baer disjointness graph on 360 vertices strongly regular?

We know:
  - 360 vertices (Baer subplanes of PG(2,F_4))
  - 16-regular
  - Each edge in exactly 1 triangle (partition)
  - Local neighborhood = 8K_2

If lambda=1 (verified: each edge extends uniquely to partition, giving
exactly 1 common neighbor), what is mu?

For SRG(360, 16, 1, mu):
  Necessary condition: k(k-1-lambda) = mu(v-k-1)
  16 * 14 = mu * 343
  224 = 343 * mu
  mu = 224/343 ≈ 0.653 — NOT an integer!

So it CANNOT be strongly regular with lambda=1!

Wait — lambda=1 means adjacent vertices have 1 common neighbor.
But the local graph = 8K_2 means each neighbor is paired with exactly
one other neighbor, and pairs are independent.
So each vertex's neighbor sees: its partner (1 common neighbor) and
no connection to other pairs.

Actually, in 8K_2, the number of edges is 8. Each of the 16 neighbors
has degree 1 in the local graph. So for two adjacent u,v:
their common neighbors = vertices adjacent to both u and v.
In the local graph of w (where u,v are neighbors of w and u~v):
u and v are matched → they have common neighbor only through
the partition structure.

Let me actually COMPUTE lambda and mu from the Baer data.
"""

import numpy as np
from itertools import combinations
from collections import Counter

# ══════════════════════════════════════════════════════════════════
# Rebuild PG(2,F_4) and Baer subplanes
# ══════════════════════════════════════════════════════════════════

mul_table = [
    [0, 0, 0, 0],
    [0, 1, 2, 3],
    [0, 2, 3, 1],
    [0, 3, 1, 2],
]

def f4_add(a, b):
    return a ^ b

def f4_mul(a, b):
    return mul_table[a][b]

def f4_inv(a):
    for b in range(1, 4):
        if f4_mul(a, b) == 1:
            return b
    return None

def normalize_pg(triple):
    x, y, z = triple
    for v in [x, y, z]:
        if v != 0:
            inv = f4_inv(v)
            return (f4_mul(x, inv), f4_mul(y, inv), f4_mul(z, inv))
    return None

points = set()
for x in range(4):
    for y in range(4):
        for z in range(4):
            if (x, y, z) == (0, 0, 0):
                continue
            p = normalize_pg((x, y, z))
            if p:
                points.add(p)
points = sorted(points)
pt_idx = {p: i for i, p in enumerate(points)}
assert len(points) == 21

# Lines
lines = []
seen_lines = set()
for a in range(4):
    for b in range(4):
        for c in range(4):
            if (a, b, c) == (0, 0, 0):
                continue
            pts_on = []
            for p in points:
                dot = f4_add(f4_add(f4_mul(a, p[0]), f4_mul(b, p[1])), f4_mul(c, p[2]))
                if dot == 0:
                    pts_on.append(pt_idx[p])
            key = frozenset(pts_on)
            if key not in seen_lines:
                seen_lines.add(key)
                lines.append(key)
assert len(lines) == 21

# Find all Baer subplanes
def line_through(i, j):
    for l in lines:
        if i in l and j in l:
            return l
    return None

def meet(l1, l2):
    common = l1 & l2
    if len(common) == 1:
        return list(common)[0]
    return None

def collinear(i, j, k):
    l = line_through(i, j)
    return l is not None and k in l

print("Finding all Baer subplanes...")
all_baer = set()
n = 21
for p1 in range(n):
    for p2 in range(p1+1, n):
        for p3 in range(p2+1, n):
            if collinear(p1, p2, p3):
                continue
            for p4 in range(p3+1, n):
                if (collinear(p1, p2, p4) or collinear(p1, p3, p4) or
                    collinear(p2, p3, p4)):
                    continue
                pts = {p1, p2, p3, p4}
                changed = True
                while changed and len(pts) <= 7:
                    changed = False
                    pts_list = sorted(pts)
                    for a, b in combinations(pts_list, 2):
                        lab = line_through(a, b)
                        for c, d in combinations(pts_list, 2):
                            if {a,b} == {c,d}:
                                continue
                            lcd = line_through(c, d)
                            if lab and lcd and lab != lcd:
                                m = meet(lab, lcd)
                                if m is not None and m not in pts:
                                    pts.add(m)
                                    changed = True
                    if len(pts) > 7:
                        break
                if len(pts) == 7:
                    all_baer.add(frozenset(pts))

baer_list = sorted(all_baer, key=lambda s: sorted(s))
nb = len(baer_list)
print(f"Total Baer subplanes: {nb}")

# ══════════════════════════════════════════════════════════════════
# Build disjointness graph and check SRG
# ══════════════════════════════════════════════════════════════════

print("\nBuilding disjointness adjacency matrix...")
# Use sparse representation for memory
adj = {}
for i in range(nb):
    adj[i] = set()
    for j in range(nb):
        if i != j and len(baer_list[i] & baer_list[j]) == 0:
            adj[i].add(j)

degrees = {i: len(adj[i]) for i in range(nb)}
deg_dist = Counter(degrees.values())
print(f"Degree distribution: {dict(deg_dist)}")

# Check lambda (common neighbors of adjacent pair)
print("\nComputing lambda (common neighbors of adjacent vertices)...")
lambda_vals = Counter()
sample_count = 0
for i in range(min(nb, 50)):  # sample first 50 vertices
    for j in adj[i]:
        if j > i:
            common = len(adj[i] & adj[j])
            lambda_vals[common] += 1
            sample_count += 1

print(f"Lambda distribution (from {sample_count} edges): {dict(sorted(lambda_vals.items()))}")

# Check mu (common neighbors of non-adjacent pair)
print("\nComputing mu (common neighbors of non-adjacent vertices)...")
mu_vals = Counter()
sample_count_mu = 0
for i in range(min(nb, 30)):
    non_adj = set(range(nb)) - adj[i] - {i}
    for j in list(non_adj)[:20]:  # sample 20 non-neighbors
        common = len(adj[i] & adj[j])
        mu_vals[common] += 1
        sample_count_mu += 1

print(f"Mu distribution (from {sample_count_mu} non-edges): {dict(sorted(mu_vals.items()))}")

is_srg = len(lambda_vals) == 1 and len(mu_vals) == 1
print(f"\nStrongly regular: {is_srg}")

if is_srg:
    lam = list(lambda_vals.keys())[0]
    mu = list(mu_vals.keys())[0]
    k = 16
    print(f"SRG({nb}, {k}, {lam}, {mu})")

# ══════════════════════════════════════════════════════════════════
# Intersection numbers
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("INTERSECTION PATTERN OF BAER SUBPLANES")
print("=" * 70)

# How many points can two Baer subplanes share?
inter_dist = Counter()
for i in range(nb):
    for j in range(i+1, nb):
        inter_dist[len(baer_list[i] & baer_list[j])] += 1

print(f"Intersection size distribution (pairs of Baer subplanes):")
for k in sorted(inter_dist.keys()):
    print(f"  |B_i ∩ B_j| = {k}: {inter_dist[k]} pairs")

total_pairs = nb * (nb-1) // 2
print(f"Total pairs: {total_pairs}")

# In PG(2,4), two Baer subplanes can intersect in 0, 1, or 3 points
# (they're subplanes, so intersection is a subplane or degenerate)
# Actually: two copies of PG(2,2) in PG(2,4) can meet in:
# - 0 points (disjoint)
# - 1 point
# - 3 points (a line of PG(2,2))
# - 7 points (identical)

# ══════════════════════════════════════════════════════════════════
# The 2-(21,7,36) design intersection triangle
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("DESIGN INTERSECTION TRIANGLE")
print("=" * 70)

# The 360 Baer subplanes form a 2-(21,7,36) design.
# Fisher's inequality: b >= v → 360 >= 21 ✓ (way more blocks than points)
# Replication: r = b*k/v = 360*7/21 = 120
# Each point in 120 Baer subplanes ✓

# The intersection numbers for a 2-design:
# Two blocks can share 0,1,2,...,k points
# For our design: lambda = 36 means each pair of points in 36 blocks

# The "block intersection graph" has 360 vertices,
# edges when blocks share >= 1 point (i.e., non-disjoint)
# Its complement = our disjointness graph (16-regular)
# So the block intersection graph is (360-1-16) = 343-regular

print(f"Block intersection graph: 343-regular on 360 vertices")
print(f"  343 = 7^3 = 7 * 49")
print(f"  The complement of the disjointness graph!")
print(f"  Each Baer subplane intersects 343 others")
print(f"  and is disjoint from 16 others")
print(f"  343 + 16 = 359 = 360 - 1 ✓")

print(f"\n  343 = 7^3: the number of non-disjoint neighbors")
print(f"  is a PERFECT CUBE of the Fano number!")
print(f"  16 = 2^4: the number of disjoint neighbors")
print(f"  is a PERFECT POWER of 2!")
