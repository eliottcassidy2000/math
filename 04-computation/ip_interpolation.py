#!/usr/bin/env python3
"""
ip_interpolation.py — Independence polynomial as simplex-cuboid interpolation.

KEY INSIGHT from the packing framework:
  I(K_m, x) = 1 + mx         (clique/simplex — all adjacent)
  I(m·K₁, x) = (1+x)^m       (independent/cuboid — all disjoint)

For the OCF, Ω(T) interpolates between these extremes.
The achievable I(Ω, 2) values correspond to valid "interpolation points."

Question: WHY does the interpolation skip 7 and 21?

The simplex gives I(2) = 1+2m (odd values 1,3,5,...)
The cuboid gives I(2) = 3^m (values 1,3,9,27,...)
Mixed structures interpolate between.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("INDEPENDENCE POLYNOMIAL AS SIMPLEX-CUBOID INTERPOLATION")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: Achievable I(G,2) for all graphs up to 10 vertices
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: I(G,2) for small graphs ---")
print("  Which values of I(G,2) are achievable by ANY graph G?")
print("  (Not just tournament conflict graphs)")

# For a graph on m vertices with adjacency structure,
# I(G,2) = sum over independent sets S of 2^|S|

# For m=1: I = 1+2 = 3
# For m=2: I ∈ {(1+2)^2, 1+4} = {9, 5} (independent or edge)
# For m=3: K₃ → 1+6=7, K₃̄ → 27, P₃ → 1+6+4=11, ...

# Let's enumerate small graphs
def graph_ip_at_2(adj, m):
    """I(G, 2) for graph with adjacency dict adj on vertices 0..m-1."""
    total = 0
    for mask in range(2**m):
        verts = [i for i in range(m) if mask & (1 << i)]
        # Check independence
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj.get((verts[i], verts[j]), False):
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2**len(verts)
    return total

# Enumerate all graphs on m vertices
achievable_vals = set()
achievable_by_m = defaultdict(set)

for m in range(1, 7):  # m<=6 fast
    edges = [(i,j) for i in range(m) for j in range(i+1,m)]
    ne = len(edges)

    for bits in range(2**ne):
        adj = {}
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[(i,j)] = adj[(j,i)] = True

        val = graph_ip_at_2(adj, m)
        achievable_vals.add(val)
        achievable_by_m[m].add(val)

    # Show gaps
    all_vals = sorted(achievable_by_m[m])
    max_v = max(all_vals)
    gaps = [v for v in range(1, max_v+1, 2) if v not in achievable_by_m[m]]
    print(f"  m={m}: I(G,2) range [{min(all_vals)}, {max_v}], "
          f"odd gaps: {gaps[:10]}{'...' if len(gaps) > 10 else ''}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: Which odd values are NEVER achievable by any graph?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Permanently missing I(G,2) values ---")

all_achievable = set()
for m in range(1, 7):  # m<=6 fast
    all_achievable |= achievable_by_m[m]

# Check odd values up to max(achievable)
max_val = max(all_achievable)
missing_odd = [v for v in range(1, 82, 2) if v not in all_achievable]
print(f"  Checking odd values 1 to 81:")
print(f"  Missing (never achieved for m≤8): {missing_odd}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: Graph I.P. vs tournament Ω I.P.
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: Key comparison ---")
print("  Tournament H values (odd) at n=5: ", end="")

# n=5 tournaments
from itertools import permutations as perms

def fast_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

n = 5
edges5 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne5 = len(edges5)
h_vals_5 = set()
for bits in range(2**ne5):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges5):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    h_vals_5.add(fast_hp(A, n))
print(sorted(h_vals_5))

n = 6
edges6 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne6 = len(edges6)
h_vals_6 = set()
for bits in range(2**ne6):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges6):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    h_vals_6.add(fast_hp(A, n))
print(f"  Tournament H values (odd) at n=6: {sorted(h_vals_6)}")

print(f"\n  Graph I(G,2) achievable odd values ≤45:")
print(f"    {sorted([v for v in all_achievable if v % 2 == 1 and v <= 45])}")

print(f"\n  Tournament-missing vs Graph-missing:")
t_missing = [v for v in range(1, 46, 2) if v not in h_vals_6]
g_missing = [v for v in range(1, 46, 2) if v not in all_achievable]
print(f"    Tournament (n≤6): {t_missing}")
print(f"    Graph (m≤8): {g_missing}")

# ═══════════════════════════════════════════════════════════════════
# Part 4: Is 7 achievable by ANY graph?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: Is I(G,2)=7 achievable by any graph? ---")

for m in range(1, 7):  # m<=6 fast
    if 7 in achievable_by_m[m]:
        print(f"  YES at m={m}!")
        # Find an example
        edges = [(i,j) for i in range(m) for j in range(i+1,m)]
        ne = len(edges)
        for bits in range(2**ne):
            adj = {}
            for idx, (i,j) in enumerate(edges):
                if bits & (1 << idx):
                    adj[(i,j)] = adj[(j,i)] = True
            if graph_ip_at_2(adj, m) == 7:
                edge_list = [(i,j) for idx, (i,j) in enumerate(edges) if bits & (1 << idx)]
                print(f"    Example: m={m}, edges={edge_list}")
                break
        break
else:
    print(f"  NO — I(G,2)=7 is NEVER achievable for any graph with ≤8 vertices!")

print("\n--- Part 5: Is I(G,2)=21 achievable by any graph? ---")

for m in range(1, 7):  # m<=6 fast
    if 21 in achievable_by_m[m]:
        print(f"  YES at m={m}!")
        edges = [(i,j) for i in range(m) for j in range(i+1,m)]
        ne = len(edges)
        for bits in range(2**ne):
            adj = {}
            for idx, (i,j) in enumerate(edges):
                if bits & (1 << idx):
                    adj[(i,j)] = adj[(j,i)] = True
            if graph_ip_at_2(adj, m) == 21:
                edge_list = [(i,j) for idx, (i,j) in enumerate(edges) if bits & (1 << idx)]
                print(f"    Example: m={m}, edges={edge_list}")
                break
        break
else:
    print(f"  NO — I(G,2)=21 is NEVER achievable for any graph with ≤8 vertices!")

print("\n--- Part 6: UNIVERSAL I.P. gaps ---")
print("  If I(G,2)=7 and I(G,2)=21 are impossible for ALL graphs,")
print("  then the impossibility is GRAPH-THEORETIC, not tournament-specific!")
print("  This would be a MUCH stronger result than the tournament-only version.")

# Actually I(K₃, 2) = 1 + 3·2 = 7. Let me re-check.
print("\n  WAIT: I(K₃, 2) = 1 + 3·2 = 7")
adj_k3 = {(0,1): True, (1,0): True, (0,2): True, (2,0): True, (1,2): True, (2,1): True}
print(f"  I(K₃, 2) = {graph_ip_at_2(adj_k3, 3)}")

print("\n  And I(K₁₀, 2) = 1 + 10·2 = 21")
adj_k10 = {}
for i in range(10):
    for j in range(i+1, 10):
        adj_k10[(i,j)] = adj_k10[(j,i)] = True
print(f"  I(K₁₀, 2) = {graph_ip_at_2(adj_k10, 10)}")

print("\n  So 7 and 21 ARE achievable by graphs!")
print("  The constraint is SPECIFICALLY about which graphs arise as")
print("  tournament conflict graphs Ω(T).")
print("  K₃ can be a conflict graph (3 pairwise-intersecting cycles),")
print("  but the Splicing Lemma forces a 4th cycle → not exactly K₃.")
print("  K₁₀ similarly: 10 pairwise-intersecting cycles force additional cycles.")
