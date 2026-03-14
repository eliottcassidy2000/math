#!/usr/bin/env python3
"""
h63_forbidden_analysis.py — opus-2026-03-14-S71g
Investigate whether H=63 is permanently forbidden.

KEY INSIGHT: 63 = 9 × 7 = 3^2 × 7.
All multiplicative decompositions of 63 = I(G,2) involve:
  63 = 7 × 9 = I(C_3, 2) × I(E_2, 2) → C_3 component in Ω
  63 = 21 × 3 = I(P_4, 2) × I(E_1, 2) → P_4 component in Ω
  63 = 63 × 1 → need a CONNECTED graph G with I(G,2)=63

If C_3 and P_4 are both impossible as Ω components (THM-201, THM-202),
then H=63 requires I(G,2)=63 with G connected and avoiding both.

We need to enumerate ALL graphs G (connected) with I(G,2) = 63
and check whether any can be Ω(T).
"""

import math
from itertools import combinations, product
from collections import defaultdict

print("=" * 70)
print("H=63 FORBIDDEN ANALYSIS")
print("opus-2026-03-14-S71g")
print("=" * 70)

# ============================================================
# Part 1: Independence polynomial at x=2
# ============================================================
print("\n" + "=" * 70)
print("PART 1: ALL FACTORIZATIONS OF 63")
print("=" * 70)

print("\n63 = 2^6 - 1 = 9 × 7 = 3^2 × 7")
print()

# I(G1 ⊔ G2, 2) = I(G1, 2) × I(G2, 2)
# So we need all factorizations of 63 into products of achievable I(G, 2) values.

# First: what values can I(G, 2) take for connected graphs G?
# I(single vertex, 2) = 3
# I(K_2, 2) = 1 + 2·2 = 5
# I(K_3, 2) = 1 + 3·2 = 7
# I(P_3, 2) = 1 + 3·2 + 2·4 = ... wait, let me compute properly.

# I(G, x) = Σ_S (independent set) x^|S|
# For P_3 (path 1-2-3):
#   Independent sets: {}, {1}, {2}, {3}, {1,3}
#   I(P_3, 2) = 1 + 2 + 2 + 2 + 4 = 11... wait
#   I(P_3, 2) should be 5 from Jacobsthal: I(P_0)=1, I(P_1)=3, I(P_2)=5
#   Path P_k has k+1 vertices and k edges.

# Let me be careful about notation.
# P_k = path with k VERTICES (not k edges).
# I(P_1, x) = 1 + x (single vertex)
# I(P_2, x) = 1 + 2x (edge: {}, {1}, {2})
# I(P_3, x) = 1 + 3x + x^2 (path 1-2-3: {}, {1},{2},{3}, {1,3})
# I(P_3, 2) = 1 + 6 + 4 = 11

# Hmm, but earlier we had I(P_k, 2) for the recurrence with k being the
# number of vertices (or edges?). Let me align.

# Using the recurrence: a(0)=1, a(1)=1+x, a(k) = a(k-1) + x·a(k-2)
# a(0) = 1 = I of empty graph
# a(1) = 1+x = I of single vertex  [= I(P_1, x)]
# a(2) = (1+x) + x·1 = 1+2x      [= I(P_2, x) = I(edge, x)]
# a(3) = (1+2x) + x(1+x) = 1+3x+x^2  [= I(P_3, x)]
# a(4) = (1+3x+x²) + x(1+2x) = 1+4x+3x²  [= I(P_4, x)]
# a(5) = ... + x(1+3x+x²) = 1+5x+6x²+x³  [= I(P_5, x)]

# At x=2:
# a(0) = 1
# a(1) = 3
# a(2) = 5
# a(3) = 11
# a(4) = 21
# a(5) = 43
# a(6) = 85

# OK so these use vertex count. I(P_k, 2) with k vertices.

def indpoly_at_2(adj, n):
    """Compute I(G, 2) for graph G on n vertices given as adjacency matrix."""
    total = 0
    for mask in range(1 << n):
        # Check if mask is an independent set
        is_indep = True
        for i in range(n):
            if not (mask & (1 << i)):
                continue
            for j in range(i+1, n):
                if not (mask & (1 << j)):
                    continue
                if adj[i][j]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            bits = bin(mask).count('1')
            total += 2**bits
    return total

# Check some known values
# K_3 (triangle): I = 1 + 3·2 = 7
adj_k3 = [[0,1,1],[1,0,1],[1,1,0]]
print(f"I(K_3, 2) = {indpoly_at_2(adj_k3, 3)} (expected 7)")

# C_3 = K_3: same thing
# P_4: path 1-2-3-4
adj_p4 = [[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]]
print(f"I(P_4, 2) = {indpoly_at_2(adj_p4, 4)} (expected 21)")

# C_4: cycle 1-2-3-4-1
adj_c4 = [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]
print(f"I(C_4, 2) = {indpoly_at_2(adj_c4, 4)} (expected 17)")

# C_5: cycle on 5 vertices
adj_c5 = [[0]*5 for _ in range(5)]
for i in range(5):
    adj_c5[i][(i+1)%5] = adj_c5[(i+1)%5][i] = 1
print(f"I(C_5, 2) = {indpoly_at_2(adj_c5, 5)} (expected 31)")

print(f"\n63 = 9 × 7 = 3² × 7")
print(f"  3 = I(single vertex) → component = isolated vertex")
print(f"  7 = I(K_3) → component = K_3 = triangle")
print(f"  9 = I(2 isolated vertices) = 3² → two isolated vertices")
print(f"  21 = I(P_4) → component = P_4")

# ============================================================
# Part 2: Enumerate all connected graphs with I(G,2) = 63
# ============================================================
print("\n" + "=" * 70)
print("PART 2: CONNECTED GRAPHS WITH I(G,2) = 63")
print("=" * 70)

# We need to find all connected graphs G with I(G,2) = 63.
# 63 = I(G,2) means: Σ_{independent S} 2^|S| = 63
# Max |S| ≤ α(G) (independence number)
# If α(G) = k, then I(G,2) ≥ 1 + n·2 ≥ 1+2k (from singletons)
# and I(G,2) ≤ 3^k (from k independent vertices contributing at most 3^k)
# Wait: I(G,2) ≤ 3^n for n vertices (all subsets).
#
# 63 ≤ 3^n → n ≥ 4 (since 3^3=27 < 63, 3^4=81 ≥ 63)
# Also I(G,2) ≥ 1 + 2n, so 63 ≥ 1+2n → n ≤ 31.
# But connected graphs with many vertices and few edges have high I.
# If G is connected on n vertices, I(G,2) ≥ I(T_n, 2) where T_n is a tree.
# For a path on n vertices: I(P_n, 2) ~ (2^{n+2}+...)/3
# I(P_n, 2) > 63 when n ≥ 7 (I(P_7,2) = 171).
# Actually I(P_6,2) = 85 > 63. I(P_5,2) = 43 < 63.
# But adding edges to P_n can only DECREASE I (more constraints).
# So connected graphs with n ≥ 7 on a spanning tree P_7 have I ≥ 171?
# No! Adding edges decreases I. The tree itself gives the maximum for trees.
# Hmm, actually denser graphs have SMALLER I.
# For K_n: I(K_n, 2) = 1+2n, which is small.
# For sparse: I is large. For dense: I is small.
# I(K_31, 2) = 63. So K_31 achieves 63!
# Also I(K_n minus edges, 2) varies.

# Let's be systematic. For connected graph on n vertices:
# I(K_n, 2) = 1+2n → to get 63: n=31 (but huge)
# I(K_n minus one edge, 2) = 1+2n-1+4 = 2n+4 for n≥3
# Wait: K_n minus edge e={u,v}: independent sets are same as K_n plus {u,v}
# I(K_n\e, 2) = I(K_n, 2) + 2^2 = 1+2n+4 = 2n+5
# Hmm: 2n+5 = 63 → n=29. Still large.

# For small n, enumerate:
print("Enumerating connected graphs on n=4,5,6 vertices with I(G,2) = 63:\n")

def generate_connected_graphs(n):
    """Generate all non-isomorphic connected graphs on n vertices (brute force)."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    results = []
    seen = set()

    for mask in range(1, 1 << m):
        # Build adjacency matrix
        adj = [[0]*n for _ in range(n)]
        edge_list = []
        for k in range(m):
            if mask & (1 << k):
                i, j = edges[k]
                adj[i][j] = adj[j][i] = 1
                edge_list.append((i, j))

        # Check connectivity (BFS)
        visited = {0}
        queue = [0]
        while queue:
            v = queue.pop()
            for u in range(n):
                if adj[v][u] and u not in visited:
                    visited.add(u)
                    queue.append(u)

        if len(visited) < n:
            continue

        # Compute I(G, 2)
        I_val = indpoly_at_2(adj, n)

        if I_val == 63:
            # Compute degree sequence for basic isomorphism check
            degs = tuple(sorted([sum(adj[i]) for i in range(n)]))
            edge_count = sum(sum(row) for row in adj) // 2

            # Compute independence number
            alpha = 0
            for s in range(1 << n):
                is_indep = True
                for i in range(n):
                    if not (s & (1 << i)): continue
                    for j in range(i+1, n):
                        if not (s & (1 << j)): continue
                        if adj[i][j]:
                            is_indep = False
                            break
                    if not is_indep: break
                if is_indep:
                    alpha = max(alpha, bin(s).count('1'))

            key = (degs, edge_count, alpha)
            if key not in seen:
                seen.add(key)
                results.append({
                    'n': n, 'edges': edge_count, 'degs': degs,
                    'alpha': alpha, 'adj': [row[:] for row in adj]
                })

    return results

for n in range(4, 8):
    print(f"n={n}:", end=" ")
    if n <= 6:
        graphs = generate_connected_graphs(n)
        print(f"Found {len(graphs)} connected graphs with I(G,2)=63")
        for g in graphs:
            print(f"    edges={g['edges']}, degrees={g['degs']}, α={g['alpha']}")
    else:
        # n=7 has too many graphs for brute force with all edges
        # But we can check specific structures
        print("(checking specific structures only)")

        # K_7 minus edges
        # I(K_7, 2) = 1+14 = 15. Need to add 48 to get 63.
        # Each removed edge adds a 2-element independent set: +4
        # But removing edges can enable larger independent sets too.
        # Removing e edges from K_7: I ≈ 15 + 4e + higher order
        # For 63: need ~48/4 = 12 edges removed, leaving C(7,2)-12 = 9 edges

        # Check specific:
        # Star K_{1,6}: I = 1 + 2 + 6·2 + C(6,2)·4 + C(6,3)·8 + ... (center vs leaves)
        # Center is adjacent to all leaves. Leaves are independent.
        # Independent sets: subsets of leaves ∪ {center alone}
        # = all subsets of {leaves} + ∅ union {center}
        # I = (1+2)^6 + 2 - 1 = 729 + 1 = 730. Way too high.
        # Wait: I(star, 2) = Σ_{S ⊆ leaves} 2^|S| + 2^1 = 3^6 + 2 = 731.
        # No: empty set + center = 1 + 2. Subsets of 6 leaves: 3^6 = 729 (including empty).
        # But empty set is counted once: I = 729 + 2 = 731 - 1 = 730. Hmm.
        # Actually: independent sets of star K_{1,6}:
        #   - Empty: 1
        #   - {center}: 1
        #   - subsets of leaves (non-empty): 2^6 - 1 = 63
        #   - {center} is not adjacent to any leaf? No! Center IS adjacent to all leaves.
        #   So center cannot be in any set with a leaf.
        #   Independent sets: all subsets of leaves (2^6 = 64) + {center} (1) = 65.
        #   Wait no: empty is a subset of leaves AND already counted.
        #   I(star, x) = (1+x)^6 + x = Σ_S x^|S| for S ⊆ leaves, plus x for {center}
        #   I(K_{1,6}, 2) = 3^6 + 2 = 729 + 2 = 731. Way too high for 63.

        # Dense graphs: K_7 minus few edges
        # I(K_7, 2) = 15
        # Removing 1 edge: +4 → 19
        # Removing 2 disjoint edges: +8 → 23, or +4+4 with possible 4-set = 23 or 27
        # This is getting complex. Let's just check for n≤6.
        pass

# ============================================================
# Part 3: Graph structures that achieve I(G,2) = 63
# ============================================================
print("\n" + "=" * 70)
print("PART 3: DISCONNECTED DECOMPOSITIONS")
print("=" * 70)

print("""
Since I(G1 ⊔ G2, 2) = I(G1) · I(G2), we need factorizations of 63:

63 = 1 × 63   → G = connected graph with I=63
63 = 3 × 21   → G = vertex ⊔ (graph with I=21)
63 = 7 × 9    → G = K_3 ⊔ (graph with I=9)
63 = 9 × 7    → same as above
63 = 21 × 3   → G = (graph with I=21) ⊔ vertex

Factorizations with more factors:
63 = 3 × 3 × 7   → vertex ⊔ vertex ⊔ K_3
63 = 3 × 21      → vertex ⊔ (I=21 graph)
63 = 7 × 3 × 3   → K_3 ⊔ vertex ⊔ vertex

So the DISCONNECTED options are:
1. E_2 ⊔ K_3 (two isolated vertices + triangle) [I = 9·7 = 63]
2. E_1 ⊔ G_{21} (one isolated vertex + I=21 component) [I = 3·21 = 63]

For option 1: Ω(T) contains K_3 as component → IMPOSSIBLE by THM-201 ✓
For option 2: need I(G,2) = 21 as a connected component of Ω.
  The five graphs with I=21 were enumerated: P_4, K_3+K_1, K_6-2e, K_8-e, K_10.
  P_4 is ruled out by THM-202.
  K_3+K_1 is disconnected (ruled out since we need connected component with I=21).

  Actually: K_3+K_1 is disconnected, so as a connected component it doesn't work.
  We'd need: Ω has components {vertex, component_with_I=21}
  The component with I=21 must be connected: K_6-2e, K_8-e, K_10, or P_4.
  P_4 is ruled out by THM-202.
  The three remaining are what we need to investigate.

For the connected case (I=63 as single connected component):
This requires finding connected graphs with I(G,2) = 63.
""")

# ============================================================
# Part 4: Quick check — is 63 achievable at specific n?
# ============================================================
print("=" * 70)
print("PART 4: IS H=63 EVER ACHIEVED?")
print("=" * 70)

# Quick check at small n
from itertools import permutations
import random

def count_hp(adj_matrix, n):
    """Count Hamiltonian paths by permutation check."""
    count = 0
    # For small n, direct count
    if n <= 8:
        for perm in permutations(range(n)):
            valid = True
            for i in range(n-1):
                if adj_matrix[perm[i]][perm[i+1]] == 0:
                    valid = False
                    break
            if valid:
                count += 1
    return count

def random_tournament(n):
    """Generate random tournament on n vertices."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

# Check at n=5: H values are {1,3,5,9,11,13,15}
# 63 > 15, so not achievable at n=5.

# Check at n=7: H range is roughly 1 to ~5040
# Need to check if 63 appears
print("Checking if H=63 appears at n=7 (sampling)...")

random.seed(42)
h63_count = 0
total_checked = 0
for _ in range(100000):
    T = random_tournament(7)
    h = count_hp(T, 7)
    if h == 63:
        h63_count += 1
    total_checked += 1

print(f"  n=7: {h63_count}/{total_checked} tournaments have H=63")

if h63_count == 0:
    print("  H=63 NOT found at n=7 in 100K samples!")
    print("  (This is consistent with 63 being forbidden)")
else:
    print(f"  H=63 found! Frequency: {h63_count/total_checked:.6f}")

# Also check n=6
print("\nChecking at n=6 exhaustively...")
h_values_6 = set()
for mask in range(1 << 15):  # C(6,2) = 15 arcs
    adj = [[0]*6 for _ in range(6)]
    idx = 0
    for i in range(6):
        for j in range(i+1, 6):
            if mask & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    h = count_hp(adj, 6)
    h_values_6.add(h)

print(f"  n=6: H values = {sorted(h_values_6)}")
print(f"  63 in H-spectrum at n=6: {63 in h_values_6}")

# ============================================================
# Part 5: Structural analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: WHY 63 MIGHT BE FORBIDDEN")
print("=" * 70)

print("""
ARGUMENT FOR H=63 BEING FORBIDDEN:

1. 63 = 3² × 7. The ONLY multiplicative decompositions are:
   1 × 63, 3 × 21, 7 × 9, 3 × 3 × 7

2. Via OCF H = I(Ω, 2), we need Ω(T) with I(Ω, 2) = 63.

3. DISCONNECTED cases (I factors):
   a) 3 × 3 × 7: Ω has components {vertex, vertex, K_3}.
      K_3 as Ω component → IMPOSSIBLE (THM-201).
   b) 3 × 21: Ω has components {vertex, G_21}.
      G_21 must be connected with I=21.
      Options: P_4 (THM-202 rules out), K_6-2e, K_8-e, K_10.
      If these are all impossible as Ω components → this case falls.
   c) 7 × 9: Ω has components {K_3, 2 vertices} or {K_3, E_2}.
      K_3 as component → IMPOSSIBLE (THM-201).

4. CONNECTED case (I(G,2) = 63 for connected G):
   Need to enumerate all connected graphs with this property.
   From n≤6 search above: check what we found.

CONCLUSION: H=63 is forbidden IF:
   (a) THM-201 (K_3 component impossible) — PROVED ✓
   (b) All connected graphs with I=21 are impossible as Ω components
   (c) All connected graphs with I=63 are impossible as Ω
""")

print("\n" + "=" * 70)
print("DONE — H=63 ANALYSIS")
print("=" * 70)
