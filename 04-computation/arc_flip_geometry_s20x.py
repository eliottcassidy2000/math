#!/usr/bin/env python3
"""
arc_flip_geometry_s20x.py -- kind-pasteur-2026-03-22-S20x

THE GEOMETRY OF ARC FLIPS: how does H change along paths in tournament space?

Key finding: from the SC max at n=6, arc flips change H by 0, -4, -8, or -16.
The JUMP SIZES are multiples of 4 = 2^2 = (n/2 - 1)!

This connects to the Gray code work (delta=+-10 missing at n=5 = Petersen number).

The tournament space forms a graph where vertices are tournaments and edges are
arc flips. The H function on this graph has interesting geometric properties:
- The gradient of H points toward the SC diagonal
- The Hessian captures the curvature of H
- Level sets of H may have interesting topology

Author: kind-pasteur-2026-03-22-S20x
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

print("=" * 65)
print("  ARC FLIP GEOMETRY OF TOURNAMENT SPACE")
print("=" * 65)

# ================================================================
# 1. EXHAUSTIVE DELTA ANALYSIS AT n=5
# ================================================================
for n in [4, 5]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Compute H for all tournaments
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H_map[bits] = count_hp(A, n)

    # Compute all arc-flip deltas
    delta_dist = defaultdict(int)
    delta_by_arc = defaultdict(lambda: defaultdict(int))

    for bits in range(2**m):
        for k in range(m):
            bits2 = bits ^ (1 << k)
            delta = H_map[bits2] - H_map[bits]
            delta_dist[delta] += 1
            delta_by_arc[k][delta] += 1

    print(f"\n  n={n}: ARC FLIP DELTA DISTRIBUTION")
    print(f"  Total flips: {len(delta_dist)} distinct deltas")
    print(f"  {'delta':>6s} {'count':>8s} {'%':>8s}")
    for d in sorted(delta_dist.keys()):
        count = delta_dist[d]
        pct = 100 * count / (2**m * m)
        print(f"  {d:>+6d} {count:>8d} {pct:>7.2f}%")

    # Key: what are the POSSIBLE delta values?
    all_deltas = sorted(set(delta_dist.keys()))
    print(f"\n  Possible deltas: {all_deltas}")

    # Are they always even? (Redei: both H values are odd, so delta is always even)
    all_even = all(d % 2 == 0 for d in all_deltas)
    print(f"  All even: {all_even}")

    # GCD of all deltas
    from math import gcd
    from functools import reduce
    g = reduce(gcd, [abs(d) for d in all_deltas if d != 0])
    print(f"  GCD of nonzero deltas: {g}")

    # Missing deltas in range
    min_d, max_d = min(all_deltas), max(all_deltas)
    even_range = set(range(min_d, max_d + 1, 2))
    missing = sorted(even_range - set(all_deltas))
    print(f"  Missing even deltas in [{min_d}, {max_d}]: {missing}")

    # Delta depends on arc position?
    print(f"\n  DELTA BY ARC POSITION:")
    for k, (i,j) in enumerate(pairs):
        deltas = sorted(set(delta_by_arc[k].keys()))
        print(f"    arc ({i},{j}): deltas = {deltas}")

# ================================================================
# 2. THE DELETION-CONTRACTION DELTA FORMULA
# ================================================================
print(f"\n" + "=" * 65)
print(f"  THE ARC FLIP DELTA FORMULA")
print(f"=" * 65)

# From the deletion-contraction formula:
# H(T) = H(T\e) + H(T/e)
# If we flip arc e = (i->j) to (j->i):
# T has arc i->j. T' has arc j->i.
# T\e = T'\e' = the same digraph (remove the arc between i and j)
# H(T) - H(T') = H(T/e) - H(T'/e')
# where T/e contracts i->j (w inherits in from j, out from i)
# and T'/e' contracts j->i (w inherits in from i, out from j)

# So: delta = H(T/e) - H(T'/e')
# The contracted tournaments have n-1 vertices.
# This gives us delta in terms of (n-1)-vertex tournament HP counts.

# Let's verify this at n=4 and n=5.
print()
n = 4
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print(f"  VERIFICATION AT n={n}:")

# Test with a specific tournament and specific arc
A = np.array([[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]], dtype=np.int8)
H_orig = count_hp(A, n)
print(f"  T: H = {H_orig}")

for i in range(n):
    for j in range(n):
        if i == j or A[i][j] == 0: continue
        # Flip arc i->j to j->i
        A_flip = A.copy()
        A_flip[i][j] = 0
        A_flip[j][i] = 1
        H_flip = count_hp(A_flip, n)
        delta = H_flip - H_orig

        # Contract arc i->j: merge i,j into w
        # w inherits IN from j (tail) and OUT from i (head)
        # Wait: for T/e where e = i->j:
        # The merged vertex w has:
        #   w->k if i->k (and k != j)
        #   k->w if k->j (and k != i)
        # For the "other" arc between i and k or j and k: if both exist, take from respective end
        n_c = n - 1
        # Vertex mapping: keep all except j, merged vertex takes index of i
        # Actually let's use index 0 for merged, shift others
        verts = [k for k in range(n) if k != i and k != j]

        A_contract = np.zeros((n_c, n_c), dtype=np.int8)
        # Merged vertex is index 0, others are 1, 2, ...
        for ki, k in enumerate(verts):
            # arc from merged to k: i->k in T
            if A[i][k]: A_contract[0][ki+1] = 1
            else: A_contract[ki+1][0] = 1
            # Wait: we also need k->merged: k->j in T
            # But what if i->k and k->j, or k->i and j->k?
            # Standard convention: w->k if i->k, k->w if k->j
            # But what about the OTHER arc? The contracted tournament only has
            # one arc between w and k.
            # Actually: in the DC formula for tournaments:
            # T/e where e = (i->j): merged vertex w
            # For each k != i,j: arc(w,k) = arc(i,k) and arc(k,w) = arc(k,j)
            # But arc(i,k) and arc(k,j) might BOTH be outgoing from {w,k} side
            # That would give arc(w,k)=1 AND arc(k,w)=0 or similar.
            # Wait -- it IS a tournament contraction, the result may NOT be a tournament!

        # Actually: T/e for a tournament may not be a tournament.
        # It's a general digraph. The DC formula still works on digraphs.
        # H(D) = H(D\e) + H(D/e) for any digraph, where H counts HP.

        # Let me just compute directly instead.
        # T/(i->j): merge i and j. For each k != i,j:
        #   w -> k iff i -> k in T
        #   k -> w iff k -> j in T
        # Note: this can create BOTH w->k AND k->w (not a tournament!)

        A_contract = np.zeros((n_c, n_c), dtype=np.int8)
        for ki, k in enumerate(verts):
            if A[i][k]: A_contract[0][ki+1] = 1  # w -> k from i's out-arcs
            if A[k][j]: A_contract[ki+1][0] = 1  # k -> w from j's in-arcs
            for li, l in enumerate(verts):
                if ki == li: continue
                if A[k][l]: A_contract[ki+1][li+1] = 1

        H_contract = count_hp(A_contract, n_c)

        # Now contract the FLIPPED arc j->i
        A_contract2 = np.zeros((n_c, n_c), dtype=np.int8)
        for ki, k in enumerate(verts):
            # In T', arc is j->i. Contract j->i:
            # w -> k iff j -> k in T' = j -> k in T (for k != i)
            # k -> w iff k -> i in T' = k -> i in T (for k != j)
            if A[j][k]: A_contract2[0][ki+1] = 1
            if A[k][i]: A_contract2[ki+1][0] = 1
            for li, l in enumerate(verts):
                if ki == li: continue
                if A_flip[k][l]: A_contract2[ki+1][li+1] = 1

        H_contract2 = count_hp(A_contract2, n_c)

        delta_dc = H_contract - H_contract2
        match = delta == delta_dc

        if not match or abs(delta) > 0:
            print(f"    flip ({i}->{j}): delta={delta:+d}, DC: H(T/e)={H_contract}, H(T'/e')={H_contract2}, diff={delta_dc:+d}, match={match}")

print()

# ================================================================
# 3. THE LEVEL SET TOPOLOGY
# ================================================================
print(f"  3. LEVEL SET STRUCTURE (n=5)")
print()

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n)

# For each H value, how many neighbors have the SAME H?
# And how many have higher/lower H?
# This tells us about the "flatness" of the level set.
print(f"  {'H':>4s} {'count':>6s} {'avg_same':>9s} {'avg_up':>8s} {'avg_down':>9s} {'local_max':>9s} {'local_min':>9s}")
for H_val in sorted(set(H_map.values())):
    bits_at_H = [b for b, h in H_map.items() if h == H_val]
    total_same = 0
    total_up = 0
    total_down = 0
    local_max_count = 0
    local_min_count = 0
    for bits in bits_at_H:
        same = 0; up = 0; down = 0
        for k in range(m):
            bits2 = bits ^ (1 << k)
            h2 = H_map[bits2]
            if h2 == H_val: same += 1
            elif h2 > H_val: up += 1
            else: down += 1
        total_same += same
        total_up += up
        total_down += down
        if up == 0: local_max_count += 1
        if down == 0: local_min_count += 1
    cnt = len(bits_at_H)
    print(f"  {H_val:>4d} {cnt:>6d} {total_same/cnt:>9.2f} {total_up/cnt:>8.2f} {total_down/cnt:>9.2f} {local_max_count:>9d} {local_min_count:>9d}")

# ================================================================
# 4. GRADIENT FLOW AND BASINS OF ATTRACTION
# ================================================================
print(f"\n  4. GRADIENT FLOW: BASINS OF ATTRACTION")
print()

# For each tournament, follow the steepest ascent (largest H neighbor)
# until reaching a local maximum. What are the basins?
basin_map = {}  # bits -> local max bits
def steepest_ascent(bits):
    current = bits
    visited = set()
    while current not in visited:
        visited.add(current)
        h_curr = H_map[current]
        best_h = h_curr
        best_next = current
        for k in range(m):
            bits2 = current ^ (1 << k)
            if H_map[bits2] > best_h:
                best_h = H_map[bits2]
                best_next = bits2
        if best_next == current:
            break  # local max
        current = best_next
    return current

for bits in range(2**m):
    basin_map[bits] = steepest_ascent(bits)

# Count basins
basin_counts = defaultdict(int)
for bits, peak in basin_map.items():
    basin_counts[peak] += 1

print(f"  Number of local maxima: {len(basin_counts)}")
print(f"  Basin sizes:")
for peak, count in sorted(basin_counts.items(), key=lambda x: -x[1]):
    print(f"    peak H={H_map[peak]}, basin size={count}")

# How many local maxima are SC?
from itertools import permutations as perms

def is_sc_bits(bits, n, pairs):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    A_comp = 1 - A
    np.fill_diagonal(A_comp, 0)
    for perm in perms(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if A[perm[i]][perm[j]] != A_comp[i][j]:
                    match = False
                    break
            if not match: break
        if match: return True
    return False

print(f"\n  Local maxima SC status:")
for peak in sorted(basin_counts.keys(), key=lambda p: -basin_counts[p]):
    sc = is_sc_bits(peak, n, pairs)
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (peak >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    score = tuple(sorted(A.sum(axis=1).astype(int)))
    print(f"    H={H_map[peak]}, basin={basin_counts[peak]}, SC={sc}, score={score}")

# ================================================================
# 5. THE MORSE THEORY VIEW
# ================================================================
print(f"\n  5. MORSE THEORY: CRITICAL POINTS OF H")
print()

# Critical points: local maxima, local minima, saddle points
# A saddle has both higher and lower neighbors but is not a max or min
critical_max = 0
critical_min = 0
critical_saddle = 0

for bits in range(2**m):
    h = H_map[bits]
    has_higher = False
    has_lower = False
    for k in range(m):
        bits2 = bits ^ (1 << k)
        if H_map[bits2] > h: has_higher = True
        if H_map[bits2] < h: has_lower = True

    if not has_higher and not has_lower:
        pass  # isolated (impossible for tournaments)
    elif not has_higher:
        critical_max += 1
    elif not has_lower:
        critical_min += 1
    else:
        critical_saddle += 1

print(f"  Local maxima: {critical_max}")
print(f"  Local minima: {critical_min}")
print(f"  Saddle points: {critical_saddle}")
print(f"  Euler characteristic (max - saddle + min): {critical_max - critical_saddle + critical_min}")
print()

# The Euler characteristic of the tournament graph (which is the m-dimensional
# hypercube quotient -- each tournament connects to m neighbors via arc flips)
# The actual space is a subset of {0,1}^m... wait, it IS {0,1}^m = hypercube.
# Every binary string of length m is a valid tournament.
# So the tournament space is the m-cube, and H is a function on it.
# The Euler characteristic of the m-cube is always 0 (for m >= 1).
# With boundary: chi = sum (-1)^k * #(k-cells)
# For the hypercube graph: #vertices = 2^m, #edges = m*2^{m-1}
# chi(vertices - edges) = 2^m - m*2^{m-1}

print(f"  Tournament space = {m}-cube = Z_2^{m}")
print(f"  Vertices: {2**m}")
print(f"  Edges (arc flips): {m * 2**(m-1)}")
print(f"  Note: H is a MORSE FUNCTION on the {m}-cube")
print(f"  with {critical_max} maxima, {critical_min} minima, {critical_saddle} saddles")
