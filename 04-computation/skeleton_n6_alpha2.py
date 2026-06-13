#!/usr/bin/env python3
"""
Tiling skeleton at n=6: the first n where alpha_2 > 0 is possible.

At n=6, vertex-disjoint 3-cycles can exist (e.g., {0,1,2} and {3,4,5}).
This means I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
and H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

The tiling model at n=6 has C(6,2) - 5 = 10 non-path arcs,
giving 2^10 = 1024 tilings.

KEY QUESTIONS:
1. How many isomorphism classes are there?
2. What is the distribution of (alpha_1, alpha_2, ...)?
3. Do tilings with alpha_2 > 0 form a connected subgraph?
4. How does the transfer matrix behave at even n=6?
   (tr(M) = 0, H = Sigma/2 = off-diag sum / 2)
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np
import sys

# DP-based Hamiltonian path count (much faster than brute force)
def ham_path_count_dp(A):
    """Count Hamiltonian paths using DP (bitmask)."""
    n = len(A)
    # dp[mask][v] = number of Hamiltonian paths through vertices in mask ending at v
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]

    full_mask = (1 << n) - 1
    return sum(dp[full_mask][v] for v in range(n))

def ham_path_count_dp_endpoints(A):
    """Count Hamiltonian paths by start and end vertex."""
    n = len(A)
    dp = [[0] * n for _ in range(1 << n)]
    # dp[mask][v] = paths through mask ending at v, starting at the first set bit
    # Better: dp_start[s][mask][v] = paths starting at s through mask ending at v
    # More efficient: just track (mask, end) and sum over starts at the end

    for s in range(n):
        dp[1 << s][s] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]

    # Now dp[(1<<n)-1][v] = total paths ending at v (from any start)
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

def tournament_canonical(A):
    """Canonical form via score sequence + partial comparison.
    Full canonicalization via all n! permutations is too slow for n=6.
    Use score sequence + sorted adjacency hash."""
    n = len(A)
    # Simple hash: sorted row sums + sorted adjacency fingerprint
    scores = tuple(sorted([sum(row) for row in A]))
    # For a more precise canonical, hash the adjacency matrix
    # under the sorted-score ordering
    return (scores, tuple(tuple(row) for row in A))

def canonical_full(A):
    """Full canonical form (expensive). Only for small batches."""
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def directed_odd_cycles_fast(A):
    """Find directed odd cycles using efficient enumeration."""
    n = len(A)
    cycles = []

    # 3-cycles
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if A[i][j] != 1: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] == 1 and A[k][i] == 1:
                    # i -> j -> k -> i
                    canonical = min((i,j,k), (j,k,i), (k,i,j))
                    if canonical not in cycles:
                        cycles.append(canonical)

    # 5-cycles
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or A[v0][v1] != 1: continue
            for v2 in range(n):
                if v2 in (v0,v1) or A[v1][v2] != 1: continue
                for v3 in range(n):
                    if v3 in (v0,v1,v2) or A[v2][v3] != 1: continue
                    for v4 in range(n):
                        if v4 in (v0,v1,v2,v3) or A[v3][v4] != 1: continue
                        if A[v4][v0] == 1:
                            cycle = (v0,v1,v2,v3,v4)
                            canonical = min(cycle[i:]+cycle[:i] for i in range(5))
                            if canonical not in cycles:
                                cycles.append(canonical)

    return cycles

def independence_poly(cycles):
    """Compute independence polynomial coefficients of cycle conflict graph."""
    m = len(cycles)
    if m == 0:
        return {0: 1}

    # Build adjacency
    adj = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if set(cycles[i]) & set(cycles[j]):
                adj[i].add(j)
                adj[j].add(i)

    # Count independent sets by size (for small m, brute force is fine)
    if m > 20:
        # Too many cycles, just compute alpha_0 and alpha_1
        return {0: 1, 1: m}

    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if verts[j] in adj[verts[i]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(verts)] += 1

    return dict(alpha)

n = 6
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)
num_tilings = 2**m

print(f"n={n}: {m} tiles, {num_tilings} tilings")
print(f"Tiles: {tiles}")

# Compute data for all tilings
print("\nComputing H and odd cycles for all tilings...")
tiling_data = {}

for bits in range(num_tilings):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)

    # Verify OCF
    I_at_2 = sum(alphas.get(k, 0) * 2**k for k in alphas)
    if H != I_at_2:
        print(f"  WARNING: OCF mismatch at bits={format(bits, f'0{m}b')}: H={H}, I(2)={I_at_2}")

    tiling_data[bits] = {
        'H': H, 'cycles': cycles, 'alphas': alphas,
        'bits': format(bits, f'0{m}b')
    }

print("Done.\n")

# =====================================================================
# H distribution
# =====================================================================
print("=" * 70)
print("H DISTRIBUTION")
print("=" * 70)

H_dist = defaultdict(int)
for bits in range(num_tilings):
    H_dist[tiling_data[bits]['H']] += 1

print(f"\n  H   Count")
print("  " + "-" * 20)
for H in sorted(H_dist.keys()):
    print(f"  {H:3d}  {H_dist[H]:4d}")

# =====================================================================
# Alpha distribution
# =====================================================================
print()
print("=" * 70)
print("ALPHA DISTRIBUTION (alpha_k for k >= 1)")
print("=" * 70)

alpha_dist = defaultdict(int)
max_alpha = 0
for bits in range(num_tilings):
    d = tiling_data[bits]
    # Key by (alpha_1, alpha_2, ...)
    a1 = d['alphas'].get(1, 0)
    a2 = d['alphas'].get(2, 0)
    a3 = d['alphas'].get(3, 0)
    alpha_dist[(a1, a2, a3)] += 1
    if a2 > 0:
        max_alpha = max(max_alpha, a2)

print(f"\n  (a1, a2, a3)    Count   H")
print("  " + "-" * 40)
for key in sorted(alpha_dist.keys()):
    a1, a2, a3 = key
    H = 1 + 2*a1 + 4*a2 + 8*a3
    count = alpha_dist[key]
    print(f"  ({a1:2d}, {a2:2d}, {a3:2d})     {count:4d}    H={H}")

print(f"\n  Max alpha_2 = {max_alpha}")
has_alpha2 = sum(c for k, c in alpha_dist.items() if k[1] > 0)
print(f"  Tilings with alpha_2 > 0: {has_alpha2} / {num_tilings}")

# =====================================================================
# The tilings with alpha_2 > 0: examine their structure
# =====================================================================
print()
print("=" * 70)
print("TILINGS WITH VERTEX-DISJOINT ODD CYCLES (alpha_2 > 0)")
print("=" * 70)

count_shown = 0
for bits in range(num_tilings):
    d = tiling_data[bits]
    a2 = d['alphas'].get(2, 0)
    if a2 > 0 and count_shown < 5:
        count_shown += 1
        a1 = d['alphas'].get(1, 0)
        print(f"\n  Tiling {d['bits']} (H={d['H']}, a1={a1}, a2={a2}):")
        print(f"    Cycles: {[list(c) for c in d['cycles']]}")

        # Find pairs of vertex-disjoint cycles
        cycles = d['cycles']
        nc = len(cycles)
        disjoint_pairs = []
        for i in range(nc):
            for j in range(i+1, nc):
                if not (set(cycles[i]) & set(cycles[j])):
                    disjoint_pairs.append((i, j))
        print(f"    Vertex-disjoint pairs:")
        for i, j in disjoint_pairs:
            print(f"      {list(cycles[i])} and {list(cycles[j])}")

# =====================================================================
# Self-flip pairs at even n=6 (blueself possible!)
# =====================================================================
print()
print("=" * 70)
print("SELF-FLIP PAIRS AT EVEN n=6")
print("=" * 70)

# At even n, blueself pairs (T isomorphic to flip(T)) are possible
# We can't compute full canonical forms for all 1024 tilings at n=6
# (6! = 720 permutations per tiling, 1024 tilings = 737,280 comparisons)
# But we can check score sequences and H values as a fast filter

self_flip_same_H = 0
for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits > bits:
        H1 = tiling_data[bits]['H']
        H2 = tiling_data[flip_bits]['H']
        if H1 == H2:
            self_flip_same_H += 1

print(f"\n  Blue pairs (T, flip(T)) with same H: {self_flip_same_H}")

# Show the H-sum distribution for blue pairs
H_sum_dist = defaultdict(int)
for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits >= bits:
        H_sum = tiling_data[bits]['H'] + tiling_data[flip_bits]['H']
        H_sum_dist[H_sum] += 1

print(f"\n  H(T) + H(flip(T)) distribution:")
for s in sorted(H_sum_dist.keys()):
    print(f"    sum={s}: {H_sum_dist[s]} pairs")

# =====================================================================
# Skeleton edges: delta_H distribution
# =====================================================================
print()
print("=" * 70)
print("DELTA_H DISTRIBUTION UNDER SINGLE TILE FLIPS")
print("=" * 70)

delta_H_dist = defaultdict(int)
for bits in range(num_tilings):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            dH = tiling_data[neighbor]['H'] - tiling_data[bits]['H']
            delta_H_dist[dH] += 1

print(f"\n  delta_H  Count")
print("  " + "-" * 20)
for dH in sorted(delta_H_dist.keys()):
    print(f"  {dH:+5d}  {delta_H_dist[dH]:5d}")

# Check: is delta_H always even? (Redei)
all_even = all(dH % 2 == 0 for dH in delta_H_dist.keys())
print(f"\n  All delta_H even? {all_even}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
n=6 tiling skeleton:
  - {num_tilings} tilings, {m} tiles per tiling
  - H ranges from {min(H_dist.keys())} to {max(H_dist.keys())}
  - alpha_2 > 0 for {has_alpha2}/{num_tilings} tilings
  - Max alpha_2 = {max_alpha}
  - delta_H always even: {all_even}

At n=6, the independence polynomial I(Omega, x) can have degree >= 2,
meaning vertex-disjoint odd cycles exist. This is the first n where
the FULL independence polynomial structure is relevant for H.
""")
