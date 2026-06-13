#!/usr/bin/env python3
"""
h21_structural.py — Why H=21 is impossible: structural analysis.

H=21 = 1 + 2·10 = 3·7 = Φ₂(2)·Φ₃(2).

The independence polynomial decomposition:
  I(Ω,2) = 21 requires exactly one of:
  (a) α₁=10, α₂=0  — 10 pairwise-intersecting cycles
  (b) α₁=8,  α₂=1  — 8 cycles, 1 disjoint pair
  (c) α₁=6,  α₂=2  — 6 cycles, 2 disjoint pairs
  (d) α₁=4,  α₂=3  — 4 cycles, 3 disjoint pairs
  (e) α₁=6,  α₂=0, α₃=1 — 6 cycles, 1 independent triple, no pairs

For each: can this (α₁,α₂) pattern actually occur in a tournament?
If NO pattern can, then H=21 is impossible.

Key insight: the problem reduces to checking whether certain
independence polynomial values are achievable.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ═══════════════════════════════════════════════════════════════════
# Part 1: What (α₁,α₂) values occur at each n?
# ═══════════════════════════════════════════════════════════════════
print("=" * 70)
print("H=21 STRUCTURAL ANALYSIS")
print("=" * 70)

def fast_hp(A, n):
    """DP Hamiltonian path count."""
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

def get_cycle_sets(A, n):
    """Get odd-cycle vertex sets."""
    cs = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    cs.add(frozenset(verts))
                    break
    return cs

def compute_ip_coeffs(cycle_sets):
    """Compute independence polynomial coefficients of intersection graph."""
    cs_list = list(cycle_sets)
    n_cs = len(cs_list)
    if n_cs == 0:
        return {0: 1}

    # Build adjacency
    adj = [[False]*n_cs for _ in range(n_cs)]
    for i in range(n_cs):
        for j in range(i+1, n_cs):
            if cs_list[i] & cs_list[j]:
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size
    coeffs = defaultdict(int)
    for mask in range(2**n_cs):
        verts = [i for i in range(n_cs) if mask & (1<<i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(verts)] += 1
    return dict(coeffs)

# ═══════════════════════════════════════════════════════════════════
# Part 2: Map (α₁,α₂) → H at n=5,6
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: (α₁,α₂) landscape ---")

for n in [5, 6]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    ne = len(edges)

    a1a2_to_H = defaultdict(set)
    H_to_a1a2 = defaultdict(set)

    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        cs = get_cycle_sets(A, n)
        a1 = len(cs)

        # Count disjoint pairs
        cs_list = list(cs)
        a2 = sum(1 for i in range(len(cs_list)) for j in range(i+1, len(cs_list))
                 if not (cs_list[i] & cs_list[j]))

        H = fast_hp(A, n)
        a1a2_to_H[(a1,a2)].add(H)
        H_to_a1a2[H].add((a1,a2))

    print(f"\n  n={n}:")

    # Show H=21 decomps and whether their (α₁,α₂) are achievable
    decomps_21 = [(10,0), (8,1), (6,2), (4,3)]
    for a1, a2 in decomps_21:
        h_vals = sorted(a1a2_to_H.get((a1,a2), set()))
        achieved = (a1,a2) in a1a2_to_H
        print(f"    (α₁,α₂)=({a1},{a2}): {'ACHIEVED' if achieved else 'NOT ACHIEVED'}, H values: {h_vals}")

    # Show what H values close to 21 map to
    for h in [17, 19, 21, 23, 25]:
        patterns = sorted(H_to_a1a2.get(h, set()))
        if patterns:
            print(f"    H={h} comes from (α₁,α₂) ∈ {patterns}")
        else:
            print(f"    H={h}: NOT ACHIEVED")

# ═══════════════════════════════════════════════════════════════════
# Part 3: What (α₁,α₂) actually exist at n=7? Focus on target range.
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: (α₁,α₂) at n=7 (exhaustive, no HP needed) ---")

n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

a1a2_count = defaultdict(int)
target_decomps = {(10,0), (8,1), (6,2), (4,3), (6,0), (2,0)}

print(f"  Scanning all {2**ne} tournaments at n=7 for (α₁,α₂) near H=21 targets...")

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Only count 3-cycles for speed (α₁ ≥ dc3)
    dc3 = 0
    for v0, v1, v2 in combinations(range(n), 3):
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            dc3 += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            dc3 += 1

    # Need α₁ in {10,8,6,4} for H=21 decomps.
    # dc3 alone gives lower bound. At n=7, we also have 5-cycles and 7-cycles.
    # Skip if dc3 > 12 (too many 3-cycles, H will be >> 21)
    if dc3 > 12:
        continue
    # Skip if dc3 < 2 (too few cycles for α₁≥4)
    if dc3 < 2:
        continue

    cs = get_cycle_sets(A, n)
    a1 = len(cs)

    if a1 < 3 or a1 > 12:
        continue

    # Count disjoint pairs
    cs_list = list(cs)
    a2 = sum(1 for i in range(len(cs_list)) for j in range(i+1, len(cs_list))
             if not (cs_list[i] & cs_list[j]))

    a1a2_count[(a1, a2)] += 1

    if bits % 500000 == 0 and bits > 0:
        print(f"  Progress: {bits}/{2**ne}")

print(f"\n  (α₁,α₂) counts for α₁ ∈ [3,12]:")
for key in sorted(a1a2_count.keys()):
    if key in target_decomps or (key[0] in [4,6,8,10] and key[1] <= 5):
        marker = " ← H=21 TARGET" if key in target_decomps else ""
        print(f"    ({key[0]:2d},{key[1]:2d}): {a1a2_count[key]:6d}{marker}")

# Check which target decomps exist
print("\n  H=21 target decomposition status:")
for a1, a2 in [(10,0), (8,1), (6,2), (4,3)]:
    count = a1a2_count.get((a1,a2), 0)
    if count > 0:
        print(f"    (α₁,α₂)=({a1},{a2}): EXISTS ({count} tournaments)")
        print(f"      But does any of them give H=21? Need to check HP.")
    else:
        print(f"    (α₁,α₂)=({a1},{a2}): DOES NOT EXIST at n=7!")
