#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Fast claw check at n=9 using targeted construction.

Instead of brute-force cycle enumeration, construct a specific tournament
where we know the claw structure, then verify with H computation.
"""

import numpy as np
from itertools import combinations
from collections import Counter
import random

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def count_3cycles(T):
    n = len(T)
    count = 0
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            count += 1
        elif T[a][c] and T[c][b] and T[b][a]:
            count += 1
    return count

def count_3cycles_through_triple(T, a, b, c):
    """Check if {a,b,c} forms a 3-cycle."""
    if T[a][b] and T[b][c] and T[c][a]:
        return 1
    if T[a][c] and T[c][b] and T[b][a]:
        return 1
    return 0

def find_3cycles(T):
    """Find all directed 3-cycles."""
    n = len(T)
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            cycles.append((a, c, b))
    return cycles

def count_5cycles_on_verts(T, verts):
    """Count directed 5-cycles on given 5 vertices."""
    assert len(verts) == 5
    count = 0
    v = list(verts)
    from itertools import permutations
    for perm in permutations(v[1:]):
        path = [v[0]] + list(perm)
        ok = True
        for i in range(5):
            if T[path[i]][path[(i+1) % 5]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count  # each directed cycle found once (v[0] fixed)

# ═══════════════════════════════════════════
# CONSTRUCT: Tournament on 9 vertices with 4 specific 3-cycles
# ═══════════════════════════════════════════
print("="*60)
print("CONSTRUCTION: Tournament with claw in Ω(T) at n=9")
print("="*60)

# Groups: A={0,1,2}, B={3,4,5}, C={6,7,8}
# Target cycles:
# C₁ = 0→1→2→0 (within A)
# C₂ = 3→4→5→3 (within B)
# C₃ = 6→7→8→6 (within C)
# C₀ = 0→3→6→0 (cross-group, hitting one from each)

# Strategy to minimize extra 3-cycles:
# Within each group: the cycle determines all arcs (cyclic tournament)
# Between groups: make it "almost transitive" to avoid extra 3-cycles

# Cross-group: A→B, B→C, C→A ... but C₀ has 0→3, 3→6, 6→0
# which IS the C→A→B→C pattern for these specific vertices

# Key: avoid cross-group 3-cycles.
# A cross-group triple like (0, 4, 7) forms a 3-cycle if
# T[0][4]=1, T[4][7]=1, T[7][0]=1  (or the reverse direction)

# Let's try: for cross-group arcs not already determined by C₀,
# use A→B→C→A pattern consistently.

def construct_tournament():
    n = 9
    T = [[0]*n for _ in range(n)]

    # Intra-group: cyclic tournaments
    # A: 0→1→2→0
    T[0][1] = T[1][2] = T[2][0] = 1
    # B: 3→4→5→3
    T[3][4] = T[4][5] = T[5][3] = 1
    # C: 6→7→8→6
    T[6][7] = T[7][8] = T[8][6] = 1

    # C₀: 0→3→6→0
    T[0][3] = T[3][6] = T[6][0] = 1

    # Remaining cross-group arcs:
    # A-B (not yet set): (0,4), (0,5), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5)
    # A-C (not yet set): (0,7), (0,8), (1,6), (1,7), (1,8), (2,6), (2,7), (2,8)
    # B-C (not yet set): (3,7), (3,8), (4,6), (4,7), (4,8), (5,6), (5,7), (5,8)

    # Strategy: A→B (all remaining), B→C (all remaining), C→A (all remaining)
    # A→B:
    for a in [0, 1, 2]:
        for b in [3, 4, 5]:
            if T[a][b] == 0 and T[b][a] == 0:
                T[a][b] = 1
    # B→C:
    for b in [3, 4, 5]:
        for c in [6, 7, 8]:
            if T[b][c] == 0 and T[c][b] == 0:
                T[b][c] = 1
    # C→A:
    for c in [6, 7, 8]:
        for a in [0, 1, 2]:
            if T[c][a] == 0 and T[a][c] == 0:
                T[c][a] = 1

    # Verify tournament property
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j] + T[j][i] != 1:
                print(f"ERROR: ({i},{j}): T={T[i][j]}, T'={T[j][i]}")
                return None
    return T

T = construct_tournament()
print("Tournament constructed.")

# Verify the 4 target cycles exist
print("\nVerifying target cycles:")
for name, verts in [("C₁", (0,1,2)), ("C₂", (3,4,5)), ("C₃", (6,7,8)), ("C₀", (0,3,6))]:
    a, b, c = verts
    fwd = T[a][b] and T[b][c] and T[c][a]
    rev = T[a][c] and T[c][b] and T[b][a]
    print(f"  {name} ({a},{b},{c}): {'forward ✓' if fwd else 'reverse ✓' if rev else '✗'}")

# Count ALL 3-cycles
all_3cycles = find_3cycles(T)
print(f"\nTotal 3-cycles: {len(all_3cycles)}")
for cyc in all_3cycles:
    # Classify: intra-group or cross-group?
    groups = [0 if v < 3 else (1 if v < 6 else 2) for v in cyc]
    if len(set(groups)) == 1:
        gname = "intra-" + "ABC"[groups[0]]
    elif len(set(groups)) == 3:
        gname = "cross-ABC"
    else:
        gname = f"mixed-{groups}"
    print(f"  {cyc} [{gname}]")

# The claw structure: C₁, C₂, C₃ are pairwise disjoint (✓ by construction)
# C₀ shares vertex with each: 0∈C₁, 3∈C₂, 6∈C₃ (✓)
# Check if C₁, C₂, C₃ are pairwise disjoint in the CONFLICT GRAPH
# (no shared vertices): {0,1,2} ∩ {3,4,5} = ∅, etc.

# Build conflict graph for the 3-cycles found
vsets = [set(c) for c in all_3cycles]
k = len(all_3cycles)
adj = [[False]*k for _ in range(k)]
for i in range(k):
    for j in range(i+1, k):
        if vsets[i] & vsets[j]:
            adj[i][j] = adj[j][i] = True

# Check for claws
from itertools import combinations as combs
has_claw = False
for center in range(k):
    nbrs = [j for j in range(k) if adj[center][j]]
    if len(nbrs) < 3:
        continue
    for triple in combs(nbrs, 3):
        a, b, c = triple
        if not adj[a][b] and not adj[a][c] and not adj[b][c]:
            has_claw = True
            print(f"\nCLAW FOUND!")
            print(f"  Center: {all_3cycles[center]}")
            print(f"  Leaf 1: {all_3cycles[a]}")
            print(f"  Leaf 2: {all_3cycles[b]}")
            print(f"  Leaf 3: {all_3cycles[c]}")
            break
    if has_claw:
        break

if not has_claw:
    print("\nNo claw found among 3-cycles. Checking with more cycles...")
    # If no claw in 3-cycles only, the structure might need 5-cycles too.
    # But for the claw, we only need the 3 disjoint cycles + center.
    print("NOTE: The claw needs 3 pairwise disjoint cycles as leaves.")
    print(f"Max independent set among {k} 3-cycles:")

    # Find max independent set in conflict graph of 3-cycles
    max_indep = 0
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            if len(verts) > max_indep:
                max_indep = len(verts)
                if len(verts) >= 3:
                    print(f"  Found independent set of size {len(verts)}: {[all_3cycles[i] for i in verts]}")
    print(f"  Max independent set: {max_indep}")

# Independence polynomial (using only 3-cycles as vertices)
coeffs = [0] * (k + 1)
for mask in range(2**k):
    verts = [i for i in range(k) if (mask >> i) & 1]
    ok = True
    for a in range(len(verts)):
        for b in range(a+1, len(verts)):
            if adj[verts[a]][verts[b]]:
                ok = False
                break
        if not ok:
            break
    if ok:
        coeffs[len(verts)] += 1
while len(coeffs) > 1 and coeffs[-1] == 0:
    coeffs.pop()

print(f"\nIndependence polynomial of 3-cycle conflict graph: {coeffs}")

# Compute H
H = hamiltonian_path_count(T)
print(f"H(T) = {H}")

# Check if polynomial at x=2 gives H (it won't — need ALL odd cycles, not just 3-cycles)
H_3only = sum(c * 2**k for k, c in enumerate(coeffs))
print(f"I(Ω₃, 2) using only 3-cycles = {H_3only} (vs H={H})")
print(f"Difference = {H - H_3only} (from 5-cycles, 7-cycles, 9-cycle)")

# ═══════════════════════════════════════════
# Vary the cross-group arcs to minimize total 3-cycles
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("Varying cross-group arcs to minimize total c3")
print("="*60)

# The cross-group arcs that are free (not determined by C₀, C₁, C₂, C₃)
# We need to find a configuration that minimizes c3

# Fixed arcs from cycles:
# Intra: all 9 arcs (3 per group)
# C₀: T[0][3]=1, T[3][6]=1, T[6][0]=1

# Free arcs: cross-group pairs not in C₀
free_pairs = []
for a in range(9):
    for b in range(a+1, 9):
        ga, gb = a//3, b//3
        if ga == gb:
            continue  # intra-group, already set
        if (a,b) in [(0,3), (3,6)] or (b,a) in [(6,0)]:
            continue  # C₀ arcs
        # Check if this arc is determined by any cycle
        if T[a][b] != 0 or T[b][a] != 0:
            # Already set by the A→B→C→A pattern
            pass
        free_pairs.append((a, b))

# Actually, all cross-group arcs are set by the construction.
# Let's count the extra 3-cycles with this specific configuration.
c3 = count_3cycles(T)
print(f"Total 3-cycles with A→B→C→A configuration: {c3}")

# Try randomizing some cross-group arcs to reduce c3
best_c3 = c3
best_extra_cycles = all_3cycles

# The mandatory arcs: the 4 target cycles
mandatory = {
    (0,1): 1, (1,2): 1, (2,0): 1,  # C₁
    (3,4): 1, (4,5): 1, (5,3): 1,  # C₂
    (6,7): 1, (7,8): 1, (8,6): 1,  # C₃
    (0,3): 1, (3,6): 1, (6,0): 1,  # C₀
}

# Try 10000 random cross-group configurations
for trial in range(10000):
    T2 = [[0]*9 for _ in range(9)]
    # Set mandatory arcs
    for (a,b), val in mandatory.items():
        T2[a][b] = val
        T2[b][a] = 1 - val

    # Randomize remaining cross-group arcs
    for a in range(9):
        for b in range(a+1, 9):
            if a//3 == b//3:
                continue  # intra-group already set
            if T2[a][b] != 0 or T2[b][a] != 0:
                continue  # already set by mandatory
            if random.random() < 0.5:
                T2[a][b] = 1
            else:
                T2[b][a] = 1

    c3_2 = count_3cycles(T2)
    if c3_2 < best_c3:
        best_c3 = c3_2
        # Check that the 4 target cycles still exist and nothing else
        cycles_2 = find_3cycles(T2)
        best_extra_cycles = cycles_2
        H2 = hamiltonian_path_count(T2)
        print(f"  Trial {trial}: c3={c3_2}, H={H2}, "
              f"3-cycles: {[(c, 'target' if set(c).issubset({0,1,2}) or set(c).issubset({3,4,5}) or set(c).issubset({6,7,8}) or set(c)=={0,3,6} else 'EXTRA') for c in cycles_2]}")

        # Is it possible to get EXACTLY c3=4 (just the 4 targets)?
        if c3_2 == 4:
            print(f"  *** PERFECT: exactly 4 three-cycles! ***")
            print(f"  H = {H2}")
            # Check claw
            vsets2 = [set(c) for c in cycles_2]
            k2 = len(cycles_2)
            adj2 = [[False]*k2 for _ in range(k2)]
            for i in range(k2):
                for j in range(i+1, k2):
                    if vsets2[i] & vsets2[j]:
                        adj2[i][j] = adj2[j][i] = True
            # Independence polynomial
            coeffs2 = [0] * (k2 + 1)
            for mask in range(2**k2):
                vs = [i for i in range(k2) if (mask >> i) & 1]
                ok = True
                for a in range(len(vs)):
                    for b in range(a+1, len(vs)):
                        if adj2[vs[a]][vs[b]]:
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    coeffs2[len(vs)] += 1
            while len(coeffs2) > 1 and coeffs2[-1] == 0:
                coeffs2.pop()
            print(f"  I(Ω₃, x) = {coeffs2}")
            print(f"  I(Ω₃, 2) = {sum(c*2**kk for kk,c in enumerate(coeffs2))} (but H={H2}, so 5+cycles contribute)")
            break

print(f"\nBest c3 found: {best_c3}")

print("\n\nDONE.")
