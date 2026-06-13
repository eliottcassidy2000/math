#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Can Ω(T) have a claw at n=9?

The claw K_{1,3} requires:
- Center C₀: an odd cycle sharing ≥1 vertex with each of C₁, C₂, C₃
- Leaves C₁, C₂, C₃: pairwise vertex-disjoint odd cycles

Minimum configuration: C₁, C₂, C₃ = three disjoint 3-cycles (9 vertices).
C₀ = a 3-cycle using one vertex from each C_i.

Construct such a tournament and check:
1. Does the claw actually appear?
2. Is the independence polynomial still real-rooted?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import random

def find_all_directed_odd_cycles(T):
    n = len(T)
    cycles = set()
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(L):
                    if T[path[i]][path[(i+1) % L]] != 1:
                        ok = False
                        break
                if ok:
                    mi = path.index(min(path))
                    canon = tuple(path[(mi + j) % L] for j in range(L))
                    cycles.add(canon)
                    break
    return list(cycles)

def build_conflict_graph(cycles):
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vsets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def has_claw(adj):
    k = len(adj)
    for center in range(k):
        neighbors = [j for j in range(k) if adj[center][j]]
        if len(neighbors) < 3:
            continue
        for triple in combinations(neighbors, 3):
            a, b, c = triple
            if not adj[a][b] and not adj[a][c] and not adj[b][c]:
                return True, center, triple
    return False, None, None

def independence_polynomial_coeffs(adj):
    k = len(adj)
    if k == 0: return [1]
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
    return coeffs

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

# ═══════════════════════════════════════════
# CONSTRUCT a tournament on 9 vertices with a claw in Ω(T)
# ═══════════════════════════════════════════
print("="*60)
print("Constructing tournament with claw in Ω(T) at n=9")
print("="*60)

# Target cycles:
# C₁ = (0,1,2): 0→1→2→0
# C₂ = (3,4,5): 3→4→5→3
# C₃ = (6,7,8): 6→7→8→6
# C₀ = (0,3,6): 0→3→6→0
#
# For a claw: C₁,C₂,C₃ pairwise disjoint (they are!),
# C₀ shares vertex 0 with C₁, vertex 3 with C₂, vertex 6 with C₃.

n = 9
T = [[0]*n for _ in range(n)]

# Set cycle arcs
# C₁: 0→1→2→0
T[0][1] = T[1][2] = T[2][0] = 1
# C₂: 3→4→5→3
T[3][4] = T[4][5] = T[5][3] = 1
# C₃: 6→7→8→6
T[6][7] = T[7][8] = T[8][6] = 1
# C₀: 0→3→6→0
T[0][3] = T[3][6] = T[6][0] = 1

# Now set remaining arcs to minimize extra odd cycles
# Strategy: make the remaining structure as "transitive" as possible
# Cross-group arcs: group A = {0,1,2}, B = {3,4,5}, C = {6,7,8}
# Want: A→B→C→A... wait, that creates more cycles.
# Better: A beats B and C, B beats C (transitive between groups).
# Except for the arcs already set.

# Set remaining intra-group arcs (complete each triple)
# In C₁ (0,1,2): 0→1✓, 1→2✓, 2→0✓. Missing: 0→2? 2→0 is set, so 0→2=0.
# Actually 0↔2: 2→0=1 (from cycle), so 0→2=0 (T[0][2]=0, T[2][0]=1) ✓
# Missing arc: 1↔0: already set (T[0][1]=1, T[1][0]=0). OK.
# All intra-C₁ arcs are determined by the cycle.

# Wait: in a tournament on {0,1,2}, cycle 0→1→2→0 determines ALL 3 arcs:
# T[0][1]=1, T[1][2]=1, T[2][0]=1. The reverses: T[1][0]=0, T[2][1]=0, T[0][2]=0.

# Similarly for C₂ and C₃.

# C₀ arcs: T[0][3]=1, T[3][6]=1, T[6][0]=1.
# Reverses: T[3][0]=0, T[6][3]=0, T[0][6]=0.

# Remaining: all arcs between the groups that aren't yet set.
# Between A={0,1,2} and B={3,4,5}:
# 0↔3: T[0][3]=1 ✓
# 0↔4, 0↔5, 1↔3, 1↔4, 1↔5, 2↔3, 2↔4, 2↔5
# Between A and C: 0↔6: T[6][0]=1 (so T[0][6]=0) ✓
# 0↔7, 0↔8, 1↔6, 1↔7, 1↔8, 2↔6, 2↔7, 2↔8
# Between B and C: 3↔6: T[3][6]=1 ✓
# 3↔7, 3↔8, 4↔6, 4↔7, 4↔8, 5↔6, 5↔7, 5↔8

# Try: make A beat B (except 0→3 already set, and 3→0=0)
# A beats B: {1,2} beat {3,4,5}
remaining_AB = [(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(0,4),(0,5)]
remaining_AC = [(0,7),(0,8),(1,6),(1,7),(1,8),(2,6),(2,7),(2,8)]
remaining_BC = [(3,7),(3,8),(4,6),(4,7),(4,8),(5,6),(5,7),(5,8)]

# Try all A→B, B→C
for i,j in remaining_AB:
    T[i][j] = 1  # A beats B
for i,j in remaining_BC:
    T[i][j] = 1  # B beats C
for i,j in remaining_AC:
    T[j][i] = 1  # C beats A (to avoid creating cycles? Actually let's try A→C)
    # Wait, C₀ has 6→0 (C beats A at vertex 0). Let's try C→A for consistency.

# Check tournament property
for i in range(n):
    for j in range(i+1, n):
        assert T[i][j] + T[j][i] == 1, f"Not tournament: ({i},{j}): {T[i][j]}, {T[j][i]}"

print("Tournament constructed. Checking cycles...")

cycles = find_all_directed_odd_cycles(T)
print(f"Total odd cycles: {len(cycles)}")
for c in sorted(cycles, key=lambda x: (len(x), x)):
    print(f"  {c} (length {len(c)})")

adj = build_conflict_graph(cycles)
result = has_claw(adj)
print(f"\nClaw found: {result[0]}")
if result[0]:
    print(f"  Center: cycle {result[1]} = {cycles[result[1]]}")
    print(f"  Leaves: cycles {result[2]} = {[cycles[i] for i in result[2]]}")

# Check real-rootedness
coeffs = independence_polynomial_coeffs(adj)
print(f"\nIndependence polynomial: {coeffs}")
H = sum(c * 2**k for k, c in enumerate(coeffs))
H_dp = hamiltonian_path_count(T)
print(f"H = {H} (OCF), H = {H_dp} (DP)")

if len(coeffs) > 1:
    np_coeffs = list(reversed(coeffs))
    zeros = np.roots(np_coeffs)
    real_zeros = [z for z in zeros if abs(z.imag) < 1e-8]
    complex_zeros = [z for z in zeros if abs(z.imag) >= 1e-8]
    print(f"Zeros: {len(real_zeros)} real, {len(complex_zeros)} complex")
    for z in sorted(zeros, key=lambda x: (abs(x.imag) > 1e-8, x.real)):
        print(f"  {z.real:.6f} + {z.imag:.6f}i")
    all_real = len(complex_zeros) == 0
    print(f"ALL REAL: {all_real}")

# ═══════════════════════════════════════════
# Try many random n=9 tournaments to find claws and check real-rootedness
# ═══════════════════════════════════════════
print(f"\n{'='*60}")
print("Random n=9 tournament sampling")
print("="*60)

def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

n_claw = 0
n_not_real = 0
n_sampled = 0
alpha3_dist = Counter()
target_samples = 500

for trial in range(target_samples):
    T = random_tournament(9)
    cycles = find_all_directed_odd_cycles(T)
    if len(cycles) < 4:
        n_sampled += 1
        continue

    adj = build_conflict_graph(cycles)
    result = has_claw(adj)
    if result[0]:
        n_claw += 1

    # Check real-rootedness
    coeffs = independence_polynomial_coeffs(adj)
    alpha3 = coeffs[3] if len(coeffs) > 3 else 0
    alpha3_dist[alpha3 > 0] += 1

    if len(coeffs) > 1:
        np_coeffs = list(reversed(coeffs))
        zeros = np.roots(np_coeffs)
        all_real = all(abs(z.imag) < 1e-6 for z in zeros)
        if not all_real:
            n_not_real += 1
            if n_not_real <= 3:
                print(f"  NON-REAL at trial {trial}! coeffs={coeffs[:8]}...")
                for z in zeros:
                    if abs(z.imag) >= 1e-6:
                        print(f"    complex zero: {z}")

    n_sampled += 1
    if trial % 100 == 0 and trial > 0:
        print(f"  ...{trial}/{target_samples}, claws={n_claw}, non-real={n_not_real}")

print(f"\nn=9 results ({n_sampled} sampled):")
print(f"  Claws in Ω(T): {n_claw} ({100*n_claw/n_sampled:.1f}%)")
print(f"  Non-real zeros: {n_not_real} ({100*n_not_real/n_sampled:.1f}%)")
print(f"  α₃ > 0: {alpha3_dist}")

print("\n\nDONE.")
