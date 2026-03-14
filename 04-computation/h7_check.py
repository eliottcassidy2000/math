#!/usr/bin/env python3
"""
h7_check.py — opus-2026-03-14-S71e

CHECK: Is H=7 actually achievable at n=7?

HYP-992 claims H=7 is never achievable (sampling only).
H=7 requires α₁+2α₂=3.

Possible: α₁=3, α₂=0 (3 odd cycles, none disjoint)
          α₁=1, α₂=1 (1 cycle + 1 disjoint pair... wait, that's 2 cycles)

Actually: I(x) = 1 + α₁x + α₂x² + α₃x³
H = I(2) = 1 + 2α₁ + 4α₂ + 8α₃ = 7
→ 2α₁ + 4α₂ + 8α₃ = 6
→ α₁ + 2α₂ + 4α₃ = 3

Only non-negative integer solutions:
  (α₁, α₂, α₃) = (3, 0, 0) or (1, 1, 0)

For (3,0,0): 3 odd cycles, no two disjoint.
  Need: dc3+dc5+dc7=3, all sharing a common vertex.

For (1,1,0): 1 odd cycle, and also a disjoint pair somehow?
  Wait: α₁=1 means 1 odd cycle total. α₂=1 means 1 disjoint pair.
  But you need AT LEAST 2 cycles to have a disjoint pair.
  If α₁=1 and α₂=1: impossible! α₂ ≤ C(α₁,2) = 0.
  So (1,1,0) is impossible.

Only (3,0,0) could give H=7.

Question: can 3 directed 3-cycles share a common vertex at n=7?
"""

import sys
from itertools import combinations, permutations
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def count_dc3(adj, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_dc5(adj, n):
    count = 0
    for verts in combinations(range(n), 5):
        v = list(verts)
        for perm in permutations(v):
            if all(adj[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                count += 1
    return count // 5

def count_dc7(adj, n):
    """Count directed 7-cycles (Hamiltonian cycles for n=7)."""
    count = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1) % n]] for i in range(n)):
            count += 1
    return count // n

def count_hp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
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
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def find_3cycles(adj, n):
    """Find all directed 3-cycles as vertex triples."""
    cycles = []
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycle = tuple(sorted([i,j,k]))
                    if cycle not in [tuple(sorted(c)) for c in cycles]:
                        cycles.append((i,j,k))
    # Deduplicate by vertex set
    seen = set()
    unique = []
    for c in cycles:
        key = tuple(sorted(c))
        if key not in seen:
            seen.add(key)
            unique.append(c)
    return unique

def alpha2_from_cycles(cycles):
    """Count vertex-disjoint pairs of cycles."""
    count = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            verts_i = set(cycles[i])
            verts_j = set(cycles[j])
            if verts_i.isdisjoint(verts_j):
                count += 1
    return count

print("=" * 70)
print("CAN 3 DIRECTED 3-CYCLES SHARE A COMMON VERTEX AT n=7?")
print("=" * 70)

# Try to construct such a tournament
# All 3 cycles contain vertex 0:
# Cycle 1: 0→1→2→0 (need adj[0][1]=1, adj[1][2]=1, adj[2][0]=1)
# Cycle 2: 0→3→4→0 (need adj[0][3]=1, adj[3][4]=1, adj[4][0]=1)
# Cycle 3: 0→5→6→0 (need adj[0][5]=1, adj[5][6]=1, adj[6][0]=1)

# Required arcs:
# 0→1, 1→2, 2→0
# 0→3, 3→4, 4→0
# 0→5, 5→6, 6→0

# This means: 0 beats {1,3,5}, and {2,4,6} beat 0.
# Out-degree of 0 = 3 (beats 1,3,5).

# Remaining arcs to determine: all pairs not involving these fixed arcs.
# Fixed: 0→1, 0→3, 0→5, 2→0, 4→0, 6→0, 1→2, 3→4, 5→6
# Remaining pairs: {1,3}, {1,4}, {1,5}, {1,6}, {2,3}, {2,4}, {2,5}, {2,6},
#                  {3,5}, {3,6}, {4,5}, {4,6}
# That's 12 remaining pairs.

n = 7
print("\n  Constructing tournaments with 3 cycles all containing vertex 0...")
print("  Fixed arcs: 0→1, 0→3, 0→5, 2→0, 4→0, 6→0, 1→2, 3→4, 5→6")

remaining_pairs = [(1,3),(1,4),(1,5),(1,6),(2,3),(2,4),(2,5),(2,6),(3,5),(3,6),(4,5),(4,6)]

h7_found = False
h7_count = 0
h_dist = {}

for bits in range(2**12):
    adj = [[0]*n for _ in range(n)]

    # Fixed arcs
    adj[0][1] = 1  # 0→1
    adj[0][3] = 1  # 0→3
    adj[0][5] = 1  # 0→5
    adj[2][0] = 1  # 2→0
    adj[4][0] = 1  # 4→0
    adj[6][0] = 1  # 6→0
    adj[1][2] = 1  # 1→2
    adj[3][4] = 1  # 3→4
    adj[5][6] = 1  # 5→6

    # Set remaining arcs
    for idx, (i, j) in enumerate(remaining_pairs):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    # Count 3-cycles
    dc3 = count_dc3(adj, n)

    # We want EXACTLY dc3=3 (our 3 constructed cycles, no more)
    if dc3 != 3:
        continue

    # Count dc5, dc7
    dc5 = count_dc5(adj, n)
    dc7 = count_dc7(adj, n)

    alpha1 = dc3 + dc5 + dc7

    # Count HP
    H = count_hp(adj, n)

    h_dist[H] = h_dist.get(H, 0) + 1

    # Find the 3-cycles and check if they share vertex 0
    cycles = find_3cycles(adj, n)
    a2 = alpha2_from_cycles(cycles)

    # Also count disjoint pairs including dc5 and dc7
    # For simplicity, just use the OCF formula
    # H = 1 + 2*alpha1 + 4*alpha2 + 8*alpha3
    # alpha2 from OCF = (H - 1 - 2*alpha1) / 4 - 2*alpha3
    # But we need alpha3 too... at n=7, alpha3 max = floor(7/3) = 2
    # Actually let's just use the direct count for now

    if H == 7:
        h7_found = True
        h7_count += 1
        print(f"\n  *** H=7 FOUND! ***")
        print(f"  dc3={dc3}, dc5={dc5}, dc7={dc7}, α₁={alpha1}")
        print(f"  3-cycles: {cycles}")
        print(f"  α₂(from 3-cycles only) = {a2}")
        print(f"  Scores: {sorted([sum(adj[i]) for i in range(n)])}")

        # Print adjacency
        for i in range(n):
            print(f"    {i}: beats {[j for j in range(n) if adj[i][j]]}")

    elif dc5 == 0 and dc7 == 0 and a2 == 0:
        print(f"  dc3=3, dc5=0, dc7=0, a2=0: H={H}, α₁={alpha1}")
        # This would give H = 1 + 2*3 = 7 if OCF holds
        if H != 7:
            print(f"  *** OCF CHECK: H={H} but expected 1+2*3=7. α₂={a2} ***")

print(f"\n  Total tournaments with dc3=3 (from 3-through-vertex-0 construction): {sum(h_dist.values())}")
print(f"  H distribution: {sorted(h_dist.items())}")
print(f"  H=7 found: {h7_found} ({h7_count} tournaments)")

print("\n" + "=" * 70)
print("PART 2: EXHAUSTIVE CHECK FOR dc3=3, α₂=0 AT n=7")
print("=" * 70)

# Even broader: enumerate all tournaments with dc3=3 at n=7
# This is 2^21 = 2M tournaments, checking dc3 for each is feasible if fast

# Actually, the construction above constrains the 3-cycles to go through
# vertex 0. But there could be other 3-cycles among {1,...,6} too!
# We filtered for dc3==3, so there are no extra 3-cycles.

# The key insight: with dc3=3 and all 3 cycles sharing vertex 0,
# the 3 cycles use vertex sets {0,1,2}, {0,3,4}, {0,5,6}.
# These are vertex-disjoint except for 0.
# No pair is fully vertex-disjoint.
# So α₂(3-cycles) = 0.

# But there might be dc5 > 0 or dc7 > 0 cycles!
# If dc5 > 0, then α₁ > 3 and H > 7.

# So we need dc3=3, dc5=0, dc7=0, AND all 3 cycles sharing a vertex.
# Let's check how many such tournaments have H=7.

print("\n  Checking dc5 and dc7 for tournaments with dc3=3 and all cycles through v0:")

dc3_3_a2_0 = []
for bits in range(2**12):
    adj = [[0]*n for _ in range(n)]
    adj[0][1]=1; adj[0][3]=1; adj[0][5]=1
    adj[2][0]=1; adj[4][0]=1; adj[6][0]=1
    adj[1][2]=1; adj[3][4]=1; adj[5][6]=1
    for idx, (i,j) in enumerate(remaining_pairs):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    dc3 = count_dc3(adj, n)
    if dc3 != 3:
        continue

    dc5 = count_dc5(adj, n)
    if dc5 > 0:
        continue

    dc7 = count_dc7(adj, n)
    if dc7 > 0:
        continue

    H = count_hp(adj, n)
    cycles = find_3cycles(adj, n)
    a2 = alpha2_from_cycles(cycles)

    dc3_3_a2_0.append((H, a2, bits, cycles))

    if len(dc3_3_a2_0) <= 10:
        print(f"  bits={bits:4d}: dc3=3, dc5=0, dc7=0, α₁=3, α₂(3cy)={a2}, H={H}")
        if H == 7:
            print(f"  *** H=7 FOUND ***")

print(f"\n  Total with dc3=3, dc5=0, dc7=0: {len(dc3_3_a2_0)}")
if dc3_3_a2_0:
    h_vals = [x[0] for x in dc3_3_a2_0]
    a2_vals = [x[1] for x in dc3_3_a2_0]
    print(f"  H values: {sorted(set(h_vals))}")
    print(f"  α₂ values: {sorted(set(a2_vals))}")

    if 7 in h_vals:
        print(f"\n  *** H=7 IS ACHIEVABLE at n=7 ***")
        print(f"  This CONTRADICTS HYP-992!")
    else:
        print(f"\n  H=7 NOT found in this construction.")
        print(f"  This is consistent with HYP-992 (but not proof).")

print("\nDone.")
