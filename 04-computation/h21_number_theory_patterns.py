"""
h21_number_theory_patterns.py — Look for number-theoretic patterns in H-spectrum.

Known facts:
- H is always odd (Rédei's theorem)
- H=1 is achievable (transitive tournament)
- H=3 is achievable at n>=3
- H=5 is achievable at n>=4
- H=7 is NEVER achievable (permanent gap)
- H=9, 11, ..., 19 are achievable at various n
- H=21 is NEVER achievable (permanent gap)
- H=23, 25, ... mostly achievable

Questions:
1. Is 21 = 3 * 7 significant? Both 3 and 7 are H-values with special roles.
2. Are there other permanent gaps beyond H=7 and H=21?
3. Is there a number-theoretic characterization of achievable H?

Author: opus-2026-03-07-S43
"""

# Load known H-spectrum data
# From h21_achievable_w_n7.c: all achievable w values at n=7
# From h21_w_n8_sample.c: sampled w values at n=8

# Known achievable H values at each n (from exhaustive + sampling):
# n=1: {1}
# n=2: {1}
# n=3: {1, 3}
# n=4: {1, 3, 5}
# n=5: {1, 3, 5, 9, 11, 15}
# n=6: {1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, ...}
# n=7: full set known from exhaustive computation

# Let's analyze the binary structure of w values
# w = (H-1)/2, H odd
# H = 2w + 1

# Known permanent gaps: w=3 (H=7), w=10 (H=21)
# In binary: w=3 = 011, w=10 = 1010

# Check if there's a pattern in the "gap structure"
print("=== PERMANENT GAP ANALYSIS ===\n")

# w values that are permanent gaps (never achievable at any n)
permanent_gaps_w = [3, 10]  # known

# Analyze these in various number systems
for w in permanent_gaps_w:
    h = 2*w + 1
    print(f"w={w} (H={h}):")
    print(f"  Binary: {bin(w)}")
    print(f"  w mod 3 = {w % 3}")
    print(f"  w mod 4 = {w % 4}")
    print(f"  w mod 7 = {w % 7}")
    print(f"  Factorization of H={h}: ", end="")
    n = h
    factors = []
    for p in range(2, n+1):
        while n % p == 0:
            factors.append(p)
            n //= p
    print(" × ".join(map(str, factors)))
    print(f"  H = {h} = sum of two achievable H values? ", end="")
    # H=7: 3+? nope (not even)
    # H=21: 3+18? 5+16? Not meaningful for odd values.
    print()

# Look at the decomposition w = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...
# This is a "base-2 signed digit" representation
# Each alpha_k counts k-element independent sets in Omega(T)

# Observation: w=3 has the property that it's the SMALLEST value where
# EVERY decomposition fails. w=10 is the NEXT such value.
# Is this a recognizable sequence?

print("\n=== CHECKING: 3, 10, ... in OEIS-style ===")
print("3, 10 — differences: 7")
print("If pattern is a_n = a_{n-1} + (next prime)? 3+7=10, 10+?")
print("If pattern involves triangular numbers: T(2)=3, T(4)=10")
print("  T(n) = n(n+1)/2. T(2)=3, T(4)=10. Skip T(3)=6.")
print("  So gaps at T(2) and T(4)? Next would be T(6)=21 or T(8)=36?")
print()

# More careful: maybe the gap is at specific values related to
# the obstruction structure
print("=== ANALYSIS: Why these specific values? ===\n")

# For w=3: decompositions (alpha_1, alpha_2) with alpha_1 + 2*alpha_2 = 3, alpha_3=0:
# (1,1): alpha_2=1 but only 1 cycle => C(1,2)=0 pairs. Infeasible.
# (3,0): 3 mutual conflicts. THM-029 blocks.
# Key constraint: alpha_2 <= C(alpha_1, 2) AND tournament-specific.

# For w=10: (4,3), (6,2), (8,1), (10,0) all blocked.
# Key: at each decomposition, either:
#   - The Omega graph is too dense (few independent pairs) for the alpha_2 needed
#   - OR the Omega graph is too constrained by tournament structure

# Let's compute: for each w, how many decompositions survive basic feasibility?
from math import comb

print("w | H  | #decomps | #graph-feasible | #need_tourn_obstruction")
print("-" * 65)

for w in range(1, 35):
    h = 2*w + 1
    decomps = 0
    feasible = 0
    for a1 in range(w + 1):
        rem = w - a1
        if rem < 0 or rem % 2 != 0:
            continue
        a2 = rem // 2
        decomps += 1
        if a2 <= comb(a1, 2):
            # Check Turán bound
            turan = a1 * a1 // 4
            if a2 <= turan:
                feasible += 1
    gap = "***GAP***" if w in permanent_gaps_w else ""
    print(f" {w:2d} | {h:3d} | {decomps:8d} | {feasible:15d} | {gap}")

# KEY INSIGHT: Look at w values with VERY FEW graph-feasible decompositions
# These are the candidates for permanent gaps (all decompositions can be blocked)

print("\n=== w VALUES WITH FEW FEASIBLE DECOMPOSITIONS ===\n")
print("w | H  | feasible | decompositions")
print("-" * 60)
for w in range(1, 50):
    h = 2*w + 1
    feasible_decomps = []
    for a1 in range(w + 1):
        rem = w - a1
        if rem < 0 or rem % 2 != 0:
            continue
        a2 = rem // 2
        if a2 <= comb(a1, 2) and a2 <= a1 * a1 // 4:
            feasible_decomps.append((a1, a2))
    if len(feasible_decomps) <= 5:
        is_gap = " ***GAP***" if w in permanent_gaps_w else ""
        print(f" {w:2d} | {h:3d} | {len(feasible_decomps):8d} | {feasible_decomps}{is_gap}")
