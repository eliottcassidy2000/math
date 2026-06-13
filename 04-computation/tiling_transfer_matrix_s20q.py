#!/usr/bin/env python3
"""
tiling_transfer_matrix_s20q.py — kind-pasteur-2026-03-22-S20q

Build the TRANSFER MATRIX for the staircase tiling recursion.

The staircase delta_k is a binary tiling on C(k+1,2) positions.
The hypotenuse has k tiles. The n -> n+1 operation adds a new strip
of k+1 tiles along the hypotenuse.

The transfer matrix T(k) maps:
  (hypotenuse state of delta_k) -> (H contribution from adding strip k+1)

T(k) is a 2^k x 2^{k+1} matrix (old hypotenuse -> new hypotenuse).

Actually simpler: fix the base-path convention and think strip-by-strip.

For the staircase delta_k, strips are:
  Strip 2: 1 tile (the innermost)
  Strip 3: 2 tiles
  ...
  Strip k+1: k tiles (the hypotenuse)

A tournament is built by filling strips 2 through k+1.
H can be computed STRIP BY STRIP using a transfer matrix that
acts on the "state" visible at the strip boundary.

What is the "state"? For the transfer matrix method on lattice models,
the state is the configuration on the boundary between filled and unfilled.

For the staircase: if we fill strips from inside to outside,
after filling strip j, the boundary is strip j itself (j-1 tiles).
Adding strip j+1 (j tiles) creates interactions between the j-1 old
boundary tiles and the j new tiles.

Author: kind-pasteur-2026-03-22-S20q
"""

import sys
import numpy as np
from math import comb
from collections import defaultdict
from itertools import product as iterproduct

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

print("=" * 72)
print("  THE STRIP TRANSFER MATRIX FOR STAIRCASE TILINGS")
print("=" * 72)

# The idea: build the tournament strip by strip, tracking how H
# depends on the boundary configuration.
#
# But H is not a "local" quantity — it depends on the ENTIRE tournament.
# So the transfer matrix approach needs to track enough state
# to compute H incrementally.
#
# ALTERNATIVE: Instead of computing H directly, compute the
# H-GENERATING FUNCTION strip by strip.
#
# Define: for a partial tiling of strips 2 through j,
# let f_j(boundary) = sum over all completions of strips j+1 through k+1
#                      of H(complete tiling).
#
# Then f_{k+1}(boundary) = H(T) for the complete tiling.
# And f_j is computed from f_{j+1} by summing over all choices for strip j+1.
#
# But this is BACKWARD recursion, not forward.
#
# SIMPLER APPROACH: Just compute the H distribution by brute force
# for each strip configuration and see the pattern.

# For n=5 (delta_3): strips are 2, 3, 4.
# Strip 2: 1 tile at (1,1) = arc (3,1)
# Strip 3: 2 tiles at (1,2),(2,1) = arcs (4,2),(4,1)
# Strip 4: 3 tiles at (1,3),(2,2),(3,1) = arcs (5,3),(5,2),(5,1)

# The STATE after strip j = the configuration of strip j (j-1 bits).

# Let me compute: for each possible state of strips 2 and 3 (inner),
# what is the DISTRIBUTION of H over completions of strip 4 (outer)?

print("\n  n=5: STRIP-BY-STRIP ANALYSIS")
print()

n = 5
# Enumerate all tilings and group by inner configuration
inner_to_H = defaultdict(list)  # inner bits -> list of H values

for strip2_bit in range(2):  # 1 tile
    for strip3_bits in range(4):  # 2 tiles
        for strip4_bits in range(8):  # 3 tiles
            # Build the tournament
            # Strip 2, tile (1,1) = arc (3,1): 1-indexed a=3, b=1
            # Strip 3, tiles (1,2)=arc(4,2) and (2,1)=arc(4,1)
            # Strip 4, tiles (1,3)=arc(5,3), (2,2)=arc(5,2), (3,1)=arc(5,1)

            A = np.zeros((n, n), dtype=np.int8)
            # Base path: 5->4->3->2->1 (0-indexed: 4->3->2->1->0)
            for i in range(1, n):
                A[i][i-1] = 1  # i beats i-1

            # Strip 2: arc (3,1) -> 0-indexed (2, 0)
            if strip2_bit == 0:
                A[2][0] = 1  # forward: 3->1
            else:
                A[0][2] = 1  # backward: 1->3

            # Strip 3: arc (4,2) -> (3, 1) and arc (4,1) -> (3, 0)
            if (strip3_bits >> 0) & 1 == 0:
                A[3][1] = 1  # forward: 4->2
            else:
                A[1][3] = 1  # backward: 2->4

            if (strip3_bits >> 1) & 1 == 0:
                A[3][0] = 1  # forward: 4->1
            else:
                A[0][3] = 1  # backward: 1->4

            # Strip 4: arcs (5,3)->(4,2), (5,2)->(4,1), (5,1)->(4,0)
            if (strip4_bits >> 0) & 1 == 0:
                A[4][2] = 1  # forward: 5->3
            else:
                A[2][4] = 1  # backward: 3->5

            if (strip4_bits >> 1) & 1 == 0:
                A[4][1] = 1  # forward: 5->2
            else:
                A[1][4] = 1  # backward: 2->5

            if (strip4_bits >> 2) & 1 == 0:
                A[4][0] = 1  # forward: 5->1
            else:
                A[0][4] = 1  # backward: 1->5

            H = count_hp(A, n)

            inner_key = (strip2_bit, strip3_bits)
            inner_to_H[inner_key].append(H)

print(f"  Inner configurations (strip 2 + strip 3): {len(inner_to_H)}")
print()
print(f"  {'Inner (s2, s3)':>20s} {'H values over strip 4':>40s} {'avg H':>8s}")

for key in sorted(inner_to_H.keys()):
    Hs = inner_to_H[key]
    avg = sum(Hs) / len(Hs)
    H_set = sorted(set(Hs))
    # Count each
    H_counts = {h: Hs.count(h) for h in H_set}
    desc = " ".join(f"{h}x{c}" for h, c in sorted(H_counts.items()))
    print(f"  {str(key):>20s} {desc:>40s} {avg:>8.2f}")

print()

# The TRANSFER MATRIX: for each inner state, what is the weight vector
# over strip-4 states that gives each H value?
print("  TRANSFER MATRIX (inner state -> H distribution over outer strip):")
print()

# For each inner state, the strip-4 choices are 8 possibilities (3 bits).
# We can build a matrix: rows = inner states (8), cols = H values.
# But more useful: a matrix that maps inner to the H-SUM.

# Actually: the transfer matrix T should map:
# T[inner_state] = sum over outer_states of H(inner_state, outer_state)
# This gives E[H | inner].

# Or: T[inner_state, H_value] = count of outer states giving that H.

all_H_vals = sorted(set(h for hs in inner_to_H.values() for h in hs))
print(f"  H values: {all_H_vals}")
print()

# Build matrix: rows = inner states, cols = H values, entries = counts
inner_states = sorted(inner_to_H.keys())
matrix = np.zeros((len(inner_states), len(all_H_vals)), dtype=int)
for i, key in enumerate(inner_states):
    for h in inner_to_H[key]:
        j = all_H_vals.index(h)
        matrix[i][j] += 1

print(f"  {'Inner':>12s}", end="")
for h in all_H_vals:
    print(f" H={h:>2d}", end="")
print(f" {'E[H]':>6s}")

for i, key in enumerate(inner_states):
    print(f"  {str(key):>12s}", end="")
    for j in range(len(all_H_vals)):
        print(f" {matrix[i][j]:>4d}", end="")
    avg = sum(inner_to_H[key]) / len(inner_to_H[key])
    print(f" {avg:>6.1f}")

print()

# THE KEY: the y=x symmetry maps (s2, s3) -> (complement of s2, complement of s3 reversed?)
# Actually y=x on the staircase:
# (r,c) -> (c,r). For strip 2: (1,1) -> (1,1). Self-symmetric.
# For strip 3: (1,2) -> (2,1) and (2,1) -> (1,2). The two tiles swap.
# For strip 4: (1,3) -> (3,1) and (3,1) -> (1,3). The outer tiles swap.
#              (2,2) -> (2,2). Self-symmetric.

# y=x symmetry on bits: strip 2 bit is self-symmetric.
# Strip 3: bit 0 (tile (1,2)) swaps with bit 1 (tile (2,1)).
# But also: y=x maps the tiling bit b to 1-b (complement).
# So: inner state (s2, s3) maps to (1-s2, complement_and_swap(s3)).

# complement_and_swap of s3: swap bits 0 and 1, then complement each.
# For s3 = 0 (00): swap -> 00, complement -> 11 = 3.
# For s3 = 1 (01): swap -> 10 = 2, complement -> 01 = 1.
# For s3 = 2 (10): swap -> 01 = 1, complement -> 10 = 2.
# For s3 = 3 (11): swap -> 11 = 3, complement -> 00 = 0.

# So the y=x map is: (s2, s3) -> (1-s2, 3-s3).
# This should map H -> H (since H(T) = H(T^op) for self-complementary).

print("  y=x SYMMETRY CHECK:")
for key in sorted(inner_to_H.keys()):
    s2, s3 = key
    partner = (1 - s2, 3 - s3)
    Hs_orig = sorted(inner_to_H[key])
    Hs_partner = sorted(inner_to_H[partner])
    match = Hs_orig == Hs_partner
    print(f"    {key} <-> {partner}: H match = {match}")

print()
print("  The y=x symmetric inner states (self-complementary inner):")
for key in sorted(inner_to_H.keys()):
    s2, s3 = key
    if (1 - s2, 3 - s3) == (s2, s3):
        Hs = sorted(set(inner_to_H[key]))
        print(f"    {key}: H values = {Hs}")

# Fixed points of (s2, s3) -> (1-s2, 3-s3): need s2 = 1-s2 (impossible for binary)
# So NO inner state is y=x self-symmetric at this level.
# This is because the y=x symmetry acts as COMPLEMENT, which swaps 0 and 1.
# No binary string is its own complement (at length > 0).

print()
print("  No inner state is y=x self-symmetric (complement has no fixed points).")
print("  Self-complementary tournaments require PAIRS of symmetric inner states.")
print("  The SC condition involves the FULL tiling, not just the inner part.")
