#!/usr/bin/env python3
"""
Order-2 coefficients for DISJOINT tile pairs: finding the pattern.
opus-2026-04-03-S28

For disjoint tiles (no shared vertex), the coefficient c₂ is NOT just
a function of the two skip values. At n=6:
  AE disjoint skip=(4,2): c₂=4  but  AE has specific geometry
  BE disjoint skip=(4,4): c₂=16
  BG disjoint skip=(4,2): c₂=4

But some skip-(2,2) pairs have c₂=0 while others have c₂=4.
  DG (3,1)×(4,2) skip=(2,2): c₂=0
  DJ (3,1)×(6,4) skip=(2,2): c₂=4
  GI (4,2)×(5,3) skip=(2,2): c₂=0
  IJ (5,3)×(6,4) skip=(2,2): c₂=0

The difference: DJ = (3,1)×(6,4) — the two tiles are at OPPOSITE corners
of the staircase! They are the grid-reflection pair. Meanwhile DG and GI
are "near each other" in the staircase.

HYPOTHESIS: For disjoint tiles (x1,y1) and (x2,y2):
  c₂ = number of Hamiltonian paths in T(flip both) that use BOTH new arcs
  This counts the "cooperation" between the two flips.

Two tiles cooperate when they can both appear in the same Hamiltonian path.
This requires: the two flipped arcs y1→x1 and y2→x2 can be concatenated
in a path (either directly or with between-vertices).

When can y1→x1 and y2→x2 appear in the same HP?
  - If x1 < y2 or x2 < y1: one arc goes "below" the other, and they can
    be in the same path by visiting the lower arc first.
  - If the ranges [y1,x1] and [y2,x2] overlap: they compete for the same
    vertices and may not coexist.

Let me compute the c₂ exactly and look for the pattern in terms of
how the arcs' ranges relate.
"""

from itertools import permutations

def tournament_matrix(n, flips=None):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                T[i][j] = 1
    if flips:
        for (x, y) in flips:
            T[x][y] = 0
            T[y][x] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def list_hps(T, n):
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

for n in range(4, 8):
    print(f"\n{'='*70}")
    print(f"DISJOINT PAIR COEFFICIENTS: n={n}")
    print(f"{'='*70}")

    # Get all tile pairs
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    m = len(tiles)

    for i in range(m):
        for j in range(i+1, m):
            x1, y1 = tiles[i]
            x2, y2 = tiles[j]

            # Check disjoint
            shared = set(tiles[i]) & set(tiles[j])
            if shared:
                continue

            skip1, skip2 = x1-y1, x2-y2

            # Compute c₂
            H0 = 1  # transitive
            H_i = count_hp(tournament_matrix(n, [(x1-1, y1-1)]), n)
            H_j = count_hp(tournament_matrix(n, [(x2-1, y2-1)]), n)
            H_ij = count_hp(tournament_matrix(n, [(x1-1, y1-1), (x2-1, y2-1)]), n)
            c2 = H_ij - H_i - H_j + H0

            # Compute the RANGE of each tile
            # Tile (x,y) "covers" vertices y, y+1, ..., x
            range1 = set(range(y1, x1+1))
            range2 = set(range(y2, x2+1))
            overlap = range1 & range2
            range_relation = "NESTED" if range1 <= range2 or range2 <= range1 else \
                           "DISJOINT" if not overlap else "OVERLAP"

            # Is one range entirely above the other?
            above = y2 > x1 or y1 > x2  # no overlap at all in vertex range
            nested = range1 <= range2 or range2 <= range1

            # The BETWEEN count: vertices strictly between the two arcs
            # If ranges are disjoint: gap = min(y-values of higher arc) - max(x-values of lower arc) - 1
            if y2 > x1:
                gap = y2 - x1 - 1
                between = gap
            elif y1 > x2:
                gap = y1 - x2 - 1
                between = gap
            else:
                between = -1  # overlapping

            label_i = chr(65+i) if i < 26 else f'T{i}'
            label_j = chr(65+j) if j < 26 else f'T{j}'

            print(f"  {label_i}{label_j} ({x1},{y1})×({x2},{y2}) "
                  f"skip=({skip1},{skip2}) {range_relation} "
                  f"ranges=[{y1}..{x1}]×[{y2}..{x2}] "
                  f"between={between}: c₂={c2}")

            # Count paths using both new arcs
            if n <= 7:
                T2 = tournament_matrix(n, [(x1-1, y1-1), (x2-1, y2-1)])
                hps = list_hps(T2, n)
                both_count = 0
                for p in hps:
                    uses_1 = any(p[k] == y1-1 and p[k+1] == x1-1 for k in range(n-1))
                    uses_2 = any(p[k] == y2-1 and p[k+1] == x2-1 for k in range(n-1))
                    if uses_1 and uses_2:
                        both_count += 1
                if both_count != c2:
                    print(f"    *** MISMATCH: paths using both arcs = {both_count}, c₂ = {c2}")
                else:
                    # c₂ = number of paths using BOTH new arcs. Good!
                    pass

    # NESTED pairs (one range inside the other): what's the formula?
    print(f"\n  NESTED DISJOINT pairs:")
    for i in range(m):
        for j in range(i+1, m):
            x1, y1 = tiles[i]
            x2, y2 = tiles[j]
            if set(tiles[i]) & set(tiles[j]):
                continue

            range1 = set(range(y1, x1+1))
            range2 = set(range(y2, x2+1))
            if not (range1 < range2 or range2 < range1):
                continue

            # One is inside the other
            if range1 < range2:
                inner, outer = tiles[i], tiles[j]
                s_inner, s_outer = x1-y1, x2-y2
            else:
                inner, outer = tiles[j], tiles[i]
                s_inner, s_outer = x2-y2, x1-y1

            H_ij = count_hp(tournament_matrix(n, [(x1-1, y1-1), (x2-1, y2-1)]), n)
            H_i = count_hp(tournament_matrix(n, [(x1-1, y1-1)]), n)
            H_j = count_hp(tournament_matrix(n, [(x2-1, y2-1)]), n)
            c2 = H_ij - H_i - H_j + 1

            print(f"    inner=({inner[0]},{inner[1]}) s={s_inner}, "
                  f"outer=({outer[0]},{outer[1]}) s={s_outer}: c₂={c2}")
