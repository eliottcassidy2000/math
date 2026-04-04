#!/usr/bin/env python3
"""
PROVING: Linear coefficient of tile (x,y) in H = 2^(x-y-1).
opus-2026-04-03-S28

Observation: In the multilinear polynomial H(x_1,...,x_m), the coefficient
of x_j (where tile j = (x,y)) is 2^(skip-1) where skip = x-y.

This has a beautiful combinatorial interpretation:
  - The coefficient of x_j equals H(x_j=1, all others=0) - H(x_j=0, all others=0)
  - But that's just two specific tournaments differing in one arc.
  - H(all zeros) = 1 (transitive tournament, unique HP).
  - H(flip tile j only) = 1 + (number of NEW Hamiltonian paths created by
    flipping one arc from the transitive tournament).

Actually, the multilinear coefficient c_j is NOT simply H(e_j) - H(0).
The multilinear polynomial is:
  H(x) = sum_S c_S * prod_{i in S} x_i

so c_j = H(e_j) - H(0) only if there are no higher-order terms involving j.
In general, c_j = sum_{S containing j} c_S * (evaluated at x=0 for non-j bits).

Wait, the coefficient IS H(e_j) - H(0) by Möbius inversion for order-1 terms:
  c_{j} = sum_{T subset {j}} (-1)^{1-|T|} H(e_T) = H(e_j) - H(0)

So: linear coefficient of tile j = H(transitive + flip j) - H(transitive)

The transitive tournament has score sequence (0, 1, 2, ..., n-1) and H=1.
Flipping the arc between vertices x and y (where x > y+1) creates a
non-transitive tournament. How many Hamiltonian paths does it have?

CLAIM: H(transitive + flip (x,y)) = 1 + 2^(x-y-1)
=> linear coefficient = 2^(x-y-1)

PROOF SKETCH:
The transitive tournament has the unique HP: n -> n-1 -> ... -> 1.
After flipping arc (x,y) to go y -> x instead of x -> y:
- The original HP n -> ... -> x -> ... -> y -> ... -> 1 is still valid
  (it uses the arc x -> x-1, ..., y+1 -> y, not x -> y directly)
  WAIT: The transitive tournament has arcs i -> j for ALL i > j.
  So the HP uses arcs n->n-1, n-1->n-2, ..., 2->1 (the base path).
  The arc x->y (with x-y >= 2) is NOT used by the unique HP.

  After flipping x->y to y->x, the base path n->n-1->...->1 still exists!
  So the original HP is preserved. But now there might be NEW HPs.

  A new HP must use the arc y->x at some point. Since all other arcs
  still go from larger to smaller vertex, the new HP must:
  - Start from some vertex, descend to y, then jump to x, then descend again.
  - The vertices between y+1 and x-1 must be visited either before y or after x.

  The number of ways to interleave the "below y" vertices, the y->x jump,
  and the "above x" vertices gives exactly 2^(x-y-1) new paths.

Let me verify this.
"""

from itertools import permutations
from collections import defaultdict

def tournament_matrix(n, flip_arc=None):
    """Create transitive tournament, optionally flip one arc."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                T[i][j] = 1  # i -> j for i > j (transitive)
    if flip_arc:
        x, y = flip_arc  # x > y, flip from x->y to y->x
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
    """Enumerate all Hamiltonian paths."""
    paths = []
    for perm in permutations(range(n)):
        valid = all(T[perm[i]][perm[i+1]] for i in range(n-1))
        if valid:
            paths.append(perm)
    return paths

print("="*70)
print("LINEAR COEFFICIENT = H(trans + flip) - H(trans) = 2^(skip-1)")
print("="*70)

for n in range(3, 9):
    print(f"\nn={n}:")
    T_trans = tournament_matrix(n)
    H_trans = count_hp(T_trans, n)
    assert H_trans == 1, f"Transitive should have H=1, got {H_trans}"

    # For each possible arc flip (x,y) with x-y >= 2 (0-indexed: x > y+1)
    for x in range(n):
        for y in range(x):
            skip = x - y
            if skip < 2:
                continue
            T_flip = tournament_matrix(n, flip_arc=(x, y))
            H_flip = count_hp(T_flip, n)
            linear_coeff = H_flip - H_trans
            predicted = 2**(skip - 1)
            match = "✓" if linear_coeff == predicted else "✗"
            if linear_coeff != predicted or n <= 6:
                print(f"  arc ({x+1},{y+1}) skip={skip}: H={H_flip}, "
                      f"coeff={linear_coeff}, 2^(skip-1)={predicted} {match}")

print(f"\n{'='*70}")
print("PROOF: Why 2^(skip-1)?")
print(f"{'='*70}")

# Let me enumerate the actual new paths at n=5, flip (5,1) -> (1,5)
# In 0-indexed: flip (4,0) arc from 4->0 to 0->4
n = 5
print(f"\nn={n}, flipping arc (5,1) [0-indexed: (4,0)]:")
T_trans = tournament_matrix(n)
T_flip = tournament_matrix(n, flip_arc=(4, 0))

paths_before = list_hps(T_trans, n)
paths_after = list_hps(T_flip, n)

print(f"Paths before (transitive): {len(paths_before)}")
for p in paths_before:
    print(f"  {' -> '.join(str(v+1) for v in p)}")

print(f"\nPaths after (flip arc 5->1 to 1->5): {len(paths_after)}")
for p in paths_after:
    uses_flip = any(p[i] == 0 and p[i+1] == 4 for i in range(n-1))
    print(f"  {' -> '.join(str(v+1) for v in p)} {'[uses 1->5]' if uses_flip else '[old path]'}")

new_paths = [p for p in paths_after if p not in paths_before]
print(f"\nNew paths: {len(new_paths)} = 2^(5-1-1) = {2**(5-1-1)}?")
print(f"Actually: 2^(skip-1) = 2^(4-1) = {2**3}")

# For n=6, flip (6,1) -> (1,6)
n = 6
print(f"\nn={n}, flipping arc (6,1) [0-indexed: (5,0)]:")
T_flip = tournament_matrix(n, flip_arc=(5, 0))
paths_after = list_hps(T_flip, n)
new_paths = [p for p in paths_after if not all(T_trans[p[i]][p[i+1]] for i in range(n-1))]
print(f"New paths: {len(new_paths)}")
print(f"2^(skip-1) = 2^(5-1) = {2**4}")
for p in new_paths:
    print(f"  {' -> '.join(str(v+1) for v in p)}")

# The STRUCTURE of new paths:
# After flipping x->y to y->x (x > y by skip >= 2):
# The new HP must use the arc y->x.
# Path: ... -> y -> x -> ...
# Vertices below y: {0, ..., y-1} — must be visited
# Vertices between y and x: {y+1, ..., x-1} — must be visited
# Vertices above x: {x+1, ..., n-1} — must be visited
#
# Since all arcs (other than y->x) go from larger to smaller vertex:
# After reaching x, we must descend. Before reaching y, we must descend.
# The "between" vertices {y+1, ..., x-1} can be visited either:
#   - Before y (while descending to y from above), OR
#   - After x (but wait, we need to go UP from y to x, and we can't go back up after x)
# Actually, we CAN visit them after x! From x, we descend. The between-vertices
# are between y and x, so they're BELOW x. We can visit them on the way down from x.
#
# The key insight: after flipping y->x, a new HP has the form:
#   [descend from top to y through some subset S of between-vertices]
#   -> y -> x ->
#   [descend from x through the remaining between-vertices and down to bottom]
#
# Each between-vertex can go in either part (before y or after x).
# There are (x-y-1) between-vertices, giving 2^(x-y-1) choices.

print(f"\n{'='*70}")
print("COMBINATORIAL PROOF:")
print("After flipping arc x->y to y->x:")
print("  Between-vertices: {y+1, ..., x-1} have (x-y-1) = (skip-1) elements")
print("  Each can be visited before y or after x: 2^(skip-1) new paths")
print(f"{'='*70}")

# WAIT: Let me also verify that the QUADRATIC coefficients have a pattern.
# At n=5, the order-2 coefficient of tiles (i,j) is:
# For shared-vertex pairs: always -2
# For disjoint pairs: can be +2, +4, or 0

# Let me check: for shared-vertex pairs, is the coefficient always -2?
print(f"\n{'='*70}")
print("ORDER-2 COEFFICIENTS: Testing c_{{ij}} pattern")
print(f"{'='*70}")

def get_tiles_1indexed(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

for n in range(4, 7):
    m = (n-1)*(n-2)//2
    tiles = get_tiles_1indexed(n)
    print(f"\nn={n}, m={m}:")

    # H values for computing Möbius
    T_trans = tournament_matrix(n)
    H_trans = count_hp(T_trans, n)  # = 1

    # Compute H for single and double flips
    H_single = {}
    for i in range(m):
        x, y = tiles[i]
        T = tournament_matrix(n, flip_arc=(x-1, y-1))
        H_single[i] = count_hp(T, n)

    H_double = {}
    for i in range(m):
        for j in range(i+1, m):
            T = [[0]*n for _ in range(n)]
            for a in range(n):
                for b in range(n):
                    if a > b:
                        T[a][b] = 1
            # Flip both arcs
            x1, y1 = tiles[i]
            T[x1-1][y1-1] = 0
            T[y1-1][x1-1] = 1
            x2, y2 = tiles[j]
            T[x2-1][y2-1] = 0
            T[y2-1][x2-1] = 1

            H_double[(i,j)] = count_hp(T, n)

    # Order-2 Möbius coefficient:
    # c_{ij} = H(e_i + e_j) - H(e_i) - H(e_j) + H(0)
    for i in range(m):
        for j in range(i+1, m):
            c2 = H_double[(i,j)] - H_single[i] - H_single[j] + H_trans
            t_i, t_j = tiles[i], tiles[j]
            shared = set(t_i) & set(t_j)
            skip_i = t_i[0] - t_i[1]
            skip_j = t_j[0] - t_j[1]

            # Is there a pattern?
            if abs(c2) > 0:
                print(f"  ({t_i[0]},{t_i[1]})×({t_j[0]},{t_j[1]}) "
                      f"skip=({skip_i},{skip_j}) "
                      f"shared={shared if shared else 'none'}: "
                      f"c2={c2}")
