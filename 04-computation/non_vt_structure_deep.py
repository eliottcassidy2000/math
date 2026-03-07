#!/usr/bin/env python3
"""
Deep structural analysis of the non-VT position-uniform tournament at n=5.

The tournament has:
  Vertices: 0, 1, 2, 3, 4
  Scores: (1, 2, 2, 2, 3)
  Aut: Z/3Z on {1,2,3}, fixing {0,4}
  Anti-Aut: 3 maps, each swapping 0<->4

Adjacency:
  0->4 (only out-edge of 0)
  4->1, 4->2, 4->3 (out-edges of 4)
  1->0, 1->3
  2->0, 2->1
  3->0, 3->2

So the structure is:
  - 0 is a "near-source" (only beats 4)
  - 4 is a "near-sink" (beats 1,2,3)
  - {1,2,3} form a directed 3-cycle: 1->3->2->1
  - All of {1,2,3} beat 0
  - 4 beats all of {1,2,3}

This looks like: take a 3-cycle C_3 on {1,2,3}, add a "source" 4 above
and a "sink" 0 below: 4 beats {1,2,3}, {1,2,3} beat 0, and 0 beats 4.
The result is a "3-cycle with poles" — a directed cone/cocone structure.

Key observation: This tournament is the COMPOSITION T_2[C_3, T_1]
where T_2 is the 2-vertex tournament, C_3 is the 3-cycle, and T_1 is
a single vertex. Actually, it's more like a "substitution" construction.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

# The non-VT scalar-M tournament
A = [
    [0, 0, 0, 0, 1],  # vertex 0: beats only 4
    [1, 0, 0, 1, 0],  # vertex 1: beats 0, 3
    [1, 1, 0, 0, 0],  # vertex 2: beats 0, 1
    [1, 0, 1, 0, 0],  # vertex 3: beats 0, 2
    [0, 1, 1, 1, 0],  # vertex 4: beats 1, 2, 3
]
n = 5

print("=" * 70)
print("STRUCTURAL ANALYSIS OF NON-VT POSITION-UNIFORM TOURNAMENT (n=5)")
print("=" * 70)

# Verify it's a tournament
for i in range(n):
    for j in range(i+1, n):
        assert A[i][j] + A[j][i] == 1, f"Not a tournament at ({i},{j})"
    assert A[i][i] == 0
print("\n  Verified: valid tournament")

# Edge structure
print("\n  Edge structure:")
print("    0 -> 4")
print("    4 -> 1, 2, 3")
print("    1 -> 0, 3  (3-cycle: 1->3)")
print("    2 -> 0, 1  (3-cycle: 2->1)")
print("    3 -> 0, 2  (3-cycle: 3->2)")
print()
print("    Pattern: 'Directed cone over C_3 with tip reversal'")
print("    Take 3-cycle {1,2,3}: 1->3->2->1")
print("    Add vertex 4 beating all of {1,2,3}")
print("    Add vertex 0 beaten by all of {1,2,3}")
print("    Edge between poles: 0->4 (surprising — the 'weak' beats the 'strong'!)")

# This is actually a well-known construction:
# "Tournament with a dominating vertex and a dominated vertex"
# 4 dominates {1,2,3}, 0 is dominated by {1,2,3}

# Key: the reversal 0<->4 is the anti-automorphism
# sigma = (0 4)(1 1)(2 3)(3 2) or similar

# Check all anti-automorphisms explicitly
print("\n  Anti-automorphisms (T -> T^op):")
antis = []
for perm in permutations(range(n)):
    if all(A[perm[i]][perm[j]] == A[j][i] for i in range(n) for j in range(n)):
        antis.append(perm)
        print(f"    {perm}")

# =====================================================================
# Is this construction generalizable?
# At n = 2k+1, take C_k on {1,...,k}, add poles 0 and n-1
# with n-1 dominating middle, 0 dominated by middle, 0->n-1
# =====================================================================
print("\n" + "=" * 70)
print("GENERALIZATION: CONE OVER C_k WITH POLE REVERSAL")
print("=" * 70)

def cone_tournament(k):
    """Build tournament on n = 2k+1 (wait, that's not right for k=1 giving n=3)
    Actually: n = k + 2 vertices (k middle + 2 poles)
    """
    n = k + 2
    A = [[0]*n for _ in range(n)]

    # Middle vertices: {1, ..., k} form cycle 1->k->...->2->1
    for i in range(k):
        v = 1 + i
        w = 1 + (i + 1) % k  # next in cycle... wait, need to match the n=5 pattern

    # Actually, for k=3: the 3-cycle is 1->3->2->1
    # More generally, a circulant C_k on {1,...,k}
    # For the construction to have position-uniformity, we probably need
    # the middle to be vertex-transitive (i.e., a circulant/regular tournament on k vertices)

    # For k=3: only tournament on 3 vertices (up to iso) is the 3-cycle
    # For k=5: Paley on 5 vertices, or C_5

    # Build: poles 0 and n-1
    # n-1 beats all middle, all middle beat 0, 0 beats n-1
    for v in range(1, k+1):
        A[n-1][v] = 1  # pole n-1 beats middle
        A[v][0] = 1     # middle beats pole 0
    A[0][n-1] = 1        # pole 0 beats pole n-1

    # Middle: need a tournament on {1,...,k}
    return A, n

# Test at k=1: n=3
A3, n3 = cone_tournament(1)
print(f"\n  k=1 (n=3): Cone over 1 vertex")
for i in range(n3):
    print(f"    {A3[i]}")

# With 1 middle vertex, the tournament is 0->2, 1->0, 2->1
# That's the 3-cycle! Which IS position-uniform.
ham3 = [p for p in permutations(range(n3)) if all(A3[p[i]][p[i+1]] == 1 for i in range(n3-1))]
print(f"  H = {len(ham3)}, position-uniform = {len(ham3) == 3}")

# Test at k=3: n=5
# Need middle tournament on {1,2,3}
# The 3-cycle: 1->3->2->1 (which is what our non-VT tournament has)

def cone_with_middle(middle_A, k):
    """Build cone tournament with given middle adjacency."""
    n = k + 2
    A = [[0]*n for _ in range(n)]

    # Copy middle tournament
    for i in range(k):
        for j in range(k):
            A[1+i][1+j] = middle_A[i][j]

    # Poles
    for v in range(1, k+1):
        A[n-1][v] = 1  # top pole beats middle
        A[v][0] = 1     # middle beats bottom pole
    A[0][n-1] = 1        # bottom beats top (reversal!)

    return A

# 3-cycle on {0,1,2}: 0->2->1->0
middle3 = [[0,0,1],[1,0,0],[0,1,0]]
A5_cone = cone_with_middle(middle3, 3)
print(f"\n  k=3 (n=5): Cone over C_3")
for i in range(5):
    print(f"    {A5_cone[i]}")

# Check if this matches our non-VT tournament
match = all(A5_cone[i][j] == A[i][j] for i in range(5) for j in range(5))
print(f"  Matches non-VT tournament: {match}")

# Check position-uniformity
ham5 = [p for p in permutations(range(5)) if all(A5_cone[p[i]][p[i+1]] == 1 for i in range(4))]
N5 = np.zeros((5,5), dtype=int)
for p in ham5:
    for j, v in enumerate(p):
        N5[v][j] += 1
pu5 = np.all(N5 == len(ham5)//5)
print(f"  H = {len(ham5)}, position-uniform = {pu5}")

# =====================================================================
# Generalize to k=5: cone over Paley_5
# =====================================================================
print(f"\n  k=5 (n=7): Cone over Paley_5")

# Paley on {0,...,4}: gen = {1,2}
middle5_paley = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i) % 5 in {1, 2}:
            middle5_paley[i][j] = 1

A7_cone = cone_with_middle(middle5_paley, 5)
print(f"  Adjacency:")
for i in range(7):
    print(f"    {A7_cone[i]}")

ham7 = [p for p in permutations(range(7)) if all(A7_cone[p[i]][p[i+1]] == 1 for i in range(6))]
H7 = len(ham7)
print(f"  H = {H7}, H/n = {H7/7:.4f}, H%7 = {H7%7}")

if H7 > 0:
    N7 = np.zeros((7,7), dtype=int)
    for p in ham7:
        for j, v in enumerate(p):
            N7[v][j] += 1
    pu7 = np.all(N7 == H7//7) if H7 % 7 == 0 else False
    print(f"  Position-uniform: {pu7}")
    print(f"  N[v,j] matrix:")
    for v in range(7):
        print(f"    vertex {v} (deg={sum(A7_cone[v])}): {list(N7[v])}")

# =====================================================================
# Also try cone over C_3^op (reverse the 3-cycle direction)
# =====================================================================
print(f"\n  k=3 (n=5): Cone over C_3^op")
middle3_op = [[0,1,0],[0,0,1],[1,0,0]]
A5_cone_op = cone_with_middle(middle3_op, 3)
ham5_op = [p for p in permutations(range(5)) if all(A5_cone_op[p[i]][p[i+1]] == 1 for i in range(4))]
H5_op = len(ham5_op)
if H5_op > 0 and H5_op % 5 == 0:
    N5_op = np.zeros((5,5), dtype=int)
    for p in ham5_op:
        for j, v in enumerate(p):
            N5_op[v][j] += 1
    pu5_op = np.all(N5_op == H5_op//5)
    print(f"  H = {H5_op}, position-uniform = {pu5_op}")
else:
    print(f"  H = {H5_op}, H%5 = {H5_op%5}")

# =====================================================================
# Does the ORIENTATION of 0->n-1 matter? Try n-1->0 instead.
# =====================================================================
print(f"\n  k=3 (n=5): Cone over C_3 with REVERSED pole edge (4->0)")
A5_rev = [row[:] for row in A5_cone]
A5_rev[0][4] = 0
A5_rev[4][0] = 1
ham5_rev = [p for p in permutations(range(5)) if all(A5_rev[p[i]][p[i+1]] == 1 for i in range(4))]
H5_rev = len(ham5_rev)
print(f"  H = {H5_rev}")
if H5_rev > 0 and H5_rev % 5 == 0:
    N5_rev = np.zeros((5,5), dtype=int)
    for p in ham5_rev:
        for j, v in enumerate(p):
            N5_rev[v][j] += 1
    pu5_rev = np.all(N5_rev == H5_rev//5)
    print(f"  Position-uniform = {pu5_rev}")
    print(f"  N matrix:")
    for v in range(5):
        print(f"    {list(N5_rev[v])}")
else:
    print(f"  H%5 = {H5_rev%5}, not position-uniform")

# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
The non-VT position-uniform tournament at n=5 is the "CONE CONSTRUCTION":
  Take a regular tournament T_mid on k vertices (here C_3)
  Add two "pole" vertices: top (dominator) and bottom (dominated)
  Top beats all of T_mid
  All of T_mid beats bottom
  Bottom beats top (the "reversal" edge)

This gives a tournament on n = k + 2 vertices with:
  - |Aut| = |Aut(T_mid)| (acts on the middle)
  - |Anti| = |Aut(T_mid)| (swap poles + middle anti-aut)
  - Score sequence: (1, (k-1)/2, ..., (k-1)/2, k) (poles have extreme scores)
  - NOT vertex-transitive (poles are in different orbits)
  - But POSITION-UNIFORM (each vertex at each position equally often)

QUESTION: Does this generalize to n=7 (cone over Paley_5)?
""")
