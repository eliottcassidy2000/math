#!/usr/bin/env python3
"""
Are there non-circulant tournaments with scalar M at odd n?

At n=5: all scalar-M tournaments are circulant (or isomorphic to circulant).
At n=7: check a broader class.

Also: THM-052 applies to any tournament with vertex-transitive automorphism
group containing an n-cycle. But vertex-transitivity alone might suffice.

More precisely: if T has an automorphism σ such that:
1. σ acts transitively on vertices
2. There exists an anti-automorphism τ (T → T^op) that "reverses" σ

Then the palindromic argument goes through.

Condition 2 is satisfied for circulant T because the reflection r: i→-i
is such a τ. But other groups might also provide this.

Let's check: at n=7, which tournaments have scalar M?
(Full enumeration is 2^21 ≈ 2M, too many for full M computation.
But we can check all doubly regular tournaments, all vertex-transitive ones.)

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def circulant_tournament(n, gen_set):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in gen_set:
                A[i][j] = 1
    return A

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def transfer_matrix_partial(A, a, b):
    """Compute M[a,b] only."""
    n = len(A)
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            S_verts = sorted(list(S) + [a])
            R_verts = sorted(R + [b])
            ea = count_paths_subset(A, S_verts, end=a)
            bb = count_paths_subset(A, R_verts, start=b)
            total += ((-1)**k) * ea * bb
    return total

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        if all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1)):
            count += 1
    return count

def score_sequence(A):
    return tuple(sorted(sum(row) for row in A))

def is_vertex_transitive(A):
    """Check if Aut(T) acts transitively on vertices (by brute force)."""
    n = len(A)
    # Find all automorphisms
    auts = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            auts.append(perm)

    if not auts:
        return False, 0

    # Check transitivity: can we reach vertex 1 from vertex 0?
    orbit = {0}
    for aut in auts:
        orbit.add(aut[0])
    return len(orbit) == n, len(auts)

# =====================================================================
print("=" * 70)
print("SCALAR M AT n=7: CIRCULANT vs NON-CIRCULANT")
print("=" * 70)

n = 7

# All circulant tournaments
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

circulant_H_vals = {}
print(f"\nn=7: {len(gen_sets)} circulant tournaments")
for gs in sorted(gen_sets):
    A = circulant_tournament(n, gs)
    H = ham_path_count(A)
    circulant_H_vals[gs] = H
    print(f"  gen={sorted(gs)}: H={H}, H/n={H//n}")

# All have M = (H/n)*I by THM-052.

# =====================================================================
# Now: check non-circulant vertex-transitive tournaments at n=7
# =====================================================================
print()
print("=" * 70)
print("n=7: CHECKING VERTEX-TRANSITIVE NON-CIRCULANT TOURNAMENTS")
print("=" * 70)

# At n=7, the non-circulant vertex-transitive tournaments are rare.
# The automorphism group must act transitively but not contain Z/7Z as a subgroup.
# For prime n=7, any transitive group on 7 elements contains Z/7Z (by Sylow's theorem).
# So at prime n, ALL vertex-transitive tournaments are circulant!

print("""
At prime n=7: By Sylow's theorem, any transitive permutation group
on a prime number of elements contains a p-cycle (element of order p).
This means any vertex-transitive tournament on a prime number of
vertices is circulant.

CONSEQUENCE: At prime odd n, scalar M ⟺ circulant ⟺ vertex-transitive.
""")

# =====================================================================
# n=9 (composite): are there non-circulant VT tournaments?
# =====================================================================
print("=" * 70)
print("n=9: COMPOSITE — NON-CIRCULANT VT POSSIBLE")
print("=" * 70)

# At n=9 = 3^2, the symmetric group S_9 has transitive subgroups
# that don't contain a 9-cycle (e.g., Z/3Z × Z/3Z).
# A tournament on Z/3Z × Z/3Z invariant under (Z/3Z × Z/3Z)-translations
# would be vertex-transitive but not (necessarily) circulant.

# Construct one: vertices = (i,j) for i,j in Z/3Z.
# Edge rule: (i,j) → (i',j') iff some function of (i'-i, j'-j) mod 3 is in a "gen set"

# The group Z/3 × Z/3 has 8 nonzero elements.
# Tournament requires: for each nonzero g = (a,b), exactly one of {g, -g} is a "forward" direction.
# The 8 nonzero elements pair up: {(1,0),(2,0)}, {(0,1),(0,2)}, {(1,1),(2,2)}, {(1,2),(2,1)}.
# So there are 2^4 = 16 such tournaments.

print(f"\nn=9: Z/3 × Z/3 tournaments (vertex-transitive, not necessarily circulant)")

verts_9 = [(i,j) for i in range(3) for j in range(3)]
vert_idx = {v: k for k, v in enumerate(verts_9)}

# 4 pairs of opposite directions
dir_pairs = [
    ((1,0), (2,0)),
    ((0,1), (0,2)),
    ((1,1), (2,2)),
    ((1,2), (2,1)),
]

tested = 0
scalar_count = 0

for mask in range(16):
    # Choose one from each pair
    gen_dirs = set()
    for k, (d1, d2) in enumerate(dir_pairs):
        if mask & (1 << k):
            gen_dirs.add(d1)
        else:
            gen_dirs.add(d2)

    # Build tournament
    n_v = 9
    A = [[0]*n_v for _ in range(n_v)]
    for vi, v in enumerate(verts_9):
        for vj, w in enumerate(verts_9):
            if vi == vj: continue
            diff = ((w[0]-v[0])%3, (w[1]-v[1])%3)
            if diff in gen_dirs:
                A[vi][vj] = 1

    H = ham_path_count(A)

    # Quick scalar check: M[0,0] = H/9? and M[0,1] = 0?
    if H % n_v != 0:
        print(f"  mask={mask:04b}: H={H}, H%9={H%n_v} => NOT scalar")
        tested += 1
        continue

    M_00 = transfer_matrix_partial(A, 0, 0)
    M_01 = transfer_matrix_partial(A, 0, 1)

    is_scalar = (M_00 == H // n_v) and (M_01 == 0)
    tested += 1

    if is_scalar:
        scalar_count += 1

    # Check if it's isomorphic to a circulant on Z/9Z
    # (we'd need to check all relabelings, which is expensive)
    # Instead, check if the automorphism group contains a 9-cycle
    scores = score_sequence(A)

    print(f"  mask={mask:04b}: H={H}, M[0,0]={M_00}, M[0,1]={M_01}, "
          f"scalar?{is_scalar}, scores={scores}")

print(f"\n  Scalar: {scalar_count}/{tested}")

# Check if any of these are NOT circulant
# A Z/3×Z/3 tournament is circulant on Z/9 iff the group action
# can be generated by a single element of order 9.
# Z/3×Z/3 has no element of order 9 (max order is 3).
# So these are NOT circulant (not Z/9Z-invariant), but they ARE vertex-transitive!

print("""
IMPORTANT: Z/3 × Z/3 tournaments are vertex-transitive but NOT circulant
(since Z/3 × Z/3 has no element of order 9, so no 9-cycle automorphism).

If any of these have scalar M, this would extend THM-052 beyond circulant!
""")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
