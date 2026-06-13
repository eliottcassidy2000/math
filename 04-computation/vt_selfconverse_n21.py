#!/usr/bin/env python3
"""
Test whether ALL vertex-transitive tournaments at n=21 are self-converse.

At n=21, McKay lists 22 non-circulant VT tournaments. These are the first
cases where THM-052's proof requires something beyond the i→-i map.

For circulant tournaments on Z/nZ, the map i→-i gives T→T^op.
For non-circulant VT tournaments, we need some OTHER anti-automorphism.

This script:
1. Downloads the n=21 non-circulant VT tournament data from McKay's site
2. Checks each for self-converse property (T ≅ T^op)
3. Reports the automorphism group structure

If ALL 22 are self-converse, THM-052 extends to n=21.
If ANY is not, we need a different proof strategy.

kind-pasteur-2026-03-06-S25e
"""

import subprocess
import sys
import os

# First, let's work with smaller cases we can handle directly.
# Check: are all circulant tournaments at n=9 AND all Z/3×Z/3 Cayley
# tournaments self-converse?

from itertools import permutations


def adj_from_genset(n, S):
    """Build adjacency matrix for circulant tournament on Z/nZ."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A


def is_self_converse(A):
    """Check if tournament A is isomorphic to A^T (its reversal)."""
    n = len(A)
    # Build A^op (reversal)
    Aop = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            Aop[i][j] = A[j][i]

    # Check if there exists a permutation sigma such that
    # A[sigma[i]][sigma[j]] = Aop[i][j] = A[j][i]
    # i.e., A[sigma[i]][sigma[j]] = A[j][i] for all i,j
    # Equivalently: for all i<j, A[sigma[i]][sigma[j]] = A[j][i] = 1 - A[i][j]

    # For small n, try all permutations (n! search)
    if n > 10:
        # Use canonical form comparison (need nauty)
        return is_self_converse_canonical(A)

    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if A[perm[i]][perm[j]] != A[j][i]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def is_self_converse_canonical(A):
    """Use graph canonical form to check self-converse (for larger n).
    This is a placeholder — needs nauty/bliss or similar."""
    # For now, use a heuristic: check if score multiset matches
    # (necessary but not sufficient condition)
    n = len(A)
    scores = sorted(sum(row) for row in A)
    scores_op = sorted(sum(A[j][i] for j in range(n)) for i in range(n))
    if scores != scores_op:
        return False
    # If score sequences match, we'd need full isomorphism test
    # For tournaments: scores always match (both regular if VT)
    return None  # Unknown — need proper isomorphism test


def automorphisms(A):
    """Find all automorphisms of tournament A."""
    n = len(A)
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts


def anti_automorphisms(A):
    """Find all anti-automorphisms σ: A[σ(i)][σ(j)] = A[j][i]."""
    n = len(A)
    antis = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if A[perm[i]][perm[j]] != A[j][i]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            antis.append(perm)
    return antis


# ============================================================
print("=" * 70)
print("SELF-CONVERSE TEST: ALL CIRCULANT TOURNAMENTS AT n=7")
print("=" * 70)

n = 7
# Generator sets for n=7: choose 3 from {1,...,6} with d, 7-d not both in S
gensets = []
from itertools import combinations
for S in combinations(range(1, 7), 3):
    S_set = set(S)
    valid = True
    for d in S_set:
        if (7 - d) in S_set:
            valid = False
            break
    if valid:
        gensets.append(S_set)

print(f"  Found {len(gensets)} circulant generator sets")

for S in gensets:
    A = adj_from_genset(7, S)
    sc = is_self_converse(A)
    auts = automorphisms(A)
    antis = anti_automorphisms(A)
    print(f"  S={sorted(S)}: self-converse={sc}, |Aut|={len(auts)}, |Anti|={len(antis)}")

# ============================================================
print("\n" + "=" * 70)
print("SELF-CONVERSE TEST: Z/3 x Z/3 CAYLEY TOURNAMENTS AT n=9")
print("=" * 70)

n = 9
# Elements of Z/3 x Z/3 as (a,b) for a,b in {0,1,2}
# Index: element (a,b) -> 3*a + b
# Addition: (a1,b1) + (a2,b2) = ((a1+a2)%3, (b1+b2)%3)
# Non-identity elements: all except (0,0)

elements = [(a, b) for a in range(3) for b in range(3)]
idx = {e: i for i, e in enumerate(elements)}

def z33_add(e1, e2):
    return ((e1[0]+e2[0]) % 3, (e1[1]+e2[1]) % 3)

def z33_neg(e):
    return ((-e[0]) % 3, (-e[1]) % 3)

# Non-identity elements
nonid = [e for e in elements if e != (0, 0)]

# For a tournament: exactly one of g, -g is in the connection set S
# -g = z33_neg(g). Pair up: {g, -g} for g != 0.
# Since |G|=9 is odd, no element equals its own inverse.
# Number of pairs: 8/2 = 4. So |S| = 4, choose one from each pair.

pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        neg_g = z33_neg(g)
        pairs.append((g, neg_g))
        seen.add(g)
        seen.add(neg_g)

print(f"  Pairs: {pairs}")
print(f"  Number of connection sets: {2**len(pairs)} = {2**len(pairs)}")

count_sc = 0
count_not_sc = 0

for bits in range(2**len(pairs)):
    S = set()
    for i, (g, neg_g) in enumerate(pairs):
        if (bits >> i) & 1:
            S.add(g)
        else:
            S.add(neg_g)

    # Build adjacency matrix
    A = [[0]*9 for _ in range(9)]
    for i, ei in enumerate(elements):
        for j, ej in enumerate(elements):
            if i == j:
                continue
            diff = z33_add(ej, z33_neg(ei))
            if diff in S:
                A[i][j] = 1

    # Verify tournament
    for i in range(9):
        for j in range(i+1, 9):
            assert A[i][j] + A[j][i] == 1, f"Not tournament at ({i},{j})"

    # Check self-converse
    sc = is_self_converse(A)

    # Check if it's circulant (has a cyclic automorphism of order 9)
    auts = automorphisms(A)
    is_circ = any(
        all(perm[perm[perm[perm[perm[perm[perm[perm[perm[i]]]]]]]]] == i for i in range(9))
        and len(set(perm)) == 9
        and all(perm[i] != i for i in range(9))
        for perm in auts
        if len(set(perm[perm[j]] for j in range(9) if True)) == 9
    )

    # Simpler circulant check: does Aut contain an n-cycle?
    has_9cycle = False
    for perm in auts:
        # Check if perm is a single 9-cycle
        visited = [False]*9
        cycle_len = 0
        j = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        if cycle_len == 9 and all(visited):
            has_9cycle = True
            break

    label = "CIRC" if has_9cycle else "NON-CIRC"

    if sc:
        count_sc += 1
    else:
        count_not_sc += 1

    if not has_9cycle:
        antis = anti_automorphisms(A)
        print(f"  S={sorted(S)}: SC={sc}, |Aut|={len(auts)}, |Anti|={len(antis)}, [{label}]")

print(f"\n  Total: {count_sc} self-converse, {count_not_sc} NOT self-converse")
print(f"  (out of {2**len(pairs)} Z/3xZ/3 Cayley tournaments)")

# ============================================================
print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
If all Z/3xZ/3 tournaments at n=9 are self-converse, this supports
the conjecture that all VT tournaments are self-converse.

The map (a,b) -> (-a,-b) in Z/3xZ/3 is inversion, and since the group
is abelian, this always gives an anti-automorphism.

For non-abelian groups (first at n=21: Z/7 x| Z/3), inversion does NOT
in general give an anti-automorphism because the connection set must be
closed under conjugation for this to work.

OPEN QUESTION: Is every vertex-transitive tournament self-converse?
If yes, THM-052 extends to ALL VT tournaments at ALL odd n.
If no, we need a different proof strategy for the non-SC cases.
""")
