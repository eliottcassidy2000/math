#!/usr/bin/env python3
"""
Verify: for the TRANSITIVE tournament (all invariants = 0),
a_k(T) = A(n, k) (ordinary Eulerian numbers).

This would explain why the constant terms in the a_{n-1-m} OCF formulas
are Eulerian numbers A(n, m).

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb

def eulerian_number(n, k):
    """A(n,k) = permutations of [n] with k descents."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def transitive_tournament(n):
    """The transitive tournament: i -> j iff i < j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

for n in range(3, 8):
    A = transitive_tournament(n)
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1

    print(f"n={n}:")
    print(f"  Forward-edge dist: {dict(sorted(dist.items()))}")
    euler = {k: eulerian_number(n, k) for k in range(n)}
    print(f"  Eulerian A(n,k):   {euler}")

    # But wait: for the transitive tournament, A[i][j]=1 iff i<j.
    # So a forward edge at step i means perm[i] < perm[i+1], i.e. an ASCENT.
    # The number of perms with k ascents = A(n, k) (Eulerian number, descent version
    # depends on convention). Let's check:
    # Standard: A(n,k) = perms with k descents. Ascents = (n-1)-descents.
    # So perms with k forward edges (ascents) = A(n, n-1-k) = A(n, k) by palindromy!

    match = all(dist[k] == eulerian_number(n, k) for k in range(n))
    # Or maybe with shift:
    match_shifted = all(dist[k] == eulerian_number(n, n-1-k) for k in range(n))
    print(f"  Direct match: {match}, shifted match: {match_shifted}")

    # Actually, Eulerian number A(n,k) counts perms with k ASCENTS in some conventions,
    # or k DESCENTS in others. Let me check which convention we're using.
    # With A(n,k) = sum (-1)^j C(n+1,j)(k+1-j)^n, we have:
    # A(n,0) = 1 (identity perm: 0 descents), A(n,n-1) = 1 (reverse: n-1 descents)
    # So this is the DESCENT convention.
    # Forward edges in transitive = ascents = (n-1) - descents.
    # So a_k = A(n, n-1-k).
    print(f"  a_k = A(n, n-1-k): {all(dist[k] == eulerian_number(n, n-1-k) for k in range(n))}")
    print()
