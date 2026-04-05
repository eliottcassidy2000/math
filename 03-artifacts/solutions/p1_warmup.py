#!/usr/bin/env python3
"""
Problem Set 7 — Problem 1: Warm-up.

(a) 2^C(n,2) labeled tournaments.
(b) The four tournaments on 4 vertices.
(c) Self-complementary identification.
(d) E[H(T)] = n!/2^{n-1} by linearity of expectation.
"""

from math import factorial, comb
from itertools import permutations

def all_tournaments(n):
    """Generate all 2^C(n,2) labeled tournaments."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(pairs)):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def canonical_form(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def score_seq(T):
    return tuple(sorted(sum(row) for row in T))

def complement(T):
    n = len(T)
    return [[T[j][i] for j in range(n)] for i in range(n)]

# (a)
for n in range(1, 7):
    print(f"n={n}: 2^C({n},2) = {2**comb(n,2)} labeled tournaments")

# (b)
print("\nNon-isomorphic tournaments on n=4:")
seen = {}
for T in all_tournaments(4):
    c = canonical_form(T)
    if c not in seen:
        seen[c] = T
        print(f"  scores={score_seq(T)}")
print(f"T(4) = {len(seen)}")

# (c)
print("\nSelf-complementary check:")
for c, T in seen.items():
    Tc = complement(T)
    if canonical_form(Tc) == c:
        print(f"  Self-complementary: scores={score_seq(T)}")

# (d)
print("\nE[H(T)] = n!/2^{n-1}:")
for n in range(1, 8):
    print(f"  n={n}: E[H] = {factorial(n)}/{2**(n-1)} = {factorial(n)/2**(n-1):.4f}")
