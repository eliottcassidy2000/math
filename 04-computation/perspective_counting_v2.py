#!/usr/bin/env python3
"""
Fast perspective counting using score-sequence bucketing.
Conjecture: P(n) = sum over iso classes of (# vertex orbits) = 2*(n-1)!

Strategy: bucket by score sequence, only compare within buckets.
For n<=5: exhaustive. For n=6: exhaustive with optimized canonical.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations
from math import factorial

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def are_isomorphic(A, B, n):
    """Check if A and B are isomorphic tournaments."""
    # Quick reject by score sequence
    if score_seq(A, n) != score_seq(B, n):
        return False
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            if not match:
                break
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != B[i][j]:
                    match = False
                    break
        if match:
            return True
    return False

def vertex_orbits(A, n):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            if not is_aut:
                break
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    is_aut = False
                    break
        if is_aut:
            for v in range(n):
                union(v, perm[v])

    return len(set(find(v) for v in range(n)))

print("PERSPECTIVE COUNTING: P(n) = sum of vertex orbits across iso classes")
print("Conjecture: P(n) = 2*(n-1)!")
print("=" * 60)

for n in range(2, 7):
    m = n*(n-1)//2
    num_t = 2**m

    # Group tournaments by score sequence for faster iso testing
    score_groups = {}
    for bits in range(num_t):
        A = tournament_from_bits(n, bits)
        ss = score_seq(A, n)
        if ss not in score_groups:
            score_groups[ss] = []
        score_groups[ss].append(A)

    # Find iso class representatives
    iso_reps = []
    for ss, group in score_groups.items():
        reps = []
        for A in group:
            found = False
            for R in reps:
                if are_isomorphic(A, R, n):
                    found = True
                    break
            if not found:
                reps.append(A)
        iso_reps.extend(reps)

    total_orbits = sum(vertex_orbits(A, n) for A in iso_reps)
    target = 2 * factorial(n-1) if n >= 2 else 1
    match = "YES" if total_orbits == target else "NO"

    print(f"  n={n}: P({n}) = {total_orbits:6d}, 2*(n-1)! = {target:6d}, match={match}, {len(iso_reps)} classes")

print("\nDONE")
