#!/usr/bin/env python3
"""
Fast perspective counting: P(n) = sum over iso classes of (# vertex orbits).
Conjecture: P(n) = 2*(n-1)!

Uses nauty-style canonical form via sorted adjacency for speed.
kind-pasteur-2026-03-06-S25g
"""

from itertools import permutations, combinations

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

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def vertex_orbits_fast(A, n):
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

from math import factorial

print("PERSPECTIVE COUNTING: P(n) = sum of vertex orbits across iso classes")
print("=" * 60)

for n in range(2, 8):
    m = n*(n-1)//2
    if 2**m > 2**21:
        print(f"  n={n}: too large ({2**m} tournaments), skipping exhaustive")
        continue

    canon_to_A = {}
    for bits in range(2**m):
        A = tournament_from_bits(n, bits)
        c = canonical(A, n)
        if c not in canon_to_A:
            canon_to_A[c] = A

    total_orbits = 0
    for A in canon_to_A.values():
        total_orbits += vertex_orbits_fast(A, n)

    target = 2 * factorial(n-1) if n >= 2 else 1
    match = "YES" if total_orbits == target else "NO"
    print(f"  n={n}: P({n}) = {total_orbits}, 2*(n-1)! = {target}, match = {match}")
    print(f"         {len(canon_to_A)} iso classes")

print("\nDONE")
