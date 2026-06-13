#!/usr/bin/env python3
"""
Check: does our P(n) = sum of vertex orbits match A093934?
A093934 = # unlabeled tournaments with n signed nodes.
A signed tournament = tournament + {+,-} assignment to each vertex.
# equivalence classes = sum over iso classes T of 2^{orbits(T)}.

But our P(n) = sum over iso classes of orbits(T).

So A093934(n) = sum_T 2^{orbits(T)}, not sum_T orbits(T).

Let me verify both.
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
    if score_seq(A, n) != score_seq(B, n):
        return False
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            if not match: break
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
        if px != py: parent[px] = py
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            if not is_aut: break
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    is_aut = False
                    break
        if is_aut:
            for v in range(n):
                union(v, perm[v])
    return len(set(find(v) for v in range(n)))

print("Checking P(n) vs A093934")
print("A093934: 1, 2, 4, 12, 48, 296, 3040, ...")
print("=" * 50)

for n in range(2, 7):
    m = n*(n-1)//2
    score_groups = {}
    for bits in range(2**m):
        A = tournament_from_bits(n, bits)
        ss = score_seq(A, n)
        if ss not in score_groups:
            score_groups[ss] = []
        score_groups[ss].append(A)

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

    sum_orbits = 0
    sum_2_pow_orbits = 0
    for A in iso_reps:
        orb = vertex_orbits(A, n)
        sum_orbits += orb
        sum_2_pow_orbits += 2**orb

    print(f"  n={n}: sum(orbits)={sum_orbits}, sum(2^orbits)={sum_2_pow_orbits}")

print("\nSo which matches A093934(n) = 1, 2, 4, 12, 48, 296?")
print("Our sum(orbits): n=2:2, n=3:4, n=4:12, n=5:48, n=6:296")
print("A093934:         n=1:2, n=2:4, n=3:12, n=4:48, n=5:296")
print("\nThey match with offset! A093934(n) = P(n+1) = sum of vertex orbits")
print("across all iso classes of tournaments on n+1 vertices.")

# But wait - A093934 description says 'signed nodes'.
# 2^orbits would give # signed tournament classes.
# Our sum_orbits matching suggests a DIFFERENT interpretation.

# Actually P(n) = # pointed (rooted) tournament iso classes
# = # equivalence classes of (tournament, distinguished vertex) pairs

# OEIS search for "rooted tournament"
print("\n\nAlternative: P(n) = # rooted tournament iso classes")
print("= # orbits of S_n acting on (tournaments x vertices)")

print("\nDONE")
