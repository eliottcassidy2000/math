#!/usr/bin/env python3
"""
Detailed perspective counting at n=6 to understand why P(6)=296 != 240.
Also check: what IS the pattern?
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

def aut_group_size(A, n):
    count = 0
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            if not is_aut: break
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    is_aut = False
                    break
        if is_aut:
            count += 1
    return count

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

def is_self_converse(A, n):
    Aop = [[A[j][i] for j in range(n)] for i in range(n)]
    return are_isomorphic(A, Aop, n)

# Detailed analysis for n=3,4,5,6
for n in range(3, 7):
    m = n*(n-1)//2
    num_t = 2**m

    score_groups = {}
    for bits in range(num_t):
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

    print(f"\n{'='*60}")
    print(f"n={n}: {len(iso_reps)} iso classes")
    print(f"{'='*60}")

    orbit_dist = {}  # orbits -> count of classes with that many orbits
    total = 0
    for A in iso_reps:
        orb = vertex_orbits(A, n)
        aut = aut_group_size(A, n)
        sc = is_self_converse(A, n)
        ss = score_seq(A, n)
        total += orb
        if orb not in orbit_dist:
            orbit_dist[orb] = 0
        orbit_dist[orb] += 1
        if n <= 5:
            print(f"  orbits={orb}, |Aut|={aut}, SC={'Y' if sc else 'N'}, scores={ss}")

    print(f"\n  Orbit distribution: {dict(sorted(orbit_dist.items()))}")
    print(f"  Total orbits P({n}) = {total}")
    print(f"  2*(n-1)! = {2*factorial(n-1)}")

    # Weighted sum: sum of orbits*aut_size?
    weighted = sum(vertex_orbits(A, n) * aut_group_size(A, n) for A in iso_reps)
    print(f"  sum(orbits * |Aut|) = {weighted}")
    print(f"  n! = {factorial(n)}")

    # Another formula: sum of n/|Aut| = # labeled / n! * n = ...
    sum_n_over_aut = sum(n / aut_group_size(A, n) for A in iso_reps)
    print(f"  sum(n/|Aut|) = {sum_n_over_aut:.1f}")
    print(f"  2^C(n,2)/n! * n = {2**m / factorial(n) * n:.1f}")

print("\nDONE")
