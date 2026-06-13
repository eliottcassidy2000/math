#!/usr/bin/env python3
"""
f_poly_representation_theory.py — Representation-theoretic structure of F(T,x).

IDEA: F(T,x) = sum_P x^{fwd(P)} over n! permutations is an S_n character
evaluation in disguise.

For the TRANSITIVE tournament (total order), fwd(P) = # ascents of P.
The descent/ascent polynomial is the Eulerian polynomial A_n(x).
The Eulerian polynomial decomposes into S_n representations via
the "descent representation" (Solomon, 1976).

For a general tournament T, fwd(P) generalizes ascents.
F(T,x) = sum_{sigma in S_n} x^{des_T(sigma)} where des_T is the
T-descent statistic.

QUESTION: Does the S_n action on the coefficients F_k decompose nicely?
The k-th "level" is L_k = {P : fwd(P) = k}, with |L_k| = F_k.
S_n acts on L_k by conjugation: sigma.P = sigma ∘ P.
But fwd(sigma ∘ P) ≠ fwd(P) in general, so this doesn't preserve levels.

ALTERNATIVE: The Aut(T) action preserves F(T,x).
For sigma in Aut(T): fwd(sigma ∘ P) = fwd(P), so Aut(T) preserves levels.

NOVEL: F(T,x) as character of the REGULAR representation restricted
to Aut(T). The trace of the "x-weighted forward-edge matrix" gives F(T,x).

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
from collections import Counter

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def automorphisms(A, n):
    """Find automorphisms of tournament A."""
    auts = []
    for sigma in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if A[i][j] != A[sigma[i]][sigma[j]]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            auts.append(sigma)
    return auts

# ============================================================
# F(T,x) AND AUTOMORPHISM GROUP
# ============================================================
print("=" * 60)
print("F(T,x) AND AUTOMORPHISM STRUCTURE")
print("=" * 60)

n = 5
m = n*(n-1)//2
seen = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    auts = automorphisms(A, n)
    aut_size = len(auts)

    # Cycle types of automorphisms
    def cycle_type(sigma):
        visited = [False] * len(sigma)
        cycles = []
        for i in range(len(sigma)):
            if not visited[i]:
                cycle_len = 0
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = sigma[j]
                    cycle_len += 1
                cycles.append(cycle_len)
        return tuple(sorted(cycles, reverse=True))

    ctypes = Counter(cycle_type(a) for a in auts)

    print(f"  H={H:3d}: |Aut|={aut_size}, F={F}, cycle types={dict(ctypes)}")

# ============================================================
# NOVEL: F_k mod |Aut(T)| — does the automorphism group divide F_k?
# ============================================================
print("\n" + "=" * 60)
print("F_k mod |Aut(T)| — ORBIT COUNTING")
print("=" * 60)
# By Burnside/orbit-counting: if Aut(T) acts on L_k (perms with fwd=k),
# then |L_k| / |Aut(T)| = # orbits (if the action is free).
# In general: |L_k| = |Aut| * # orbits only if action is free.

n = 5
m = n*(n-1)//2
seen = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    aut_size = len(automorphisms(A, n))

    # Check divisibility
    divs = [F[k] % aut_size for k in range(n)]
    all_div = all(d == 0 for d in divs)

    print(f"  H={H:3d}: |Aut|={aut_size}, F mod |Aut| = {divs}, "
          f"{'ALL DIVIDE' if all_div else 'FAILS'}")

# ============================================================
# NOVEL: CHROMATIC QUASISYMMETRIC FUNCTION COMPARISON
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) AS TOURNAMENT DESCENT POLYNOMIAL — COMPARISON WITH")
print("CHROMATIC QUASISYMMETRIC FUNCTION")
print("=" * 60)
# Shareshian-Wachs: the chromatic quasisymmetric function X_G(x; q)
# has descent polynomial interpretation when G is the incomparability graph.
#
# For tournaments, the "descent" w.r.t. T is: position i is a descent
# if A[P_i][P_{i+1}] = 0 (backward edge). So fwd(P) = n-1 - des(P).
#
# The Rédei-Berge polynomial u_T(m) counts "m-colorings" of the tournament.
# The Grujic-Stojadinovic framework has u_T(m) = sum over compositions
# of n into m parts of products of H values.
#
# KEY: F(T,x) IS the "X-descent polynomial" in the Rédei-Berge Hopf algebra!
# This means F(T,x) carries the same information as the peak set,
# descent set, and shuffle statistics of a poset.

# Let me check: for n=5, how does F relate to the score sequence?
print("\nn=5 — F(T,x) vs score sequence:")
n = 5
m = n*(n-1)//2
score_to_F = {}

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    scores = sorted([sum(A[i]) for i in range(n)])
    F = compute_F(A, n)

    key_s = tuple(scores)
    if key_s not in score_to_F:
        score_to_F[key_s] = set()
    score_to_F[key_s].add(tuple(F))

for scores in sorted(score_to_F.keys()):
    F_set = score_to_F[scores]
    print(f"  scores={scores}: {len(F_set)} distinct F: {sorted(F_set)}")

# ============================================================
# NOVEL: F(T,x) FROBENIUS CHARACTERISTIC — does it decompose into
# Schur functions with nice coefficients?
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) COEFFICIENT ANALYSIS")
print("=" * 60)
# The Schur expansion of the Eulerian polynomial descent representation
# is known (Brenti, 1994). For tournaments, does F(T,x) have a similar
# expansion?
#
# Simpler question: is F_k (the coefficient) always a sum of binomial
# coefficients, or multinomial coefficients?

# At n=5, F has 8 distinct values. Let's see what determines F_k.
print("\nn=5 — What determines F_k?")
n = 5
m = n*(n-1)//2

def count_t3(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    return t3

seen = {}
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue

    H = F[n-1]
    t3 = count_t3(A, n)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    seen[key] = {'H': H, 't3': t3, 'scores': scores, 'F': F}

print(f"  {'H':>3} {'t3':>3} {'scores':>15} {'F':>25}")
for key, data in sorted(seen.items(), key=lambda x: x[1]['H']):
    print(f"  {data['H']:3d} {data['t3']:3d} {str(data['scores']):>15} {str(data['F']):>25}")

# Check: is F_1 = n*(n-1)/2 + something * t3?
print(f"\n  F_1 values: {sorted(set(d['F'][1] for d in seen.values()))}")
print(f"  Check if F_1 = a + b*t3:")
for key, data in sorted(seen.items(), key=lambda x: x[1]['t3']):
    print(f"    t3={data['t3']}, F_1={data['F'][1]}, F_1 - 2*n = {data['F'][1] - 2*n}")

# ============================================================
# NOVEL: MÖBIUS FUNCTION ON TOURNAMENT POSET
# ============================================================
print("\n" + "=" * 60)
print("TOURNAMENT WEAK ORDER AND MÖBIUS FUNCTION")
print("=" * 60)
# Define a partial order on {n} by: i ≤_T j if there's a directed path
# from i to j in T. For tournaments, this gives a total preorder
# (strong components are the equivalence classes, ordered linearly).
# The Möbius function of this poset gives inclusion-exclusion formulas.

# For a tournament, the strong components C_1 < C_2 < ... < C_r
# are linearly ordered. H(T) = prod H(T[C_i]) * multinomial coefficient.

# Actually: H(T) = n! / prod |C_i|! * prod H(T[C_i]) * multinomial...
# No. The correct formula: H(T) = (n! / (|C_1|! * |C_2|! * ... * |C_r|!))
#                                  * prod H(T[C_i])
# because paths must visit C_1 fully, then C_2 fully, etc.
# Wait: that's not right either. Within each C_i, the vertices can be
# visited in any Hamiltonian path order of T[C_i], and between components
# the order is fixed.
# So: H(T) = prod H(T[C_i]) * (n! / prod |C_i|!) ... no.
# Actually H(T) = prod H(T[C_i]) when components are arranged linearly.
# Because a Hamiltonian path must go through the components in order
# C_1, C_2, ..., C_r, and within each component, any HP of T[C_i] works.
# But the interleaving: actually no, the path must FIRST visit all of C_1,
# THEN all of C_2, etc., because all arcs go from earlier to later components.
# So H(T) = prod H(T[C_i]) exactly.

n = 5
m = n*(n-1)//2
print(f"\nn=5 — strong component decomposition:")

for bits in [0, 7, 15, 31, (1<<m)-1]:
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    H = F[n-1]

    # Find strong components via simple DFS
    # (for small n, just check all subsets)
    # A vertex v is in a strong component with u if u→...→v and v→...→u
    reach = [[False]*n for _ in range(n)]
    for i in range(n):
        reach[i][i] = True
    for _ in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][j]:
                    for k in range(n):
                        if A[j][k]:
                            reach[i][k] = True

    # Strong components: i and j in same iff reach[i][j] and reach[j][i]
    comp = [-1] * n
    comp_id = 0
    for i in range(n):
        if comp[i] == -1:
            comp[i] = comp_id
            for j in range(i+1, n):
                if reach[i][j] and reach[j][i]:
                    comp[j] = comp_id
            comp_id += 1

    comp_sizes = Counter(comp)
    print(f"  bits={bits:5d}: H={H:3d}, components={sorted(comp_sizes.values(), reverse=True)}")
