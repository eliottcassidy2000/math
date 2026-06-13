#!/usr/bin/env python3
"""
f_poly_betti_deep.py - Deep investigation of F-polynomial vs path homology.

KEY FINDING from tournament_quasisym_betti.py:
  At n=5: F=[9,30,42,30,9] has β₁=0 (120 tournaments) AND β₁=1 (120 tournaments)
  At n=6: β₃>0 has two F-polys: [45,117,198,198,117,45] (31) and [9,81,270,270,81,9] (8)

So F-polynomial does NOT determine Betti numbers!
But: some F-polynomials FORCE specific Betti numbers.

QUESTION: What ADDITIONAL invariant, combined with F, determines β?
Candidates: t₃ (3-cycle count), score sequence, Pfaffian, ...

At n=5: For F=[9,30,42,30,9], what distinguishes β₁=0 from β₁=1?

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict
import numpy as np

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def score_seq(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def f_polynomial(A, n):
    counts = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        counts[fwd] += 1
    return tuple(counts)

def is_strongly_connected(A, n):
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] and u not in visited:
                visited.add(u)
                queue.append(u)
    if len(visited) < n:
        return False
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] and u not in visited:
                visited.add(u)
                queue.append(u)
    return len(visited) == n

# === n=5 EXHAUSTIVE: What distinguishes F=[9,30,42,30,9] tournaments? ===
print("=" * 60)
print("n=5: F=[9,30,42,30,9] — what distinguishes β₁?")
print("=" * 60)

n = 5
m = n*(n-1)//2
target_F = (9, 30, 42, 30, 9)

b1_0_data = []
b1_1_data = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    F = f_polynomial(A, n)
    if F != target_F:
        continue

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    t3 = count_t3(A, n)
    sc = is_strongly_connected(A, n)
    score = score_seq(A, n)

    data = {'t3': t3, 'sc': sc, 'score': score, 'bits': bits}
    if b1 == 0:
        b1_0_data.append(data)
    else:
        b1_1_data.append(data)

print(f"\nF=[9,30,42,30,9]: {len(b1_0_data)} with β₁=0, {len(b1_1_data)} with β₁=1")

# Compare t₃
t3_b0 = Counter(d['t3'] for d in b1_0_data)
t3_b1 = Counter(d['t3'] for d in b1_1_data)
print(f"\n  β₁=0: t₃ dist = {dict(sorted(t3_b0.items()))}")
print(f"  β₁=1: t₃ dist = {dict(sorted(t3_b1.items()))}")

# Compare score sequence
sc_b0 = Counter(d['score'] for d in b1_0_data)
sc_b1 = Counter(d['score'] for d in b1_1_data)
print(f"\n  β₁=0: score dist = {dict(sorted(sc_b0.items()))}")
print(f"  β₁=1: score dist = {dict(sorted(sc_b1.items()))}")

# Compare strong connectivity
ssc_b0 = Counter(d['sc'] for d in b1_0_data)
ssc_b1 = Counter(d['sc'] for d in b1_1_data)
print(f"\n  β₁=0: SC dist = {dict(ssc_b0)}")
print(f"  β₁=1: SC dist = {dict(ssc_b1)}")

# AH-HA: score (2,2,2,2,2) is the regular tournament at n=5
# Regular tournaments exist only at odd n
# n=5 has C₅ (the Paley tournament) and its complement as the only regular tournaments
# C₅: score (2,2,2,2,2), t₃=5, β₁=1
# So score determines things within this F-class?

# But wait, F = [9,30,42,30,9]. Let me check: H = F_4 = 9.
# Regular tournaments at n=5 have H = 15 (max), not 9.
# So this F doesn't correspond to regular tournaments.

# Let me also check: what F-poly do regular 5-tournaments have?
print("\n\nRegular n=5 tournament (C₅) check:")
A_c5 = [[0]*5 for _ in range(5)]
for i in range(5):
    A_c5[i][(i+1)%5] = 1
    A_c5[i][(i+2)%5] = 1
F_c5 = f_polynomial(A_c5, 5)
beta_c5 = path_betti_numbers(A_c5, 5)
print(f"  F = {list(F_c5)}")
print(f"  β = {[int(b) for b in beta_c5]}")
print(f"  t₃ = {count_t3(A_c5, 5)}")

# === What's the FULL classification at n=5? ===
print("\n" + "=" * 60)
print("n=5: FULL (F-poly, t₃, SC, score) → β₁ classification")
print("=" * 60)

n = 5
full_class = defaultdict(Counter)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    t3 = count_t3(A, n)
    sc = is_strongly_connected(A, n)
    score = score_seq(A, n)

    key = (t3, sc, score)
    full_class[key][b1] += 1

print("\n(t₃, SC, score) → β₁ distribution:")
for key in sorted(full_class.keys()):
    dist = dict(full_class[key])
    if len(dist) > 1 or 1 in dist:  # Mixed or has β₁=1
        print(f"  t₃={key[0]}, SC={key[1]}, score={key[2]}: {dist}")

# KEY: Is (t₃, score) sufficient to determine β₁?
print("\n\nIs (t₃, score) a COMPLETE β₁ determinant at n=5?")
ambiguous = 0
for key, dist in full_class.items():
    if len(dist) > 1:
        ambiguous += 1
        print(f"  AMBIGUOUS: t₃={key[0]}, SC={key[1]}, score={key[2]}: {dict(dist)}")
if ambiguous == 0:
    print("  YES! (t₃, SC, score) completely determines β₁ at n=5.")
