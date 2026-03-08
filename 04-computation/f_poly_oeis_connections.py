#!/usr/bin/env python3
"""
f_poly_oeis_connections.py — Search for OEIS connections in F(T,x) data.

Compute various sequences derived from F(T,x) and check against known sequences.

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
from functools import reduce
from math import gcd

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

# ============================================================
# SEQUENCE 1: F_k for the REGULAR tournament
# ============================================================
print("=" * 60)
print("F_k FOR SPECIAL TOURNAMENTS")
print("=" * 60)

# Transitive tournament (all forward)
print("\nTransitive tournament F(T,x) coefficients:")
for n in [3, 4, 5, 6, 7]:
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F(A, n)
    print(f"  n={n}: F = {F}")
    # These should be reversed Eulerian numbers

# Regular tournament (when it exists at odd n)
print("\nRegular tournament (Paley-like) F(T,x) coefficients:")
for n in [3, 5, 7]:
    m = n*(n-1)//2
    # Find a regular tournament
    best_F = None
    best_H = 0
    for bits in (range(1<<m) if n <= 5 else [random.getrandbits(m) for _ in range(500)]):
        A = tournament_from_bits(bits, n)
        scores = [sum(A[i]) for i in range(n)]
        if all(s == (n-1)//2 for s in scores):
            F = compute_F(A, n)
            if F[n-1] > best_H:
                best_H = F[n-1]
                best_F = F
    if best_F:
        print(f"  n={n}: F = {best_F} (H={best_F[n-1]})")

# ============================================================
# SEQUENCE 2: F(T, 2) for transitive tournament
# ============================================================
print("\n" + "=" * 60)
print("F(T, 2) FOR TRANSITIVE TOURNAMENT")
print("=" * 60)
# For transitive T, fwd(P) = # ascents of P = n-1-des(P)
# F_k = A(n, n-1-k) (reversed Eulerian)
# F(T, 2) = sum A(n, k) * 2^{n-1-k} = ??

for n in range(2, 9):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F(A, n)
    F2 = sum(F[k] * 2**k for k in range(n))
    F3 = sum(F[k] * 3**k for k in range(n))
    Fhalf = sum(F[k] * (0.5)**k for k in range(n))
    print(f"  n={n}: F(2)={F2}, F(3)={F3}, F(1/2)*2^(n-1)={Fhalf * 2**(n-1):.0f}")

print("\n  F(transitive, 2) sequence: ", end="")
for n in range(2, 9):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F(A, n)
    F2 = sum(F[k] * 2**k for k in range(n))
    print(f"{F2}, ", end="")
print()
# This should be: sum A(n,k) * 2^{n-1-k} = 2^{n-1} * sum A(n,k) / 2^k
# = 2^{n-1} * A_n(1/2) where A_n(x) = Eulerian polynomial

# ============================================================
# SEQUENCE 3: Maximum F_1 (coefficient of x) across all tournaments
# ============================================================
print("\n" + "=" * 60)
print("MAX/MIN F_k ACROSS ALL TOURNAMENTS")
print("=" * 60)

for n in [3, 4, 5, 6, 7]:
    m = n*(n-1)//2

    if n <= 6:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(500)]

    seen = set()
    all_F = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    for k in range(n):
        min_fk = min(F[k] for F in all_F)
        max_fk = max(F[k] for F in all_F)
        print(f"  n={n}, k={k}: F_k in [{min_fk}, {max_fk}]", end="")
        if k == 0 or k == n-1:
            print(f"  (= H range)")
        else:
            print()

# ============================================================
# SEQUENCE 4: F(T, -1) = S(T)*(-1)^{n-1} distribution
# ============================================================
print("\n" + "=" * 60)
print("S(T) = SIGNED HP PERMANENT DISTRIBUTION")
print("=" * 60)

for n in [3, 4, 5, 6, 7]:
    m = n*(n-1)//2

    if n <= 6:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(500)]

    seen = set()
    S_vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        S = (-1)**(n-1) * sum(F[k] * (-1)**k for k in range(n))
        S_vals.append(S)

    S_set = sorted(set(S_vals))
    if n <= 4:
        print(f"  n={n}: S values = {S_set}")
    else:
        print(f"  n={n}: S range = [{min(S_set)}, {max(S_set)}], "
              f"# distinct = {len(S_set)}, "
              f"GCD = {reduce(gcd, [abs(s) for s in S_set if s != 0]) if any(s != 0 for s in S_set) else 0}")

# ============================================================
# SEQUENCE 5: Number of distinct F polynomials
# ============================================================
print("\n" + "=" * 60)
print("NUMBER OF DISTINCT F(T,x) POLYNOMIALS")
print("=" * 60)

for n in range(2, 8):
    m = n*(n-1)//2
    seen = set()
    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        seen.add(tuple(F))
    print(f"  n={n}: {len(seen)} distinct F polynomials (out of {1<<m} tournaments)")

# ============================================================
# SEQUENCE 6: F_1 (number of perms with exactly 1 forward edge)
# ============================================================
print("\n" + "=" * 60)
print("F_1 SEQUENCE FOR TRANSITIVE TOURNAMENT")
print("=" * 60)
# For transitive T: F_1 = A(n, n-2) = Eulerian number
for n in range(2, 9):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F(A, n)
    print(f"  n={n}: F_1 = {F[1] if n > 1 else 'N/A'}")

# ============================================================
# NOVEL: F(T,x) GENERATING FUNCTION IN n
# ============================================================
print("\n" + "=" * 60)
print("SUM OF F_k OVER ALL TOURNAMENTS AT n (generating function in n)")
print("=" * 60)
# For each n, compute sum_{all T} F_k(T). This is a sequence in n.
# By symmetry over T, each permutation P has the same total fwd count
# (since choosing T randomly means each edge is forward with prob 1/2).
# So sum_{T} fwd(P for fixed P) = C(n-1, 2) * 2^{m-1} if edges are independent.
# Actually: for a fixed permutation P, each of the n-1 consecutive pairs
# has its direction independently chosen across tournaments. The non-consecutive
# pairs don't affect fwd. Wait: the tournament is defined on ALL C(n,2) pairs,
# but fwd(P) only depends on the n-1 consecutive pairs in P.
# For a fixed P, fwd(P) depends on n-1 edges, each of which is "forward"
# in exactly half the tournaments (by symmetry of bit encoding).
# So E[fwd(P)] = (n-1)/2 for each P, and sum_{T} fwd(P) = (n-1)/2 * 2^m.
# And sum_{T} sum_P x^{fwd(P)} = n! * sum_{T} average_P(x^{fwd}).

for n in range(2, 7):
    m = n*(n-1)//2
    total_by_k = [0] * n
    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        for k in range(n):
            total_by_k[k] += F[k]

    # Normalize by 2^m
    avg_by_k = [total_by_k[k] / (1 << m) for k in range(n)]
    print(f"  n={n}: total F_k = {total_by_k}")
    print(f"         avg F_k = {[f'{a:.1f}' for a in avg_by_k]}")
    print(f"         avg = n!*C(n-1,k)/2^(n-1)? = {[f'{math.factorial(n)*math.comb(n-1,k)/2**(n-1):.1f}' for k in range(n)]}")
    # Check: avg F_k = n! * C(n-1, k) / 2^{n-1}
    # This would mean the average F polynomial is n!/2^{n-1} * (1+x)^{n-1}
    # = n!/2^{n-1} * sum C(n-1,k) x^k
    expected = [math.factorial(n) * math.comb(n-1, k) / 2**(n-1) for k in range(n)]
    match = all(abs(avg_by_k[k] - expected[k]) < 0.01 for k in range(n))
    print(f"         match: {match}")
