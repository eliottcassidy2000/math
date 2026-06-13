#!/usr/bin/env python3
"""
Detailed analysis of M=3*I for non-VT n=5 tournaments.

At n=5:
  scores=(2,2,2,2,2) H=15 |Aut|=5 VT, M=3*I  (24 such tournaments)
  scores=(1,2,2,2,3) H=15 |Aut|=3 NOT VT, M=3*I ???

Verify and understand WHY non-VT tournaments can have scalar M.

kind-pasteur-2026-03-06-S25b (continuation)
"""

from itertools import permutations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_entry(T, n, a, a)
        for b in range(a+1, n):
            M[a, b] = compute_M_entry(T, n, a, b)
            M[b, a] = M[a, b]
    return M

def check_automorphism(T, n, sigma):
    for i in range(n):
        for j in range(n):
            if i != j:
                if T.get((i,j), 0) != T.get((sigma[i], sigma[j]), 0):
                    return False
    return True

def get_automorphisms(T, n):
    auts = []
    for perm in permutations(range(n)):
        if check_automorphism(T, n, perm):
            auts.append(perm)
    return auts


# Find ALL n=5 tournaments with M=3*I
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

scalar_tournaments = []

for bits in range(1 << len(pairs)):
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = {}
    idx = 0
    for (i,j) in pairs:
        if b_list[idx]:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
        idx += 1

    H = count_H(T, n)
    if H != 15:
        continue

    # Quick check
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)
    if M00 != 3 or M01 != 0:
        continue

    M = compute_M(T, n)
    if np.array_equal(M, 3 * np.eye(n, dtype=int)):
        scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
        scalar_tournaments.append((scores, T.copy()))

print(f"Total n=5 tournaments with M=3*I: {len(scalar_tournaments)}")

# Classify by score sequence
from collections import Counter
score_counts = Counter(s for s, t in scalar_tournaments)
print(f"Score sequences: {dict(score_counts)}")

# Analyze the non-regular ones
for scores, T in scalar_tournaments:
    if scores != (2,2,2,2,2):
        out_degrees = tuple(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n))
        auts = get_automorphisms(T, n)

        print(f"\n  Non-regular M=3*I tournament:")
        print(f"    Out-degrees: {out_degrees}")
        print(f"    |Aut| = {len(auts)}")
        print(f"    Automorphisms: {auts}")

        # Show adjacency
        print(f"    Adjacency (i->j = 1):")
        for i in range(n):
            row = [T.get((i,j), 0) for j in range(n)]
            print(f"      {row}")

        # Check position distribution
        pos_count = np.zeros((n, n), dtype=int)
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                for k in range(n):
                    pos_count[perm[k], k] += 1

        print(f"    Position distribution P[v,k] (v=vertex, k=position):")
        for v in range(n):
            print(f"      v={v}: {list(pos_count[v,:])}")

        uniform = all(pos_count[v,k] == 3 for v in range(n) for k in range(n))
        print(f"    Uniform position: {uniform}")

        # Show the full M matrix
        M = compute_M(T, n)
        print(f"    M = {M.tolist()}")

        break  # Just show one example


# ============================================================
# KEY QUESTION: At n=5, is M=3*I the same as "ALL H=15"?
# ============================================================
print("\n" + "=" * 70)
print("n=5: Are ALL H=15 tournaments M=3*I?")
print("=" * 70)

total_H15 = 0
scalar_H15 = 0

for bits in range(1 << len(pairs)):
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = {}
    idx = 0
    for (i,j) in pairs:
        if b_list[idx]:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
        idx += 1

    H = count_H(T, n)
    if H != 15:
        continue
    total_H15 += 1

    M = compute_M(T, n)
    if np.array_equal(M, 3 * np.eye(n, dtype=int)):
        scalar_H15 += 1

print(f"  H=15 tournaments: {total_H15}")
print(f"  M=3*I tournaments: {scalar_H15}")
print(f"  ALL H=15 have M=3*I: {total_H15 == scalar_H15}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
