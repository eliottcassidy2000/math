#!/usr/bin/env python3
"""
KEY TEST: Are the consec(a,b,j) values ALWAYS palindromic
for position-uniform tournaments?

consec(a,b,j) = #{Ham paths in T' = T + {a->b} where a at position j, b at j+1}

If consec(a,b,j) = consec(a,b,n-2-j) (palindrome) and n is odd,
then sum_j (-1)^j consec(a,b,j) = 0, proving M[a,b] = 0.

NOTE: This is different from the e_s test! e_s aggregates by SUBSET SIZE,
while consec aggregates by POSITION. They give different distributions.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T


def compute_consec(T, n, a, b):
    """Compute consec(a,b,j) for T' = T + {a->b}."""
    Tp = dict(T)
    Tp[(a,b)] = 1

    consec = [0] * (n - 1)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if Tp.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        if prod > 0:
            for j in range(n-1):
                if perm[j] == a and perm[j+1] == b:
                    consec[j] += 1
    return consec


# ============================================================
# n=5: ALL position-uniform tournaments
# ============================================================
print("=" * 70)
print("n=5: Consec palindrome test for ALL position-uniform tournaments")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
palindrome_count = 0
non_palindrome_count = 0
total_pairs = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        continue

    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T[(a,b)] == 1: continue

            consec = compute_consec(T, n, a, b)
            total_pairs += 1

            is_pal = all(consec[j] == consec[n-2-j] for j in range(n-1))
            if is_pal:
                palindrome_count += 1
            else:
                non_palindrome_count += 1
                if non_palindrome_count <= 5:
                    alt = sum((-1)**j * consec[j] for j in range(n-1))
                    print(f"  bits={bits}: ({a},{b}) consec={consec}, alt={alt}")

print(f"\n  Total non-edge pairs: {total_pairs}")
print(f"  Palindromic consec: {palindrome_count}")
print(f"  Non-palindromic consec: {non_palindrome_count}")

if non_palindrome_count == 0:
    print("  ALL consec vectors are palindromic!")


# ============================================================
# Also check: non-uniform tournaments for comparison
# ============================================================
print("\n" + "=" * 70)
print("n=5: Consec palindrome for NON-uniform tournaments (sample)")
print("=" * 70)

non_pal_non_uniform = 0
pal_non_uniform = 0
checked = 0
import random
random.seed(42)

for _ in range(100):
    bits = random.randint(0, (1 << len(pairs)) - 1)
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if is_uniform:
        continue

    checked += 1
    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T[(a,b)] == 1: continue

            consec = compute_consec(T, n, a, b)
            is_pal = all(consec[j] == consec[n-2-j] for j in range(n-1))
            if is_pal:
                pal_non_uniform += 1
            else:
                non_pal_non_uniform += 1
                if non_pal_non_uniform <= 3:
                    alt = sum((-1)**j * consec[j] for j in range(n-1))
                    print(f"  bits={bits}: ({a},{b}) consec={consec}, alt={alt}")

    if checked >= 20:
        break

print(f"  Non-uniform checked: {checked}")
print(f"  Palindromic: {pal_non_uniform}, Non-palindromic: {non_pal_non_uniform}")


# ============================================================
# n=7: Circulant tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=7: Consec palindrome for circulant tournaments")
print("=" * 70)

n = 7
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

palindrome_count = 0
non_palindrome_count = 0

for gs in sorted(gen_sets):
    T = circulant_tournament(n, gs)

    # By circulant symmetry, check (0, b) for non-edges
    for b in range(1, n):
        if T[(0,b)] == 1: continue

        consec = compute_consec(T, n, 0, b)
        is_pal = all(consec[j] == consec[n-2-j] for j in range(n-1))

        if is_pal:
            palindrome_count += 1
        else:
            non_palindrome_count += 1
            alt = sum((-1)**j * consec[j] for j in range(n-1))
            print(f"  gs={sorted(gs)}: (0,{b}) consec={consec}, alt={alt}")

print(f"\n  Palindromic: {palindrome_count}, Non-palindromic: {non_palindrome_count}")

if non_palindrome_count == 0:
    print("""
  THEOREM CANDIDATE (CONSEC PALINDROME):
  For ALL position-uniform tournaments at odd n,
  the consecutive position counts are palindromic:

    consec(a,b,j) = consec(a,b,n-2-j)

  where consec(a,b,j) = #{Ham paths in T + {a->b} with a at pos j, b at pos j+1}.

  Since n-1 is even (odd n), the palindrome implies:
    sum_j (-1)^j consec(a,b,j) = 0

  This proves M[a,b] = 0 for all non-edges, hence M = (H/n)*I.
""")


print("=" * 70)
print("DONE")
print("=" * 70)
