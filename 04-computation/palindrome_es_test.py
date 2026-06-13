#!/usr/bin/env python3
"""
TEST: Is e_s palindromic (e_s = e_{|U|-s}) for all non-edge pairs
in position-uniform tournaments?

If yes at odd n (where |U| = n-2 is odd), this directly proves NONHAM=0:
  sum (-1)^s e_s = sum pairs [(-1)^s + (-1)^{|U|-s}] e_s = 0 since |U| odd.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

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

def E_paths(T, verts, a):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == a else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != a: continue
        if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            count += 1
    return count

def B_paths(T, verts, b):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == b else 0
    count = 0
    for p in permutations(verts):
        if p[0] != b: continue
        if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            count += 1
    return count

def compute_es(T, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    e_s = [0] * (len(U) + 1)
    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        s = len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})
        ea = E_paths(T, S_set, a)
        bb = B_paths(T, R_set, b)
        e_s[s] += ea * bb
    return e_s


# ============================================================
# n=5: ALL position-uniform tournaments
# ============================================================
print("=" * 70)
print("n=5: Palindrome test for ALL position-uniform tournaments")
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

            e_s = compute_es(T, n, a, b)
            U_len = n - 2
            total_pairs += 1

            is_palindrome = all(e_s[s] == e_s[U_len - s] for s in range(U_len + 1))
            if is_palindrome:
                palindrome_count += 1
            else:
                non_palindrome_count += 1
                if non_palindrome_count <= 5:
                    print(f"  bits={bits}: ({a},{b}) e_s={e_s} NOT palindromic")

print(f"\n  Total non-edge pairs: {total_pairs}")
print(f"  Palindromic: {palindrome_count}")
print(f"  Non-palindromic: {non_palindrome_count}")


# ============================================================
# n=7: ALL circulant tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=7: Palindrome test for ALL circulant tournaments")
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
total_pairs = 0

for gs in sorted(gen_sets):
    T = circulant_tournament(n, gs)

    # By circulant symmetry, check only a=0
    for b in range(1, n):
        if T[(0,b)] == 1: continue

        e_s = compute_es(T, n, 0, b)
        U_len = n - 2
        total_pairs += 1

        is_palindrome = all(e_s[s] == e_s[U_len - s] for s in range(U_len + 1))
        if is_palindrome:
            palindrome_count += 1
        else:
            non_palindrome_count += 1
            if non_palindrome_count <= 5:
                print(f"  gs={sorted(gs)}: (0,{b}) e_s={e_s} NOT palindromic")

    # Also check: does palindrome hold for NON-edges only, or also edges?
    for b in range(1, n):
        if T[(0,b)] == 0: continue  # Only edges

        e_s = compute_es(T, n, 0, b)
        U_len = n - 2
        is_palindrome = all(e_s[s] == e_s[U_len - s] for s in range(U_len + 1))
        if not is_palindrome:
            if non_palindrome_count <= 10:
                alt_sum = sum((-1)**s * e_s[s] for s in range(len(e_s)))
                print(f"  gs={sorted(gs)}: (0,{b}) T[0,{b}]=1: e_s={e_s} NOT palindromic, alt_sum={alt_sum}")

print(f"\n  Total non-edge pairs: {total_pairs}")
print(f"  Palindromic: {palindrome_count}")
print(f"  Non-palindromic: {non_palindrome_count}")

if non_palindrome_count == 0:
    print("""
  THEOREM CANDIDATE: For position-uniform tournaments at odd n,
  the size-stratified E*B sums are palindromic:

    e_s = e_{|U|-s}   where e_s = sum_{|S|=s} E_a(S+a) * B_b(R+b)

  Since |U| = n-2 is odd, palindrome symmetry directly gives:
    M[a,b] = sum_s (-1)^s e_s = 0

  This is the missing step in the proof chain!
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
