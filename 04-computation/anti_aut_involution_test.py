#!/usr/bin/env python3
"""
Test: Is every anti-automorphism of an SC tournament an involution?
===================================================================
Theory: sigma anti-aut => sigma^2 is an automorphism.
So sigma is involution iff sigma^2 = id, which holds when Aut(T) = {id}.
But Aut(T) can be non-trivial, so sigma^2 could be a non-trivial automorphism.

This script tests exhaustively at n=4,5,6,7 (7 is slow).

kind-pasteur-2026-03-06-S18e
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import tournament_from_bits, opposite_tournament
from itertools import permutations
from collections import defaultdict
import time

def canonical_form(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_self_converse(T):
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

def find_all_anti_auts(T):
    n = len(T)
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def find_all_auts(T):
    n = len(T)
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != T[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def compose(sigma, tau):
    return tuple(sigma[tau[i]] for i in range(len(sigma)))

def order_of(sigma):
    n = len(sigma)
    x = list(range(n))
    for k in range(1, 2 * n + 1):
        x = [sigma[xi] for xi in x]
        if all(x[i] == i for i in range(n)):
            return k
    return None

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    print(f"\n{'='*70}")
    print(f"n={n}: Exhaustive anti-automorphism analysis")
    print("=" * 70)

    seen = {}
    non_involution_count = 0
    involution_count = 0
    total_sc = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cf = canonical_form(T)
        if cf in seen:
            continue
        seen[cf] = True

        if not is_self_converse(T):
            continue
        total_sc += 1

        anti_auts = find_all_anti_auts(T)
        auts = find_all_auts(T)

        all_involutions = all(compose(s, s) == tuple(range(n)) for s in anti_auts)
        non_inv = [s for s in anti_auts if compose(s, s) != tuple(range(n))]

        if all_involutions:
            involution_count += 1
        else:
            non_involution_count += 1
            print(f"  NON-INVOLUTION FOUND! bits={bits}")
            print(f"    |Aut(T)| = {len(auts)}, |AntiAut(T)| = {len(anti_auts)}")
            for s in non_inv:
                s2 = compose(s, s)
                print(f"    sigma={s}, sigma^2={s2}, order(sigma)={order_of(s)}")

    print(f"\nn={n}: {total_sc} SC iso classes, {involution_count} all-involution, {non_involution_count} with non-involutions")

# n=7: sample random SC tournaments
import random
random.seed(42)

print(f"\n{'='*70}")
print("n=7: Random SC tournament anti-automorphism check")
print("=" * 70)

n = 7
non_inv_7 = 0
inv_7 = 0

# Find SC tournaments by random sampling
for trial in range(100):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1

    if not is_self_converse(T):
        continue

    anti_auts = find_all_anti_auts(T)
    all_inv = all(compose(s, s) == tuple(range(n)) for s in anti_auts)

    if all_inv:
        inv_7 += 1
    else:
        non_inv_7 += 1
        print(f"  NON-INV at trial {trial}")
        auts = find_all_auts(T)
        print(f"    |Aut|={len(auts)}, |AntiAut|={len(anti_auts)}")
        for s in anti_auts:
            s2 = compose(s, s)
            if s2 != tuple(range(n)):
                print(f"    sigma={s}, sigma^2={s2}, ord={order_of(s)}")
                break

print(f"\nn=7 random: {inv_7} all-involution, {non_inv_7} with non-involutions (out of {inv_7+non_inv_7} SC)")
print("\nDone.")
