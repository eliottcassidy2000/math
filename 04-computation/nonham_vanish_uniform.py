#!/usr/bin/env python3
"""
Test NONHAM=0 for position-uniform tournaments at n=7.

Since NONHAM != 0 for general tournaments (confirmed at n=3,4,5),
we now test whether NONHAM=0 specifically for POSITION-UNIFORM tournaments
at n=7 (extending the n=5 result).

Strategy: enumerate circulant tournaments (guaranteed vertex-transitive,
hence position-uniform at odd n), test NONHAM for each.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np
from collections import defaultdict
import time

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
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def circulant_tournament(n, gen_set):
    """Build circulant tournament on Z/nZ with edges i->j iff (j-i)%n in gen_set."""
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def check_nonham_for_pair(T, n, a, b):
    """Check if NONHAM(a,b) = 0 for given tournament T."""
    U = [v for v in range(n) if v != a and v != b]
    nonham_sum = 0

    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        for p in permutations(S_set):
            if p[-1] != a: continue
            p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
            if not p_valid: continue

            for q in permutations(R_set):
                if q[0] != b: continue
                q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                if not q_valid: continue

                concat = list(p) + list(q)
                is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                if not is_ham:
                    nonham_sum += sign

    return nonham_sum


# ============================================================
# n=5: ALL position-uniform (confirmation)
# ============================================================
print("=" * 70)
print("n=5: NONHAM=0 for ALL position-uniform tournaments?")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
uniform_count = 0
fail_count = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        continue

    uniform_count += 1
    for a in range(n):
        for b in range(n):
            if a == b: continue
            nh = check_nonham_for_pair(T, n, a, b)
            if nh != 0:
                fail_count += 1
                print(f"  FAIL: bits={bits}, ({a},{b}): NONHAM={nh}")
                break
        if fail_count > 0:
            break

print(f"  Position-uniform: {uniform_count}, fails: {fail_count}")
if fail_count == 0:
    print("  CONFIRMED: NONHAM=0 for ALL position-uniform n=5 tournaments")


# ============================================================
# n=7: Circulant tournaments (vertex-transitive => position-uniform at odd n)
# ============================================================
print("\n" + "=" * 70)
print("n=7: NONHAM=0 for circulant tournaments?")
print("=" * 70)

n = 7
# Generate all valid circulant generating sets
# For n=7, gen_set subset of {1,2,3} with complement {4,5,6} = {n-3,n-2,n-1}
# Must have |gen_set| = 3 (half of n-1=6 non-zero elements)
# And if d in gen_set then n-d not in gen_set (tournament property)
half = list(range(1, (n+1)//2))  # {1, 2, 3}
# Each subset S of half determines gen_set: for each d in half,
# either d in gen_set or n-d in gen_set
gen_sets = []
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.append(frozenset(gs))

# Remove duplicates
gen_sets = list(set(gen_sets))
print(f"  {len(gen_sets)} distinct circulant tournaments at n={n}")

t0 = time.time()
fail_count = 0
checked = 0

for gs in gen_sets:
    T = circulant_tournament(n, gs)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))

    checked += 1
    print(f"  gen_set={sorted(gs)}, H={H}, uniform={is_uniform}", end="")

    if not is_uniform:
        print(" (skipping)")
        continue

    # For circulant, by symmetry we only need to check pairs (0,b)
    found_fail = False
    for b in range(1, n):
        nh = check_nonham_for_pair(T, n, 0, b)
        if nh != 0:
            found_fail = True
            fail_count += 1
            print(f" FAIL: (0,{b}): NONHAM={nh}")
            break

    if not found_fail:
        print(" OK: NONHAM=0 for all pairs")

    elapsed = time.time() - t0
    print(f"    ({elapsed:.1f}s elapsed)")

print(f"\n  Checked {checked} circulant tournaments, {fail_count} failures")
if fail_count == 0:
    print("  CONFIRMED: NONHAM=0 for ALL circulant n=7 tournaments")


# ============================================================
# n=7: Also check non-circulant position-uniform (sample)
# ============================================================
print("\n" + "=" * 70)
print("n=7: Sampling non-circulant tournaments for position-uniformity")
print("=" * 70)

import random
random.seed(42)

n = 7
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
sampled = 0
uniform_found = 0
nonham_fails = 0
max_samples = 5000

for _ in range(max_samples):
    bits = random.randint(0, (1 << len(pairs)) - 1)
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    sampled += 1

    if not is_uniform:
        continue

    uniform_found += 1
    # Check a few pairs
    found_fail = False
    for a in range(n):
        for b in range(a+1, n):
            nh = check_nonham_for_pair(T, n, a, b)
            if nh != 0:
                found_fail = True
                nonham_fails += 1
                print(f"  FAIL: bits={bits}, ({a},{b}): NONHAM={nh}, H={H}")
                break
        if found_fail:
            break

    if uniform_found % 5 == 0:
        elapsed = time.time() - t0
        print(f"  ... {uniform_found} uniform found, {nonham_fails} fails ({elapsed:.1f}s)")

    if uniform_found >= 20:
        break

print(f"\n  Sampled {sampled}, found {uniform_found} uniform, {nonham_fails} NONHAM failures")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
