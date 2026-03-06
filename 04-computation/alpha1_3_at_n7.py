#!/usr/bin/env python3
"""
Check if alpha_1=3 is still impossible at n=7, even though the common-vertex
property fails.

At n=7, c3=3 tournaments CAN have triples spanning 7 vertices without a
common vertex. But do they still always have c5 >= 1 (or c7 >= 1)?

kind-pasteur-2026-03-06-S22
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations
import random

def count_directed_odd_cycles(T, n):
    total = 0
    by_length = {}
    for k in range(3, n+1, 2):
        count = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    count += c
        by_length[k] = count
        total += count
    return total, by_length

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

print("=" * 60)
print("n=7: Does c3=3 still force alpha_1 >= 4?")
print("=" * 60)

random.seed(42)
n = 7
m = n*(n-1)//2

tested = 0
alpha1_dist = {}

# Sample tournaments with c3=3
for trial in range(2000000):
    bits = random.randint(0, (1 << m) - 1)
    T = tournament_from_bits(n, bits)
    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    if sum(s*s for s in scores) != 85:  # c3=3 constraint
        continue

    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue

    tested += 1
    alpha1, by_length = count_directed_odd_cycles(T, n)
    alpha1_dist[alpha1] = alpha1_dist.get(alpha1, 0) + 1

    if alpha1 == 3:
        print(f"  ALPHA_1=3 FOUND! bits={bits}, by_length={by_length}")
        print(f"  triples={[sorted(t) for t in triples]}")
        print(f"  This would DISPROVE our theorem!")

    if tested % 100 == 0 and tested <= 500:
        print(f"  tested={tested}, current dist={alpha1_dist}", flush=True)

    if tested >= 2000:
        break

print(f"\nTested {tested} tournaments with c3=3 at n=7")
print(f"alpha_1 distribution: {dict(sorted(alpha1_dist.items()))}")
print(f"alpha_1=3 found? {3 in alpha1_dist}")

# Also check: what is the min alpha_1 among c3=3 tournaments?
if alpha1_dist:
    print(f"Min alpha_1: {min(alpha1_dist.keys())}")

print("\nDone.")
