#!/usr/bin/env python3
"""
alpha3_structure.py -- kind-pasteur-2026-03-13-S60

At Paley p=11: alpha_3 = 1155. These are triples of directed cycles with
pairwise disjoint vertex sets. 

Possible size combinations (k1,k2,k3) with k1+k2+k3 <= 11, all odd:
- (3,3,3): 9 vertices used, 2 left over
- (3,3,5): 11 vertices used (exact partition of Z_11)
- No (3,5,5) since 13 > 11
- No (5,5,5) since 15 > 11

So alpha_3 decomposes into (3,3,3) and (3,3,5) contributions.

The (3,3,5) triples form EXACT PARTITIONS of Z_11! 
This is a packing/covering structure.
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


p = 11
S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
A = build_adj(p, S_qr)

print(f"Paley p={p}")

# Get all active vertex sets
active = {}
for k in range(3, p+1, 2):
    for subset in combinations(range(p), k):
        fs = frozenset(subset)
        nc = count_ham_cycles(A, list(subset))
        if nc > 0:
            active[fs] = (k, nc)

# Separate by size
active_3 = [(fs, nc) for fs, (k, nc) in active.items() if k == 3]
active_5 = [(fs, nc) for fs, (k, nc) in active.items() if k == 5]

print(f"  Active 3-sets: {len(active_3)}, active 5-sets: {len(active_5)}")

# Count (3,3,3) triples: three pairwise disjoint active 3-sets
alpha3_333 = 0
triple_count_333 = 0
for i in range(len(active_3)):
    for j in range(i+1, len(active_3)):
        if active_3[i][0] & active_3[j][0]:
            continue
        for l in range(j+1, len(active_3)):
            if (active_3[i][0] & active_3[l][0]) or (active_3[j][0] & active_3[l][0]):
                continue
            n1, n2, n3 = active_3[i][1], active_3[j][1], active_3[l][1]
            alpha3_333 += n1 * n2 * n3
            triple_count_333 += 1

print(f"\n  (3,3,3) vertex-set triples: {triple_count_333}")
print(f"  alpha3_333 = {alpha3_333} (weighted by cycle counts)")

# What vertices are LEFT OVER in (3,3,3) triples?
leftover = defaultdict(int)
for i in range(len(active_3)):
    for j in range(i+1, len(active_3)):
        if active_3[i][0] & active_3[j][0]:
            continue
        for l in range(j+1, len(active_3)):
            if (active_3[i][0] & active_3[l][0]) or (active_3[j][0] & active_3[l][0]):
                continue
            used = active_3[i][0] | active_3[j][0] | active_3[l][0]
            left = frozenset(range(p)) - used
            leftover[left] += 1

print(f"  Leftover vertex sets: {dict((tuple(sorted(k)), v) for k, v in leftover.items())}")

# Count (3,3,5) triples: two disjoint active 3-sets + one active 5-set covering the rest
alpha3_335 = 0
triple_count_335 = 0
for i in range(len(active_3)):
    for j in range(i+1, len(active_3)):
        if active_3[i][0] & active_3[j][0]:
            continue
        # The 5-set must be the complement of the two 3-sets
        remaining = frozenset(range(p)) - active_3[i][0] - active_3[j][0]
        if len(remaining) == 5 and remaining in dict(active_5 + [(frozenset(), 0)]):
            # Check if remaining is active
            for fs, nc in active_5:
                if fs == remaining:
                    n1, n2 = active_3[i][1], active_3[j][1]
                    alpha3_335 += n1 * n2 * nc
                    triple_count_335 += 1
                    break

# Actually, the 5-set doesn't have to be the exact complement. 
# It just needs to be disjoint from both 3-sets. 
# With 3+3+5=11 = p, the 5-set IS the complement.
# But what about 3+3+5 < 11? That's impossible since 3+3+5=11=p.

# Also need (3,5,3) and (5,3,3) orderings? No: we're counting UNORDERED triples.
# But the cycle types differ: two 3-cycles and one 5-cycle.
# The formula counts: choose unordered pair of 3-sets, then one 5-set.
# For alpha_3 as unordered triples of cycles: {C1, C2, C3} with pairwise disjoint vertex sets.
# When two are from size-3 and one from size-5: number of such triples = 
# (# disjoint 3-set pairs) * (# of active 5-sets = complement) * product of cycle counts.

print(f"\n  (3,3,5) vertex-set triples (exact Z_11 partitions): {triple_count_335}")
print(f"  alpha3_335 = {alpha3_335}")

print(f"\n  Total alpha_3 = {alpha3_333 + alpha3_335}")
print(f"  Expected: 1155")

# Decomposition of the H contribution
H_333 = 8 * alpha3_333
H_335 = 8 * alpha3_335
print(f"\n  8 * alpha3_333 = {H_333}")
print(f"  8 * alpha3_335 = {H_335}")
print(f"  Total alpha_3 contribution to H: {H_333 + H_335}")

# Now let's understand the (3,3,5) partitions: how many ways can we partition
# Z_11 into two 3-sets and one 5-set, all active?
print(f"\n\n=== Z_11 PARTITIONS INTO (3,3,5) ===")

# Count disjoint 3-set pairs
disj_3_pairs = []
for i in range(len(active_3)):
    for j in range(i+1, len(active_3)):
        if not (active_3[i][0] & active_3[j][0]):
            disj_3_pairs.append((i, j))

print(f"  Disjoint active 3-set pairs: {len(disj_3_pairs)}")

# For each pair, check if complement is an active 5-set
active_5_dict = {fs: nc for fs, nc in active_5}
partitions = []
for i, j in disj_3_pairs:
    comp = frozenset(range(p)) - active_3[i][0] - active_3[j][0]
    if comp in active_5_dict:
        partitions.append((active_3[i], active_3[j], (comp, active_5_dict[comp])))

print(f"  Partitions where complement is active 5-set: {len(partitions)}")
print(f"  Partitions where complement is NOT active: {len(disj_3_pairs) - len(partitions)}")

# Show a few partitions
for p_data in partitions[:5]:
    (v1, n1), (v2, n2), (v3, n3) = p_data
    print(f"    {sorted(v1)} + {sorted(v2)} + {sorted(v3)}: n=({n1},{n2},{n3})")

# By Z_p symmetry: how many orbits of partitions?
partition_orbits = set()
for (v1, n1), (v2, n2), (v3, n3) in partitions:
    # Canonical: apply all Z_p translations, take min of sorted tuple of sorted sets
    for t in range(p):
        s1 = frozenset((v+t)%p for v in v1)
        s2 = frozenset((v+t)%p for v in v2)
        s3 = frozenset((v+t)%p for v in v3)
        canon = tuple(sorted([tuple(sorted(s1)), tuple(sorted(s2)), tuple(sorted(s3))]))
    partition_orbits.add(canon)

print(f"\n  Partition orbit types: {len(partition_orbits)}")
# Each orbit should have p=11 members
print(f"  Expected: {len(partitions)}/{p} = {len(partitions)/p} orbits")

print("\nDONE.")
