#!/usr/bin/env python3
"""
partition_bibd_connection.py -- kind-pasteur-2026-03-13-S60

The (3,3,5) partitions of Z_11 connect to BIBD theory.

At Paley p=7: the 14 directed 3-cycles form a (7,3,1)-BIBD
(each pair of vertices appears in exactly 1 cyclic triple).
This maximizes alpha_1 = 80 and minimizes alpha_2 = 7.

At Paley p=11: 55 directed 3-cycles. Does the arrangement have
BIBD-like properties? In a (11,3,lambda)-BIBD, each pair appears
in lambda blocks. With 55 blocks and C(11,2)=55 pairs, lambda = 55*3/55 = 3.
Wait: for a (v,k,lambda)-BIBD, b*k*(k-1) = v*(v-1)*lambda.
55 * 3 * 2 = 330 = 11 * 10 * lambda => lambda = 3.

So at Paley p=11: every pair of vertices appears in exactly 3 cyclic triples!
This is a (11, 3, 3)-BIBD (also called a triple system).

CONNECTION TO ALPHA_3:
For the (3,3,3) triples: we need THREE pairwise disjoint blocks from the BIBD.
This is a PARALLEL CLASS question: does the BIBD have a resolvable structure?

For (3,3,5) partitions: we need two disjoint blocks whose union leaves an
active 5-set. The BIBD structure constrains which pairs are disjoint.

Test: cross-orientation comparison of partition counts.
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


# BIBD analysis at p=11
p = 11
m = (p - 1) // 2
N = 1 << m

print(f"p={p}, testing all {N} orientations")
print(f"{'='*70}")

for bits in range(N):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)

    # Count cyclic triples (active 3-sets)
    active_3 = []
    for subset in combinations(range(p), 3):
        nc = count_ham_cycles(A, list(subset))
        if nc > 0:
            active_3.append(frozenset(subset))

    c3 = len(active_3)

    # Check pair coverage (BIBD property)
    pair_count = defaultdict(int)
    for fs in active_3:
        verts = sorted(fs)
        for i in range(3):
            for j in range(i+1, 3):
                pair_count[(verts[i], verts[j])] += 1

    lambdas = list(pair_count.values())
    lambda_set = set(lambdas)

    # Count disjoint 3-set pairs
    disj_pairs = 0
    for i in range(len(active_3)):
        for j in range(i+1, len(active_3)):
            if not (active_3[i] & active_3[j]):
                disj_pairs += 1

    S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    is_paley = sorted(S) == S_qr

    tag = " <-- PALEY" if is_paley else ""
    if len(lambda_set) == 1:
        lam = list(lambda_set)[0]
        print(f"  bits={bits}: c3={c3}, lambda={lam} (BIBD!), "
              f"disj_pairs={disj_pairs}{tag}")
    else:
        print(f"  bits={bits}: c3={c3}, lambda={sorted(lambda_set)}, "
              f"disj_pairs={disj_pairs}{tag}")

# Detailed Paley analysis
print(f"\n\n{'='*70}")
print(f"DETAILED PALEY p=11 BIBD ANALYSIS")
print(f"{'='*70}")

S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
A = build_adj(p, S_qr)

active_3 = []
for subset in combinations(range(p), 3):
    nc = count_ham_cycles(A, list(subset))
    if nc > 0:
        active_3.append(frozenset(subset))

print(f"  Active 3-sets (blocks): {len(active_3)}")
print(f"  This is a ({p}, 3, lambda)-BIBD with lambda = {len(active_3)*3*2//(p*(p-1))}")

# Resolvability: can we partition the blocks into parallel classes?
# A parallel class = set of blocks covering every vertex exactly once.
# With p=11 and k=3: 11/3 is not integer, so no EXACT parallel class.
# But 3*3=9 < 11, leaving 2 vertices uncovered.
# A NEAR-parallel class covers 9 of 11 vertices.

# How many (3,3,3) near-parallel classes exist?
npc = 0
for i in range(len(active_3)):
    for j in range(i+1, len(active_3)):
        if active_3[i] & active_3[j]:
            continue
        for l in range(j+1, len(active_3)):
            if (active_3[i] & active_3[l]) or (active_3[j] & active_3[l]):
                continue
            npc += 1

print(f"  Near-parallel classes (3,3,3): {npc}")
print(f"  Per leftover pair: {npc / 55}")

# Cross-orientation: how does the number of near-parallel classes vary?
print(f"\n\n{'='*70}")
print(f"NEAR-PARALLEL CLASSES ACROSS ORIENTATIONS")
print(f"{'='*70}")

for bits in range(N):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)
    active_3 = []
    for subset in combinations(range(p), 3):
        if count_ham_cycles(A, list(subset)) > 0:
            active_3.append(frozenset(subset))

    c3 = len(active_3)

    # Count (3,3,3) near-parallel classes
    npc = 0
    for i in range(len(active_3)):
        for j in range(i+1, len(active_3)):
            if active_3[i] & active_3[j]:
                continue
            for l in range(j+1, len(active_3)):
                if (active_3[i] & active_3[l]) or (active_3[j] & active_3[l]):
                    continue
                npc += 1

    # Count disjoint pairs
    dp_count = sum(1 for i in range(len(active_3))
                   for j in range(i+1, len(active_3))
                   if not (active_3[i] & active_3[j]))

    is_paley = sorted(S) == S_qr
    tag = " <-- PALEY" if is_paley else ""
    print(f"  bits={bits}: c3={c3}, disj_pairs={dp_count}, npc_333={npc}{tag}")

print("\nDONE.")
