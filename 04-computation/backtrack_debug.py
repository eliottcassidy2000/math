#!/usr/bin/env python3
"""
backtrack_debug.py -- kind-pasteur-2026-03-13-S60

Debug: backtracking gives alpha_2=13 for Interval p=7 but direct counting gives 14.
Both have 59 cycles. Run both methods side by side with detailed comparison.
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


p = 7
S = [1, 2, 3]  # Interval
A = build_adj(p, S)

# Enumerate cycles
cycles = []
for k in range(3, p + 1, 2):
    for subset in combinations(range(p), k):
        verts = list(subset)
        n_cyc = count_ham_cycles(A, verts)
        for _ in range(n_cyc):
            cycles.append(frozenset(subset))

n = len(cycles)
print(f"Interval p=7: {n} directed cycles")
print(f"  k=3: {sum(1 for c in cycles if len(c)==3)}")
print(f"  k=5: {sum(1 for c in cycles if len(c)==5)}")
print(f"  k=7: {sum(1 for c in cycles if len(c)==7)}")

# Method 1: Direct counting of disjoint pairs
alpha_2_direct = 0
disjoint_pairs = []
for i in range(n):
    for j in range(i + 1, n):
        if not (cycles[i] & cycles[j]):
            alpha_2_direct += 1
            disjoint_pairs.append((i, j))

print(f"\nDirect alpha_2 = {alpha_2_direct}")
print(f"Disjoint pairs:")
for i, j in disjoint_pairs:
    print(f"  cycles[{i}]={set(cycles[i])}, cycles[{j}]={set(cycles[j])}")

# Method 2: Backtracking (from alpha3_p7_only.py)
nbr = [0] * n
for i in range(n):
    for j in range(i + 1, n):
        if cycles[i] & cycles[j]:
            nbr[i] |= (1 << j)
            nbr[j] |= (1 << i)

alpha_bt = [0] * (n + 1)

def backtrack(v, mask, size):
    alpha_bt[size] += 1
    for w in range(v + 1, n):
        if not (mask & (1 << w)):
            backtrack(w, mask | nbr[w], size + 1)

import sys
sys.setrecursionlimit(10000)
backtrack(-1, 0, 0)

max_j = max(j for j in range(len(alpha_bt)) if alpha_bt[j] > 0)
print(f"\nBacktrack alpha: {[alpha_bt[j] for j in range(max_j+1)]}")
print(f"Backtrack alpha_2 = {alpha_bt[2]}")

H_bt = sum(alpha_bt[j] * (2**j) for j in range(len(alpha_bt)))
print(f"H(backtrack) = {H_bt}")

# Verify: which pairs does backtracking miss?
if alpha_bt[2] != alpha_2_direct:
    print(f"\n*** MISMATCH: direct={alpha_2_direct}, bt={alpha_bt[2]} ***")
    print(f"Difference: {alpha_2_direct - alpha_bt[2]}")

    # Try to find which pairs are missed by examining the conflict graph
    # Check: are any pairs (i,j) with both i,j having ALL neighbors covered by some
    # other cycle in the independent set order?
    for i, j in disjoint_pairs:
        # Check if pair (i,j) can be found by the backtracking
        # The backtracking finds pair (i,j) only if i < j and neither is blocked
        # at the time they're considered
        print(f"  Pair ({i},{j}): verts {set(cycles[i])} U {set(cycles[j])}")
        # Check how many cycles conflict with i
        conflicts_i = bin(nbr[i]).count('1')
        conflicts_j = bin(nbr[j]).count('1')
        print(f"    conflicts_i={conflicts_i}, conflicts_j={conflicts_j}")
else:
    print(f"\nBacktrack matches direct: alpha_2 = {alpha_2_direct}")

# Also check: number of edges in conflict graph
n_edges = sum(bin(nbr[i]).count('1') for i in range(n)) // 2
print(f"\nConflict graph: {n} vertices, {n_edges} edges")
print(f"Independent number alpha_2(direct) = {alpha_2_direct}")
