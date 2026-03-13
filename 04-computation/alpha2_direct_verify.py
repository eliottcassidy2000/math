#!/usr/bin/env python3
"""
alpha2_direct_verify.py -- kind-pasteur-2026-03-13-S60

Direct verification of alpha_2 (disjoint cycle pairs) at p=7.
The backtracking gives alpha_2=13 for Interval, but overlap analysis says 7 disjoint pairs.
H from Held-Karp = 175, OCF = 171. Discrepancy = 4 = 4*1, suggesting alpha_2 off by 1.
"""

import cmath
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


def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


p = 7
for name, S in [("Paley", [1, 2, 4]), ("Interval", [1, 2, 3])]:
    A = build_adj(p, S)
    H_hk = compute_H_heldkarp(A, p)

    print(f"\n{'='*70}")
    print(f"{name}, S={S}, H(Held-Karp)={H_hk}")
    print(f"{'='*70}")

    # Enumerate ALL directed cycles with full detail
    cycles = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            if n_cyc > 0:
                print(f"  k={k}, verts={verts}: {n_cyc} directed cycles")
            for _ in range(n_cyc):
                cycles.append(frozenset(subset))

    n = len(cycles)
    print(f"\n  Total directed cycles: {n}")

    # Count by length
    by_len = defaultdict(int)
    for c in cycles:
        by_len[len(c)] += 1
    for k in sorted(by_len):
        print(f"    k={k}: {by_len[k]} directed cycles")

    # Direct alpha_2 computation: count ALL disjoint pairs
    alpha_2_direct = 0
    disjoint_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if not (cycles[i] & cycles[j]):
                alpha_2_direct += 1
                disjoint_pairs.append((i, j, len(cycles[i]), len(cycles[j])))

    print(f"\n  Direct alpha_2 = {alpha_2_direct}")
    print(f"  Disjoint pairs by type:")
    type_counts = defaultdict(int)
    for i, j, ki, kj in disjoint_pairs:
        type_counts[(ki, kj)] += 1
    for (ki, kj), cnt in sorted(type_counts.items()):
        print(f"    ({ki},{kj})-disjoint pairs: {cnt}")

    # Direct alpha_3: count ALL mutually disjoint triples
    alpha_3_direct = 0
    for i in range(n):
        for j in range(i + 1, n):
            if cycles[i] & cycles[j]:
                continue
            for k in range(j + 1, n):
                if not (cycles[i] & cycles[k]) and not (cycles[j] & cycles[k]):
                    alpha_3_direct += 1

    print(f"  Direct alpha_3 = {alpha_3_direct}")

    # Verify H = 1 + 2*n + 4*alpha_2 + 8*alpha_3
    H_ocf = 1 + 2*n + 4*alpha_2_direct + 8*alpha_3_direct
    print(f"\n  H(OCF) = 1 + 2*{n} + 4*{alpha_2_direct} + 8*{alpha_3_direct} = {H_ocf}")
    print(f"  H(Held-Karp) = {H_hk}")
    print(f"  Match: {H_ocf == H_hk}")

    if H_ocf != H_hk:
        print(f"  MISMATCH! Difference = {H_hk - H_ocf}")

    # Also check backtracking alpha
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

    if n <= 25:
        backtrack(-1, 0, 0)
        max_j = max((j for j in range(len(alpha_bt)) if alpha_bt[j] > 0), default=0)
        print(f"\n  Backtrack alpha: {[alpha_bt[j] for j in range(max_j+1)]}")
        H_bt = sum(alpha_bt[j] * (2**j) for j in range(len(alpha_bt)))
        print(f"  H(backtrack) = {H_bt}")
        print(f"  alpha_2(direct) = {alpha_2_direct}, alpha_2(bt) = {alpha_bt[2]}")
        if alpha_bt[2] != alpha_2_direct:
            print(f"  *** BACKTRACKING BUG: alpha_2 differs! ***")

print("\nDONE.")
