#!/usr/bin/env python3
"""
Does t3 parity ALWAYS flip under tiling flip for GS tilings specifically?

If so, this would explain the bipartiteness of the skeleton.

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import Counter

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

for n in [3, 5, 7, 9]:
    m = num_tiling_bits(n)
    pairs, fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, pairs, fixed)

    print(f"n={n}: {len(gs_tilings)} GS tilings, m={m}")

    parity_flips = 0
    parity_same = 0
    diffs = []

    for bits in gs_tilings:
        A = tournament_from_tiling(n, bits)
        t3_before = count_3cycles(A, n)

        flipped = flip_tiling(bits, m)
        A_f = tournament_from_tiling(n, flipped)
        t3_after = count_3cycles(A_f, n)

        if t3_before % 2 != t3_after % 2:
            parity_flips += 1
        else:
            parity_same += 1
        diffs.append(t3_after - t3_before)

    print(f"  Parity flips: {parity_flips}, Parity same: {parity_same}")
    print(f"  ALWAYS flips: {parity_same == 0}")
    print(f"  Diff distribution: {dict(sorted(Counter(diffs).items()))}")

    # Also check: (t3_before + t3_after) mod 2 for GS tilings
    sum_mod2 = Counter((count_3cycles(tournament_from_tiling(n, b), n) +
                         count_3cycles(tournament_from_tiling(n, flip_tiling(b, m)), n)) % 2
                        for b in gs_tilings)
    print(f"  (t3_before + t3_after) mod 2: {dict(sum_mod2)}")
    print()
