#!/usr/bin/env python3
"""
Prove: for GS tilings, t3(T) + t3(flip(T)) is always odd.

Decompose into triple types:
- Type 1 (consecutive): {i, i+1, i+2} — contributes exactly 1 to sum (proved)
- Type 2+ (non-consecutive): need to show these contribute even total for GS

Also check: the diff distribution has a pattern. At n=5, diffs are odd in
{-3,-1,1,3} with counts {2,6,6,2} = 2*{1,3,3,1} = 2*C(3,k).
At n=7, diffs in {-5,-3,-1,1,3,5} with counts {16,80,160,160,80,16}
  = 16*{1,5,10,10,5,1} = 16*C(5,k).
At n=9, diffs in {-7,-5,-3,-1,1,3,5,7} with counts {512,...,512}
  = 512*{1,7,21,35,35,21,7,1} = 512*C(7,k).

Pattern: diffs are odd, symmetric, with counts = 2^(n-3) * C(n-2, k)!
This is a BINOMIAL distribution on n-2 consecutive triples!

Each consecutive triple independently contributes +1 or -1 to the diff,
and the non-consecutive triples have NO net effect on t3 under GS flip!

kind-pasteur-2026-03-06-S25h
"""
from itertools import combinations
from collections import Counter
from math import comb

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

for n in [5, 7, 9]:
    m = num_tiling_bits(n)
    pairs, fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, pairs, fixed)

    print(f"\nn={n}: {len(gs_tilings)} GS tilings")

    # For each GS tiling, decompose t3_diff into consecutive and non-consecutive
    consec_diffs = []
    nonconsec_diffs = []

    for bits in gs_tilings:
        A = tournament_from_tiling(n, bits)
        Af = tournament_from_tiling(n, flip_tiling(bits, m))

        # Consecutive triple contribution
        consec_before = 0
        consec_after = 0
        for i in range(n-2):
            a, b, c = i, i+1, i+2
            for x, y, z in [(a,b,c), (a,c,b)]:
                if A[x][y] and A[y][z] and A[z][x]:
                    consec_before += 1
                if Af[x][y] and Af[y][z] and Af[z][x]:
                    consec_after += 1

        # Total
        total_before = sum(1 for i,j,k in combinations(range(n),3)
                          for x,y,z in [(i,j,k),(i,k,j)]
                          if A[x][y] and A[y][z] and A[z][x])
        total_after = sum(1 for i,j,k in combinations(range(n),3)
                         for x,y,z in [(i,j,k),(i,k,j)]
                         if Af[x][y] and Af[y][z] and Af[z][x])

        nonconsec_before = total_before - consec_before
        nonconsec_after = total_after - consec_after

        consec_diffs.append(consec_after - consec_before)
        nonconsec_diffs.append(nonconsec_after - nonconsec_before)

    print(f"  Consecutive diff distribution: {dict(sorted(Counter(consec_diffs).items()))}")
    print(f"  Non-consecutive diff distribution: {dict(sorted(Counter(nonconsec_diffs).items()))}")
    print(f"  Non-consec sum parities: {set(d % 2 for d in nonconsec_diffs)}")

    # Check: are consecutive diffs always ODD?
    print(f"  Consecutive diffs always odd: {all(d % 2 == 1 for d in consec_diffs)}")

    # Check binomial pattern
    total_diff_counts = Counter(d for d in (cd + nd for cd, nd in zip(consec_diffs, nonconsec_diffs)))
    print(f"  Total diff: {dict(sorted(total_diff_counts.items()))}")

    # Expected: 2^(n-3) * C(n-2, k) for diff = n-2 - 2k
    print(f"  Expected (binomial):")
    for k in range(n-1):
        diff_val = (n-2) - 2*k
        expected = 2**(len(pairs) + len(fixed) - (n-2)) * comb(n-2, k)
        actual = total_diff_counts.get(diff_val, 0)
        print(f"    diff={diff_val:+d}: expected={expected}, actual={actual}, match={'YES' if expected==actual else 'NO'}")
