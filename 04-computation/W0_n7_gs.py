#!/usr/bin/env python3
"""
W(0) for GS tilings at n=7 — does the vanishing pattern generalize?

At n=5: W(0) = 1 - t3 + 2*t5, vanishes for ALL odd-t3 GS tilings.
At n=7: W(0) = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc (fractional!)

Question: Is W(0) = 0 mod something for one side of the bipartition?
Or is there a different pattern?

Also check: W(0) as function of t3 parity.

Optimized: compute W(0) directly as product of s_e values.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction

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

def compute_W0_fast(A, n):
    """W(0) = sum_perm prod(A[p_i,p_{i+1}] - 1/2)."""
    half = Fraction(1, 2)
    total = Fraction(0)
    for perm in permutations(range(n)):
        prod = Fraction(1)
        for i in range(n-1):
            prod *= (Fraction(A[perm[i]][perm[i+1]]) - half)
        total += prod
    return total

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs, fixed, seen = [], [], set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen: continue
        ti, tj = n-1-j, n-1-i
        if ti > tj: ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx); seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx)); seen.add(idx); seen.add(tidx)
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

# n=7: 512 GS tilings, each requires 7! = 5040 permutations
# Total: 2.5M operations — should be feasible but slow
# Let's do a subset first

n = 7
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)
gs_dof = len(pairs) + len(fixed)

print(f"W(0) for GS tilings at n={n}")
print(f"  {len(gs_tilings)} GS tilings, m={m}, gs_dof={gs_dof}")
print(f"  Computing W(0) for first 32 GS tilings...")

# Check first 32 pairs
results = []
for idx in range(min(64, len(gs_tilings))):
    bits = gs_tilings[idx]
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    W0 = compute_W0_fast(A, n)
    results.append((idx, bits, t3, W0))

    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)
    if idx <= flip_idx and idx < 32:
        Af = tournament_from_tiling(n, flipped)
        t3f = count_3cycles(Af, n)
        W0f = compute_W0_fast(Af, n)
        print(f"  pair {idx}<->{flip_idx}: t3={t3}->{t3f}, W(0)={float(W0):.4f}->{float(W0f):.4f}, sum={float(W0+W0f):.4f}, prod={float(W0*W0f):.4f}")

# Analyze by t3 parity
print(f"\n  By t3 parity:")
odd_w0 = [W0 for _, _, t3, W0 in results if t3 % 2 == 1]
even_w0 = [W0 for _, _, t3, W0 in results if t3 % 2 == 0]
print(f"    Odd t3 ({len(odd_w0)}): W(0) values = {sorted(set(float(x) for x in odd_w0))}")
print(f"    Even t3 ({len(even_w0)}): W(0) values = {sorted(set(float(x) for x in even_w0))}")

# Check: does W(0) have a common denominator?
print(f"\n  All W(0) denominators: {sorted(set(W0.denominator for _, _, _, W0 in results))}")

# Check: 4*W(0) pattern
print(f"\n  4*W(0) values:")
for _, _, t3, W0 in results[:32]:
    print(f"    t3={t3}, 4*W(0)={4*W0} (t3%2={t3%2})")
