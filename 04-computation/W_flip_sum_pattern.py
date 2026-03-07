#!/usr/bin/env python3
"""
W(0) + W_flip(0) pattern — what determines the sum?

At n=5: W(0)+Wf(0) took values {-1, 1, 3}.
  W(0) = 1 - t3 + 2*t5
  W(0)+Wf(0) = 2 - (t3+t3f) + 2*(t5+t5f)

For GS tilings at n=5, t3+t3f has only a few values.
THM-060 tells us t3+t3f is always odd.

Actually the key formula is:
  W_sum(0) = 2*F_{n-1}(0) + 2*F_{n-3}(0)*(t3+t3f) + 2*F_{n-5}(0)*(t5+t5f) + ...

And from the t3_gs_parity analysis:
  t3+t3f follows distribution: C(n-2, k) * 2^{gs_dof-(n-2)}

Let's explore: what determines W(0)+Wf(0)? Is it the Hamming weight of the
GS tiling vector?

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

def compute_W0(A, n):
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

def hamming_weight(x):
    return bin(x).count('1')

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)
gs_dof = len(pairs) + len(fixed)

print(f"n={n}, m={m}, gs_dof={gs_dof}")
print(f"{'='*70}")

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)
    if idx > flip_idx:
        continue

    A = tournament_from_tiling(n, bits)
    Af = tournament_from_tiling(n, flipped)
    t3 = count_3cycles(A, n)
    t3f = count_3cycles(Af, n)
    W0 = compute_W0(A, n)
    W0f = compute_W0(Af, n)

    # GS free bits
    hw = hamming_weight(idx)  # Hamming weight of GS free vector
    hwf = hamming_weight(flip_idx)

    print(f"  gs_free={idx:04b} (HW={hw}) <-> {flip_idx:04b} (HW={hwf}): "
          f"t3={t3}+{t3f}={t3+t3f}, W(0)+Wf(0)={W0+W0f}")

# H(T) + H(flip) pattern
print(f"\n{'='*70}")
print(f"H(T) + H(flip) pattern at n={n}")
print(f"{'='*70}")

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)
    if idx > flip_idx:
        continue

    A = tournament_from_tiling(n, bits)
    Af = tournament_from_tiling(n, flipped)
    t3 = count_3cycles(A, n)
    t3f = count_3cycles(Af, n)

    # H = 1 + 2*t3 + 4*... = I(Omega,2) via OCF
    # At n=5: H = 1 + 2*t3 + 2*t5 (since bc terms only start at n=6)
    def count_5cycles(A, n):
        t5 = 0
        for combo in combinations(range(n), 5):
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    t5 += 1
        return t5 // 5

    t5 = count_5cycles(A, n)
    t5f = count_5cycles(Af, n)
    H = 1 + 2*t3 + 2*t5  # Actually H = I(Omega,2) which is exact
    Hf = 1 + 2*t3f + 2*t5f

    hw = hamming_weight(idx)
    print(f"  gs_free={idx:04b} (HW={hw}): t3={t3}+{t3f}={t3+t3f}, t5={t5}+{t5f}={t5+t5f}, "
          f"H={H}+{Hf}={H+Hf}")

print(f"\n{'='*70}")
print(f"DEEPER: t3+t3f vs GS Hamming weight at n={n}")
print(f"{'='*70}")

# t3 + t3f should depend on the GS free bits somehow
from collections import Counter
hw_to_t3sum = {}
for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped)
    if idx > flip_idx:
        continue

    A = tournament_from_tiling(n, bits)
    Af = tournament_from_tiling(n, flipped)
    t3 = count_3cycles(A, n)
    t3f = count_3cycles(Af, n)
    hw = hamming_weight(idx)
    print(f"  HW({idx:04b})={hw}: t3+t3f={t3+t3f}, dt3={t3-t3f}")
