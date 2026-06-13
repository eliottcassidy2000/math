#!/usr/bin/env python3
"""
W(0) vanishing pattern under GS flip.

DISCOVERY: At n=5, every GS flip pair has W(0)*W_flip(0) = 0.
One member always has W(0) = 0. Which side?

W(0) = 1 - t3 + 2*t5 at n=5
W(0) = 0 when t3 is odd (t3 in {1,3,5,...}) and t5 has the right value.

Actually: W(0) = F_4(0) + 2*F_2(0)*t3 + 2*F_0(0)*t5
        = 1 - t3 + 2*t5

For W(0) = 0: t3 = 1 + 2*t5, i.e. t3 must be odd.
For W(0) != 0: t3 is even, OR t3 is odd but t5 doesn't satisfy.

Since t3 parity flips under GS flip (THM-060), exactly one member of each
pair has odd t3. And the data shows W(0)=0 whenever t3 is odd!

Hmm wait, that's not quite right. Let me check more carefully.

Also: check at n=7 (much richer).

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

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
    """W(0) = sum over perms of prod(s_i) where s_i = A[p_i,p_{i+1}] - 1/2."""
    total = Fraction(0)
    for perm in permutations(range(n)):
        prod = Fraction(1)
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            prod *= s
        total += prod
    return total

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_5cycles(A, n):
    t5 = 0
    for combo in combinations(range(n), 5):
        sub = list(combo)
        for perm in permutations(sub):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
    return t5 // 5

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

n = 5
print(f"W(0) VANISHING at n={n}")
print("="*70)
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Check ALL GS tilings
odd_t3_vanish = 0
even_t3_vanish = 0
odd_t3_total = 0
even_t3_total = 0

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    t5 = count_5cycles(A, n)
    W0 = compute_W0(A, n)

    if t3 % 2 == 1:
        odd_t3_total += 1
        if W0 == 0: odd_t3_vanish += 1
    else:
        even_t3_total += 1
        if W0 == 0: even_t3_vanish += 1

    print(f"  tiling {idx:04b}: t3={t3}, t5={t5}, W(0)={W0}, t3%2={t3%2}")

print(f"\n  Odd t3: {odd_t3_vanish}/{odd_t3_total} have W(0)=0")
print(f"  Even t3: {even_t3_vanish}/{even_t3_total} have W(0)=0")

# Now: at n=5, W(0) = 1 - t3 + 2*t5
# t3 odd => W(0) even (since 1-odd+2t5 = even)
# t3 even => W(0) odd (since 1-even+2t5 = odd)
# So W(0) is NEVER zero when t3 is even!
# And W(0) CAN be zero when t3 is odd!

print(f"\n  PARITY: W(0) mod 2 when t3 odd: always even")
print(f"  PARITY: W(0) mod 2 when t3 even: always odd (so NEVER 0!)")
print(f"  This means W(0)=0 can ONLY happen on the odd-t3 side of bipartition!")

# Check n=3
print(f"\n{'='*70}")
print(f"W(0) at n=3")
print(f"{'='*70}")
n = 3
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

for idx in range(len(gs_tilings)):
    bits = gs_tilings[idx]
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    W0 = compute_W0(A, n)
    print(f"  tiling {idx}: t3={t3}, W(0)={W0}")

# General formula for W(0) parity
print(f"\n{'='*70}")
print(f"GENERAL W(0) PARITY ANALYSIS")
print(f"{'='*70}")
print("  W(0) = F_{n-1}(0) + sum 2^p * F_f(0) * I(T)")
print("  F_{2k}(0) = (-1)^k * T_{k+1} / 4^k where T_k are tangent numbers")
print("  T_1=1, T_2=2, T_3=16, T_4=272, T_5=7936")
print()
print("  At n=5: F_4(0) = 1, F_2(0) = -1/2, F_0(0) = 1")
print("  W(0) = 1 + 2*(-1/2)*t3 + 2*1*t5 = 1 - t3 + 2*t5")
print("  W(0) mod 2 = 1 - t3 mod 2 = (1+t3) mod 2")
print("  If t3 even: W(0) is odd => NEVER zero")
print("  If t3 odd: W(0) is even => CAN be zero")
print()
print("  At n=7: F_6(0) = -17/4, F_4(0) = 1, F_2(0) = -1/2, F_0(0) = 1")
print("  W(0) = -17/4 + 2*1*t3 + 2*(-1/2)*t5 + 2*1*t7 + 4*(-1/2)*bc")
print("       = -17/4 + 2*t3 - t5 + 2*t7 - 2*bc")
print("  This is NOT integer in general! Let's check...")

n = 7
m = num_tiling_bits(n)
print(f"\n  At n=7 with a few tilings:")
for bits in [0, 1, (1 << m) - 1, 42]:
    A = tournament_from_tiling(n, bits)
    W0 = compute_W0(A, n)
    t3 = count_3cycles(A, n)
    print(f"    bits={bits}: t3={t3}, W(0)={W0} (float: {float(W0):.4f})")
