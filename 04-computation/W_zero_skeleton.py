#!/usr/bin/env python3
"""
W(0) under GS flip — connecting tangent numbers to skeleton bipartiteness.

FINDINGS SO FAR:
- F_f(-1/2) = (-1)^f  (!!!)
- W(-1/2) = W(1/2) = H(T) at odd n (W has even parity)
- W(0) = 1 - t3 + 2*t5 at n=5

Key insight: W(0) is the "tangent evaluation" of the tournament.
Via THM-059: W(0) = sum_I 2^{parts(I)} * F_f(0) * I(T)
where F_{2k}(0) = (-1)^k * T_{k+1} / 4^k (tangent numbers!)

So W(0) is a SIGNED weighted independence polynomial!

Question: What happens to W(0) under GS flip?
W(0) - W_flip(0) = sum_I 2^{parts(I)} * F_f(0) * (I(T) - I(flip(T)))

For the bipartite skeleton (THM-060), t3 parity flips.
So W(0) changes sign modulo something?

Let's check at n=5,7.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict
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

def compute_W_poly(A, n):
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(n)):
        poly = {0: Fraction(1)}
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            new_poly = {}
            for power, coeff in poly.items():
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly
        for power, coeff in poly.items():
            result[power] += coeff
    return dict(result)

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

def poly_eval(poly, r):
    return sum(c * r**p for p, c in poly.items())

for n in [3, 5]:
    print(f"\n{'='*70}")
    print(f"W(0) UNDER GS FLIP at n={n}")
    print(f"{'='*70}")

    m = num_tiling_bits(n)
    pairs, fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, pairs, fixed)

    print(f"  {len(gs_tilings)} GS tilings")

    w0_values = []
    for idx in range(len(gs_tilings)):
        bits = gs_tilings[idx]
        flipped = flip_tiling(bits, m)
        flip_idx = gs_tilings.index(flipped)

        A = tournament_from_tiling(n, bits)
        W = compute_W_poly(A, n)
        t3 = count_3cycles(A, n)
        W0 = poly_eval(W, Fraction(0))
        H = poly_eval(W, Fraction(1,2))

        Af = tournament_from_tiling(n, flipped)
        Wf = compute_W_poly(Af, n)
        t3f = count_3cycles(Af, n)
        W0f = poly_eval(Wf, Fraction(0))
        Hf = poly_eval(Wf, Fraction(1,2))

        if idx <= flip_idx:
            print(f"  pair {idx}<->{flip_idx}: t3={t3}->{t3f}")
            print(f"    W(0)={W0}, W_flip(0)={W0f}, sum={W0+W0f}, diff={W0-W0f}")
            print(f"    H={H}, H_flip={Hf}")
            w0_values.append((W0, W0f))

    # Check: is W(0) + W_flip(0) constant?
    sums = [a+b for a,b in w0_values]
    print(f"\n  W(0)+W_flip(0) values: {sums}")
    print(f"  All equal? {len(set(sums)) == 1}")
    if len(set(sums)) == 1:
        print(f"  Constant = {sums[0]}")

    # Check: W(0) * W_flip(0)?
    prods = [a*b for a,b in w0_values]
    print(f"  W(0)*W_flip(0) values: {prods}")

print()
print("="*70)
print("F_f(-1/2) PATTERN")
print("="*70)
print("  This is a key identity: F_f(-1/2) = (-1)^f")
print("  Proof: F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k")
print("  At r=-1/2: (r+1/2) = 0, (r-1/2) = -1")
print("  Only k=f survives: A(f+1,f) * 0^0 * (-1)^f = 1 * (-1)^f")
print("  Since A(n, n-1) = 1 (only the reverse permutation has n-1 descents)")
print("  PROVED: F_f(-1/2) = (-1)^f")

print()
print("="*70)
print("CONSEQUENCE FOR W(-1/2)")
print("="*70)
print("  W(-1/2) = F_{n-1}(-1/2) + sum 2^{parts} * F_f(-1/2) * I(T)")
print("  = (-1)^{n-1} + sum 2^{parts} * (-1)^f * I(T)")
print("  At odd n (n-1 even): (-1)^{n-1} = 1")
print("  f = (n-1) - 2*|pi|, so (-1)^f = (-1)^{n-1} * (-1)^{-2|pi|} = 1 * 1 = 1")
print("  So W(-1/2) = 1 + sum 2^{parts} * I(T) = I(Omega,2) = H(T)")
print("  This recovers W(-1/2) = W(1/2) = H(T) at odd n!")
print()
print("  At even n (n-1 odd): (-1)^{n-1} = -1")
print("  f = (n-1) - 2*|pi|, (-1)^f = (-1)^{n-1} * 1 = -1")
print("  So W(-1/2) = -1 - sum 2^{parts} * I(T) = -H(T)")
print("  W(-1/2) = -H(T) at even n! Let's verify.")

for n in [4, 6]:
    print(f"\n  n={n}: checking W(-1/2) = -H(T)")
    m = num_tiling_bits(n)
    # Just check a few random tilings
    for bits in [0, 1, (1 << m) - 1, 5]:
        A = tournament_from_tiling(n, bits)
        W = compute_W_poly(A, n)
        H = poly_eval(W, Fraction(1,2))
        anti = poly_eval(W, Fraction(-1,2))
        print(f"    bits={bits}: H={H}, W(-1/2)={anti}, -H={-H}, match={anti == -H}")
