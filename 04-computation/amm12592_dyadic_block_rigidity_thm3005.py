#!/usr/bin/env python3
"""Classification of composition-balanced blocks for AMM 12592.

A BLOCK [N, N+l) is a rule that decides every critical value n in [N, N+l)
by flip N+l.  Its terminal words are

    0^N y   (y in {0,1}^l, y != 0^l)      and     1^N y  (y != 1^l),

and the block is BALANCED if every fixed-composition class of that set
splits evenly between heads and tails -- the condition that makes the
Exactness Lemma apply block by block.

CLAIM.  The class-size generating function is
    (1+u^N)(1+u)^l - 1 - u^(N+l),
so the block is balanceable iff, over F_2,
    (1+u^N)(1+u)^l = 1 + u^(N+l).
CLASSIFICATION.  That holds iff N = 2^a and N+l = 2^g with g > a.
Hence every balanced block is a dyadic interval [2^a, 2^g), and in
particular l >= N: no balanced block has ratio (N+l)/N < 2.
"""
from itertools import product
from math import comb


def class_sizes(N, l):
    """Exact composition-class sizes of the block's terminal set."""
    sizes = {}
    for j in range(0, l + 1):                    # branch 0: weight j, y != 0
        if j == 0:
            continue
        sizes[j] = sizes.get(j, 0) + comb(l, j)
    for j in range(0, l + 1):                    # branch 1: weight N+j, y != 1^l
        if j == l:
            continue
        sizes[N + j] = sizes.get(N + j, 0) + comb(l, j)
    return sizes


def balanceable(N, l):
    return all(v % 2 == 0 for v in class_sizes(N, l).values())


def f2_criterion(N, l):
    """(1+u^N)(1+u)^l == 1 + u^(N+l) over F_2 ?"""
    lhs = [0] * (N + l + 1)
    for k in range(l + 1):
        c = comb(l, k) % 2
        if c:
            lhs[k] ^= 1
            lhs[k + N] ^= 1
    rhs = [0] * (N + l + 1)
    rhs[0] ^= 1
    rhs[N + l] ^= 1
    return lhs == rhs


def is_pow2(x):
    return x >= 1 and (x & (x - 1)) == 0


def classification(N, l):
    return is_pow2(N) and is_pow2(N + l)


def brute_balanced(N, l):
    """Direct check: does some heads/tails labelling bisect every class?
       (Possible iff every class has even size.)"""
    words = []
    for y in product((0, 1), repeat=l):
        if any(y):
            words.append((0,) * N + y)
        if not all(y):
            words.append((1,) * N + y)
    cnt = {}
    for w in words:
        cnt[sum(w)] = cnt.get(sum(w), 0) + 1
    return all(v % 2 == 0 for v in cnt.values())


if __name__ == "__main__":
    print("N   l   balanceable  F2-criterion  [2^a,2^g)  brute")
    bad = 0
    for N in range(1, 17):
        for l in range(1, 17):
            b = balanceable(N, l)
            f = f2_criterion(N, l)
            c = classification(N, l)
            br = brute_balanced(N, l) if l <= 12 else b
            if not (b == f == c == br):
                print(f"MISMATCH N={N} l={l}: {b} {f} {c} {br}")
                bad += 1
            if b:
                print(f"{N:2d} {l:3d}   {str(b):5s}        {str(f):5s}"
                      f"        {str(c):5s}     {str(br):5s}"
                      f"   ratio={(N+l)/N:g}")
    print(f"\nmismatches: {bad}")

    # wider sweep of the three criteria only (no enumeration)
    bad2 = 0
    for N in range(1, 300):
        for l in range(1, 300):
            if not (balanceable(N, l) == f2_criterion(N, l) == classification(N, l)):
                bad2 += 1
                if bad2 < 5:
                    print("WIDE MISMATCH", N, l)
    print(f"wide sweep N,l <= 299 mismatches: {bad2}")

    # the two corollaries
    viol = [(N, l) for N in range(1, 300) for l in range(1, 300)
            if balanceable(N, l) and l < N]
    print("balanced blocks with ratio < 2 (l < N):", viol)
    viol2 = [(N, l) for N in range(1, 300) for l in range(1, 300)
             if balanceable(N, l) and not is_pow2(N)]
    print("balanced blocks with non-dyadic start:", viol2)
