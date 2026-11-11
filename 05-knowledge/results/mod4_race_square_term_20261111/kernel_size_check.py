#!/usr/bin/env python3
"""Expected kernel size of the adjacency matrix of a uniformly random labelled tournament over F_2:
E|ker A| = 1 + 2^(1-m) (N01(m) - 1) = 2 + 2^((1-m)/2) cos((m-1) pi/4) - 2^(1-m), N01(m) = #{U : |U| = 0, 1 (mod 4)}.
Exhaustive check over all labelled tournaments of order m <= MMAX (default 6)."""
import sys
from fractions import Fraction as Fr
from math import comb, cos, pi
def rank2(rows, m):
    rows = rows[:]; r = 0
    for c in range(m):
        p = None
        for i in range(r, m):
            if (rows[i] >> c) & 1: p = i; break
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
        r += 1
    return r
mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for m in range(1, mmax + 1):
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]; tot = 0
    for bits in range(1 << len(pairs)):
        rows = [0] * m
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
        tot += 2 ** (m - rank2(rows, m))
    E = Fr(tot, 2 ** len(pairs)); N01 = sum(comb(m, k) for k in range(m + 1) if k % 4 in (0, 1))
    pred = 1 + Fr(N01 - 1, 2 ** (m - 1)); trig = 2 + 2 ** ((1 - m) / 2) * cos((m - 1) * pi / 4) - 2.0 ** (1 - m)
    print("m = %d: E|ker A| = %s; formula %s (%s); trigonometric form %.10f" % (m, E, pred, "equal" if E == pred else "DIFFERENT", trig))
