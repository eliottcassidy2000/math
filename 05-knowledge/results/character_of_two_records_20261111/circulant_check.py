#!/usr/bin/env python3
"""Theorem Q on circulant tournaments: for every connection set S of Z_n (n odd, S and -S partition the non-zero
residues) the tournament x -> x + s (s in S) has the anti-involution x -> -x fixing 0.  Compares H mod 4 (subset DP)
with 1 + 2 * (number of directed paths of the folded tournament ending in D_1).  Usage: circulant_check.py NMAX"""
import sys, itertools
from fold import fold, paths_ending_in, ham
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 13
for n in range(3, nmax + 1, 2):
    half = list(range(1, (n + 1) // 2)); bad = 0; cnt = 0; res = {1: 0, 3: 0}
    for signs in itertools.product((0, 1), repeat=len(half)):
        S = [h if s == 0 else n - h for h, s in zip(half, signs)]
        out = [0] * n
        for x in range(n):
            for s in S: out[x] |= 1 << ((x + s) % n)
        psi = [(-x) % n for x in range(n)]
        Z, D1, plus = fold(out, n, psi)
        H = ham(out, n); pred = (1 + 2 * paths_ending_in(Z, D1)) % 4
        cnt += 1; res[H % 4] += 1
        if H % 4 != pred: bad += 1
    print("n = %2d: %4d connection sets; H = 1 (mod 4) for %d, H = 3 (mod 4) for %d; disagreements with Theorem Q: %d" % (n, cnt, res[1], res[3], bad)); sys.stdout.flush()
