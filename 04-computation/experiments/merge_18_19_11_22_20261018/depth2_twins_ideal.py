#!/usr/bin/env python3
"""Ideal (N -> infinity) value of the depth-2 equal-depth meet sums by clock lag (mac-mini-2026-10-18-sixcolour).
For a node t = r (mod 27) at scale x = t/N, its depth-2 preimage words (p, s) have clock tau(p) + tau(s), source scale
y = s/N and source residue s mod 3; the ideal source coefficient is h(s mod 3, y) (THM-4690(4)). Pairs first meeting
at t (different children p != p') contribute (-1)^(clock1 + clock2) h h' at lag |clock1 - clock2|. Bands: node scale
9/2^(b+1) < x <= 9/2^b, density 1/72 per class mod 27, so weight 1/(16 * 2^b); the lag-resolved sums are not periodic
in b, so they are summed directly to b < BMAX (float, then rational recognition).
Check: the total over lags must equal (E_2 - E_1)/2 of the exact ideal ladder (S41 audit a4)."""
import sys
from fractions import Fraction as Fr
from collections import defaultdict
kappa = (0, -3, 3); lam = (-1, -1, 2)
def h(r, y):
    tot, k = 0, 0
    while 4 ** k * y <= 2:
        tot += kappa[(r + k) % 3]; k += 1
    if 4 ** k * y <= 4:
        tot += lam[(r + k) % 3]
    return tot
def preimages(rmod, mod, x, limit):
    """(a, j, residue of s mod (mod/3), y) for t = rmod (mod mod) at t/N = x, preimages s_j = (2^a B^j(t) - 1)/3 with
    scale y <= limit. Intermediate nodes may lie above 4N (paths climb): a child at scale y1 has grandchildren at scale
    >= 2 y1/3, so children up to scale 6 matter."""
    out = []; rho = rmod % mod
    for j in range(0, 80):
        ys = 4.0 ** j * x / 3
        if 2 * ys > limit: break
        if rho % 3 == 1: a = 2
        elif rho % 3 == 2: a = 1
        else: a = None
        if a is not None:
            out.append((a, j, ((2 ** a * rho - 1) % mod) // 3, 2 ** a * ys))
        rho = (4 * rho + 1) % mod
    return out
BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 110
lag = defaultdict(float)
for b in range(BMAX):
    x = 0.75 * 9 / 2 ** b                              # interior point of the band 9/2^(b+1) < x <= 9/2^b
    w = 1.0 / (16 * 2 ** b)
    for r in range(27):
        words = []                                     # (child index, clock, signed value) of depth-2 words
        for ci, (a1, j1, rp, y1) in enumerate(preimages(r, 27, x, 6.0)):       # child p: residue mod 9, scale y1
            for (a2, j2, rs, y2) in preimages(rp, 9, y1, 4.0):                 # grandchild s: residue mod 3
                v = h(rs, y2)
                if v: words.append((ci, a1 + 2 * j1 + a2 + 2 * j2, (-1) ** (a1 + a2) * v))
        bycl = defaultdict(list)
        for ci, cl, v in words: bycl[cl].append((ci, v))
        cls = sorted(bycl)
        for i, c1 in enumerate(cls):
            for c2 in cls[i:]:
                D = c2 - c1
                s = 0.0
                for (ci, v) in bycl[c1]:
                    for (cj, u) in bycl[c2]:
                        if ci != cj: s += v * u
                if D == 0: s /= 2                          # unordered pairs within one clock
                if s: lag[D] += w * s
tot = sum(lag.values())
E1, E2 = Fr(32898, 9709), Fr(572764898734162, 222399981598543)
print('BMAX = %d; total depth-2 meet sum = %.12f ; (E_2 - E_1)/2 (exact ideal ladder) = %.12f' % (BMAX, tot, float((E2 - E1) / 2)))
for D in sorted(lag):
    if abs(lag[D]) > 1e-7:
        print('  lag %2d: %+.12f   %s' % (D, lag[D], Fr(lag[D]).limit_denominator(100000)))
