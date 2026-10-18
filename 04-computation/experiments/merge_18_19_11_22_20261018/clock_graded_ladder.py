#!/usr/bin/env python3
"""Clock-graded ideal ladder (mac-mini-2026-10-18-sixcolour): THM-4667's synchronous cross terms E_d - E_(d-1), resolved by
literal clock lag, in the ideal N -> infinity band model of clock_duality's R7 / 3c.

F(e, r, b): signed clock profile {clock: value} of all backward words of length e from a node of class r mod 3^(e+1)
whose scale x = t/N lies in band b, i.e. M_e/2^(b+1) < x <= M_e/2^b with M_e = 4 (3/2)^e. A child s_j of such a node
(a_j = 2 if B^j(t) = 1 mod 3, a_j = 1 if = 2 mod 3, none if 0) has clock a_j + 2j, sign (-1)^(a_j), class
(2^(a_j) B^j(t) - 1)/3 mod 3^e, and lies exactly in band b + 1 - a_j - 2j of M_(e-1). Bands b < 0 carry nothing (all
sources beyond 4N). F(0, r, b) = h(r, y) for y in (2^(1-b), 2^(2-b)].
Node weight per class mod 3^(d+1) and band b: (3/8)/3^(d+1) * M_d/2^(b+1).
Cross term at depth d = 2 * sum over nodes of sum_(j<j') <V_j, V_j'> split by |clock - clock'|.
Checks against the exact ideal energies E_1 = 32898/9709, E_2 = 572764898734162/222399981598543 (S41 audit a4).
Usage: python3 clock_graded_ladder.py DMAX BMAX"""
import sys
from functools import lru_cache
from collections import defaultdict
from fractions import Fraction as Fr
kappa = (0, -3, 3); lam = (-1, -1, 2)
def h(r, y):
    tot, k = 0, 0
    while 4 ** k * y <= 2:
        tot += kappa[(r + k) % 3]; k += 1
    if 4 ** k * y <= 4:
        tot += lam[(r + k) % 3]
    return tot
DMAX, BMAX = int(sys.argv[1]), int(sys.argv[2])

def children(e, r, b):
    """(sign, clock, child class mod 3^e, child band) for a node of class r mod 3^(e+1) in band b."""
    mod = 3 ** (e + 1); out = []; rho = r % mod; j = 0
    while True:
        cb = b + 1 - 2 * j - 1            # largest possible child band index is b + 1 - a - 2j with a >= 1
        if cb < 0: break
        if rho % 3 == 1: a = 2
        elif rho % 3 == 2: a = 1
        else: a = None
        if a is not None:
            band = b + 1 - a - 2 * j
            if band >= 0:
                out.append(((-1) ** a, a + 2 * j, ((2 ** a * rho - 1) % mod) // 3, band))
        rho = (4 * rho + 1) % mod; j += 1
    return out

@lru_cache(maxsize=None)
def F(e, r, b):
    if b < 0: return ()
    if e == 0:
        v = h(r % 3, 3.0 * 2.0 ** (-b))
        return ((0, v),) if v else ()
    acc = defaultdict(int)
    for sg, cl, rc, bc in children(e, r, b):
        for c, v in F(e - 1, rc, bc):
            acc[c + cl] += sg * v
    return tuple(sorted((c, v) for c, v in acc.items() if v))

for d in range(1, DMAX + 1):
    lag = defaultdict(float); Ed = 0.0
    M = 4 * 1.5 ** d
    for b in range(BMAX):
        w = (3 / 8) / 3 ** (d + 1) * M / 2 ** (b + 1)
        for r in range(3 ** (d + 1)):
            kids = [(cl, sg, F(d - 1, rc, bc)) for sg, cl, rc, bc in children(d, r, b)]
            prof = [{c + cl: sg * v for c, v in Fk} for cl, sg, Fk in kids]
            tot = sum(sum(p.values()) for p in prof); Ed += w * tot * tot
            for i in range(len(prof)):
                for k in range(i + 1, len(prof)):
                    for c1, v1 in prof[i].items():
                        for c2, v2 in prof[k].items():
                            lag[abs(c1 - c2)] += w * v1 * v2
    cross = 2 * sum(lag.values())
    print('d = %d: E_d = %.12f ; cross term E_d - E_(d-1) = %+.12f ; lag 0: %+.8f lag 1: %+.8f lag 2: %+.8f lag 3: %+.8f (x2 = THM-4667 units) ; F-cache %d'
          % (d, Ed, cross, 2 * lag[0], 2 * lag[1], 2 * lag[2], 2 * lag[3], F.cache_info().currsize), flush=True)
    if d == 1: print('   exact E_1 = 32898/9709 = %.12f' % (32898 / 9709))
    if d == 2: print('   exact E_2 = %.12f' % float(Fr(572764898734162, 222399981598543)))
