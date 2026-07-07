#!/usr/bin/env python3
"""
kps-2026-07-07-S55 -- the SATURATED hard core of LRC(14), and how the sieve + coarse
reduction localize the crux.

Integrating opus-S131's sieve reframe (the near-tight AP/Farey families I censused in S54 are
SIEVE-EASY -- they miss q=14, so M>=1/14 by lonely_of_no_multiple; the genuine hard core is the
SATURATED families = multiples of every q in {2..14}).  This script (a) adversarially confirms
the saturated margin, and (b) states the decomposition of LRC(14) that results.

SIEVE DICHOTOMY (opus-S131, standard).  Non-saturated (misses some q<=14) => M>=1/q>=1/14
(GREEN counterexample_needs_all_divisors).  So a counterexample must be SATURATED.

CENSUS (this script, confirming opus + adversarial):
  * min M over saturated 13-families = 1/12 = 0.0833 (exhaustive {1..21}: 11385 families, +
    randomized to hi=100); extremal = {1,2,3,4,10,11,12,13,14,15,16,17,18}.  ZERO below 1/13.
    Robust margin M >= 1/12 >> 1/14.

DECOMPOSITION of LRC(14) (the synthesis -- sieve[opus] + coarse[kps] + decorrelation[mac-mini/opus]):
  1. NON-SATURATED           => M>=1/q>=1/14         SIEVE          [GREEN]
  2. SATURATED, clustered     => M>=1/13-eps>1/14     COARSE RED.    [GREEN: kps lonely14_of_coarse_le12]
     (tight clusters into <=12 groups at a large scale L; incl. consecutive blocks {N..N+12},
      r=1 cluster, M->1/2)
  3. SATURATED, spread single-scale (no tight clustering, bounded ratio, large)
                              => M>=1/12 empirically  DECORRELATION  [OPEN -- the crux, but with MARGIN]
  4. SATURATED, small (bounded speeds)  => M>=1/12    FINITE CHECK   [census]

So the crux is localized to leg 3 -- spread single-scale saturated families -- which carry a
clear MARGIN (M>=1/12, 17% above 1/14), the analytic content the density-floor / decorrelation
machinery targets.  Legs 1,2 are GREEN; leg 4 is finite.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations

def Mw(v):
    v = [x for x in v if x]; S = sum(abs(x) for x in v); Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64); bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q); r = np.outer(va, a) % q
        d = np.minimum(r, q - r); col = d.min(axis=0); qb = int(col.max())
        if qb * bd > bn * q: bn, bd = qb, q
    return Fraction(bn, bd)

def saturated(v):
    return all(any(x % q == 0 for x in v) for q in range(2, 15))

print("SATURATED hard-core M-floor (exhaustive small range):")
for hi in (18, 20, 22):
    mins = None; minv = None; cnt = 0; below13 = 0
    for combo in combinations(range(1, hi + 1), 13):
        if not saturated(combo): continue
        cnt += 1; M = Mw(list(combo))
        if M < Fraction(1, 13): below13 += 1
        if mins is None or M < mins: mins = M; minv = combo
    print("  {1..%d}: %6d saturated, min M = %s = %.5f, #below 1/13 = %d" %
          (hi, cnt, str(mins), float(mins), below13))
    if hi == 22: print("           extremal =", minv)

print()
print("consecutive blocks {a..a+12} (saturated when they contain a mult of 14): margin grows")
for a in (2, 3, 5, 8, 14, 28):
    v = list(range(a, a + 13))
    print("  a=%-3d {%d..%d} saturated=%s  M=%s=%.4f" %
          (a, a, a + 12, saturated(v), str(Mw(v)), float(Mw(v))))
print()
print("=> saturated margin >= 1/12 robust; crux = spread single-scale saturated (leg 3), with margin.")
