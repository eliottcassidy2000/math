#!/usr/bin/env python3
"""
klein-2026-07-05-S132 (HYP-4095) - THE ONE-WINDOW PEEL LEMMA: verify + map the corner.

LEMMA: V = W u {v*}, all speeds positive integers. Suppose t* has
  min_{w in W} dist(w t*) >= beta > 1/14.
Window J = [t* - d, t* + d], d = (beta - 1/14)/maxW: every w in W stays >= 1/14 on J
(Lipschitz, slope <= w <= maxW). The killer's 1/14-bad set is a union of arcs of length
exactly (1/7)/v* separated by good runs; an interval longer than one bad arc contains a
1/14-good point. So if
    2*(beta - 1/14)/maxW > 1/(7 v*)    <=>    v* > maxW / (14*(beta - 1/14))
then some t in J is 1/14-good for ALL of V => Lonely 14.

Thresholds: beta=1/13 -> exactly 13*maxW (the dominant/compressed line!);
            beta=2/25 -> (25/3)*maxW; beta>=1/7 -> threshold <= maxW (any killer).

This script: (a) numeric sanity of the lemma on random/structured families;
(b) MAP THE CORNER: compressed covering primitive families where the largest-peel base
is loose (M(base)>2/25 say) yet the killer is BELOW threshold -- what are their M values?
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def dist(x):  # distance to nearest integer, x Fraction
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def exact_M(V, denbound=None):
    """max over breakpoint-grid candidates of min_i dist(v_i t) (exact; grid exhaustive
    for the optimum: crossings a/(u+v), a/(v-u), peaks a/(2u))."""
    qs = set()
    Vl = list(V)
    for i, u in enumerate(Vl):
        qs.add(2*u)
        for v in Vl[i+1:]:
            qs.add(u+v); qs.add(v-u)
    qs.discard(0)
    best = F(0); bt = None
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            t = F(a, q)
            m = min(dist(v*t) for v in Vl)
            if m > best: best, bt = m, t
    return best, bt

def covering(V):
    return all(any(v % q == 0 for v in V) for q in range(2, 15))

# (a) lemma sanity: for families where threshold holds, confirm Lonely (M(V) >= 1/14).
random.seed(5)
viol = 0; checked = 0
for _ in range(400):
    W = sorted(random.sample(range(1, 60), 12))
    MW, tstar = exact_M(W)
    if MW <= F(1,14): continue
    maxW = max(W)
    thresh = F(maxW, 1) / (14 * (MW - F(1,14)))
    vstar = int(thresh) + random.randint(1, 50)
    V = W + [vstar]
    checked += 1
    # lemma says Lonely; verify directly at the constructed witness:
    d = (MW - F(1,14)) / maxW
    # find good point of v* in [t*-d, t*+d]
    lo, hi = tstar - d, tstar + d
    s1 = vstar * lo
    # if s1 good take lo else m = round(s1), witness s = m + 1/14
    if dist(s1) >= F(1,14):
        t = lo
    else:
        m = int(s1 + F(1,2))
        t = F(m + F(1,14), 1) / vstar
        t = (F(m) + F(1,14)) / vstar
    ok = (lo <= t <= hi) and all(dist(w*t) >= F(1,14) for w in W) and dist(vstar*t) >= F(1,14)
    if not ok:
        viol += 1
        print("VIOLATION", W, vstar, float(MW))
print(f"(a) lemma constructive check: {checked} cases, {viol} violations")

# (b) THE CORNER: compressed covering primitive, loose base, killer below threshold.
print("\n(b) corner map: compressed covering primitive, base loose, killer sub-threshold")
corner = []
seeds = []
# structured seeds: base near-AP-but-loose x killer in (maxW, thresh)
for base in [list(range(1,12))+[24], list(range(1,12))+[13], [1,2,3,4,5,6,7,8,9,10,11,14],
             [1,2,3,4,5,6,7,8,9,10,12,13], [2,3,4,5,6,7,8,9,10,11,12,13]]:
    MW, _ = exact_M(base)
    maxW = max(base)
    if MW <= F(1,14): continue
    thresh = F(maxW,1)/(14*(MW - F(1,14)))
    for vstar in range(maxW+1, min(13*maxW, int(thresh))+1):
        V = base + [vstar]
        if not covering(V): continue
        if gcd(*V) != 1: continue
        MV, _ = exact_M(V)
        seeds.append((float(MV), V, float(MW), float(thresh)))
        if MV < F(1,13):
            corner.append((MV, V))
seeds.sort()
print(f"corner families checked: {len(seeds)}; below 1/13: {len(corner)}")
for mv, V in corner[:10]:
    print("  BELOW 1/13:", V, mv, float(mv))
if seeds:
    print("lowest five corner M values:")
    for mv, V, mw, th in seeds[:5]:
        print(f"  M={mv:.6f} V={V} (base M={mw:.5f}, thresh={th:.1f})")
