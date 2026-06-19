#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of the GAP dimension penalty claim.
Written from scratch with fractions.Fraction. Does NOT import the engine.

Claim under test (the angle): "every proper d>=2 generalized AP E with |E|=k
(k=8,9,10) has L_y(E) < cap_k, exploitable by a crude bound."

Adversarial finding to verify: the crude-bound premise is FALSE because the
binding worst non-AP set is a MINIMAL single-defect near-AP with TINY margin
(k=9 margin 0.00697), not a spread GAP. We re-derive everything exactly.
"""
from fractions import Fraction as F
import itertools
from math import gcd
from functools import reduce

# ---- exact distribution of N(x) ----
def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)  # frac part in [0,1)
        s = (v.numerator * 7) // v.denominator   # floor(7v)
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        p[N_at(E, mid)] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            g.append(F((t - 1) * (t - 2) * (t - 4) * (t - 5), 40))
        elif k in (9, 10):
            g.append(F(-(t - 2) * (t - 3) * (t - 6), 36))
        else:
            g.append(F((t - 3) * (t - 4), 12))
    return g

def L_y(E, k):
    p = dist_p(E)
    g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7)), p

def primitive(E):
    E = sorted(set(E))
    return reduce(gcd, E) == 1

# caps (from MEMORY, exact comparison done as float for reporting only)
CAP = {8: 0.38153, 9: 0.49426, 10: 0.6044}

def consec(k):
    return list(range(k))

print("=== Independent verification: GAP dimension penalty / single-defect maximizer ===\n")

# 1. Reproduce the claimed worst-set values exactly
claimed = {
    8:  ([0,2,3,4,5,6,7,8], F(30306,100000)),   # value reported ~0.30306
    9:  ([0,1,2,3,4,5,6,7,9], None),
    10: ([0,1,2,3,4,5,6,7,8,10], None),
}
for k in [8, 9, 10]:
    E, _ = claimed[k]
    Lc, pc = L_y(consec(k), k)
    Le, pe = L_y(E, k)
    print(f"k={k}: consec L_y = {Lc} = {float(Lc):.6f}")
    print(f"      worst-cand {E} L_y = {Le} = {float(Le):.6f}")
    print(f"      margin to cap_{k}({CAP[k]}) = {CAP[k]-float(Le):.6f}   beats consec? {Le>Lc}")
    print(f"      p(worst)={[str(x) for x in pe]}")
    print()
