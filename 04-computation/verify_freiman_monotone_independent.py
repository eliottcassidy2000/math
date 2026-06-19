#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of the Freiman-partition / monotone replacement claim.
Own implementation of dist_p, L_y (exact Fraction). No imports from the project engine.

Claims under test:
(M) max L_y over PRIMITIVE NON-AP sets equals: k=8 -> 0.30306, k=9 -> 0.48729, k=10 -> 0.56241,
    each attained at minimal AP perturbation [0,2,3,4,5,6,7,8] / [0,1..7,9] / [0,1..8,10], sigma=2.0.
(C) caps: cap_8=0.38153, cap_9=0.49426, cap_10=0.6044, all > the non-AP max (margins listed).
(D) the binary (a)|(b) dichotomy: every 8-10 elt set sits in SOME 2-dim GAP; "stranger" type nearly vacuous.
(MONO) max L_y non-increasing in doubling sigma=|E+E|/k.
Hunt: any primitive non-AP set with L_y exceeding the stated non-AP max, or any set with L_y >= cap.
"""
import sys, itertools
from fractions import Fraction
from math import gcd
from functools import reduce

def frac_part(v: Fraction) -> Fraction:
    return v - (v.numerator // v.denominator)

def N_at(E, x: Fraction) -> int:
    hit = set()
    for e in E:
        v = frac_part(e * x)
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
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
            g.append(Fraction((t-1)*(t-2)*(t-4)*(t-5), 40))
        elif k in (9, 10):
            g.append(Fraction(-(t-2)*(t-3)*(t-6), 36))
        else:
            g.append(Fraction((t-3)*(t-4), 12))
    return g

def L_y(E, k):
    p = dist_p(E)
    g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7))

def is_primitive(E):
    d = reduce(gcd, [e for e in E if e != 0])
    return d == 1

def is_AP(E):
    s = sorted(set(E))
    if len(s) <= 2:
        return True
    d = s[1] - s[0]
    return all(s[i+1]-s[i] == d for i in range(len(s)-1))

def sigma_doubling(E):
    s = sorted(set(E))
    sums = {a+b for a in s for b in s}
    return Fraction(len(sums), len(s))

# Sanity: reproduce stated values for the claimed maximizers and caps reference
checks = {
    8: [0,2,3,4,5,6,7,8],
    9: [0,1,2,3,4,5,6,7,9],
    10:[0,1,2,3,4,5,6,7,8,10],
}
print("=== sanity: claimed non-AP maximizers ===")
for k, E in checks.items():
    L = L_y(E, k)
    print(f"k={k} E={E} primitive={is_primitive(E)} AP={is_AP(E)} sigma={float(sigma_doubling(E))} L_y={L} = {float(L):.6f}")

print("\n=== sanity: consec (AP) L_y per k ===")
for k in (8,9,10):
    C = list(range(k))
    print(f"k={k} consec L_y={float(L_y(C,k)):.6f}")
