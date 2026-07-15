#!/usr/bin/env python3
"""ramanujan_fourier_profile_kps_S128c25.py -- kind-pasteur S128 cont.25.
THE RAMANUJAN-FOURIER EXPANSION OF THE INTERVAL-CORE GOOD SET (THM-873):
for lambda < 1/(k+1), with overlap corrections on closed Farey gaps,
   ghat(h) = -sum_{l<=k} c_l(h) sin(2 pi h lambda / l)/(pi h)  [+ closed-gap corrections],
c_l = Ramanujan sums. Referee vs exact interval transforms; disc_v vs THM-732 exact form."""
import sys, cmath
from fractions import Fraction as F
from math import comb, gcd, pi, sin
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
def ram(l, h):
    return sum(cmath.exp(2j*pi*h*c/l).real for c in range(l) if gcd(c, l) == 1)
def farey(k):
    return sorted(set(F(a, i) for i in range(1, k + 1) for a in range(i)))
def good_intervals(k, lam):
    pieces = []
    for w in range(1, k + 1):
        r = lam / w
        for a in range(w):
            c = F(a, w); lo, hi = c - r, c + r
            if lo < 0: pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1: pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else: pieces.append((lo, hi))
    pieces.sort(); out = []; cur = F(0)
    for lo, hi in pieces:
        if lo > cur: out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1: out.append((cur, F(1)))
    return out
def ghat_direct(k, lam, h):
    tot = 0j
    for a, b in good_intervals(k, lam):
        tot += (cmath.exp(-2j*pi*h*float(b)) - cmath.exp(-2j*pi*h*float(a))) / (-2j*pi*h)
    return tot
def ghat_formula(k, lam, h):
    tot = 0.0
    for l in range(1, k + 1):
        tot += ram(l, h) * sin(2*pi*h*float(lam)/l) / (pi*h)
    # closed-gap overlap corrections
    fr = farey(k)
    corr = 0j
    for x, y in zip(fr, fr[1:] + [F(1)]):
        i = x.denominator; j = (y.denominator if y != 1 else 1)
        if lam * (i + j) > 1/F(i*j)*F(i*j):  # placeholder never true
            pass
        ov = float(lam)/i + float(lam)/j - float(y - x)
        if ov > 0:
            lo = float(y) - float(lam)/j
            hi = float(x) + float(lam)/i
            corr += (cmath.exp(-2j*pi*h*hi) - cmath.exp(-2j*pi*h*lo)) / (-2j*pi*h)
    return -(tot - corr.real - 1j*corr.imag)  # ghat = -Bhat; Bhat = arcs - overlaps
print("REFEREE: ghat formula vs direct (max |diff| over h=1..40)")
for k in [5, 8, 12]:
    for lam in [F(1, 3*k), F(1, k+2)]:
        if lam >= F(1, k+1): continue
        mx = max(abs(ghat_formula(k, lam, h) - ghat_direct(k, lam, h)) for h in range(1, 41))
        print("  k=%2d lam=%s: max|formula-direct| = %.2e %s" % (k, lam, mx, "OK" if mx < 1e-10 else "** FAIL"))
print()
print("disc_v at the deep well {1..12}, lam=1/14 via the Ramanujan form (h-sum, H=4000) vs |G'|^2 scale:")
k = 12; lam = F(1, 14)
for v in [182, 84, 13]:
    s = sum(abs(ghat_formula(k, lam, h*v))**2 for h in range(1, 4000//v + 1)) * 2  # +-h
    print("  v=%3d: disc_v ~ %.6e  (crude bound r^2/(3v^2) with r=12: %.3e)" % (v, s, 144/(3*v*v)))
print()
print("h-profile concentration at v=182 (|ghat(hv)|^2, h=1..12):")
print("  ", ["%.2e" % abs(ghat_formula(12, F(1,14), h*182))**2 for h in range(1, 13)])
print("DONE")
