#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): close the window. The worst non-AP k-set is a single-defect near-AP.
A single-defect near-AP (primitive, step 1) is consec_k with ONE internal/edge gap
enlarged: E = {0,1,...,p-1, p+g, p+g+1, ..., } i.e. insert a jump of size (1+g) at
position p, g>=1.  Equivalently E = consec block of length L on the left + consec block
of length k-L on the right, separated by a gap of size s>=2:
    E = {0..L-1} U {L-1+s .. L-1+s+(k-L)-1},  s>=2, 1<=L<=k-1.
These are exactly the 2-dim-GAP sets of "near-AP" type (excess small).  Scan L and s
over a LARGE range; confirm the supremum over this family is the s=2 minimal defect found
earlier and is < cap_k.  This shows no large-position defect beats the in-window worst.

kind-pasteur-2026-06-19.
"""
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        p[N_at(E, (lo+hi)/2)] += (hi-lo)
    return p

def g_poly(k):
    if k == 8:
        return [Fraction((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)]
    return [Fraction(-(t-2)*(t-3)*(t-6), 36) for t in range(7)]

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

CAP = {8: 0.38153, 9: 0.49426, 10: 0.6044}

if __name__ == "__main__":
    SMAX = 60
    for k in [8, 9, 10]:
        cap = CAP[k]
        best = (Fraction(-1), None)
        for L in range(1, k):                  # left block length
            for s in range(2, SMAX+1):         # gap size (>=2 -> non-AP)
                left = list(range(L))
                right = [L-1+s + i for i in range(k-L)]
                E = left + right
                if reduce(gcd, [e for e in E if e]) != 1:
                    continue
                Lv = L_y(E, k)
                if Lv > best[0]:
                    best = (Lv, (L, s, E))
        Lv, (L, s, E) = best
        print(f"k={k} cap={cap}")
        print(f"  worst single-defect near-AP: left len={L}, gap s={s}, E={E}")
        print(f"    L_y={float(Lv):.6f}  margin to cap = {cap-float(Lv):.6f}")
        print()
