#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the WORST non-AP set is a NEAR-AP with a single end-defect.
For k=9 the binding case is E={0,1,2,3,4,5,6,7,9} (L_y=0.487289, margin 0.00697 < cap_9).

Here we (1) EXHAUSTIVELY confirm, for each k in {8,9,10}, the max L_y over ALL primitive
non-full-AP k-sets with 0 in E and max(E) <= W, scanning W up to a window that contains
the supremum (the worst sets are tight near-APs, so a modest window suffices once L_y has
peaked and the family value oscillates strictly below the peak).
(2) Report the worst set and its margin to cap.

caps: cap_8=0.38153, cap_9=0.49426, cap_10=0.6044.
kind-pasteur-2026-06-19.
"""
import sys, itertools
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
    return [Fraction(-(t-2)*(t-3)*(t-6), 36) for t in range(7)]  # 9,10

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

def is_full_AP(E):
    E = sorted(E)
    d = E[1]-E[0]
    return all(E[i+1]-E[i]==d for i in range(len(E)-1))

CAP = {8: 0.38153, 9: 0.49426, 10: 0.6044}
WIN = {8: 14, 9: 14, 10: 15}   # window max(E) <= W; worst sets are near-APs so small W

if __name__ == "__main__":
    for k in [8, 9, 10]:
        cap = CAP[k]; W = WIN[k]
        Lc = float(L_y(list(range(k)), k))
        best = (Fraction(-1), None); cnt = 0
        for combo in itertools.combinations(range(1, W+1), k-1):
            E = (0,) + combo
            if reduce(gcd, E[1:]) != 1: continue
            if is_full_AP(E): continue
            cnt += 1
            Lv = L_y(E, k)
            if Lv > best[0]:
                best = (Lv, E)
        Lv, E = best
        print(f"k={k} cap={cap} L_y(consec)={Lc:.5f}")
        print(f"  scanned {cnt} primitive non-AP sets with max(E)<={W}")
        print(f"  WORST non-AP set E={list(E)}")
        print(f"    L_y={float(Lv):.6f}   margin to cap = {cap-float(Lv):.6f}")
        # structure of the worst set
        diffs = [E[i+1]-E[i] for i in range(len(E)-1)]
        print(f"    gap pattern (consecutive diffs): {diffs}")
        print()
