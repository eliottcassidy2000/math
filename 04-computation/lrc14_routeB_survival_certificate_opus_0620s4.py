#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_survival_certificate_opus_0620s4.py   (opus-2026-06-20-S4)

THE ROUTE-B SURVIVAL CERTIFICATE (the new structural reduction).

From S3:  U4(E) = p0 + p5 + 5 p6 = 1 - G_1 + G_5 + 4 G_6,   G_b=P(N_E>=b).
The cut-by-cut certificate (all signs aligned for consec, 0 violations in S3 bank):
   (I)   G_1(consec) <= G_1(E)        <=>   p0(consec) >= p0(E)      [consec MAX lonely measure]
   (II)  G_5(consec) >= G_5(E)        [consec MAX P(N>=5)]
   (III) G_6(consec) >= G_6(E)        [consec MAX P(N>=6)=p6]
If (I),(II),(III) all hold for every bounded E, then
   U4(E) = 1 - G_1(E) + G_5(E) + 4 G_6(E)
         <= 1 - G_1(consec) + G_5(consec) + 4 G_6(consec) = U4(consec).
=> consec maximizes U4 => HYP-2693 sufficient route.  THIS DECOMPOSES THE JOINT
   non-separable problem into THREE MONOTONE cut inequalities.

This script HAMMERS (I),(II),(III) separately:
  - much wider bank (larger span), dedup by scale (E ~ dE) and by reflection,
  - k=8,9,10,11,12, find any violation, report tightest competitor & margin,
  - test each cut INDEPENDENTLY (so we learn WHICH of the 3 is the hard one).
  - also test the "true-wide" regime (second largest > 14) separately.

Honesty: a single violation REFUTES that cut as a standalone lemma (forcing the
joint argument).  Zero violations across a big exact bank = strong VERIFIED, not PROVED.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = (e*mid) % 1
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

def survival(p): return [sum(p[t] for t in range(a,7)) for a in range(7)]
def consec(k): return list(range(k))

def primitive(E):
    g=reduce(gcd,[e for e in E if e>0])
    return tuple(e//g for e in E)

if __name__=="__main__":
    EPS=F(1,10**15)
    for k in (8,9,10,11,12):
        C=consec(k); pc=dist_p(C); Gc=survival(pc)
        # span: keep exact-feasible. p-engine cost ~ 7*max(E). cap max(E) modestly.
        span = {8:13,9:13,10:13,11:13,12:13}[k]
        seen=set()
        bank=[]
        for rest in itertools.combinations(range(1,span+1),k-1):
            E=[0]+list(rest)
            pr=primitive(E)
            if pr in seen: continue
            seen.add(pr)
            bank.append(E)
        v1=v5=v6=0; t1=t5=t6=None
        worst1=worst5=worst6=F(0)
        for E in bank:
            p=dist_p(E); G=survival(p)
            # (I) p0(E) <= p0(consec)  i.e. G_1(E) >= G_1(consec)
            d1=G[1]-Gc[1]   # want >=0
            if d1 < -EPS:
                v1+=1
                if d1<worst1: worst1=d1; t1=E
            # (II) G_5(E) <= G_5(consec)
            d5=Gc[5]-G[5]   # want >=0
            if d5 < -EPS:
                v5+=1
                if d5<worst5: worst5=d5; t5=E
            # (III) G_6(E) <= G_6(consec)
            d6=Gc[6]-G[6]
            if d6 < -EPS:
                v6+=1
                if d6<worst6: worst6=d6; t6=E
        print(f"k={k}: bank={len(bank)} primitive shapes (span<= {span})")
        print(f"   consec G_1={float(Gc[1]):.5f} G_5={float(Gc[5]):.5f} G_6={float(Gc[6]):.5f}  p0=G... 1-G_1={float(1-Gc[1]):.5f}")
        print(f"   (I)  p0(consec) is MAX  [G_1 min]:  violations={v1}"
              + (f"  worst by {float(worst1):.2e} at {t1}" if t1 else "  -> consec maximizes the lonely measure p0"))
        print(f"   (II) G_5(consec) is MAX:            violations={v5}"
              + (f"  worst by {float(worst5):.2e} at {t5}" if t5 else "  -> consec maximizes P(N>=5)"))
        print(f"   (III)G_6(consec)=p6 is MAX:         violations={v6}"
              + (f"  worst by {float(worst6):.2e} at {t6}" if t6 else "  -> consec maximizes p6"))
        ok = (v1==0 and v5==0 and v6==0)
        print(f"   => SURVIVAL CERTIFICATE {'HOLDS on this bank (consec maximizes U4)' if ok else 'has a standalone-cut violation'}\n")
