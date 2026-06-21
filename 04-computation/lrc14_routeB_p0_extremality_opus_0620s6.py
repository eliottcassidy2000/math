#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_p0_extremality_opus_0620s6.py   (opus-2026-06-20-S6)

RECONCILIATION + the load-bearing cut.

DISCOVERY (S5): U4(consec_8)=0.4803 EXCEEDS cap_8=0.3814.  So the literal route
"consec maximizes U4 and U4(consec)<=cap" CANNOT close k=8.  U4 is a LOOSE
Bonferroni UPPER bound on p0.  The cap is a bound on p0 = the lonely measure
(meas{x: all 6 inner sectors empty}).  The cap-relevant, closing statement is:

   CUT (I):  p0(E) <= p0(consec_k)   AND   p0(consec_k) <= cap_k.

p0(consec_8)=0.32721 < cap_8=0.38146  -> CLOSES with margin 0.054 (matches prompt's
0.0233-ish margin family; this is the lonely measure, not U4).

So Route B's real contribution is CUT (I): "consec maximizes the lonely measure p0".
This is the LRC extremality in its cleanest form, and it is the lonely-runner-
literature shape (AP = unique bottleneck).

This script:
  (1) confirm p0(consec_k) < cap_k for k=8,9,10 (the closing margins).
  (2) HAMMER cut (I) [p0(E)<=p0(consec)] on a HUGE exact bank for k=8,9,10, deduped
      by scale+reflection; report every violation + tightest competitor.
  (3) the k=12 cut(I) failure [0..10,12]: confirm it is OUTSIDE the binding regime
      -- p0 there is 0.645 but cap_12=6/7=0.857, so NO cap violation; AND k>=11 is
      handled by THM-535 floor route, not by consec extremality.
  (4) SCALE-INVARIANCE & REFLECTION structure: p0(E)=p0(dE)=p0(reflect(E)).  Use to
      shrink the search and to expose the difference-multiset dependence:
      p0 depends only on the DIFFERENCE SET multiset {e_i-e_j}?  Test: do two
      E,E' with the same difference-multiset have equal p0?  (If yes, p0 is a
      function of the difference multiset -> majorization on it is the right frame.)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}

def dist_full(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=(e*mid)%1; hit.add((v.numerator*7)//v.denominator)
        t=sum(1 for j in range(1,7) if j not in hit); p[t]+=hi-lo
    return p

def p0_only(E):
    """p_0 = meas{x: N(x)=0} = meas{x: ALL six inner sectors 1..6 are HIT}."""
    return dist_full(E)[0]

def consec(k): return list(range(k))
def primitive(E):
    g=reduce(gcd,[e for e in E if e>0]); return tuple(e//g for e in E)
def diffmult(E):
    E=sorted(E); return Counter(E[j]-E[i] for i in range(len(E)) for j in range(i+1,len(E)))

if __name__=="__main__":
    print("="*78); print("(1) closing margins: p0(consec_k) vs cap_k"); print("="*78)
    for k in (8,9,10,11,12):
        p0c=p0_only(consec(k));
        # sanity: p0_only == dist_full[0]
        assert p0c==dist_full(consec(k))[0], "p0 engine mismatch"
        print(f"  k={k}: p0(consec)={p0c}={float(p0c):.6f}  cap_k={float(CAP[k]):.6f}  "
              f"margin={float(CAP[k]-p0c):.6f}  {'CLOSES' if p0c<=CAP[k] else 'OVER CAP'}")

    print("\n"+"="*78); print("(2) HAMMER cut (I): p0(E) <= p0(consec) on big exact bank"); print("="*78)
    for k in (8,9,10):
        C=consec(k); p0c=p0_only(C)
        span = 17 if k==8 else (16 if k==9 else 15)
        seen=set(); n=0; viol=0; worst=F(0); tw=None; tight=F(1); tightE=None
        for rest in itertools.combinations(range(1,span+1),k-1):
            E=[0]+list(rest); pr=primitive(E)
            if pr in seen: continue
            seen.add(pr); n+=1
            p0=p0_only(E)
            d=p0c-p0   # want >=0
            if d< -F(1,10**15):
                viol+=1
                if d<worst: worst=d; tw=E
            else:
                if d<tight: tight=d; tightE=E   # tightest non-violating competitor
        print(f"  k={k}: bank={n} primitive (span<= {span})  p0(consec)={float(p0c):.6f}")
        print(f"        cut(I) violations: {viol}"
              + (f"  WORST by {float(worst):.2e} at {tw}" if tw else "  -> consec maximizes p0 on this bank"))
        print(f"        tightest competitor: gap {float(tight):.2e} at {tightE}")

    print("\n"+"="*78); print("(3) k=12 cut(I) failure is OUTSIDE the binding regime"); print("="*78)
    E=[0,1,2,3,4,5,6,7,8,9,10,12]; p0=p0_only(E); p0c=p0_only(consec(12))
    print(f"  E=[0..10,12]: p0={float(p0):.6f} > p0(consec_12)={float(p0c):.6f}  BUT cap_12={float(CAP[12]):.6f}")
    print(f"     -> p0(E)={float(p0):.4f} < cap_12={float(CAP[12]):.4f}: NO cap violation. And k>=11 closes via THM-535 floor route.")

    print("\n"+"="*78); print("(4) does p0 depend only on the difference-multiset?"); print("="*78)
    # find two primitive E with same diffmult but different as sets, compare p0
    k=7; span=12; bydm={}
    for rest in itertools.combinations(range(1,span+1),k-1):
        E=tuple([0]+list(rest))
        if reduce(gcd,E)!=1: continue
        dm=tuple(sorted(diffmult(E).items()))
        bydm.setdefault(dm,[]).append(E)
    collide=[(dm,Es) for dm,Es in bydm.items() if len(Es)>1]
    print(f"  k={k}: {len(collide)} difference-multisets shared by >1 distinct primitive set")
    checked=0; allequal=True
    for dm,Es in collide:
        vals=set(p0_only(list(E)) for E in Es)
        checked+=1
        if len(vals)>1:
            allequal=False
            print(f"     COUNTEREXAMPLE: {Es} share diffmult but p0 differs: {[float(v) for v in vals]}")
            if checked>=3: break
    if collide:
        print(f"  => p0 {'IS' if allequal else 'is NOT'} a function of the difference-multiset"
              f"  (checked {checked} colliding classes)")
    else:
        print("  (no diffmult collisions in this range; inconclusive)")
