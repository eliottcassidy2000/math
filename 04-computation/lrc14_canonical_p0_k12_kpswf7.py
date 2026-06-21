#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- reconcile my W_a-sum finding with the CANONICAL p0=measS7 over [0,1].

There are TWO measS7's in the repo:
  (A) route3 W_a-sum: sum_{a=1..6} W_a(E), coverage over the 6 ACTIVE cells only
      (excludes cell_0, the always-covered "short"). This is what mac-mini ROUTE-3
      and my lrc14_*_kpswf7 scripts used.
  (B) canonical (threadB / Delsarte): measS7(E) = meas{y in [0,1] : all 7 sectors hit}
      over the WHOLE circle. THIS is p0, the wall object: the wall is p0(E) <= cap_k.

These are DIFFERENT functionals. My k=12 result (E*=(0..10,12) beats consec on the
W_a-sum) must be re-tested on the CANONICAL p0. Critically, the wall is
  p0(E) <= cap_k,   cap_k - Q(k-1) ~ 0.185,0.132,0.157,0.194,0.255 (k=8..12).
So even if consec is not the exact argmax of p0, the wall can still hold with large
margin. Determine BOTH:
  (1) Does consec maximize the CANONICAL p0 on full-residue/primitive/span<=14, k=8..12?
  (2) For any rival beating consec, does it still satisfy p0 <= cap_k (margin)?
  (3) Print cap_k = Q(k-1)+margin and Q(k-1) for reference.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
P = 7

def measS7_full(E):
    """CANONICAL: coverage over whole [0,1] (= p0)."""
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        d = P*abs(e)
        for m in range(0, d+1): bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        hit=set()
        for e in E:
            v=(e*mid)%1
            hit.add(int(v*P))
        if len(hit)==P: tot += hi-lo
    return tot

def stirling2(n,k): return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1))//math.factorial(k)
def iid_k(k):
    if k<P: return F(0)
    return F(math.factorial(P)*stirling2(k,P), P**k)

def Q(k):  # Q(k-1) target = consec base value; compute as measS7_full(consec(k))? NO.
    # Q(k-1) in canon = decorrelated plateau = measS7_full(consec(k))? Check empirically.
    return measS7_full(tuple(range(k)))

def full_residue(E): return set(e % P for e in E) == set(range(P))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))

# cap_k - Q(k-1) margins from canon (codex S46): k=8..12
MARGIN = {8:F(185,1000),9:F(132,1000),10:F(157,1000),11:F(194,1000),12:F(255,1000)}

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*88)
    print("CANONICAL p0=measS7_full over [0,1]: does consec maximize? does wall hold?")
    print("="*88)
    # First: relation between W_a-sum and full. consec_8:
    for k in [8,9,10,11,12]:
        C=consec(k)
        print(f"k={k}: p0_full(consec) = {float(measS7_full(C)):.6f}")

    print("\n--- E* k=12 canonical check ---")
    Es=(0,1,2,3,4,5,6,7,8,9,10,12); C=consec(12)
    pE=measS7_full(Es); pC=measS7_full(C)
    print(f"  p0_full(E*={Es}) = {pE} = {float(pE):.8f}")
    print(f"  p0_full(consec)  = {pC} = {float(pC):.8f}")
    print(f"  E* beats consec on CANONICAL p0? {pE>pC}  delta={float(pE-pC):+.8f}")

    print("\n--- global argmax + wall check, CANONICAL p0, full-residue/primitive ---")
    for k in [8,9,10,11,12]:
        for box in [14]:
            S=stratum(k,box)
            pC=measS7_full(consec(k))
            cap = pC + MARGIN[k]   # cap_k ~ Q(k-1)+margin; consec value used as Q proxy
            best=None; bv=F(-1); nbeat=0; ncapviol=0; worst_over_consec=F(0)
            for E in S:
                v=measS7_full(E)
                if v>bv: bv=v; best=E
                if v>pC+F(1,10**9):
                    nbeat+=1
                    if v-pC>worst_over_consec: worst_over_consec=v-pC
                if v>cap: ncapviol+=1
            print(f"  k={k} box{box}: |strat|={len(S):4d} | consec p0={float(pC):.5f} | argmax={'CONSEC' if best==consec(k) else best}")
            print(f"        max p0={float(bv):.5f} | #rivals>consec={nbeat} (worst +{float(worst_over_consec):.5f}) | cap_k~{float(cap):.5f} | #cap-violations={ncapviol}")
