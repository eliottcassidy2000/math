#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_consec_vs_L7_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D crux)

CRUX TEST: does the CONSECUTIVE-MAXIMIZER claim subsume / contradict the L7
balanced-multi-cluster worst case?

The exhaustive search only covered BOUNDED N (e.g. k=8 within {1..16}).  L7 is about
sets with LARGE spread (far clusters f1,f2 >> base).  If consecutive [1..k] really is
the global measS7-maximizer over ALL positive-integer k-sets, then:
   - the worst case for the LRC cover bound is consec (a single family);
   - the balanced-multi-cluster regime (L7) has STRICTLY SMALLER p0 and is DOMINATED;
   - L7 would be unnecessary -- a much simpler closure.

This is a strong claim; test it hard against the very configurations L7 fears:
  (a) wide-spread random k-sets (elements up to 2000): can any beat consec?
  (b) the balanced two-cluster family {1..t} u {f, f+1, ..., f+(k-t-1)} for growing f
      (the L7 'two far clusters' shape) -- track p0 as f -> inf;
  (c) r-cluster balanced families (r=2,3) with large gaps;
  (d) arithmetic-progression sets with various common differences d.
If NONE beats consec, flag CONSECUTIVE-MAXIMIZER as a strong HYP that potentially
SUBSUMES L7.  If something beats it, that is itself the new worst case to chase.
"""
import itertools, random
from fractions import Fraction as Fr
from math import gcd

P=7
def sector(yf): return int(P*yf)
def breakpoints(E):
    bp={Fr(0),Fr(1)}
    for e in E:
        if e==0: continue
        for t in range(0,P*e): bp.add(Fr(t,P*e))
    return sorted(bp)
def measS7(E):
    E=[int(e) for e in E if int(e)!=0]
    xs=breakpoints(E); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        if len(set(sector((e*mid)%1) for e in E))==P: tot+=(b-a)
    return tot

CAP={8:Fr(2243,5880),9:Fr(1979,4004),10:Fr(55,91)}

def main():
    random.seed(20260621)
    print("#"*80)
    print("# CRUX: does consecutive [1..k] dominate the L7 spread/cluster families?")
    print("#"*80)
    for k in [8,9]:
        consec=measS7(list(range(1,k+1)))
        cap=CAP[k]
        print(f"\n===== k={k}:  consec [1..{k}] p0 = {float(consec):.5f}   cap={float(cap):.5f} =====")

        # (a) random wide-spread sets
        maxr=Fr(0); maxE=None; nbeat=0
        for _ in range(4000):
            E=sorted(random.sample(range(1,2001),k))
            v=measS7(E)
            if v>maxr: maxr=v; maxE=E
            if v>consec+Fr(1,10**9): nbeat+=1
        print(f"  (a) 4000 random wide sets in [1,2000]: max p0={float(maxr):.5f}"
              f"  beats-consec={nbeat}   {'<= consec' if maxr<=consec else 'EXCEEDS consec'}")
        if maxE: print(f"      argmax sample: {maxE[:12]}{'...' if len(maxE)>12 else ''}")

        # (b) balanced two-cluster {1..t} u {f..f+(k-t-1)}, growing f
        print(f"  (b) two-cluster {{1..t}} u {{f..}}  (L7 shape), t=k//2:")
        t=k//2
        for f in [20,50,100,300,1000,3000]:
            E=list(range(1,t+1))+list(range(f,f+(k-t)))
            v=measS7(E)
            print(f"      f={f:5d}: p0={float(v):.5f}  {'(>consec!)' if v>consec else ''}")

        # (c) r=3 balanced clusters
        print(f"  (c) three-cluster balanced (k split 3 ways), large gaps:")
        for f1,f2 in [(50,2500),(100,10000),(40,1600)]:
            sz=k//3
            E=list(range(1,sz+1))+list(range(f1,f1+sz))+list(range(f2,f2+(k-2*sz)))
            v=measS7(E)
            print(f"      clusters at 1,{f1},{f2}: p0={float(v):.5f}  {'(>consec!)' if v>consec else ''}")

        # (d) AP sets with common difference d (these dilate to consec when gcd handles d)
        print(f"  (d) AP {{1,1+d,...,1+(k-1)d}} for various d (note: dilation-invariant):")
        for d in [1,2,3,5,7,11]:
            E=[1+i*d for i in range(k)]
            v=measS7(E)
            print(f"      d={d:2d}: p0={float(v):.5f}  {'(>consec!)' if v>consec+Fr(1,10**9) else ('(=consec)' if abs(v-consec)<Fr(1,10**9) else '')}")

    print("\nDONE.")

if __name__=="__main__":
    main()
