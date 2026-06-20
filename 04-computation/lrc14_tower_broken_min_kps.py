#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TARGETED: minimize meas(G_C) over TOWER-BROKEN 12-cores (miss >=1 of {1,2,4,8}).
If min over tower-broken >= 426/35035, then meas<426/35035 ==> tower-full (HYP-2661 general).
Tower-broken cores are 'larger' (missing small dense elements) so should have BIGGER measure;
we confirm the min directly.  Exhaustive over [1,cap] restricted to tower-broken.  EXACT.
kind-pasteur-2026-06-19.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
THETA=Fraction(1,14); THR2=Fraction(426,35035)
TOWER=(1,2,4,8)

def lonely_measure(C):
    segs=[]
    for d in C:
        w=THETA/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for s in (-1,0,1):
                a=max(lo+s,Fraction(0)); b=min(hi+s,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); cur=Fraction(-1); U=[]
    for a,b in segs:
        if a>cur: U.append([a,b]); cur=b
        elif b>cur: U[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in U)

def run(cap):
    best=None; cntbroken=0
    for C in itertools.combinations(range(1,cap+1),12):
        if all(b in C for b in TOWER): continue   # skip tower-full
        if reduce(gcd,C)!=1: continue
        cntbroken+=1
        L=lonely_measure(C)
        if best is None or L<best[0]: best=(L,C)
    print(f"cap={cap}: tower-BROKEN primitive cores={cntbroken}")
    L,C=best
    print(f"   MIN tower-broken meas={L}={float(L):.6f}  C={C}")
    print(f"   THR2=426/35035={float(THR2):.6f}   min>=THR2? {L>=THR2}  (if True: HYP-2661 holds on [1,{cap}])")
    return best

if __name__=="__main__":
    for cap in [16,18,20,22]:
        run(cap); print()
