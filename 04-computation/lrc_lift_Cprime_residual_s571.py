#!/usr/bin/env python3
"""Lift lemma -> C': n|v => loose. Test the partial proof (large multiples dodgeable
via one-arc-radius perturbation of the LRC(n-1) optimum) and isolate the residual
(small multiples, esp. v=n). opus-2026-06-02-S571."""
from fractions import Fraction
from math import gcd
import random
def dist(x): x%=1; return min(x,1-x)
def safe_measure(V,n):
    THR=Fraction(1,n); eps=set()
    for v in V:
        for k in range(v+1):
            for s in(-1,1): eps.add(Fraction(k*n+s,n*v)%1)
    pts=sorted(eps); meas=Fraction(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1)
        if all(dist(v*((a+ln/2)%1))>=THR for v in V): meas+=ln
    return meas
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def report(n, sel):
    m=n-1; B=n+10; rng=random.Random(n*7+sel)
    cnt=0; loose=0; minmu=None; bad=[]
    others_pool=[x for x in range(1,B+1) if x%n!=0]
    for _ in range(800):
        others=set(rng.sample(others_pool, m-1))
        if sel==0: mult=n            # v=n exactly (w=1, hardest)
        elif sel==1: mult=n*rng.randint(1,3)   # small multiples
        else:        mult=n*rng.randint(4,B//n if B//n>=4 else 4)  # large multiples
        V=tuple(sorted(others|{mult}))
        if len(V)!=m: continue
        Vp=prim(V)
        if not any(v%n==0 for v in Vp):  # primitivization may kill the multiple
            continue
        cnt+=1; mu=safe_measure(Vp,n)
        if mu>0: loose+=1
        else: bad.append(Vp)
        if minmu is None or mu<minmu: minmu=mu
    tag={0:'v=n exactly',1:'small mult (w<=3)',2:'large mult (w>=4)'}[sel]
    print(f'  n={n:2d} [{tag:18s}]: {cnt} configs; loose={loose}/{cnt}; min measure={float(minmu) if minmu is not None else 0:.5f}; tight-with-mult={len(bad)}')
if __name__=='__main__':
    print("C' test: n|v => loose, split by multiplier size (residual = small multiples)")
    for n in [6,8,10,12,14]:
        for sel in (0,1,2): report(n,sel)
