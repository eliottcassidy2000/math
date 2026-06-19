#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_subtorus17_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL, TARGETED)
Most dangerous regime for the 1/7 consecutive-min claim: SUBTORUS / common-factor / dense-core
shapes (many short relation vectors => Weyl identity predicts dips for 2/7; we test whether ANY
dip below consecutive appears for 1/7). Focus k=11,12,13 (smallest consec margin).
Exhaustive-confirm any candidate within 0.02 (float) of consec.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(31337)
ONE7=F(1,7); ONE7f=1.0/7.0

def mu17(E):
    E=sorted(set(E)); k=len(E)
    if k==1: return F(1)
    bp={F(0),F(1)}
    for i in range(k):
        for j in range(i+1,k):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bps=sorted(b for b in bp if F(0)<=b<=F(1)); tot=F(0)
    for a,b in zip(bps,bps[1:]):
        if a==b: continue
        mid=(a+b)/2; pts=[]
        for e in E:
            val=e*mid; n=val-(val%1); pts.append((val-n,e,n))
        pts.sort(key=lambda t:t[0]); m=len(pts); ivs=[]
        for i in range(m):
            (_,ei,ni)=pts[i]
            if i<m-1: (_,ej,nj)=pts[i+1]; al=F(ej-ei); be=F(ni-nj)
            else: (_,e0,n0)=pts[0]; al=F(e0-ei); be=F(ni-n0+1)
            if al==0:
                if be>ONE7: ivs.append((a,b))
            else:
                xs=(ONE7-be)/al
                if al>0: lo,hi=max(a,xs),b
                else: lo,hi=a,min(b,xs)
                if lo<hi: ivs.append((lo,hi))
        if not ivs: continue
        ivs.sort(); clo,chi=ivs[0]
        for lo,hi in ivs[1:]:
            if lo<=chi: chi=max(chi,hi)
            else: tot+=chi-clo; clo,chi=lo,hi
        tot+=chi-clo
    return tot
def mu17f(E,N=10000):
    cnt=0
    for t in range(N):
        x=(t+0.5)/N
        pts=sorted((e*x)%1.0 for e in E)
        mg=0.0
        for i in range(len(pts)-1):
            g=pts[i+1]-pts[i]
            if g>mg: mg=g
        w=pts[0]+1.0-pts[-1]
        if w>mg: mg=w
        if mg>ONE7f: cnt+=1
    return cnt/N

for k in (11,12,13):
    consec=mu17(list(range(k))); consecf=float(consec)
    cands=set()
    # subtorus: E = c*A union small fill, for many c and core sizes
    for c in (2,3,4,5,6,7,9,11,13):
        for core in range(2,k):
            A=[c*i for i in range(core)]
            fill=list(range(1,(k-core)+1))
            E=tuple(sorted(set([0]+A[1:]+fill)))
            if len(E)==k: cands.add(E)
            # also append a high tail instead of low fill
            E2=tuple(sorted(set(A+list(range(c*core+1,c*core+1+(k-core))))))
            if len(E2)==k and 0 in E2: cands.add(E2)
    # union of two APs with different steps (two short-relation cores)
    for s1 in (1,2,3):
        for s2 in (2,3,5,7):
            for n1 in range(2,k):
                n2=k-n1
                if n2<1: continue
                A=[s1*i for i in range(n1)]
                off=max(A)+s2
                B=[off+s2*i for i in range(n2)]
                E=tuple(sorted(set(A+B)))
                if len(E)==k and 0 in E: cands.add(E)
    # dense perforated cores at moderate spread
    for sp in range(k, 2*k+6):
        for _ in range(2000):
            body=sorted(random.sample(range(1,sp+1),min(k-1,sp)))
            E=tuple(sorted(set([0]+body)))
            if len(E)==k: cands.add(E)
    viol=[]; near=0
    for E in cands:
        if len(set(E))!=k or 0 not in E: continue
        mf=mu17f(list(E),8000)
        if mf<consecf+0.02:
            near+=1
            me=mu17(list(E))
            if me<consec: viol.append((me,E))
    if viol:
        viol.sort()
        print(f"  k={k}: *** VIOLATION *** consec={consecf:.5f} min={float(viol[0][0]):.5f} at {viol[0][1]}")
    else:
        print(f"  k={k}: consec={consecf:.5f}; {len(cands)} subtorus/two-AP/dense cands, {near} near-ties, NONE below consec")
print("\nsubtorus-targeted: 1/7 consecutive-min survives =", "see above (no VIOLATION => holds)")
