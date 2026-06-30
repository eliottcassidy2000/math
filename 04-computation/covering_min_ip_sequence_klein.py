#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Covering-min IP (danger-circulant) -- PIN the sequence n=7..14 (dedup rows for speed). klein-S36."""
import math, sys
import numpy as np
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds

def dist0(x,D): x%=D; return min(x,D-x)
def gap_exact(S,Dmax):
    best=F(0)
    for D in range(2,Dmax+1):
        bj=0
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m<=bj: break
            if m>bj: bj=m
        if F(bj,D)>best: best=F(bj,D)
    return best
def primes_upto(n): return [p for p in range(2,n+1) if all(p%q for q in range(2,int(p**.5)+1))]

def feasible(n,t,Smax,Dmax):
    S=list(range(1,Smax+1))
    rowset=set()
    for D in range(2,Dmax+1):
        r=math.floor(t*D)
        for a in range(1,D//2+1):                       # a and D-a give same row; a<=D/2 suffices
            row=tuple(1 if dist0(s*a,D)<=r else 0 for s in S)
            if not any(row): return None
            rowset.add(row)
    for b in range(2,n+1):
        row=tuple(1 if s%b==0 else 0 for s in S)
        if not any(row): return None
        rowset.add(row)
    for p in primes_upto(n):
        rowset.add(tuple(1 if s%p!=0 else 0 for s in S))
    A=np.array(list(rowset),dtype=float)
    res=milp(c=np.ones(len(S)),constraints=LinearConstraint(A,np.ones(A.shape[0]),np.full(A.shape[0],np.inf)),
             integrality=np.ones(len(S)),bounds=Bounds(0,1))
    if not res.success: return None
    ch=sorted(S[i] for i in range(len(S)) if res.x[i]>0.5)
    return ch if len(ch)<=n-1 else None

def covering_min(n,Smax,Dmax):
    cand=sorted({F(j,D) for D in range(2,Dmax+1) for j in range(1,D//2+1) if F(1,2*n)<=F(j,D)<=F(1,3)})
    lo,hi=0,len(cand)-1; best=None
    while lo<=hi:
        mid=(lo+hi)//2; S=feasible(n,cand[mid],Smax,Dmax)
        if S is not None: best=(cand[mid],S); hi=mid-1
        else: lo=mid+1
    return best

ns = [int(x) for x in sys.argv[1:]] or [8,10,11,12,13]
print("n  | parity | covering-min M | binding D | primitive S")
print("-"*78)
for n in ns:
    Smax=4*n; Dmax=2*Smax
    res=covering_min(n,Smax,Dmax)
    if res is None: print(f"{n:2d} | infeasible in range"); continue
    t,S=res; M=gap_exact(S,2*max(S)+2); g=math.gcd(*S)
    # binding D = where gap is achieved
    bD=next(D for D in range(2,2*max(S)+3) if any(F(min((s*a)%D,D-(s*a)%D),D)==M for a in range(1,D) for s in [0]) or
            max(min(dist0(s*a,D) for s in S) for a in range(1,D))==M.numerator and D==M.denominator) if False else M.denominator
    print(f"{n:2d} | {'even' if n%2==0 else 'odd ':4s} | {str(M):>8s}={float(M):.5f} | {M.denominator:9d} | {S}")
