#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Best-known covering-min = min(small-spread IP, construction n/Phi6); the FAREY RUNG invariant. klein-S37."""
import math
import numpy as np
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
def dist0(x,D): x%=D; return min(x,D-x)
def gap_exact(S,Dmax):
    best=F(0)
    for D in range(2,Dmax+1):
        bj=0
        for a in range(1,D):
            m=min(dist0(s*a,D) for s in S)
            if m>bj: bj=m
        if F(bj,D)>best: best=F(bj,D)
    return best
def primes_upto(n): return [p for p in range(2,n+1) if all(p%q for q in range(2,int(p**.5)+1))]
def feasible(n,t,Smax,Dmax):
    S=list(range(1,Smax+1)); rowset=set()
    for D in range(2,Dmax+1):
        r=math.floor(t*D)
        for a in range(1,D//2+1):
            row=tuple(1 if dist0(s*a,D)<=r else 0 for s in S)
            if not any(row): return None
            rowset.add(row)
    for b in range(2,n+1):
        row=tuple(1 if s%b==0 else 0 for s in S)
        if not any(row): return None
        rowset.add(row)
    for p in primes_upto(n): rowset.add(tuple(1 if s%p!=0 else 0 for s in S))
    A=np.array(list(rowset),dtype=float)
    res=milp(c=np.ones(len(S)),constraints=LinearConstraint(A,np.ones(A.shape[0]),np.full(A.shape[0],np.inf)),
             integrality=np.ones(len(S)),bounds=Bounds(0,1))
    if not res.success: return None
    ch=sorted(S[i] for i in range(len(S)) if res.x[i]>0.5)
    return ch if len(ch)<=n-1 else None
def small_spread_min(n,Smax,Dmax):
    cand=sorted({F(j,D) for D in range(2,Dmax+1) for j in range(1,D//2+1) if F(1,2*n)<=F(j,D)<=F(1,3)})
    lo,hi=0,len(cand)-1; best=None
    while lo<=hi:
        mid=(lo+hi)//2; S=feasible(n,cand[mid],Smax,Dmax)
        if S is not None: best=(cand[mid],S); hi=mid-1
        else: lo=mid+1
    return best

def rung(M,n):  # M = [0; n-1, k]; k = 1/(1/M - (n-1)); inf if M=1/(n-1)
    d=F(1,1)/M-(n-1)
    return F(1,1)/d if d!=0 else math.inf

known={7:(F(2,13),[1,2,5,6,7,8]),8:(F(2,15),None),9:(F(4,33),None),11:(F(3,31),None)}
print("n  | best-known covering-min | rung k | source        | construction n/Phi6 (rung n)")
print("-"*92)
for n in [7,8,9,10,11,12,13,14]:
    constr=F(n, n*n-n+1)                 # = n/Phi_6(n) = [0; n-1, n], rung n
    if n in known:
        M=known[n][0]; src="spread (IP)"
    else:
        Smax=4*n; Dmax=2*Smax
        ss=small_spread_min(n,Smax,Dmax)
        Mss=ss[0] if ss else F(1,1)
        M=min(Mss,constr); src="spread (IP)" if M==Mss and M<constr else ("construction" if M==constr else "spread (IP)")
    k=rung(M,n)
    print(f"{n:2d} | {str(M):>8s}={float(M):.5f}      | {str(k):>6s} | {src:13s} | {str(constr)}={float(constr):.5f}")
print("\nRUNG SEQUENCE k(n) (best-known): the Farey-ladder position of the covering-min, M(n)=[0;n-1,k]=k/(k(n-1)+1).")
print("Floor 1/n = rung 1; construction n/Phi6 = rung n; ceiling 1/(n-1) = rung inf. Covering-min rung in [2,n].")
