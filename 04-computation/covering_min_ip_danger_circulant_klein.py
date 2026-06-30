#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Covering-min as an IP over the DANGER-CIRCULANT (klein-S36, corrected).

M(S) = max over D' of gap_{D'}(S); the COVERING-MIN = min over primitive covering sets of M(S).
CORRECT IP (binary search on t): a primitive size-(n-1) covering S has honest M(S) <= t iff for EVERY
modulus D' and rotation a, some speed s has dist(s.a mod D',0) <= floor(t.D')  [the danger covers at every
scale -- the SET-COVER reframe, observer at 0]. The DUAL (owner's 'copy observer to all n points'): the n
runner-positions are an INDEPENDENT SET in the danger-circulant C_{D'}({1..floor(t.D')}) at every D'.
Constraints: + resonance-killing (THM-523: each b<=n divides some speed) + primitivity (each prime p<=n
fails to divide some speed). Smallest feasible t = the covering-min."""
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
            m=D
            for s in S:
                d=dist0(s*a,D)
                if d<m: m=d
                if m<=bj: break
            if m>bj: bj=m
        if F(bj,D)>best: best=F(bj,D)
    return best

def primes_upto(n): return [p for p in range(2,n+1) if all(p%q for q in range(2,int(p**.5)+1))]

def feasible(n, t, Smax, Dmax):
    """min #speeds in {1..Smax} with gap_{D'}<=t for all D'<=Dmax, resonance b<=n, primitive. <= n-1 ?"""
    S=list(range(1,Smax+1)); idx={s:i for i,s in enumerate(S)}
    rows=[]
    for D in range(2,Dmax+1):
        r=math.floor(t*D)
        for a in range(1,D):
            row=[1.0 if dist0(s*a,D)<=r else 0.0 for s in S]
            if sum(row)==0: return None            # uncoverable -> infeasible at this t
            rows.append(row)
    for b in range(2,n+1):
        row=[1.0 if s%b==0 else 0.0 for s in S]; 
        if sum(row)==0: return None
        rows.append(row)
    for p in primes_upto(n):                        # primitivity: some speed not divisible by p
        rows.append([1.0 if s%p!=0 else 0.0 for s in S])
    A=np.array(rows)
    res=milp(c=np.ones(len(S)), constraints=LinearConstraint(A,np.ones(len(rows)),np.full(len(rows),np.inf)),
             integrality=np.ones(len(S)), bounds=Bounds(0,1))
    if not res.success: return None
    chosen=sorted(S[i] for i in range(len(S)) if res.x[i]>0.5)
    return chosen if len(chosen)<=n-1 else None

print("="*86); print(" COVERING-MIN via the danger-circulant IP (binary search on the honest gap t)"); print("="*86)
targets={7:F(2,13), 9:F(4,33)}
for n,Smax in [(7,14),(9,36)]:
    Dmax=2*Smax+2
    cand=sorted({F(j,D) for D in range(2,Dmax+1) for j in range(1,D//2+1)})
    cand=[t for t in cand if F(1,2*n)<=t<=F(1,3)]    # covering-min is in [~1/(2n), 1/3]
    lo,hi=0,len(cand)-1; best=None
    while lo<=hi:                                     # smallest feasible t
        mid=(lo+hi)//2; S=feasible(n,cand[mid],Smax,Dmax)
        if S is not None: best=(cand[mid],S); hi=mid-1
        else: lo=mid+1
    t,S=best; Mtrue=gap_exact(S,2*max(S)+2)
    g=math.gcd(*S) if len(S)>1 else S[0]
    print(f"\n n={n}: IP covering-min t={t}={float(t):.5f}, S={S} (|S|={len(S)}, gcd={g}, primitive={g==1})")
    print(f"   honest M(S)={Mtrue}={float(Mtrue):.5f}; target {targets[n]}={float(targets[n]):.5f} -> "
          f"{'*** MATCHES ***' if Mtrue==targets[n] else 'differs'}")
