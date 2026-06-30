#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Farey-neighbor reduction + the small-depth SPREAD FAMILY: achievable rungs. klein-S38."""
import math
import numpy as np
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
def dist0(x,D): x%=D; return min(x,D-x)
def cf(fr):
    a=[]; p,q=fr.numerator,fr.denominator
    while q: a.append(p//q); p,q=q,p-(p//q)*q
    return a
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

# (1) The Farey-neighbor reduction, verified on the known covering-mins
print("="*86); print(" (1) Farey-neighbor <=> binding D = j(n-1)+1 <=> D ≡ 1 mod (n-1). Verify on covering-mins:"); print("="*86)
cmins={7:F(2,13),8:F(2,15),9:F(4,33),10:F(4,37),11:F(3,31),12:F(12,133),13:F(13,157),14:F(14,183)}
for n,M in cmins.items():
    j,D=M.numerator,M.denominator
    print(f"  n={n:2d}: M={str(M):7s} CF={str(cf(M)):12s} len={len(cf(M))}  D={D}, D mod (n-1)={D%(n-1)}, D=j(n-1)+1? {D==j*(n-1)+1}  Farey-nbr-of-1/(n-1)? {abs(j*(n-1)-D)==1}")

# (2) The small-depth SPREAD FAMILY: which rungs k are achievable (IP, small speeds Smax=4n)?
print("\n"+"="*86); print(" (2) Achievable low rungs k (IP, Smax=4n): the small-depth spread family"); print("="*86)
print(" n  | rung: 2    3    4    5    6   | min achievable rung")
for n in range(7,13):
    Smax=4*n; Dmax=2*Smax; row=[]; minr=None
    for k in range(2,7):
        t=F(k,k*(n-1)+1); S=feasible(n,t,Smax,Dmax)
        row.append("Y" if S else ".")
        if S and minr is None: minr=k
    print(f" {n:2d} |      {'    '.join(row)}   | {minr}")
