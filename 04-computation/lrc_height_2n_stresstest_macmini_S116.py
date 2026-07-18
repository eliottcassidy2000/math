#!/usr/bin/env python3
"""Stress-test max(A) <= 2n at the SPORADIC n (4,5,7) with bounds far above 2n."""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def M_exact(S):
    S=sorted(set(S)); best=F(0); dens={2*a for a in S}
    for a,b in combinations(S,2):
        dens.add(a+b); dens.add(b-a)
    dens.discard(0)
    for q in dens:
        for m in range(1,q):
            num=q
            for v in S:
                r=(v*m)%q; d=min(r,q-r)
                if d<num: num=d
                if num==0: break
            if num>0:
                c=F(num,q)
                if c>best: best=c
    return best
def census(n,N,T=8000,tol=2e-4,chunk=3000):
    L=1.0/(n+1); Lex=F(1,n+1); REQ=list(range(2,n+1))
    full=(1<<len(REQ))-1
    mask={v:sum((1<<i) for i,q in enumerate(REQ) if v%q==0) for v in range(1,N+1)}
    cands=[]
    for A in combinations(range(1,N+1),n):
        m=0
        for v in A: m|=mask[v]
        if m!=full: continue
        if reduce(gcd,A)!=1: continue
        cands.append(A)
    C=np.array(cands,dtype=np.int64)
    print(f"=== n={n}, N={N} (2n={2*n}): covering+primitive candidates {len(C):,}",flush=True)
    if not len(C): return
    t=np.arange(1,T)/T
    Cl=np.empty((N+1,T-1))
    for v in range(1,N+1):
        x=(v*t)%1.0; Cl[v]=np.minimum(x,1-x)
    keep=[]
    for s in range(0,len(C),chunk):
        B=C[s:s+chunk]; cur=Cl[B[:,0]].copy()
        for j in range(1,n): np.minimum(cur,Cl[B[:,j]],out=cur)
        for i in np.where(cur.max(axis=1)<=L+tol)[0]: keep.append(tuple(B[i]))
    tight=[list(A) for A in keep if M_exact(A)==Lex]
    print(f"    prescreen survivors {len(keep):,}; PRIMITIVE TIGHT {len(tight)}",flush=True)
    for A in sorted(tight,key=max): print(f"      {A} max={max(A)} max/n={max(A)/n:.3f}",flush=True)
    mx=max(max(A) for A in tight) if tight else 0
    print(f"    => max={mx}={mx/n:.3f}*n ; <= 2n? {mx<=2*n}",flush=True); print(flush=True)
census(4,44); census(5,50); census(7,40)
