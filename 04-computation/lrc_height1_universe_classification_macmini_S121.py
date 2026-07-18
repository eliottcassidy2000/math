#!/usr/bin/env python3
"""
Complete classification of TIGHT sets in the HEIGHT-<=1 UNIVERSE (mac-mini-S121)
================================================================================
Universe: A subset of {1..n} u {(n+1)+r : 1<=r<=n}, |A| = n  (every element is a
base speed or a height-1 lift). BOTH known mechanisms live here:
  - GW doubling: n=7 -> 12 = 4+8 ; n=13 -> 24 = 10+14
  - Lucas-type: n=4 -> 7 = 2+5 ; n=5 -> 9 = 3+6 ; n=7 -> 11,13 = 3+8, 5+8
Prunes: covering lemma (multiple of every q<=n) + sound batched numeric prescreen
(grid max <= true M, so rejection is never a false negative); exact Q confirm.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def M_exact(S, np1):
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
                if num*np1<q: break
            if num*np1>=q:
                c=F(num,q)
                if c>best: best=c
    return best

for n in range(4,14):
    np1=n+1; L=1.0/np1; Lex=F(1,np1)
    pool=list(range(1,n+1))+[np1+r for r in range(1,n+1)]
    N=max(pool)
    REQ=list(range(2,n+1)); full=(1<<len(REQ))-1
    mask={v:sum((1<<i) for i,q in enumerate(REQ) if v%q==0) for v in pool}
    cands=[]
    for A in combinations(pool,n):
        m=0
        for v in A: m|=mask[v]
        if m!=full: continue
        if reduce(gcd,A)!=1: continue
        cands.append(A)
    C=np.array(cands,dtype=np.int64) if cands else np.zeros((0,n),dtype=np.int64)
    T=6000; t=np.arange(1,T)/T
    Cl=np.empty((N+1,T-1))
    for v in range(1,N+1):
        x=(v*t)%1.0; Cl[v]=np.minimum(x,1-x)
    keep=[]
    for s in range(0,len(C),3000):
        B=C[s:s+3000]
        cur=Cl[B[:,0]].copy()
        for j in range(1,n): np.minimum(cur,Cl[B[:,j]],out=cur)
        for i in np.where(cur.max(axis=1)<=L+2e-4)[0]: keep.append(tuple(B[i]))
    tight=[list(A) for A in keep if M_exact(A,np1)==Lex]
    seg=list(range(1,n+1))
    spor=[A for A in tight if A!=seg]
    print(f"n={n:2d}: covering+primitive candidates {len(C):8d} | prescreen {len(keep):4d} | TIGHT {len(tight)} | sporadic {len(spor)}", flush=True)
    for A in spor:
        rem=[v for v in range(1,n+1) if v not in A]; add=[v for v in A if v>n]
        print(f"        {A}   removed {rem}  added {add} = (n+1)+{[v-np1 for v in add]}", flush=True)
