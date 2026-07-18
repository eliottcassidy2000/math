#!/usr/bin/env python3
"""
Batched-prescreen height census: PIN n=8 and n=9  (mac-mini-2026-07-18-S115)
===========================================================================
Goal: max(A) over PRIMITIVE TIGHT n-sets (M(A)=1/(n+1)); conjecture max <= 3n.

Two prunes:
 (P1) COVERING LEMMA (proved, S114/THM-1031): a tight n-set contains a multiple of
      every q <= n.  Required generators: q in {5,6,7,8} (n=8), {5,6,7,8,9} (n=9);
      2,3,4 are implied by 6 and 8.
 (P2) BATCHED NUMERIC PRESCREEN (sound): on a fixed grid, gridmax(A) <= M(A) always.
      So gridmax > 1/(n+1) + tol  =>  M > 1/(n+1)  =>  NOT tight.  Rejection is
      therefore never a false negative; survivors are confirmed in exact Q.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys

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

def census(n, N, REQ, T=6000, tol=2e-4, chunk=4000):
    L=1.0/(n+1); Lex=F(1,n+1)
    print(f"=== n={n}: max element <= {N} (conjecture max <= 3n = {3*n}) ===", flush=True)
    # (P1) covering filter via bitmasks
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
    print(f"  covering+primitive candidates: {len(C):,}", flush=True)
    if len(C)==0: return
    # (P2) batched prescreen
    t=np.arange(1,T)/T
    Cl=np.empty((N+1,T-1))
    for v in range(1,N+1):
        x=(v*t)%1.0; Cl[v]=np.minimum(x,1-x)
    keep=[]
    for s in range(0,len(C),chunk):
        B=C[s:s+chunk]
        cur=Cl[B[:,0]]
        for j in range(1,n):
            np.minimum(cur,Cl[B[:,j]],out=cur)
        gm=cur.max(axis=1)
        idx=np.where(gm<=L+tol)[0]
        for i in idx: keep.append(tuple(B[i]))
    print(f"  survived numeric prescreen: {len(keep):,}  ({100*len(keep)/len(C):.3f}%)", flush=True)
    # exact confirmation
    tight=[]
    for A in keep:
        if M_exact(A)==Lex: tight.append(list(A))
    print(f"  PRIMITIVE TIGHT (exact): {len(tight)}", flush=True)
    for A in tight: print(f"     {A}  max={max(A)}  max/n={max(A)/n:.3f}", flush=True)
    if tight:
        mx=max(max(A) for A in tight)
        print(f"  => n={n}: MAX over primitive tight = {mx} = {mx/n:.3f}*n ; <= 3n? {mx<=3*n}", flush=True)
    print(flush=True)

census(11,25,list(range(2,12)))
census(12,26,list(range(2,13)))
