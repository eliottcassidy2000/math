#!/usr/bin/env python3
"""
Height-bound census: pin the n=8/9 threshold  (mac-mini-2026-07-18-S114)
=======================================================================
CONJECTURE (S113): a PRIMITIVE tight n-set has max(A) <= 3n.
Data so far: max/n = 1.00 (n=3), 1.75 (n=4), 1.80 (n=5), 1.00 (n=6), 2.57 (n=7).
n=8..12 untested -> pin it.

PRUNE (proved): a tight n-set contains a MULTIPLE OF EVERY q <= n.
  At t=a/q (gcd(a,q)=1), some v has ||va/q|| <= M(A) = 1/(n+1), i.e. |va|_q <= q/(n+1) < 1
  since q <= n; so |va|_q = 0, i.e. q | va, and gcd(a,q)=1 gives q | v.  QED
(This is the corpus's "covering" condition; here it is the search prune.)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys

def M_full(S):
    S=sorted(set(S)); best=F(0); bq=None; bval=None; dens={2*a for a in S}
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
                if c>best: best=c; bq=q; bval=num
    return best,bval,bq

def covers_all(A,n):
    for q in range(2,n+1):
        if not any(v%q==0 for v in A): return False
    return True

for n,N in [(8,30),(9,32)]:
    L=F(1,n+1)
    print(f"=== n={n}, searching max element up to {N} (conjecture says max <= {3*n}) ===", flush=True)
    cand=0; tight=[]
    for A in combinations(range(1,N+1),n):
        if not covers_all(A,n): continue
        if reduce(gcd,A)!=1: continue
        cand+=1
        M,val,q=M_full(A)
        if M==L: tight.append((list(A),val,q))
    print(f"  covering+primitive candidates: {cand}", flush=True)
    print(f"  PRIMITIVE TIGHT: {len(tight)}", flush=True)
    mx=0
    for A,val,q in tight:
        mx=max(mx,max(A))
        print(f"     {A}  max={max(A)} max/n={max(A)/n:.2f} val={val} q={q}", flush=True)
    if tight:
        print(f"  => n={n}: MAX over primitive tight = {mx} = {mx/n:.3f}*n ; <= 3n? {mx<=3*n}", flush=True)
    print(flush=True)
