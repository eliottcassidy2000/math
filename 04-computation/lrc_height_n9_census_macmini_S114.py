#!/usr/bin/env python3
"""n=9 height census (pin the threshold). Covering prune: multiples of 5,6,7,8,9 required."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
n=9; N=28; L=F(1,n+1)
REQ=[5,6,7,8,9]   # 2,3,4 implied by 6,8,9
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
print(f"=== n={n}, max element <= {N} (conjecture max <= 3n = {3*n}) ===",flush=True)
cand=0; tight=[]
for A in combinations(range(1,N+1),n):
    ok=True
    for q in REQ:
        f=False
        for v in A:
            if v%q==0: f=True; break
        if not f: ok=False; break
    if not ok: continue
    if reduce(gcd,A)!=1: continue
    cand+=1
    M,val,q=M_full(A)
    if M==L: tight.append((list(A),val,q))
print(f"  covering+primitive candidates: {cand}",flush=True)
print(f"  PRIMITIVE TIGHT: {len(tight)}",flush=True)
mx=0
for A,val,q in tight:
    mx=max(mx,max(A)); print(f"     {A} max={max(A)} max/n={max(A)/n:.2f} val={val} q={q}",flush=True)
if tight: print(f"  => n={n}: MAX = {mx} = {mx/n:.3f}*n ; <= 3n? {mx<=3*n}",flush=True)
