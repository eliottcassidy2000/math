#!/usr/bin/env python3
"""
Closing the height>=2 gap at n=12  (mac-mini-2026-07-18-S122)
=============================================================
State of the n=12 gap BEFORE this run:
  - S115: ALL 12-subsets of {1..26} exhausted (covering+primitive) -> only {1..12}.
          That already includes DEEP sets (13, 26 are in range).
  - THM-770: all FULL-RESIDUE packets with lift heights <=12 -> only dilates.
  - THM-1001: single-coordinate winding at ALL heights -> excluded.
  => the genuine remaining gap is max(A) >= 27, i.e. lift height >= 2 on some
     coordinate with the set not already covered above.
This run pushes the exhaustive census from N=26 up to N=34.
Prunes: covering lemma (multiple of every q<=12) + SOUND batched prescreen
(grid max <= true M, so a rejection is never a false negative) + exact Q.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys

n=12; np1=13; L=1.0/np1; Lex=F(1,np1)

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
                if num*np1<q: break
            if num*np1>=q:
                c=F(num,q)
                if c>best: best=c
    return best

N=int(sys.argv[1]) if len(sys.argv)>1 else 34
REQ=list(range(2,13)); full=(1<<len(REQ))-1
mask=[0]*(N+1)
for v in range(1,N+1):
    mask[v]=sum((1<<i) for i,q in enumerate(REQ) if v%q==0)
T=6000; t=np.arange(1,T)/T
Cl=np.empty((N+1,T-1))
for v in range(1,N+1):
    x=(v*t)%1.0; Cl[v]=np.minimum(x,1-x)

print(f"=== n=12 exhaustive census to N={N} (S115 did N=26; gap is max>=27) ===",flush=True)
batch=[]; nb=0; ncand=0; nkeep=0; tight=[]
BS=4000
def flush(batch):
    global nkeep,tight
    if not batch: return
    B=np.array(batch,dtype=np.int64)
    cur=Cl[B[:,0]].copy()
    for j in range(1,n): np.minimum(cur,Cl[B[:,j]],out=cur)
    idx=np.where(cur.max(axis=1)<=L+2e-4)[0]
    nkeep+=len(idx)
    for i in idx:
        A=tuple(int(x) for x in B[i])
        if M_exact(A)==Lex: tight.append(list(A))
for A in combinations(range(1,N+1),n):
    m=0
    for v in A: m|=mask[v]
    if m!=full: continue
    if reduce(gcd,A)!=1: continue
    ncand+=1
    batch.append(A)
    if len(batch)>=BS:
        flush(batch); batch=[]; nb+=1
        if nb%40==0: print(f"   ... {ncand:,} candidates processed",flush=True)
flush(batch)
print(f"  covering+primitive candidates: {ncand:,}",flush=True)
print(f"  survived prescreen: {nkeep:,}",flush=True)
print(f"  PRIMITIVE TIGHT: {len(tight)}",flush=True)
for A in tight: print(f"     {A}  max={max(A)}",flush=True)
spor=[A for A in tight if A!=list(range(1,13))]
print(f"  SPORADIC (non-segment): {len(spor)} {spor}",flush=True)
