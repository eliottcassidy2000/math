#!/usr/bin/env python3
"""Independent monomial-orbit enumeration of the bounded rational ansatz."""
from itertools import permutations
from math import comb
from collections import defaultdict
import numpy as np
import json
from verify_period_ranks import rank_forward
PERMS=list(permutations(range(3)))
def comps(n):
    for a in range(n+1):
        for b in range(n-a+1):yield(a,b,n-a-b)
def make_basis(d,k,eps):
    reps=set()
    for aa in comps(d+3):
        for bb in comps(k):
            pairs=tuple(sorted(zip(aa,bb)))
            if any(b==0 and a<2 for a,b in pairs):continue
            if eps==-1 and len(set(pairs))<3:continue
            reps.add(pairs)
    basis=[]
    for pairs in sorted(reps):
        poly={}
        for p in PERMS:
            terms=[pairs[j]for j in p]
            key=tuple(a for a,b in terms)+tuple(b for a,b in terms)
            sign=eps**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))
            if key in poly and poly[key]!=sign:raise RuntimeError('bad orbit')
            poly[key]=sign
        basis.append(poly)
    return basis

def transform_monomial(t,which):
    a,b,c,A,B,C=t
    if which==1:
        for j in range(b+1):
            for J in range(B+1):
                yield (j,a+b-j,c,J,A+B-J,C),comb(b,j)*comb(B,J)*(-1)**(a+b-j+A+B-J)
    else:
        for j in range(a+1):
            for J in range(A+1):
                yield (j+b,a-j,c,J+B,A-J,C),comb(a,j)*comb(A,J)*(-1)**(j+b+J+B)
def shift(t,i):
    tt=list(t);tt[i]+=1;return tuple(tt)
def certificate_column(poly):
    out=defaultdict(int)
    for t,v in poly.items():
        out[shift(t,0)]+=v;out[shift(t,1)]-=v
        for u,c in transform_monomial(t,1):out[shift(u,0)]-=v*c
        for u,c in transform_monomial(t,2):out[shift(u,1)]+=v*c
    return {k:v for k,v in out.items()if v}
results=[];prime=100003
for k,maxd in ((1,30),(2,12),(3,3)):
    for d in range(maxd+1):
        for eps in (1,-1):
            cols=[certificate_column(b)for b in make_basis(d,k,eps)]
            rows=sorted(set().union(*(set(col)for col in cols)))if cols else[]
            ix={t:i for i,t in enumerate(rows)}
            A=np.zeros((len(rows),len(cols)),dtype=np.int64)
            for j,col in enumerate(cols):
                for t,v in col.items():A[ix[t],j]=v%prime
            rank=rank_forward(A,prime);null=len(cols)-rank
            expected=int((d,k,eps)in((12,2,1),(3,3,-1)))
            if null!=expected:raise RuntimeError((d,k,eps,null,expected))
            results.append({'d':d,'auxiliary_degree':k,'symmetry':eps,'columns':len(cols),'rows':len(rows),'rank_mod_100003':rank,'nullity_mod_100003':null})
    print('PASS ansatz auxiliary degree',k,'through original degree',maxd,flush=True)
with open('period_ansatz_receipt.json','w')as f:json.dump({'status':'PASS','cases':results},f,indent=2)
print('ALL ANSATZ CHECKS PASSED',flush=True)
