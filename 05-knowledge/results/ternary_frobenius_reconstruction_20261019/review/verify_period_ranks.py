#!/usr/bin/env python3
"""Separate dense forward-elimination audit; imports no producer code."""
from math import comb
from functools import lru_cache
from collections import defaultdict
import numpy as np
import json

@lru_cache(None)
def monomials(n):return [(i,j,n-i-j)for i in range(n+1)for j in range(n-i+1)]
def orbits(n,kind):
    seen=set();result=[]
    for t in monomials(n):
        if kind=='interior' and 0 in t:continue
        if kind=='supported' and 0 not in t:continue
        if t in seen:continue
        orbit={t,t[1:]+t[:1],t[2:]+t[:2]};seen.update(orbit);result.append(sorted(orbit))
    return result
@lru_cache(None)
def integer_columns(n,kind):
    result=[]
    for orbit in orbits(n,kind):
        d=defaultdict(int)
        for a,b,c in orbit:
            d[a,b,c]+=1
            for j in range(b+1):d[j,a+b-j,c]+=(-1)**(a+b-j)*comb(b,j)
            for j in range(a+1):d[j+b,a-j,c]+=(-1)**(b+j)*comb(a,j)
        result.append({t:v for t,v in d.items()if v})
    return result

def matrix_mod(n,kind,p):
    cols=integer_columns(n,kind)
    rows=sorted(set().union(*(set(c)for c in cols))) if cols else []
    index={t:j for j,t in enumerate(rows)}
    A=np.zeros((len(rows),len(cols)),dtype=np.int64)
    for j,col in enumerate(cols):
        for t,v in col.items():A[index[t],j]=v%p
    return A

def rank_forward(A,p):
    # Column-pivot forward elimination; never uses floating arithmetic.
    A=A.copy();height,width=A.shape;r=0
    for c in range(width):
        possible=np.nonzero(A[r:,c])[0]
        if len(possible)==0:continue
        pivot=r+int(possible[0])
        if pivot!=r:A[[r,pivot],c:]=A[[pivot,r],c:]
        inv=pow(int(A[r,c]),-1,p)
        A[r,c:]=A[r,c:]*inv%p
        for start in range(r+1,height,64):
            end=min(start+64,height)
            multipliers=A[start:end,c].copy()
            A[start:end,c:]=(A[start:end,c:]-multipliers[:,None]*A[r,c:])%p
        r+=1
        if r==height:break
    return r

if __name__=='__main__':
    thresholds={2:11,3:11,5:18,7:24,11:34,13:38,17:46,19:50}
    result={'status':'PASS','thresholds':[],'supported':[]}
    for p,delta in thresholds.items():
        vals=[]
        for n in range(delta+1):
            A=matrix_mod(n,'interior',p);rank=rank_forward(A,p);null=A.shape[1]-rank
            if null!=int(n==delta):raise RuntimeError(('threshold',p,n,A.shape,rank,null))
            vals.append({'degree':n,'columns':A.shape[1],'rank':rank,'nullity':null})
        result['thresholds'].append({'p':p,'minimum':delta,'degrees':vals})
        print('PASS complete minimum',p,delta,flush=True)
    for p in (5,7,11,13):
        vals=[]
        for n in list(range(1,40,2))+[49,51,121,125,127]:
            A=matrix_mod(n,'supported',p);rank=rank_forward(A,p);null=A.shape[1]-rank
            expected=0;power=p
            while power<=n:expected+=1;power*=p
            if null!=expected:raise RuntimeError(('supported',p,n,null,expected))
            vals.append({'degree':n,'columns':A.shape[1],'rank':rank,'nullity':null})
        result['supported'].append({'p':p,'degrees':vals})
        print('PASS supported classification',p,flush=True)
    with open('period_ranks_receipt.json','w')as f:json.dump(result,f,indent=2)
    print('ALL RANK CHECKS PASSED',flush=True)
