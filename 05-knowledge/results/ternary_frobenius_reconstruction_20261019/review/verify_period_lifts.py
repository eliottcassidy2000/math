#!/usr/bin/env python3
"""Independent cyclic lifting ranks and supported integral Smith controls."""
from verify_period_ranks import integer_columns
from math import comb
from collections import defaultdict
import numpy as np
import json
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form

def augmented_rank(A,p):
    A=A.copy();height,width=A.shape;r=0;last_pivot=False
    for c in range(width):
        possible=np.nonzero(A[r:,c])[0]
        if len(possible)==0:continue
        pivot=r+int(possible[0])
        if pivot!=r:A[[r,pivot],c:]=A[[pivot,r],c:]
        A[r,c:]=A[r,c:]*pow(int(A[r,c]),-1,p)%p
        for start in range(r+1,height,64):
            end=min(start+64,height);multipliers=A[start:end,c].copy()
            A[start:end,c:]=(A[start:end,c:]-multipliers[:,None]*A[r,c:])%p
        if c==width-1:last_pivot=True
        r+=1
        if r==height:break
    return r,last_pivot
out={'status':'PASS','full_cyclic_lifts':[],'supported_smith':[]}
for p in (5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61):
    cols=integer_columns(p,'full')
    rhs={(j,p-j,0):((-1)**(p-j)*comb(p,j)//p)%p for j in range(1,p)}
    rows=sorted(set(rhs).union(*(set(c)for c in cols)));ix={t:i for i,t in enumerate(rows)}
    A=np.zeros((len(rows),len(cols)+1),dtype=np.int64)
    for j,col in enumerate(cols):
        for t,v in col.items():A[ix[t],j]=v%p
    for t,v in rhs.items():A[ix[t],-1]=v
    rank,increased=augmented_rank(A,p)
    if not increased:raise RuntimeError(('unexpected_lift',p))
    nullity=len(cols)-(rank-1)
    if nullity!=(2 if p==59 else 1):raise RuntimeError(('nullity',p,nullity))
    out['full_cyclic_lifts'].append({'p':p,'columns':len(cols),'norm_rank':rank-1,'augmented_rank':rank,'W_dimension':nullity})
    print('PASS unrestricted cyclic nonlift',p,'dimension',nullity,flush=True)
for n in list(range(1,32,2))+[49]:
    cols=integer_columns(n,'supported')
    rows=sorted(set().union(*(set(c)for c in cols)));ix={t:i for i,t in enumerate(rows)}
    A=Matrix.zeros(len(rows),len(cols))
    for j,col in enumerate(cols):
        for t,v in col.items():A[ix[t],j]=v
    diagonal=smith_normal_form(A,domain=ZZ)
    entries=[abs(int(diagonal[j,j]))for j in range(len(cols))]
    if not all(entries):raise RuntimeError(('rational_kernel',n))
    for p in (5,7,11,13):
        vals=[]
        for x in entries:
            v=0
            while x%p==0:x//=p;v+=1
            vals.append(v)
        expected=0;power=p
        while power<=n:expected+=1;power*=p
        if max(vals,default=0)>1 or sum(v==1 for v in vals)!=expected:
            raise RuntimeError(('smith_lift',n,p,vals,expected))
    out['supported_smith'].append({'degree':n,'invariant_factors':entries})
print('PASS supported Smith controls',flush=True)
with open('period_lifts_receipt.json','w')as f:json.dump(out,f,indent=2)
