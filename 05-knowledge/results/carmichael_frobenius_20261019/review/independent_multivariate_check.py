#!/usr/bin/env python3
"""Independent finite normal-form audit of the multivariate axis clock."""
from itertools import product
from pathlib import Path
from math import prod
import json

checks=0
def require(x,why):
    global checks
    checks+=1
    if not x: raise RuntimeError(str(why))

def val(n,p):
    r=0
    while n%p==0: r+=1; n//=p
    return r

def test_ring(p,e,weights):
    terms=sorted(weights,key=lambda g:(sum(g),g)); mods=[p**weights[g] for g in terms]
    index={g:i for i,g in enumerate(terms)}; dims=len(terms[0]); zero=(0,)*dims
    require(terms[0]==zero,('constant first',terms))
    triples=[]
    for i,g in enumerate(terms):
        for j,h in enumerate(terms):
            gh=tuple(a+b for a,b in zip(g,h))
            if gh in index: triples.append((i,j,index[gh]))
    def mul(x,y):
        out=[0]*len(terms)
        for i,j,k in triples: out[k]+=x[i]*y[j]
        return tuple(z%m for z,m in zip(out,mods))
    def power(x,n):
        out=(1,)+(0,)*(len(terms)-1)
        while n:
            if n&1: out=mul(out,x)
            x=mul(x,x); n//=2
        return out
    s=val(e,p); H=0
    for i in range(dims):
        h=0
        while True:
            g=tuple(p**h if j==i else 0 for j in range(dims))
            if g not in weights: break
            H=max(H,h+weights[g]); h+=1
    R=(max(weights[zero]-1,H)+s-1)//s
    roots={x%p:x for x in range(mods[0]) if pow(x,e,mods[0])==x}
    require(len(roots)==p,('scalar roots',p,e,weights))
    observed=0; witness=None
    for x in product(*(range(m) for m in mods)):
        goal=(roots[x[0]%p],)+(0,)*(len(terms)-1)
        y=x; r=0
        while y!=goal and r<=R: y=power(y,e); r+=1
        require(y==goal and r<=R,('upper',p,e,weights,x,R,r))
        if r>observed: observed=r; witness=x
    require(observed==R,('sharp',p,e,weights,R,observed))
    return {'p':p,'e':e,'weights':{str(g):w for g,w in weights.items()},
            'R':R,'observed':observed,'ring_size':prod(mods),'witness':witness}

rows=[]
for a in range(1,4):
    for ax,ay in product(range(1,a+1),repeat=2):
        for axy in range(0,min(ax,ay)+1):
            if a+ax+ay+axy>7: continue
            ws={(0,0):a,(1,0):ax,(0,1):ay}
            if axy: ws[1,1]=axy
            for e in (3,9,21): rows.append(test_ring(3,e,ws))
for bs,e in [((4,2),3),((4,2),21),((3,3),3),((2,2,2),3)]:
    ws={g:1 for g in product(*(range(b) for b in bs))}
    rows.append(test_ring(3,e,ws))
# A torsion-sensitive mixed ring with a live cubic axis and nonlinear exponent.
rows.append(test_ring(3,21,{(0,0):2,(1,0):1,(2,0):1,(3,0):1,(0,1):1,(1,1):1}))

out={'status':'PASS','checks':checks,'rings':len(rows),'cases':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','checks','rings')}))
