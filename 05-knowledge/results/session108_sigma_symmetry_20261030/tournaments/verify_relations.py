#!/usr/bin/env python3
"""Independently enumerate tiny binary relation spaces and marked variants."""
from itertools import combinations
import json

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def spaces(n):
    return [frozenset(s) for k in [1,2,4,8] if k<=2**n for rest in combinations(range(1,2**n),k-1) for s in [(0,)+rest] if all(x^y in s for x in s for y in s)]
def dim(W):return len(W).bit_length()-1
def rotate(x,n):return ((x<<1)&((1<<n)-1))|(x>>(n-1))
def action(item,n,marked):
    if marked:
        W,w=item;return frozenset(rotate(x,n) for x in W),rotate(w,n)
    return frozenset(rotate(x,n) for x in item)
def orbits(items,n,marked=False):
    todo=set(items);result=[]
    while todo:
        x=next(iter(todo));orb={x};y=action(x,n,marked)
        while y!=x:orb.add(y);y=action(y,n,marked)
        require(orb<=set(items),'closed action');todo-=orb;result.append(orb)
    return len(result)
out={}
for n in (1,2,3):
    S=spaces(n)
    GL=[W for W in S if dim(W)<=3]
    AGL=[W for W in S if (2**n-1) in W and dim(W)<=4]
    H=[(W,w) for W in S for w in W if dim(W)<=3 and (w!=0 or dim(W)<=2)]
    out[n]={'GL_transitive':len(GL),'AGL_transitive':len(AGL),'marked_transitive':len(H)}
    if n==3:
        out[n].update(GL_cyclic=orbits(GL,n),AGL_cyclic=orbits(AGL,n),marked_cyclic=orbits(H,n,True))
require(out[2]=={'GL_transitive':5,'AGL_transitive':2,'marked_transitive':11},'degree two')
require(out[3]=={'GL_transitive':16,'AGL_transitive':5,'marked_transitive':50,'GL_cyclic':8,'AGL_cyclic':3,'marked_cyclic':20},'degree three')
print(json.dumps(out,indent=2,sort_keys=True))
