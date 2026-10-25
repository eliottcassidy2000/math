#!/usr/bin/env python3
from itertools import combinations, product
from collections import Counter
import json

def need(x,m):
    if not x: raise RuntimeError(m)

def analyse(n,es):
    deg=[0]*n
    for u,v in es: deg[u]+=1;deg[v]+=1
    scores=Counter(); forests=0
    for bits in product((0,1),repeat=len(es)):
        d=[0]*n; blocks=[{v} for v in range(n)]; acyclic=True
        for bit,(u,v) in zip(bits,es):
            d[v if bit else u]+=1
            if bit and acyclic:
                a=next(b for b in blocks if u in b);b=next(b for b in blocks if v in b)
                if a is b: acyclic=False
                else: blocks.remove(a);blocks.remove(b);blocks.append(a|b)
        scores[tuple(d)]+=1
        forests+=acyclic
    need(forests==len(scores),'forest-score mismatch')
    even=all(d%2==0 for d in deg)
    need((forests%2==1)==even,'Eulerian parity failure')
    fixed=0; mass=0;pairs=0
    for d,w in scores.items():
        rev=tuple(a-b for a,b in zip(deg,d))
        need(scores[rev]==w,'complement multiplicity')
        if rev==d: fixed+=1
        elif d<rev: pairs+=1;mass+=w
    need(fixed==int(even),'central score count')
    if not fixed:
        need(2*pairs==forests,'matching count')
        need(2*mass==2**len(es),'edge mass')
    return forests,pairs,mass,fixed
count=0
for n in range(1,6):
    complete=list(combinations(range(n),2))
    for bits in product((0,1),repeat=len(complete)):
        analyse(n,[e for e,b in zip(complete,bits) if b]);count+=1
paths=[]
for n in range(3,8):
    result=analyse(n,[(u,v) for u in range(n) for v in range(u+2,n)])
    need(result[3]==0,'path central fixed score')
    paths.append({'n':n,'fibres':result[0],'coarse_pairs':result[1],'pair_edge_mass':result[2]})
try: need(False,'sentinel')
except RuntimeError: pass
else: raise RuntimeError('negative control disabled')
print(json.dumps({'status':'PASS','all_simple_graphs_through_n5':count,'path_complements':paths,'exception_safe_control':'PASS'},indent=2))
