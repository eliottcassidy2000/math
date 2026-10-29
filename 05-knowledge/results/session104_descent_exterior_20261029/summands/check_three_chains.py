#!/usr/bin/env python3
from math import isqrt
import json
checks=0
def need(ok, msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)
def graph(n):
    return {x:{y for y in range(1,n+1) if y!=x and isqrt(x+y)**2==x+y} for x in range(1,n+1)}
def components(g):
    unseen=set(g); out=[]
    while unseen:
        stack=[min(unseen)]; c=set()
        while stack:
            x=stack.pop()
            if x in c: continue
            c.add(x); stack.extend(g[x]-c)
        unseen-=c; out.append(c)
    return out
parts=[{1,3,6,8,10},{2,7,9},{4,5,11,12}]
for n in range(4,13):
    g=graph(n)
    need(components(g)==[s&set(range(1,n+1)) for s in parts], 'three inherited components')
    need(max(map(len,g.values()))<=2,'paths degree')
    need(sum(map(len,g.values()))//2==n-3,'linear forest edge count')
need(len(components(graph(13)))==2,'first merge')
need(graph(13)[13]=={3,12},'first new neighbours')
need(len(components(graph(14)))==1,'second merge')
need(graph(14)[14]=={2,11},'second new neighbours')
founders=[]
for n in range(1,301):
    new={s*s-n for s in range(isqrt(n)+1,isqrt(2*n-1)+1)}
    if not new: founders.append(n)
    if n>=6: need(bool(new),'every later vertex attaches')
    need(all(1<=x<n for x in new),'strict earlier neighbours')
need(founders==[1,2,4],'founders')
for n in (14,15,32,100,300): need(len(components(graph(n)))==1,'connected controls')
for k in range(1,101):
    m=2*k*k
    need(isqrt(2*m)**2==2*m,'even-balance start is square-diagonal source')
need(graph(32)[32]=={4,17},'32 has two parents despite diagonal deletion')
for mask in range(1<<10):
    T={i+2 for i in range(10) if mask>>i&1}
    founder=lambda n:not(T&set(range(n+1,2*n)))
    need(not(founder(4) and founder(6)) or founder(5),'founder obstruction')
print(json.dumps({'status':'PASS','require_calls':checks,'components_12':[sorted(s) for s in parts],'component_counts_12_13_14':[3,2,1],'founders_through_300':founders,'positive_provenance_recovery':True},sort_keys=True,indent=2))
