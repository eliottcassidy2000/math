#!/usr/bin/env python3
"""Independent exact cellular rank controls for prefix-pentagon/suffix-square filling."""
import itertools as it,json
import sympy as s
results=[]
for k in range(1,5):
 ws=[w for n in range(k+1) for w in it.permutations(range(k),n)]
 edges=set()
 for w in ws:
  if w:edges.add(tuple(sorted((w,w[:-1]))))
  for j in range(len(w)-1):
   v=w[:j]+(w[j+1],w[j])+w[j+2:];edges.add(tuple(sorted((w,v))))
 edges=sorted(edges); pos={e:i for i,e in enumerate(edges)};rows=[]
 def cycle(c):
  row=[0]*len(edges)
  for u,v in zip(c,c[1:]+c[:1]):
   e=tuple(sorted((u,v)));row[pos[e]]+=1 if e==(u,v) else -1
  rows.append(row)
 for p in ws:
  unused=set(range(k))-set(p)
  for a,b in it.combinations(unused,2):cycle([p,p+(a,),p+(a,b),p+(b,a),p+(b,)])
  for j in range(len(p)-1):
   q=p[:j]+(p[j+1],p[j])+p[j+2:]
   if p<q:
    for c in unused:cycle([p,q,q+(c,),p+(c,)])
 rank=s.Matrix(rows).rank() if rows else 0
 want=len(edges)-len(ws)+1
 if rank!=want:raise RuntimeError((k,rank,want))
 results.append(dict(k=k,vertices=len(ws),edges=len(edges),cells=len(rows),cell_rank=rank,cycle_rank=want))
print(json.dumps({'status':'PASS','controls':results},sort_keys=True,indent=2))
