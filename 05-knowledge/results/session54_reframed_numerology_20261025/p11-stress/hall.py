#!/usr/bin/env python3
import itertools as it, math, collections as co, json, pathlib
ROOT=pathlib.Path(__file__).resolve().parent
def check(p,msg):
 if not p:raise RuntimeError(msg)
def census(n):
 V=range(n);Q={x*x%n for x in range(1,n)};E=[(u,v) for u in V for v in V if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)}
 subsets={k:list(it.combinations(V,k)) for k in range(n+1)}
 rows=[];allmin=1000;families={};tested=0
 for k in range(2,n+1):
  low=1000;num=0;costhist=co.Counter()
  for S in subsets[k]:
   # Edges from S, partitioned by head. T permits heads in T.
   heads=[0]*n
   for i,(u,v) in enumerate(E):
    if u in S:heads[v]|=1<<i
   allmask=0
   for h in heads:allmask|=h
   for T in subsets[k-2]:
    keep=0
    for t in T:keep|=heads[t]
    B=allmask^keep;cost=B.bit_count();tested+=1
    costhist[cost]+=1;low=min(low,cost);allmin=min(allmin,cost)
    if cost==n-1:
     families.setdefault(B,[]).append({'S':list(S),'T':list(T)})
  rows.append({'S_size':k,'minimum_cost':low,'certificate_count_at_minimum':costhist[low],'cost_histogram':dict(sorted(costhist.items()))})
 stars={sum(1<<i for i,(u,v) in enumerate(E) if u==c or v==c):c for c in V}
 check(tested==math.comb(2*n,n-2),'universe')
 check(allmin==n-1,'minimum Hall cost')
 check(not(set(families)&set(stars)),'stars are not Hall')
 # Stars have a residual one-path-cycle factor, proved by the explicit Hamilton cycle on V\{c}.
 def cycle_without(c):
  W=[v for v in V if v!=c];start=W[0]
  def dfs(path,left):
   if not left:return path if (path[-1],start) in EI else None
   for w in sorted(left):
    if (path[-1],w) in EI:
     out=dfs(path+[w],left-{w})
     if out:return out
   return None
  return dfs([start],set(W[1:]))
 cycles={c:cycle_without(c) for c in V}
 check(all(cycles.values()),'star factor witness')
 return {'n':n,'edges':E,'tested_certificates':tested,'minimum_Hall_cost':allmin,'distinct_cost_n_minus_one_Hall_blockers':len(families),'stars':len(stars),'rows':rows,'blockers':[{'mask':b,'arcs':[E[i] for i in range(len(E)) if b>>i&1],'certificates':certs} for b,certs in sorted(families.items())],'star_factor_cycles':cycles}
print(json.dumps({'status':'PASS','scope':'Hall certificates exhaustive; full HP blockers NOT classified','censuses':[census(7),census(11)]},sort_keys=True,indent=2))
