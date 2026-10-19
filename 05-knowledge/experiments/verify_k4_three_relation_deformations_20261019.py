#!/usr/bin/env python3
"""Bounded K4 deformation check with two closure algorithms and exact witnesses."""
from itertools import combinations,permutations
from collections import Counter
import hashlib,json
BASE=(0,0,0,0,0,0,44,35,26,21,54,57);N=12

def require(x,msg):
 if not x:raise RuntimeError(msg)

def closure(arcs):
 p=list(BASE)
 for a,b in arcs:p[b]|=1<<a
 for k in range(N):
  for j in range(N):
   if p[j]>>k&1:p[j]|=p[k]
 if any(p[j]>>j&1 for j in range(N)):return None
 return tuple(p)

def independent_closure(arcs):
 edges={(a,b) for b in range(N) for a in range(N) if BASE[b]>>a&1}|set(arcs)
 out=[]
 for target in range(N):
  reached=set();todo=[target]
  while todo:
   v=todo.pop()
   for a,b in edges:
    if b==v and a not in reached:reached.add(a);todo.append(a)
  if target in reached:return None
  out.append(sum(1<<a for a in reached))
 return tuple(out)

def chain(mask,d):
 return all(d[a]>>b&1 or d[b]>>a&1 for a,b in combinations([v for v in range(N) if mask>>v&1],2))

def witness(p,u):
 for dual,(d,v) in enumerate(((p,u),(u,p))):
  for a,b in combinations(range(N),2):
   if not (d[a]|v[a])>>b&1 and d[a]==d[b] and chain(v[a]&~v[b],d) and chain(v[b]&~v[a],d):return [dual,a,b]
 return None

def main():
 arcs=list(permutations(range(6),2))+list(permutations(range(6,12),2))
 seen=set();hist=Counter();records=[];by_k=[]
 for k in range(1,4):
  for aa in combinations(arcs,k):
   hist['patterns']+=1;p=closure(aa)
   if p is None:hist['cyclic_patterns']+=1;continue
   if p in seen:hist['duplicate_closures']+=1;continue
   seen.add(p);u=[sum(1<<j for j in range(N) if p[j]>>i&1) for i in range(N)]
   deg=[N-1-(p[i]|u[i]).bit_count() for i in range(N)]
   if max(deg)>7:hist['distinct_not_seven_thin']+=1;continue
   require(independent_closure(aa)==p,'Independent closure differs')
   w=witness(p,u);require(w is not None,'Seven-thin deformation lacks very good pair')
   hist['distinct_seven_thin']+=1;hist['very_good_witnesses']+=1
   records.append({'added_arcs':aa,'pred':p,'degrees':deg,'very_good_witness':w})
  by_k.append({'max_additions':k,'cumulative_counts':dict(hist)})
 require(dict(hist)=={'patterns':36050,'distinct_not_seven_thin':33336,'cyclic_patterns':1850,'duplicate_closures':240,'distinct_seven_thin':624,'very_good_witnesses':624},'Count regression')
 digest=hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest()
 print(json.dumps({'status':'ALL_BOUNDED_DEFORMATIONS_HAVE_VERY_GOOD_PAIRS','scope':'All subsets of one, two or three directed arcs within the original levels; no cross-level additions.','counts':dict(hist),'by_max_additions':by_k,'record_digest':digest,'records':records},indent=2))
if __name__=='__main__':main()
