#!/usr/bin/env python3
"""Second exact algorithm: tuple-state Bellman recursion and BFS thresholds.
Does not import check.py. Every check survives python -O.
"""
from collections import Counter,defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import permutations,combinations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
count=0

def verify(ok,why):
 global count
 count+=1
 if not ok:raise ValueError(why)

def main():
 receipt=json.loads((HERE/'receipt.json').read_text()); allrecords=json.loads((HERE/'path_records.json').read_text()); summary=[]
 for row in receipt['sources']:
  recs=allrecords[str(row['source'])];A=row['adjacency'];N=len(recs)
  paths=[tuple(r['path']) for r in recs]
  actual=[]
  # A disjoint next-vertex recurrence, rather than checking n! words.
  frontier=[(v,) for v in range(8)]
  while frontier:
   p=frontier.pop()
   if len(p)==8:actual.append(p);continue
   frontier.extend(p+(v,) for v in range(8) if v not in p and A[p[-1]][v])
  verify(sorted(actual)==paths and N==row['H'],'literal complete source paths')
  results={}
  for label in ('phase','target','path'):
   colors=[tuple(r['path']) if label=='path' else r[label] for r in recs]
   @lru_cache(None)
   def solve(S):
    if len({colors[i] for i in S})==1:return 0,0,0
    Ds=[];Es=[]
    for j in range(8):
     blocks=defaultdict(list)
     for i in S:blocks[paths[i][j]].append(i)
     if len(blocks)==1:continue
     answers=[solve(tuple(v)) for v in blocks.values()]
     Ds.append(1+max(v[0] for v in answers));Es.append((len(S)+sum(v[1] for v in answers),1+max(v[2] for v in answers)))
    E=min(Es);return min(Ds),E[0],E[1]
   @lru_cache(None)
   def constrained(S,d):
    if len({colors[i] for i in S})==1:return 0
    if d==0:return None
    choices=[]
    for j in range(8):
     blocks=defaultdict(list)
     for i in S:blocks[paths[i][j]].append(i)
     if len(blocks)==1:continue
     children=[constrained(tuple(v),d-1) for v in blocks.values()]
     if all(v is not None for v in children):choices.append(len(S)+sum(children))
    return min(choices) if choices else None
   indices=tuple(range(N));D,E,Ed=solve(indices);C=constrained(indices,D)
   expected=row['objectives'][label]
   verify((D,str(Fraction(E,N)),Ed,str(Fraction(C,N)))==(expected['worst_depth'],expected['optimal_mean'],expected['mean_policy_worst_depth'],expected['best_mean_at_optimal_worst_depth']),'independent adaptive exact optimum')
   verify(constrained(indices,D-1) is None,'independent explicit depth-infeasibility recurrence')
   valid=[]
   for mask in range(1<<8):
    observation={};ok=True
    for p,v in zip(paths,colors):
     key=tuple(p[j] for j in range(8) if mask>>j&1)
     if key in observation and observation[key]!=v:ok=False;break
     observation[key]=v
    if ok:valid.append(mask)
   nonadaptive=min(mask.bit_count() for mask in valid)
   verify(nonadaptive==6,'complete nonadaptive clock search')
   results[label]=dict(adaptive=D,nonadaptive=nonadaptive,mean=str(Fraction(E,N)),depth_capped_mean=str(Fraction(C,N)),states=solve.cache_info().currsize)
  fibres=defaultdict(list)
  for i,r in enumerate(recs):fibres[(r['path'][0],r['path'][-1],r['phase'])].append(i)
  hist=Counter();maxfibre=0
  for key,ids in fibres.items():
   for b in range(7):
    unseen=set(ids);components=[]
    while unseen:
     start=min(unseen);seen={start};unseen.remove(start);todo=[start]
     while todo:
      u=todo.pop();adj={v for v in unseen if sum(a!=b for a,b in zip(paths[u],paths[v]))<=b}
      unseen-=adj;seen|=adj;todo.extend(adj)
     components.append(seen)
    if len(components)==1:break
   hist[b]+=1;maxfibre=max(maxfibre,len(ids))
  verify(dict(hist)=={int(k):v for k,v in row['reconfiguration']['budget_histogram'].items()},'all exact reconfiguration threshold components')
  summary.append(dict(source=row['source'],objectives=results,budgets=dict(sorted(hist.items())),maximum_fibre_size=maxfibre))
 result=dict(status='PASS',checks=count,sources=summary)
 print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
