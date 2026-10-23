#!/usr/bin/env python3
"""Exact extra queries after both actual endpoint owners are already shared."""
from pathlib import Path
from functools import lru_cache
from collections import Counter,defaultdict
import json
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(p,s):
 global CHECKS
 CHECKS+=1
 if not p:raise RuntimeError(s)

def costs(paths,values):
 @lru_cache(None)
 def tuples(S):
  if len({values[i] for i in S})==1:return 0
  candidates=[]
  for j in range(1,7):
   parts=defaultdict(list)
   for i in S:parts[paths[i][j]].append(i)
   if len(parts)>1:candidates.append(1+max(tuples(tuple(v)) for v in parts.values()))
  return min(candidates)
 columns=[[sum(1<<i for i,p in enumerate(paths) if p[j]==v) for v in range(8)] for j in range(1,7)]
 colors=[sum(1<<i for i,v in enumerate(values) if v==x) for x in set(values)]
 @lru_cache(None)
 def bits(S):
  if any(S&C==S for C in colors):return 0
  options=[]
  for col in columns:
   parts=[S&C for C in col if S&C]
   if len(parts)>1:options.append(1+max(bits(C) for C in parts))
  return min(options)
 n=len(paths);a=tuples(tuple(range(n)));b=bits((1<<n)-1)
 need(a==b,'independent conditional adaptive recursions')
 fixed=[]
 for mask in range(1<<6):
  answers={};ok=True
  for p,v in zip(paths,values):
   key=tuple(p[j+1] for j in range(6) if mask>>j&1)
   if key in answers and answers[key]!=v:ok=False;break
   answers[key]=v
  if ok:fixed.append(mask.bit_count())
 return dict(adaptive=a,nonadaptive=min(fixed))

records=json.loads((HERE/'path_records.json').read_text());report=[]
for source,rs in records.items():
 fibres=defaultdict(list)
 for r in rs:fibres[(r['path'][0],r['path'][-1])].append(r)
 rows=[];hist={key:Counter() for key in ('phase','target','path')};fixedhist={key:Counter() for key in hist}
 for pins,rows0 in sorted(fibres.items()):
  ps=[r['path'] for r in rows0];answers={}
  for key in hist:
   vs=[tuple(r['path']) if key=='path' else r[key] for r in rows0];answer=costs(ps,vs);answers[key]=answer
   hist[key][answer['adaptive']]+=1;fixedhist[key][answer['nonadaptive']]+=1
  rows.append(dict(endpoints=list(pins),actual_paths=len(ps),costs=answers))
 need(sum(r['actual_paths'] for r in rows)==len(rs),'every actual path in exactly one endpoint fibre')
 need(all(max(c)<=3 for c in hist.values()),'three extra adaptive clocks suffice on each pinned fibre')
 report.append(dict(source=int(source),extra_adaptive_depth_distributions={k:dict(sorted(v.items())) for k,v in hist.items()},
                    extra_nonadaptive_depth_distributions={k:dict(sorted(v.items())) for k,v in fixedhist.items()},fibres=rows))
print(json.dumps(dict(status='PASS',checks=CHECKS,sources=report),sort_keys=True,indent=2))
