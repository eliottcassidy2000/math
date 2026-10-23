#!/usr/bin/env python3
import itertools,collections,json
from pathlib import Path

def arcs(n,m):return {(a,b) if m>>k&1 else (b,a) for k,(a,b) in enumerate(itertools.combinations(range(n),2))}
def paths(n,A):return [p for p in itertools.permutations(range(n)) if all((a,b) in A for a,b in zip(p,p[1:]))]
controls=0
for n in range(1,5):
 for m in range(1<<(n*(n-1)//2)):
  c=collections.Counter((p[0],p[-1]) for p in paths(n,arcs(n,m)))
  if max(c.values(),default=0)>1:raise RuntimeError('smaller endpoint ambiguity')
  controls+=1
A=arcs(5,4);P=(4,0,3,2,1);Q=(4,2,0,3,1);result=[]
for p in(P,Q):
 E=set(zip(p,p[1:]))
 if not E<=A:raise RuntimeError('not actual path')
 B={e if e in E else tuple(reversed(e)) for e in A}
 m=sum(1<<k for k,e in enumerate(itertools.combinations(range(5),2)) if e in B)
 result.append({'path':p,'target_mask':m,'H':len(paths(5,B)),'scores':sorted(sum(a==v for a,b in B) for v in range(5))})
if [x['H'] for x in result]!=[15,11] or len(paths(5,A))!=5:raise RuntimeError('counts')
if result[0]['scores']!=result[1]['scores']:raise RuntimeError('score match')
out={'status':'PASS','smaller_tournaments':controls,'source_mask':4,'source_H':5,'targets':result}
Path(__file__).with_name('source-endpoint-hostile-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
