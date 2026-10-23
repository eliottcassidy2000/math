from math import isqrt
from itertools import combinations
from pathlib import Path
import json
parts=[[17,19,23,31,47],[29,41,59,71],[2,3,5,7,11,13]]
P=sorted(sum(parts,[]));E=[(a,b) for a,b in combinations(P,2) if isqrt(a+b)**2==a+b]
rem=set(P);cs=[]
while rem:
 stack=[rem.pop()];c=[]
 while stack:
  a=stack.pop();c.append(a)
  for u,v in E:
   b=v if u==a else u if v==a else None
   if b in rem:rem.remove(b);stack.append(b)
 cs.append(sorted(c))
cs=sorted(cs,key=lambda c:c[0]);expected=[P]
if cs!=expected or len(E)!=14:raise RuntimeError('graph mismatch')
rows=[{'vertices':p,'edges':[e for e in E if set(e)<=set(p)]} for p in parts]
if [len(r['edges']) for r in rows]!=[2,2,3]:raise RuntimeError('blocks')
leaves=[v for v in P if sum(v in e for e in E)==1]
if leaves!=[3,11,19,31,71]:raise RuntimeError('leaves')
out={'leaves':leaves,'status':'PASS exact integer squares','all_edges':E,'components':cs,'parts':rows,'cycle_rank':len(E)-len(P)+len(cs)}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
