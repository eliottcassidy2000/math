from itertools import permutations,combinations
from pathlib import Path
import json,math
ROOT=Path(__file__).resolve().parent.parent;OUT=Path(__file__).parent

def count(out):
 n=len(out)
 h=sum(all(out[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
 cycles={}
 for l in range(3,n+1,2):
  raw=sum(all(out[p[i]]>>p[(i+1)%l]&1 for i in range(l)) for p in permutations(range(n),l))
  if raw%l:raise RuntimeError('rotation count')
  cycles[l]=raw//l
 return h,cycles
cases=[]
for row in json.loads((ROOT/'score_first_failure.json').read_text()):
 out=row['out_masks'];h,c=count(out);score=[x.bit_count() for x in out]
 if (h,c[3],c[5],score)!=(row['H'],row['c3'],row['c5'],row['scores']):raise RuntimeError('first failure')
 cases.append(dict(out_masks=out,H=h,cycles=c,scores=score))
checks=0;classes={}
for n in range(1,6):
 table={}
 for mask in range(1<<(n*(n-1)//2)):
  out=[0]*n
  for j,(u,v) in enumerate(combinations(range(n),2)):
   if mask>>j&1:out[u]|=1<<v
   else:out[v]|=1<<u
  h,c=count(out);scores=tuple(sorted(x.bit_count() for x in out));c3=math.comb(n,3)-sum(math.comb(s,2) for s in scores)
  if c.get(3,0)!=c3 or h!=1+2*sum(c.values()):raise RuntimeError('small-order identity')
  table.setdefault(scores,set()).add(h);checks+=1
 if n<=4 and any(len(v)!=1 for v in table.values()):raise RuntimeError('premature failure')
 classes[n]=len(table)
r={'verdict':'PASS','tournaments_exhaustively_checked':checks,'score_multiset_classes':classes,'first_failure_order':5,'witnesses':cases,'cycle_counter':'all ordered distinct tuples divided by length, independently of min-vertex canonicalization'}
print(json.dumps(r,indent=2));(OUT/'root_addendum_independent.json').write_text(json.dumps(r,indent=2)+'\n')
