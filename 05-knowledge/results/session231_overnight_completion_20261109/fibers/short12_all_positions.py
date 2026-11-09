import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
from collections import defaultdict
import json
from core import compose
root=Path(__file__).parent;w=json.loads((root/'SHORT12_UNIT_LOOP_WITNESS.json').read_text());target=defaultdict(set)
for n in range(3,10):
 for r in map(json.loads,(root/f'SATURATED{n}.jsonl').read_text().splitlines()):
  if r[:2]==[45,1440]:target[r[2],r[3]].add(tuple(r[:7]))
one=[(1,2,b,b,0,b,b) for b in (0,1)]
def fold(ps):
 s=ps[0]
 for p in ps[1:]:s=compose(s,p)
 return s
rows=[]
for k in range(1,13):
 X=list(map(int,w['X']))*k;Y=list(map(int,w['Y']))*k;L=len(X);pre=[fold([one[1],one[0]])];suf=[None]*(L+1);suf[L]=one[1]
 for x in X:pre.append(compose(pre[-1],one[x]))
 for i in range(L-1,-1,-1):suf[i]=compose(one[Y[i]],suf[i+1])
 a,g=w['core'][2:4];vals={};positions=set();sig=None
 for i in range(L+1):
  for s in target.get((a,g),[]):
   z=compose(compose(pre[i],s),suf[i]);sigi=z[:4]+z[5:]
   if sig is not None and sig!=sigi:raise ValueError('interface')
   sig=sigi;vals[z[4]]=(i,s);positions.add(i)
  if i<L:a,g=2*a+45*Y[i]-1440*X[i],(g+1440*Y[i]-45*X[i])//2
 vs=sorted(vals);best=1
 for t in vs:
  if t-2 in vals:continue
  length=1
  while t+2*length in vals:length+=1
  best=max(best,length)
 rows.append(dict(k=k,physical_positions=sorted(positions),fiber=len(vs),max_run=best,offsets=[v-vs[0] for v in vs]))
print(json.dumps(rows,indent=2))
