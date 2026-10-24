from pathlib import Path
from collections import defaultdict,deque
from itertools import combinations
import json
H=Path(__file__).resolve().parent;P=H.parent
opts={r['source']:r for r in json.loads((H/'inventory.json').read_text())['records']};plus=[1394,1425,1457,1480,2176,2185,2581,2835,2945];minus=[1477,1549,1691,1921,1974,2837,2905,3000,3369];base=plus+minus;outside=sorted(set(opts)-set(base));L={s:(sum(opts[s]['source_label'][j]<<j for j in range(4)),opts[s]['source_label'][4]) for s in opts}
def sig(m,ss):
 n=0;p=0;x=0
 for j,s in enumerate(ss):
  if m>>j&1:n+=1;p^=L[s][0];x+=L[s][1]
 return n,p,x%128
states=[(0,0,0)]
for s in base:states += [(n+1,p^L[s][0],(x+L[s][1])%128) for n,p,x in states]
by=defaultdict(list)
for m,t in enumerate(states):by[t].append(m)
startbase=261632;endbase=511
def fiber(buf,initial):
 target=sig(startbase|(initial<<18),base+list(buf));out=[]
 for b in range(1<<len(buf)):
  n,p,x=sig(b,buf);out.extend(m|(b<<18) for m in by[(target[0]-n,target[1]^p,(target[2]-x)%128)])
 return sorted(out)
def path(states,start,end,budget):
 prev={start:None};q=deque([start])
 while q:
  a=q.popleft()
  if a==end:
   out=[]
   while a is not None:out.append(a);a=prev[a]
   return out[::-1]
  for b in states:
   if b not in prev and (a^b).bit_count()<=budget:prev[b]=a;q.append(b)
 return None
results=[]
for k in [1]:
 for buf in combinations(outside,k):
  for initial in range(1<<k):
   states=fiber(buf,initial);start=startbase|(initial<<18);end=endbase|(initial<<18)
   for budget in range(4,19,2):
    p=path(states,start,end,budget)
    if p:break
   results.append(dict(buffers=buf,initial=initial,states=len(states),budget=budget,path=p,step_sizes=[(a^b).bit_count() for a,b in zip(p,p[1:])]))
 print('k',k,sorted([r for r in results if len(r['buffers'])==k],key=lambda r:(r['budget'],len(r['path'])))[:5],flush=True)
(H/'eighteen37-search-results.json').write_text(json.dumps(dict(base=base,results=results),indent=2)+'\n')
