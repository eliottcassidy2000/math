from pathlib import Path
from collections import defaultdict,deque
from itertools import combinations, product
import json
H=Path(__file__).resolve().parent;P=H.parent
pool=json.loads((H/'inventory.json').read_text())['records'];opts={r['source']:r for r in pool}
c=json.loads((H/'circuit16.json').read_text());base=c['plus']+c['minus'];outside=[r['source'] for r in pool if r['source'] not in base]
labels={s:(sum(opts[s]['source_label'][j]<<j for j in range(4)),opts[s]['source_label'][4]) for s in opts}
def sig(mask,sources):
 n=0;p=0;x=0
 for j,s in enumerate(sources):
  if mask>>j&1:n+=1;p^=labels[s][0];x+=labels[s][1]
 return n,p,x%128
by=defaultdict(list)
for mask in range(1<<16):by[sig(mask,base)].append(mask)
startbase=0xff00;endbase=0xff

def fiber(buffers,initial):
 target=sig(startbase|(initial<<16),base+list(buffers));states=[]
 for b in range(1<<len(buffers)):
  n,p,x=sig(b,buffers);key=(target[0]-n,target[1]^p,(target[2]-x)%128)
  states.extend(a|(b<<16) for a in by.get(key,[]))
 return sorted(states)
def route(states,start,end,threshold):
 prev={start:None};q=deque([start])
 while q:
  a=q.popleft()
  if a==end:
   path=[]
   while a is not None:path.append(a);a=prev[a]
   return path[::-1]
  for b in states:
   if b not in prev and (a^b).bit_count()<=threshold:prev[b]=a;q.append(b)
 return None
out=[]
for k in range(1,4):
 for buffers in combinations(outside,k):
  for initial in range(1<<k):
   states=fiber(buffers,initial);start=startbase|(initial<<16);end=endbase|(initial<<16)
   for budget in range(4,17,2):
    path=route(states,start,end,budget)
    if path:
     row=dict(buffers=buffers,costs=[opts[s]['cost'] for s in buffers],initial=initial,states=len(states),budget=budget,path=path,step_sizes=[(a^b).bit_count() for a,b in zip(path,path[1:])]);out.append(row);break
 print('k',k,'best',min(x['budget'] for x in out if len(x['buffers'])==k),'trials',sum(len(x['buffers'])==k for x in out),flush=True)
 print(sorted([x for x in out if len(x['buffers'])==k],key=lambda r:(r['budget'],max(r['costs']),len(r['path'])))[:5],flush=True)
(H/'search-results.json').write_text(json.dumps(out,indent=2)+'\n')
