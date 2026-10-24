from pathlib import Path
from collections import defaultdict,deque,Counter
from itertools import combinations
import json,time
H=Path(__file__).resolve().parent;P=H.parent
rows=json.loads((H/'inventory.json').read_text())['records'];opts={r['source']:r for r in rows}
c=json.loads((H/'circuit16.json').read_text());base=c['plus']+c['minus'];outside=[r['source'] for r in rows if r['source'] not in base];src=base+outside
L={r['source']:(sum(r['source_label'][j]<<j for j in range(4)),r['source_label'][4]) for r in rows};d=defaultdict(list)
for i,j in combinations(range(28),2):d[(L[src[i]][0]^L[src[j]][0],(L[src[i]][1]+L[src[j]][1])%128)].append((1<<i)|(1<<j))
rels=[(a,b,a|b) for pairs in d.values() for a,b in combinations(pairs,2)]
def route(start,end,relations):
 prev={start:None};q=deque([start])
 while q:
  u=q.popleft()
  if u==end:
   path=[]
   while u is not None:path.append(u);u=prev[u]
   return path[::-1],len(prev)
  for a,b,m in relations:
   t=u&m
   if t==a or t==b:
    v=u^m
    if v not in prev:prev[v]=u;q.append(v)
 return None,len(prev)
results=[];stats=[];found=False
for k in range(13):
 stat=Counter(k=k)
 for ix in combinations(range(12),k):
  allowed=65535|sum(1<<(16+j) for j in ix);rr=[r for r in rels if not r[2]&~allowed];stat['reservoirs']+=1
  basis={}
  for a,b,v in rr:
   while v:
    bit=v.bit_length()-1
    if bit not in basis:basis[bit]=v;break
    v^=basis[bit]
  rem=65535
  for bit in sorted(basis,reverse=True):
   if rem>>bit&1:rem^=basis[bit]
  if rem:stat['span_excluded_reservoirs']+=1;continue
  for init in range(1<<k):
   bm=sum(((init>>h)&1)<<(16+j) for h,j in enumerate(ix));start=0xff00|bm;end=0xff|bm
   path,n=route(start,end,rr);stat['initializations']+=1;stat['visited']+=n
   if path:
    buffers=[outside[j] for j in ix];row=dict(buffers=buffers,initial=init,path=path,states_discovered=n,moves=len(path)-1);results.append(row);print(row,flush=True);found=True
 stats.append(dict(stat));print(dict(stat),flush=True)
 if found:break
(H/'four-search-results.json').write_text(json.dumps(dict(source_order=src,relations=rels,stats=stats,successes=results),indent=2)+'\n')
