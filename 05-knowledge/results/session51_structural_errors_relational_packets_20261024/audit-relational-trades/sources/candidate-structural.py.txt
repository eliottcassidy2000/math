#!/usr/bin/env python3
"""Hostile controls for signed-equality synchronization and charge packet length."""
from itertools import product,combinations
from pathlib import Path
import json,argparse
P=Path(__file__).resolve().parent;checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def test():
 graphs=endpoints=0
 for n in range(1,5):
  pairs=list(combinations(range(n),2))
  for rel in product((-1,0,1),repeat=len(pairs)):
   edges=[(i,j,b) for (i,j),b in zip(pairs,rel) if b>=0];graphs+=1
   accepted=[x for x in range(1<<n) if all((((x>>i)^(x>>j))&1)==b for i,j,b in edges)]
   adj=[set() for _ in range(n)]
   for i,j,b in edges:adj[i].add(j);adj[j].add(i)
   unseen=set(range(n));comps=[]
   while unseen:
    start=min(unseen);unseen.remove(start);todo=[start];comp={start}
    while todo:
     z=todo.pop()
     for w in adj[z]&unseen:unseen.remove(w);comp.add(w);todo.append(w)
    comps.append(sum(1<<i for i in comp))
   for x,y in combinations(accepted,2):
    changed=x^y;blocks=[b for b in comps if b&changed];need(all(b&changed==b for b in blocks),'endpoint difference is union of relation components')
    margin=lambda z:sum((1+i%2)*((z>>i)&1) for i in range(n))
    charge=lambda z:sum((i+1)*((z>>i)&1) for i in range(n))%6
    if (margin(x),charge(x))!=(margin(y),charge(y)):continue
    endpoints+=1
    actual={z for z in accepted if ((z^x)&~changed)==0 and (margin(z),charge(z))==(margin(x),charge(x))}
    predicted=set()
    for bits in product((0,1),repeat=len(blocks)):
     move=0
     for bit,b in zip(bits,blocks):
      if bit:move|=b
     z=x^move
     if (margin(z),charge(z))==(margin(x),charge(x)):predicted.add(z)
    need(actual==predicted,'relation components plus neutral charge exactly recover accepted interval')
 # Long primitive sequence for C128 x C2^4, checking all subset counts via DP.
 seq=[(1,0)]*127+[(0,1<<j) for j in range(4)]+[(1,15)]
 counts={(0,0):1}
 for a,b in seq:
  nxt=dict(counts)
  for (x,y),count in counts.items():key=((x+a)%128,y^b);nxt[key]=nxt.get(key,0)+count
  counts=nxt
 need(len(seq)==132 and counts[(0,0)]==2,'only empty and full zero-sum subsequences in sharp sequence')
 # Equality path may force arbitrarily many physical changes even with constant charge.
 path=[]
 for n in range(1,10):
  states=[x for x in range(1<<n) if all(((x>>j)&1)==((x>>(j+1))&1) for j in range(n-1))]
  need(states==[0,(1<<n)-1],'pure equality path has only two states')
  path.append({'vertices':n,'forced_simultaneous_slots':n,'treewidth':0 if n==1 else 1})
 return {'status':'PASS','signed_graphs':graphs,'equal_margin_charge_endpoint_pairs':endpoints,'checks':checks,'davenport_constant':132,'minimal_sequence_zero_sum_subsets':counts[(0,0)],'equality_path_controls':path}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();out=test();target=P/'structural-certificate.json'
 if args.write:target.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(out==json.loads(target.read_text()),'frozen structural certificate')
 print(json.dumps(out,sort_keys=True))
