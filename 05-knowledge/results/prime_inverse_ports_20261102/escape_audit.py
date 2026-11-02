#!/usr/bin/env python3
"""Independent path-cube census and source/port-preserving escape audit."""
from itertools import permutations,combinations
from collections import Counter,deque
from pathlib import Path
import json
from research import need,strong
HERE=Path(__file__).resolve().parent
EDGES=list(combinations(range(6),2));IDX={e:i for i,e in enumerate(EDGES)}
def path_bits(w):
 used=value=0
 for a,b in zip(w,w[1:]):
  j=IDX[tuple(sorted((a,b)))];used|=1<<j
  if a<b:value|=1<<j
 return used,value
PATHS=[(w,*path_bits(w)) for w in permutations(range(6))]
def rows(mask):
 out=[0]*6
 for j,(a,b) in enumerate(EDGES):
  if mask>>j&1:out[a]|=1<<b
  else:out[b]|=1<<a
 return out

def main():
 h=[0]*32768
 for w,used,value in PATHS:
  free=32767^used;sub=free
  while True:
   h[value|sub]+=1
   if not sub:break
   sub=(sub-1)&free
 adj=[[] for _ in range(32768)]
 for idx,(w,used,value) in enumerate(PATHS):
  free=32767^used;sub=free
  while True:
   x=value|sub;y=x^free;adj[x].append((y,idx))
   if not sub:break
   sub=(sub-1)&free
 minima=[x for x in range(32768) if h[x]!=17 and all(abs(h[x^(1<<i)]-17)>abs(h[x]-17) for i in range(15))]
 need(len(minima)==768,'strict minima census')
 targets=[x for x in range(32768) if h[x]==17];need(len(targets)==1440,'target census')
 nd={x:0 for x in targets};nxt={};todo=deque(targets)
 while todo:
  x=todo.popleft()
  for y,idx in adj[x]:
   if y not in nd:nd[y]=nd[x]+1;nxt[y]=(x,idx);todo.append(y)
 native_hist=Counter(nd.get(x,-1) for x in minima)
 def reached(bound):
  seen=set(targets);todo=deque(targets)
  while todo:
   x=todo.popleft()
   for i in range(15):
    y=x^(1<<i)
    if y not in seen and abs(h[y]-17)<=bound:seen.add(y);todo.append(y)
  return seen
 need(not(set(minima)&reached(4)),'barrier lower bound')
 need(set(minima)<=reached(6),'barrier upper bound')
 unseen=set(minima);records=[];distance_counts=Counter()
 while unseen:
  x=min(unseen);r=rows(x);orbit=set()
  for perm in permutations(range(6)):
   # new vertex i carries old vertex perm[i]
   y=sum(1<<j for j,(a,b) in enumerate(EDGES) if r[perm[a]]>>perm[b]&1);orbit.add(y)
  need(orbit<=unseen,'disjoint orbit cover');unseen-=orbit
  oldpaths=[w for w,u,v in PATHS if x&u==v]
  nearest=None
  for d in range(1,4):
   for changes in combinations(range(15),d):
    delta=sum(1<<j for j in changes);y=x^delta
    if h[y]==17:
     common=[w for w,u,v in PATHS if x&u==v and y&u==v]
     if common:nearest=dict(distance=d,changes=changes,target=y,shared_path=common[0]);break
   if nearest:break
  need(nearest is not None,'three-bit common path escape')
  # Every closer target is excluded, without any common-path restriction.
  need(all(h[x^sum(1<<j for j in c)]!=17 for d in range(1,nearest['distance']) for c in combinations(range(15),d)),'ambient distance lower bound')
  barrier_walk=None
  for change in permutations(range(15),3):
   w=[x]
   for j in change:w.append(w[-1]^(1<<j))
   if h[w[-1]]!=17 or max(abs(h[y]-17) for y in w)>6:continue
   common=[p for p,u,v in PATHS if all(y&u==v for y in w)]
   if common:
    barrier_walk=dict(changes=change,masks=w,counts=[h[y] for y in w],shared_path=common[0]);break
  need(barrier_walk is not None,'three-step barrier-optimal path keeps port')
  # Native single path-complement is separately checked and fails.
  need(all(h[x^(32767^u)]!=17 for _,u,v in PATHS if x&u==v),'path complement does not repair')
  scores=sorted(y.bit_count() for y in r);target_scores=[1,1,2,3,4,4]
  score_lower=sum(abs(a-b) for a,b in zip(scores,target_scores))//2
  distance_counts[nearest['distance']]+=len(orbit)
  native=[];cur=x
  while cur in nxt:
   target,idx=nxt[cur];w,used,value=PATHS[idx]
   need(cur&used==value and target==cur^(32767^used),'native guarded complement')
   native.append(dict(source=cur,target=target,counts=[h[cur],h[target]],guard_path=w));cur=target
  need(h[cur]==17,'native target reached')
  records.append(dict(native_complement_walk=native,source=x,H=h[x],orbit_size=len(orbit),automorphism_order=720//len(orbit),strong=strong(r),scores=scores,degree_distance_lower_bound=score_lower,nearest=nearest,barrier_walk=barrier_walk))
 out=dict(status='PASS',method='sum the 720 actual Hamilton-path Boolean subcubes; no endpoint DP',labelled_sources=32768,strict_minima=768,target_sources=1440,barrier_absolute=6,barrier_squared=36,distance_histogram=dict(distance_counts),native_complement_distance_histogram=dict(native_hist),orbit_representatives=records)
 producer=json.loads((HERE/'local_escape.json').read_text())
 need({str(k):v for k,v in distance_counts.items()}==producer['strict_minimum_distance_histogram'],'independent distances agree')
 (HERE/'escape_audit.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
