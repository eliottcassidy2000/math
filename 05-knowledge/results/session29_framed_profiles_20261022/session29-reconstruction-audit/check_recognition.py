from itertools import combinations
from collections import deque
from pathlib import Path
import json

def req(b,m):
 if not b: raise RuntimeError(m)
def closure(a):
 r=[x|(1<<i) for i,x in enumerate(a)]
 for k in range(len(a)):
  for i in range(len(a)):
   if r[i]>>k&1:r[i]|=r[k]
 return r
def bfs(a):
 seen={0};q=deque([0]);e=set()
 while q:
  u=q.popleft()
  for v in range(len(a)):
   if a[u]>>v&1 and v not in seen:seen.add(v);q.append(v);e.add((u,v))
 return e
strong=deleted=bridges=0
for n in (3,4,5):
 ps=list(combinations(range(n),2));full=(1<<n)-1
 for mask in range(1<<len(ps)):
  a=[0]*n
  for k,(i,j) in enumerate(ps):
   u,v=(i,j) if mask>>k&1 else (j,i);a[u]|=1<<v
  if any(x!=full for x in closure(a)):continue
  strong+=1
  rev=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)]
  candidates=bfs(a)|{(v,u) for u,v in bfs(rev)}
  expected={}
  for cut in range(1,full):
   es=[(u,v) for u in range(n) if cut>>u&1 for v in range(n) if not(cut>>v&1) and a[u]>>v&1]
   if len(es)==1:expected[cut]=es[0]
  recovered={}
  for u in range(n):
   for v in range(n):
    if not(a[u]>>v&1):continue
    deleted+=1;b=a.copy();b[u]&=~(1<<v);r=closure(b)
    if all(x==full for x in r):continue
    bridges+=1;req((u,v) in candidates,'strong bridge missing from trees')
    req(all((r[x]>>y&1) or (r[y]>>x&1) for x in range(n) for y in range(n)),'non-total SCC reachability')
    suffixes={x for x in r if x!=full}
    req(suffixes=={cut for cut,e in expected.items() if e==(u,v)},'suffix cuts incorrect')
    for cut in suffixes:recovered[cut]=(u,v)
  req(recovered==expected,'cut recovery incomplete')
print(json.dumps({'status':'PASS','strong_tournaments':strong,'arc_deletions':deleted,'strong_bridge_instances':bridges,'method':'independent all-arc deletion, Boolean transitive closure, and exhaustive cuts; no producer code imported'},indent=2))
