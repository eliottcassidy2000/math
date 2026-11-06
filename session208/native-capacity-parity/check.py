from itertools import combinations,permutations
from collections import Counter
from pathlib import Path
import json
n=4;pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)};full=(1<<len(pairs))-1;N=1<<len(pairs);adj=[set() for _ in range(N)]
for t in range(N):
 for p in permutations(range(n)):
  keep=0
  for a,b in zip(p,p[1:]):
   k=index[tuple(sorted((a,b)))]
   if bool(t>>k&1)!=(a<b):break
   keep|=1<<k
  else:adj[t].add(t^(full^keep))
assert all(t in adj[u] for t in range(N) for u in adj[t])
deg=[len(a) for a in adj];assert Counter(deg)=={1:24,3:16,5:24}
low={v for v in range(N) if deg[v] in (1,3)};assert len(low)==40 and all(not(adj[v]&low) for v in low)
leaves={tuple(sorted((v,next(iter(adj[v]))))) for v in range(N) if deg[v]==1}
assert len(leaves)==24 and len({next(iter(adj[v])) for v in range(N) if deg[v]==1})==24
sources={v for v in range(N) if deg[v]==3}
choices={v:[] for v in sources}
for v in sources:
 for u in adj[v]:
  for w in adj[u]&sources:
   if w!=v:choices[v].append((w,u))
def pair(remaining,used):
 if not remaining:return []
 v=min(remaining)
 for w,u in sorted(choices[v]):
  if w in remaining and u not in used:
   tail=pair(remaining-{v,w},used|{u})
   if tail is not None:return [(v,u,w)]+tail
 return None
triples=pair(sources,set());assert triples is not None and len(triples)==8
J=leaves|{tuple(sorted(e)) for v,u,w in triples for e in ((v,u),(u,w))}
assert len(J)==40
jd=Counter(v for e in J for v in e);assert all(jd[v]%2==1 for v in range(N))
parent=list(range(N))
def root(v):
 while parent[v]!=v:v=parent[v]
 return v
def add(e):
 a,b=map(root,e)
 if a==b:return False
 parent[a]=b;return True
assert all(add(e) for e in sorted(J))
T=set(J)
for v in range(N):
 for u in sorted(adj[v]):
  if v<u and add((v,u)):T.add((v,u))
assert len(T)==N-1
loads={}
for e in sorted(T):
 seen={e[0]};todo=[e[0]]
 while todo:
  v=todo.pop()
  for a,b in T-{e}:
   w=b if a==v else a if b==v else None
   if w is not None and w not in seen:seen.add(w);todo.append(w)
 load=sum(u not in seen for v in seen for u in adj[v]);loads[str(e)]=load
 assert load%2==len(seen)%2
 assert bool(load%2)==(e in J)
weighted=Counter()
for e in T:
 for v in e:weighted[v]+=loads[str(e)]
assert all(weighted[v]>=deg[v] and (weighted[v]-deg[v])%2==0 for v in range(N))
result={'status':'PASS','degree_counts':dict(sorted(Counter(deg).items())),'independent_low_set_size':len(low),'minimum_odd_load_edges':40,'tree_edges':len(T),'odd_forest':sorted(J),'pairing_triples':triples,'tree':sorted(T),'load_histogram':dict(sorted(Counter(loads.values()).items()))}
print(json.dumps(result,indent=2))
