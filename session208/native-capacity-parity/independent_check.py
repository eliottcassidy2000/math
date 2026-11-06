import itertools,json,collections
from pathlib import Path
r=json.load(open(Path(__file__).resolve().parent/'check.out'));pairs=list(itertools.combinations(range(4),2));ix={e:k for k,e in enumerate(pairs)};E=set()
for p in itertools.permutations(range(4)):
 fixed={ix[tuple(sorted((a,b)))]:int(a<b) for a,b in zip(p,p[1:])};free=set(range(6))-fixed.keys()
 for values in itertools.product((0,1),repeat=3):
  t=sum(v<<k for k,v in fixed.items())+sum(v<<k for k,v in zip(sorted(free),values));u=t^sum(1<<k for k in free);E.add(tuple(sorted((t,u))))
deg=collections.Counter(v for e in E for v in e);assert collections.Counter(deg.values())=={1:24,3:16,5:24};low={v for v in deg if deg[v]<=3};assert len(low)==40;assert all(not(a in low and b in low) for a,b in E)
T={tuple(e) for e in r['tree']};J={tuple(e) for e in r['odd_forest']};assert len(T)==63 and T<=E and J<=T and len(J)==40
adj=collections.defaultdict(list)
for a,b in T:adj[a].append(b);adj[b].append(a)
loads=collections.Counter();transit=collections.Counter()
for a,b in E:
 q=[a];prev={a:None}
 for v in q:
  for u in adj[v]:
   if u not in prev:prev[u]=v;q.append(u)
 assert len(prev)==64
 path=[];v=b
 while v is not None:path.append(v);v=prev[v]
 for v,u in zip(path,path[1:]):loads[tuple(sorted((v,u)))]+=1
 for v in path[1:-1]:transit[v]+=1
assert {e for e in T if loads[e]%2}==J
for v in range(64):assert sum(loads[tuple(sorted((v,u)))] for u in adj[v])==deg[v]+2*transit[v]
assert dict(collections.Counter(loads.values()))=={int(k):v for k,v in r['load_histogram'].items()}
print(json.dumps({'vertices':64,'edges':len(E),'independent_set':40,'tree_edges':len(T),'odd_load_edges':len(J),'path_incidence_local_identity':'PASS','status':'PASS'},indent=2))
