import itertools,collections,json,random
n=4;pairs=list(itertools.combinations(range(n),2));ix={p:i for i,p in enumerate(pairs)};G={t:set() for t in range(64)}
for p in itertools.permutations(range(n)):
 fixed={ix[tuple(sorted((a,b)))]:int(a<b) for a,b in zip(p,p[1:])};free=[k for k in range(6) if k not in fixed]
 for bits in itertools.product((0,1),repeat=3):
  t=sum(v<<k for k,v in fixed.items())+sum(v<<k for k,v in zip(free,bits));u=t^sum(1<<k for k in free);G[t].add(u)
leaves={v for v in G if len(G[v])==1}
assert len(leaves)==24 and len({next(iter(G[v])) for v in leaves})==24
assert all(len(G[next(iter(G[v]))])==5 for v in leaves)
C={v for v in G if len(G[v])>1};E=sorted((u,v) for u in C for v in G[u]&C if u<v);A={e for e in E if len(G[e[0]])!=len(G[e[1]])};B=set(E)-A

def paths(T):
 adj={u:set() for u in C}
 for u,v in T:adj[u].add(v);adj[v].add(u)
 ps=[]
 for a,b in E:
  q=[a];par={a:None}
  for v in q:
   for u in adj[v]:
    if u not in par:par[u]=v;q.append(u)
  p=[];v=b
  while v!=a:u=par[v];p.append(tuple(sorted((u,v))));v=u
  ps.append(p)
 return ps

def stats(T):
 ps=paths(T);sa=sum(len(p) for e,p in zip(E,ps) if e in A);sb=sum(len(p) for e,p in zip(E,ps) if e in B);ca=sum(sum(f in A for f in p) for p in ps);cb=sum(sum(f in B for f in p) for p in ps)
 return sa,sb,ca,cb
from fractions import Fraction as F

def comps(adj):
 unseen=set(adj);out=[]
 while unseen:
  q=[min(unseen)];seen=set(q)
  for u in q:
   for v in adj[u]:
    if v not in seen:seen.add(v);q.append(v)
  unseen-=seen;out.append(sorted(seen))
 return out
HA={u:{v for v in G[u]&C if tuple(sorted((u,v))) in A} for u in C};ac=comps(HA)
assert [len(c) for c in ac]==[20,20]
D={u for u in C if len(G[u])==5};HB={u:G[u]&D for u in D};bc=comps(HB);assert list(map(len,bc))==[4]*6
for c in ac:
 low={u for u in c if len(G[u])==3};high=set(c)-low;assert len(low)==8 and len(high)==12
 cube={u:set() for u in low}
 for v in high:
  ns=sorted(HA[v]);assert len(ns)==2;a,b=ns;cube[a].add(b);cube[b].add(a)
 assert all(len(ns)==3 for ns in cube.values())
 color={min(low):0};q=[min(low)]
 for u in q:
  for v in cube[u]:
   if v not in color:color[v]=1-color[u];q.append(v)
   assert color[v]!=color[u]
 assert sum(color.values())==4
 # Cubic bipartite graph with4+4 vertices is K4,4 minus a perfect matching, hence cube.
# Local edge detours.
for a,b in E:
 dist={a:0};q=[a]
 for u in q:
  for v in G[u]&C:
   if {u,v}=={a,b}:continue
   if v not in dist:dist[v]=dist[u]+1;q.append(v)
 assert dist[b]==(5 if (a,b) in A else 3)
# Verify all48 relabelling/converse automorphisms and exactly two edge orbits.
maps=[]
for p in itertools.permutations(range(4)):
 for converse in (0,1):
  f={}
  for t in G:
   z=0
   for i,(a,b) in enumerate(pairs):
    u,v=p[a],p[b];bit=(t>>i)&1
    if u>v:bit=1-bit
    if converse:bit=1-bit
    z|=bit<<ix[tuple(sorted((u,v)))]
   f[t]=z
  assert all({f[v] for v in G[u]}==G[f[u]] for u in G)
  maps.append(f)
for S in (A,B):
 a,b=min(S);assert {tuple(sorted((f[a],f[b]))) for f in maps}==S

def bfs(root):
 q=[root];seen={root};T=[]
 for u in q:
  for v in sorted(G[u]&C):
   if v not in seen:seen.add(v);q.append(v);T.append(tuple(sorted((u,v))))
 assert len(T)==39
 return sorted(T)
trees={r:bfs(r) for r in (14,19,51)};st={r:stats(t) for r,t in trees.items()}
assert st=={14:(160,74,160,74),19:(128,122,196,54),51:(168,66,156,78)}
w=F(29,32);assert w*F(st[14][0],48)+(1-w)*F(st[19][0],48)==F(157,48)
assert w*F(st[14][1],24)+(1-w)*F(st[19][1],24)==F(157,48)
assert F(st[51][2],48)==F(st[51][3],24)==F(13,4)
# Contract squares: eight distinct triples, each repeated twice, the octahedron faces.
qi={v:i for i,c in enumerate(bc) for v in c};triples=collections.Counter(tuple(sorted(qi[v] for v in G[u])) for u in C if len(G[u])==3)
assert len(triples)==8 and set(triples.values())=={2}
nonpairs=[list(p) for p in itertools.combinations(range(6),2) if not any(set(p)<=set(t) for t in triples)]
assert len(nonpairs)==3 and len(set(sum(nonpairs,[])))==6
result={'status':'PASS','degree_census':dict(sorted(collections.Counter(map(len,G.values())).items())),'core_vertices':len(C),'core_edges':len(E),'A_edges':len(A),'B_edges':len(B),'subdivided_cubes':ac,'squares':bc,'octahedron_opposite_pairs':nonpairs,'edge_orbit_sizes':[48,24],'automorphism_count':len(maps),'trees':{r:{'edges':trees[r],'stretch_A_sum':st[r][0],'stretch_B_sum':st[r][1],'load_A_sum':st[r][2],'load_B_sum':st[r][3]} for r in trees},'bounds':{'stretch_lower':'19/8','stretch_upper':'157/48','capacity_lower':'79/36','capacity_upper':'13/4'},'stretch_mixture':{'root14':'29/32','root19':'3/32','symmetrize':'uniform relabelling and converse'},'capacity_distribution':{'root51':'1','symmetrize':'uniform relabelling and converse'}}
print(json.dumps(result,indent=2))
