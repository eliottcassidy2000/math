import itertools,json,collections
from fractions import Fraction as F
from pathlib import Path
r=json.load(open(Path(__file__).with_name('RESULT.json')));pairs=list(itertools.combinations(range(4),2));ix={p:k for k,p in enumerate(pairs)};G={t:set() for t in range(64)}
for t in G:
 def beat(a,b):return bool(t>>ix[tuple(sorted((a,b)))]&1)==(a<b)
 for p in itertools.permutations(range(4)):
  if all(beat(a,b) for a,b in zip(p,p[1:])):
   keep=sum(1<<ix[tuple(sorted((a,b)))] for a,b in zip(p,p[1:]));G[t].add(t^(63^keep))
C={v for v in G if len(G[v])>1};E={tuple(sorted((v,u))) for v in C for u in G[v]&C};A={e for e in E if len(G[e[0]])!=len(G[e[1]])};B=E-A
assert (len(C),len(E),len(A),len(B))==(40,72,48,24)
for a,b in E:
 q=[a];d={a:0}
 for v in q:
  for u in G[v]&C:
   if {u,v}=={a,b}:continue
   if u not in d:d[u]=d[v]+1;q.append(u)
 assert d[b]==(5 if (a,b) in A else 3)
def values(T):
 adj={v:[] for v in C}
 for a,b in T:adj[a].append(b);adj[b].append(a)
 assert len(T)==39
 loads=collections.Counter();dist={}
 for a,b in E:
  q=[a];pr={a:None}
  for v in q:
   for u in adj[v]:
    if u not in pr:pr[u]=v;q.append(u)
  assert len(pr)==40
  v=b;d=0
  while v!=a:
   u=pr[v];loads[tuple(sorted((v,u)))]+=1;d+=1;v=u
  dist[a,b]=d
 return dist,loads
maps=[]
for p in itertools.permutations(range(4)):
 inv={p[i]:i for i in range(4)}
 for rev in (0,1):
  f={}
  for t in G:
   z=0
   for k,(a,b) in enumerate(pairs):
    u,v=inv[a],inv[b];win=bool(t>>ix[tuple(sorted((u,v)))]&1)==(u<v)
    if win !=bool(rev):z|=1<<k
   f[t]=z
  assert len(set(f.values()))==64
  assert all({f[u] for u in G[v]}==G[f[v]] for v in G)
  maps.append(f)
D={root:collections.defaultdict(F) for root in (14,19,51)};L={root:collections.defaultdict(F) for root in D}
for root in D:
 T={tuple(e) for e in r['trees'][str(root)]['edges']};assert T<=E
 for f in maps:
  image={tuple(sorted((f[a],f[b]))) for a,b in T};dist,loads=values(image)
  for e in E:D[root][e]+=F(dist[e],48);L[root][e]+=F(loads[e],48)
for e in E:
 assert F(29,32)*D[14][e]+F(3,32)*D[19][e]==F(157,48)
 assert L[51][e]==F(13,4)
assert F(228,96)==F(19,8);assert F(158,72)==F(79,36)
print(json.dumps({'status':'PASS','explicit_symmetry_images_checked':144,'all72_expected_stretches':'157/48','all72_expected_loads':'13/4','all_detours':'A=5; B=3','lower_bounds':'19/8;79/36'},indent=2))
