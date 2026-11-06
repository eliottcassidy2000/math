import itertools,json,heapq,collections
from fractions import Fraction as F
from pathlib import Path
R=json.load(open(Path(__file__).with_name('WITNESS.json')))
pairs=list(itertools.combinations(range(4),2));ix={p:k for k,p in enumerate(pairs)};G={t:set() for t in range(64)}
for p in itertools.permutations(range(4)):
 fixed={ix[tuple(sorted((a,b)))]:int(a<b) for a,b in zip(p,p[1:])};free=[k for k in range(6) if k not in fixed]
 for bits in itertools.product((0,1),repeat=3):
  t=sum(v<<k for k,v in fixed.items())+sum(v<<k for k,v in zip(free,bits));u=t^sum(1<<k for k in free);G[t].add(u)
C={v for v in G if len(G[v])>1};E=sorted((a,b) for a in C for b in G[a]&C if a<b);A={e for e in E if len(G[e[0]])!=len(G[e[1]])};B=set(E)-A
D={v for v in C if len(G[v])==5};unseen=set(D);squares=[]
while unseen:
 q=[min(unseen)];seen=set(q)
 for a in q:
  for b in G[a]&D:
   if b not in seen:seen.add(b);q.append(b)
 unseen-=seen;squares.append(sorted(seen))
assert len(A)==48 and len(B)==24 and list(map(len,squares))==[4]*6

def distances(source,forbidden,wa,wb):
 dist={source:0};Q=[(0,source)]
 while Q:
  du,u=heapq.heappop(Q)
  if dist[u]!=du:continue
  for v in G[u]&C:
   e=tuple(sorted((u,v)))
   if e in forbidden:continue
   z=du+(wa if e in A else wb)
   if z<dist.get(v,100000):dist[v]=z;heapq.heappush(Q,(z,v))
 return dist
checks=[]
for wa,wb,expected,detour in [(1,1,[20,16,12,6],5),(4,5,[84,68,52,30],22)]:
 for a,b in A:assert distances(a,{(a,b)},wa,wb)[b]==detour
 tab=[set() for _ in range(4)]
 for S in squares:
  edges=sorted(e for e in B if set(e)<=set(S));assert len(edges)==4
  for mask in range(15):
   missing={e for i,e in enumerate(edges) if not(mask>>i&1)}
   total=sum(distances(a,missing,wa,wb)[b] for a,b in edges);tab[mask.bit_count()].add(total)
 assert tab==[{v} for v in expected]
 checks.append({'length_A':wa,'length_B':wb,'A_detour':detour,'square_distance_sums':expected,'subsets_tested':90})
# Forty-eight exact graph automorphisms.
maps=[]
for p in itertools.permutations(range(4)):
 for rev in (0,1):
  f={}
  for t in G:
   z=0
   for i,(a,b) in enumerate(pairs):
    u,v=p[a],p[b];bit=(t>>i)&1
    if u>v:bit=1-bit
    if rev:bit=1-bit
    z|=bit<<ix[tuple(sorted((u,v)))]
   f[t]=z
  assert all({f[v] for v in G[u]}==G[f[u]] for u in G);maps.append(f)
for S in (A,B):
 a,b=min(S);assert {tuple(sorted((f[a],f[b]))) for f in maps}==S

def evaluate(T):
 T={tuple(e) for e in T};assert len(T)==39 and T<=set(E);adj={v:set() for v in C}
 for a,b in T:adj[a].add(b);adj[b].add(a)
 paths={};loads={e:0 for e in E}
 for a,b in E:
  par={a:None};q=[a]
  for u in q:
   for v in adj[u]:
    if v not in par:par[v]=u;q.append(v)
  assert len(par)==40;p=[];v=b
  while v!=a:u=par[v];e=tuple(sorted((u,v)));p.append(e);loads[e]+=1;v=u
  paths[a,b]=len(p)
 s=[sum(paths[e] for e in A),sum(paths[e] for e in B),sum(loads[e] for e in A),sum(loads[e] for e in B)]
 return s,paths,loads
out={}
for name,packet in R.items():
 expected=F(packet['bound']);w=F(packet['weight_first']);aver={e:F(0) for e in E}
 for tag,wt in [('first',w),('second',1-w)]:
  tree=packet['trees'][tag];assert evaluate(tree['edges'])[0]==tree['sums']
  for f in maps:
   image=[sorted((f[a],f[b])) for a,b in tree['edges']];s,paths,loads=evaluate(image);values=paths if name=='stretch' else loads
   for e in E:aver[e]+=wt*values[e]/48
 assert set(aver.values())=={expected};out[name]=str(expected)
assert F(1308,480)==F(109,40) and F(858,312)==F(11,4)
print(json.dumps({'status':'PASS','square_tables':checks,'upper_bounds':out,'lower_bounds':{'stretch':'109/40','capacity':'11/4'},'image_trees_checked':192},indent=2))
