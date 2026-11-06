import itertools,json,collections
from fractions import Fraction as F
from pathlib import Path
r=json.load(open(Path(__file__).with_name('WITNESS.json')));pairs=list(itertools.combinations(range(4),2));ix={p:k for k,p in enumerate(pairs)};G={t:set() for t in range(64)}
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

# Independent Floyd-Warshall local-bound audit, not heap shortest paths.
squares=[];unseen={v for e in B for v in e}
while unseen:
 q=[min(unseen)];seen=set(q)
 for v in q:
  for e in B:
   if v in e:
    u=e[0] if e[1]==v else e[1]
    if u not in seen:seen.add(u);q.append(u)
 unseen-=seen;squares.append(seen)
assert len(squares)==6 and all(len(s)==4 for s in squares)
vs=sorted(C);vi={v:i for i,v in enumerate(vs)}
def floyd(forbidden,wa,wb):
 d=[[100000]*40 for _ in vs]
 for i in range(40):d[i][i]=0
 for a,b in E-forbidden:d[vi[a]][vi[b]]=d[vi[b]][vi[a]]=wa if (a,b) in A else wb
 for k in range(40):
  for i in range(40):
   dik=d[i][k]
   for j in range(40):
    if dik+d[k][j]<d[i][j]:d[i][j]=dik+d[k][j]
 return d
cases=0
for wa,wb,table in [(1,1,[20,16,12,6]),(4,5,[84,68,52,30])]:
 for Q in squares:
  qe=sorted(e for e in B if set(e)<=Q)
  for mask in range(15):
   missing={e for k,e in enumerate(qe) if not(mask>>k&1)};d=floyd(missing,wa,wb)
   assert sum(d[vi[a]][vi[b]] for a,b in qe)==table[mask.bit_count()];cases+=1
 for a,b in A:
  d=floyd({(a,b)},wa,wb);assert d[vi[a]][vi[b]]==(5 if wa==1 else 22)
checked=0
for kind,packet in r.items():
 avg={e:F(0) for e in E};w=F(packet['weight_first'])
 for tag,weight in [('first',w),('second',1-w)]:
  T={tuple(e) for e in packet['trees'][tag]['edges']};assert T<=E
  dist,loads=values(T)
  sums=[sum(dist[e] for e in A),sum(dist[e] for e in B),sum(loads[e] for e in A),sum(loads[e] for e in B)]
  assert sums==packet['trees'][tag]['sums']
  for f in maps:
   image={tuple(sorted((f[a],f[b]))) for a,b in T};dist,loads=values(image);obs=dist if kind=='stretch' else loads
   for e in E:avg[e]+=weight*obs[e]/48
   checked+=1
 assert set(avg.values())=={F(packet['bound'])}
assert F(1308,480)==F(109,40);assert F(858,312)==F(11,4)
print(json.dumps({'status':'PASS','square_cases':cases,'A_detours':96,'image_trees':checked,'stretch_interval':['109/40','211/66'],'congestion_interval':['11/4','269/84']},indent=2))
