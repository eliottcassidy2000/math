import json,itertools,collections,hashlib,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
B=ROOT.parent/'biclique-bridge'; P=ROOT.parent/'biclique-path-counts'
def check(ok,msg):
 if not ok: raise RuntimeError(msg)
# Independent DFS through actual bipartite vertices, with unique L-start orientation.
paths=[]
def dfs(seq):
 if len(seq)==8:
  edges=frozenset((min(a,b),max(a,b)-4) for a,b in zip(seq,seq[1:])); paths.append((tuple(seq),edges));return
 for v in (range(4,8) if seq[-1]<4 else range(4)):
  if v not in seq: dfs(seq+[v])
for v in range(4):dfs([v])
check(len(paths)==576 and len({e for _,e in paths})==576,'paths')
S=((0,0),(1,0),(0,1),(1,1)); E=list(itertools.product(range(4),repeat=2))
labels={(i,j):(S[j][0]-S[i][0],S[j][1]-S[i][1]) for i,j in E}
check(sorted(collections.Counter(labels.values()).values())==[1]*4+[2]*4+[4],'multiplicities')
pm=[sum(1<<(4*i+j) for i,j in e) for _,e in paths]
pl=[set(labels[e] for e in es) for _,es in paths]
cycles={frozenset(es|{(seq[0],seq[-1]-4)}) for seq,es in paths}
check(len(cycles)==72,'cycles')
for sets,expected in [([es for _,es in paths],(252,72,120)),(cycles,(36,12,20))]:
 check({sum(e in es for es in sets) for e in E}=={expected[0]},'edge incidence')
 for adjacent, exp in [(True,expected[1]),(False,expected[2])]:
  got={sum(e in es and f in es for es in sets) for e,f in itertools.combinations(E,2) if (e[0]==f[0] or e[1]==f[1])==adjacent}
  check(got=={exp},'pair incidence')
dist=collections.Counter(); optimal=[]
for i,j in itertools.combinations(range(576),2):
 if len(pl[i]|pl[j])==9:
  n=len(paths[i][1]|paths[j][1]);dist[n]+=1
  if n==9:optimal.append((i,j))
check(dict(dist)=={9:204,10:1852,11:5600,12:7720,13:3760,14:1488},'union distribution')
check(sum(dist.values())==20624,'pair count')
# Deliberately independent of producer incidence-bitset intersection algorithm:
# enumerate all 31,684,800 unordered triples and directly OR their edge masks.
triples=sum((a|b|c)==65535 for a,b,c in itertools.combinations(pm,3))
check(triples==675265,'triples')
check(sum(len(a|b)==16 for a,b in itertools.combinations(cycles,2))==36,'cycle covers')
fixed=[]
for a in itertools.permutations(range(4)):
 for b in itertools.permutations(range(4)):
  for swap in (False,True):
   if all(labels[i,j]==labels[(b[j],a[i]) if swap else (a[i],b[j])] for i,j in E):fixed.append((a,b,swap))
check(fixed==[((0,1,2,3),(0,1,2,3),False),((3,2,1,0),(3,2,1,0),True)],'fixed-label group')
sharp=[[0,7,1,5,2,4,3,6],[0,7,2,5,1,6,3,4]]
se=[]
for seq in sharp:
 check(len(set(seq))==8,'sharp path vertices');se.append(frozenset((min(a,b),max(a,b)-4) for a,b in zip(seq,seq[1:])))
check(len(se[0]|se[1])==9 and len(se[0]&se[1])==5 and len({labels[e] for e in se[0]|se[1]})==9,'sharp witness')
# Direct reconstruction indices and integer-coefficient line check.
indices=[(0,1),(3,0),(0,2),(1,2),(0,0),(2,1),(2,0),(0,3),(1,0)]
M=[(1,*labels[e]) for e in indices]
check(M==[(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)],'reconstruction')
lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
check(all(tuple(sum(M[k][d] for k in line) for d in range(3))==(3,0,0) for line in lines),'magic lines')
# Independent polynomial field implementation, including supplied-modulus irreducibility.
def rem(f,g,p):
 r=list(f)
 while len(r)>=len(g):
  z=r[-1];off=len(r)-len(g)
  for j,v in enumerate(g):r[off+j]=(r[off+j]-z*v)%p
  while r and r[-1]==0:r.pop()
 return r
class F:
 def __init__(self,p,n,f):
  self.p,self.n,self.q,self.f=p,n,p**n,f
  check(f[-1]==1 and len(f)==n+1,'monic modulus')
  if n>1:
   check(all(rem(f,list(g)+[1],p) for d in range(1,n//2+1) for g in itertools.product(range(p),repeat=d)),'irreducible modulus')
 def digits(self,a):return [(a//(self.p**i))%self.p for i in range(self.n)]
 def pack(self,a):return sum((v%self.p)*self.p**i for i,v in enumerate(a))
 def add(self,a,b):return self.pack([x+y for x,y in zip(self.digits(a),self.digits(b))])
 def neg(self,a):return self.pack([-x for x in self.digits(a)])
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  c=[0]*(2*self.n-1)
  for i,x in enumerate(self.digits(a)):
   for j,y in enumerate(self.digits(b)):c[i+j]=(c[i+j]+x*y)%self.p
  return self.pack(rem(c,self.f,self.p))
cert=json.loads((B/'certificate.json').read_text()); gauges=[]
for w in cert['field_witnesses']:
 f=F(w['p'],w['n'],w['modulus']);sq={f.mul(i,i) for i in range(f.q)};c,x,y=w['c'],w['x'],w['y'];ss=[0,x,y,f.add(x,y)]
 check(len(sq)==(f.q+1)//2,'square count')
 m=w['matrix'];check(len(set(m))==9 and set(m)<=sq,'field squares')
 check(len({f.add(f.add(m[a],m[b]),m[d]) for a,b,d in lines})==1,'field magic')
 forbidden=[]
 for g in range(f.q):
  L=[f.sub(f.sub(c,s),g) for s in ss];R=[f.add(s,g) for s in ss]
  edges=[f.add(l,r) for l in L for r in R]
  check(set(edges)==set(m),'cross sums')
  if len(set(L+R))<8:forbidden.append(g)
 check(len(forbidden)==9 and forbidden==w['forbidden_gauges'] and f.q-9==w['good_gauge_count'],'gauge count')
 gauges.append({'q':f.q,'good':f.q-9})
# Read-only comparison of producer stored finite enumerations.
r=json.loads((P/'results.json').read_text());extra=json.loads((P/'minimal_union.json').read_text())
check(r['number_edge_cover_triples']==triples and r['number_label_cover_pairs']==sum(dist.values()),'stored counts')
check(extra['label_cover_pair_union_edge_distribution']=={str(k):v for k,v in sorted(dist.items())},'stored distribution')
files=[B/'THEOREM.md',B/'verify.py',B/'certificate.json',P/'README.md',P/'verify.py',P/'results.json',P/'minimal_union.json']
report={'status':'PASS','paths':576,'cycles':72,'label_cover_pairs':sum(dist.values()),'edge_cover_triples':triples,'minimal_nine_edge_pairs':len(optimal),'union_distribution':dict(sorted(dist.items())),'gauges':gauges,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
print(json.dumps(report,indent=2,sort_keys=True))
