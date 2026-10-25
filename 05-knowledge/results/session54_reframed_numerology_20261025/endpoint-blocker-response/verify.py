#!/usr/bin/env python3
"""Independent response checks; only standard library; no asserts."""
import collections as C,itertools as I,json,pathlib
D=pathlib.Path(__file__).resolve().parent

def ck(x,msg):
 if not x: raise RuntimeError(msg)

def potentials(n,edges,B):
 adj=[[] for _ in range(2*n)]
 for j,(u,v) in enumerate(edges):
  val=int(bool(B>>j&1));adj[u].append((n+v,-val));adj[n+v].append((u,val))
 p=[None]*(2*n);p[0]=0;todo=[0]
 for u in todo:
  for v,k in adj[u]:
   if p[v] is None:p[v]=p[u]+k;todo.append(v)
   else:ck(p[v]==p[u]+k,'potential consistency')
 ck(None not in p,'connected support');a,d=p[:n],p[n:]
 ck(sum(a)-sum(d)==2,'potential deficiency two')
 return a,d

def fixed_start_counts(n,adj):
 dp={(1,0):1}
 for size in range(1,n):
  nxt=C.defaultdict(int)
  for (mask,u),count in dp.items():
   for v in adj[u]:
    if not mask>>v&1:nxt[mask|1<<v,v]+=count
  dp=nxt
 return [dp.get(((1<<n)-1,v),0) for v in range(n)]

out=[]
source=json.loads((D/'inputs/result.json').read_text())['classification']
for rec in source:
 n=rec['n'];Q={x*x%n for x in range(1,n)};edges=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:j for j,e in enumerate(edges)};adj=[[v for u,v in edges if u==x] for x in range(n)]
 e0=fixed_start_counts(n,adj)
 ck(e0==[0]+[(1 if n==7 else 628) if v in Q else (8 if n==7 else 1101) for v in range(1,n)],'independent endpoint DP')
 endpoint=[[e0[(t-s)%n] for t in range(n)] for s in range(n)]
 ck(all(endpoint[s][t]>0 for s in range(n) for t in range(n) if s!=t),'all ordered endpoints')
 paths=[]
 if n==7:
  def visit(p):
   if len(p)==n:paths.append(p);return
   for v in adj[p[-1]]:
    if v not in p:visit(p+[v])
  for s in range(n):visit([s])
  counts=C.Counter((p[0],p[-1]) for p in paths)
  ck(all(counts[s,t]==endpoint[s][t] for s in range(n) for t in range(n)),'full paths endpoint DP')
 histtypes=C.Counter();familycounts=C.Counter()
 rows=[sum(1<<j for j,e in enumerate(edges) if e[0]==v) for v in range(n)]
 cols=[sum(1<<j for j,e in enumerate(edges) if e[1]==v) for v in range(n)]
 basic={rows[u]|rows[v]:'outgoing_pair' for u,v in I.combinations(range(n),2)}
 basic.update({cols[u]|cols[v]:'incoming_pair' for u,v in I.combinations(range(n),2)})
 basic.update({rows[v]|cols[v]:'star' for v in range(n)})
 for B in rec['all_minimum_blockers']:
  a,d=potentials(n,edges,B);h=C.Counter()
  for s in range(n):
   for t in range(n):
    if s!=t:h[2-a[t]+d[s]]+=endpoint[s][t]
  if paths:
   actual=C.Counter(sum(bool(B>>EI[e]&1) for e in zip(p,p[1:])) for p in paths)
   ck(actual==h,'direct full-path histogram')
   for p in paths:ck(sum(bool(B>>EI[e]&1) for e in zip(p,p[1:]))==2-a[p[-1]]+d[p[0]],'pathwise endpoint identity')
  typ=basic.get(B,'exceptional_Hall');familycounts[typ]+=1;histtypes[(typ,tuple(sorted(h.items())))]+=1
  ck(dict(h)==({1:54,2:135} if n==7 and typ!='exceptional_Hall' else {1:57,2:129,3:3} if n==7 else {1:17290,2:77805}),'predicted histogram')
 out.append({'n':n,'H':sum(map(sum,endpoint)),'full_paths_enumerated':len(paths),'family_counts':dict(familycounts),'histograms':[{'family':typ,'counts':dict(h),'blockers':v} for (typ,h),v in sorted(histtypes.items())]})

def labels(n):
 return [('outgoing_pair',p) for p in I.combinations(range(n),2)]+[('incoming_pair',p) for p in I.combinations(range(n),2)]+[('star',(v,)) for v in range(n)]
def bit(label,s,t):
 typ,p=label
 return int(t in p) if typ=='outgoing_pair' else int(s in p) if typ=='incoming_pair' else int(s in p)+int(t in p)
def signature(label,L):return tuple(bit(label,0,i) for i in L),tuple(bit(label,i,0) for i in L)
def decode(sig,L,w):
 x,y=sig;m=len(L);z=(0,)*m;o=(1,)*m
 if x==y:
  if x==z:return 'star',(w,)
  if x==o:return 'star',(0,)
  ck(sum(x)==1,'equal sides must be unit');return 'star',(L[x.index(1)],)
 def pair(vec,c):
  p=[L[j] for j,b in enumerate(vec) if b]
  if c:p.append(0)
  if len(p)==1:p.append(w)
  ck(len(p)==2,'decoded pair size');return tuple(sorted(p))
 if y in (z,o) and x!=o:return 'outgoing_pair',pair(x,int(y==o))
 ck(x in (z,o),'incoming constant side');return 'incoming_pair',pair(y,int(x==o))
finite=[]
for n in range(5,32):
 labs=labels(n);ck(len(labs)==n*n,'label count')
 for w in range(1,n):
  L=[i for i in range(1,n) if i!=w];sigs=[signature(b,L) for b in labs]
  ck(len(set(sigs))==n*n,'signature injectivity')
  ck(all(decode(s,L,w)==b for s,b in zip(sigs,labs)),'explicit decoder')
 for u,v in I.combinations(range(1,n),2):
  # Even retaining every other query cannot distinguish these pairs.
  queries=[(0,i) for i in range(1,n) if i not in (u,v)]+[(i,0) for i in range(1,n)]
  b1=('outgoing_pair',(0,u));b2=('outgoing_pair',(0,v))
  ck(all(bit(b1,*q)==bit(b2,*q) for q in queries),'lower bound terminal witness')
  queries=[(t,s) for s,t in queries]
  b1=('incoming_pair',(0,u));b2=('incoming_pair',(0,v))
  ck(all(bit(b1,*q)==bit(b2,*q) for q in queries),'lower bound start witness')
 for u,v in I.permutations(range(1,n),2):
  queries=[(0,i) for i in range(1,n) if i!=u]+[(i,0) for i in range(1,n) if i!=v]
  pair=('outgoing_pair',tuple(sorted((u,v))));star=('star',(v,))
  ck(all(bit(pair,*q)==bit(star,*q) for q in queries),'mismatched omission collision')
 finite.append({'n':n,'blocker_labels':len(labs),'omissions_checked':n-1,'queries':2*n-4,'optimal_designs':n-1,'mismatched_omissions_checked':(n-1)*(n-2)})
ck(signature(('outgoing_pair',(1,2)),[1,2])==signature(('incoming_pair',(0,3)),[1,2]),'n4 collision hostile')
# At P7 exceptional Hall responses need not be bits: intersection 3 yields -1.
ck(all(3 not in dict(h) for (typ,h) in histtypes),'P11 has no third intersection')
print(json.dumps({'status':'PASS','paley_checks':out,'signature_finite_checks':finite,'hostiles':{'n4_algebraic_collision':True,'P7_exceptional_response_minus_one':True,'scalar_histogram_not_identifying_P11':True,'arbitrary_query_optimum_not_claimed':True}},indent=2,sort_keys=True))
