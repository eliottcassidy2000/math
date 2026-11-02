#!/usr/bin/env python3
"""Bounded independent native-meeting controls; not a Collatz convergence test."""
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from itertools import combinations
import json,pathlib,hashlib
N=0
def check(v,why):
 global N
 N+=1
 if not v:raise RuntimeError(why)
def C(n):return n//2 if n%2==0 else 3*n+1
def path(n,root):
 out=[n]
 while out[-1]!=root:
  n=C(out[-1]);check(n not in out,'root-stopped path has no repetition');out.append(n)
 return out
def kids(n):
 return [((2**a*n-1)//3,a) for a in range(1,7) if 2**a*n%9 in (4,7)]
def tree(root,h):
 parent={root:None};exponents={};depth={root:0};level=[root]
 for d in range(1,h+1):
  new=[]
  for q in level:
   for u,a in kids(q):
    check(u not in parent,'distinct selected odd labels');parent[u]=q;exponents[u]=a;depth[u]=d;new.append(u)
  level=new
 return parent,exponents,depth
root=5;parent,exps,depth=tree(root,6)
paths={u:path(u,root) for u in parent}
vertices=set().union(*map(set,paths.values()));edges=set().union(*({(a,b) for a,b in zip(p,p[1:])} for p in paths.values()))
groups=defaultdict(list)
for u,a in exps.items():groups[parent[u]].append(a)
check(len(vertices)==len(parent)+sum(map(max,groups.values())),'merged native vertex formula')
check(len(edges)==len(vertices)-1,'expanded graph tree count')
independent=len(parent)+sum(exps.values())
duplicates=sum(sum(a)-max(a) for a in groups.values())
check(independent-len(vertices)==duplicates,'exact duplicate formula')
check(duplicates==sum(min(a) for a in groups.values()),'binary duplicates equal min')
adj=defaultdict(set)
for a,b in edges:adj[a].add(b);adj[b].add(a)
def odd_path(u):
 out=[u]
 while parent[out[-1]] is not None:out.append(parent[out[-1]])
 return out
op={u:odd_path(u) for u in parent}
pairs=0
for u,v in combinations(parent,2):
 pu,pv=paths[u],paths[v];vp={z:i for i,z in enumerate(pv)}
 iu=next(i for i,z in enumerate(pu) if z in vp);meeting=pu[iu];iv=vp[meeting]
 q=next(z for z in op[u] if z in op[v])
 ou,ov=op[u].index(q),op[v].index(q)
 tu=sum(exps[z]+1 for z in op[u][:ou]);tv=sum(exps[z]+1 for z in op[v][:ov])
 m=0 if q in (u,v) else min(exps[op[u][ou-1]],exps[op[v][ov-1]])
 check((meeting,iu,iv)==(2**m*q,tu-m,tv-m),'odd LCA first meeting and asynchronous clocks')
 # Literal native-tree BFS, independent of the clock-derived distance.
 queue=deque([u]);dist={u:0}
 while v not in dist:
  z=queue.popleft()
  for w in adj[z]:
   if w not in dist:dist[w]=dist[z]+1;queue.append(w)
 check(dist[v]==tu+tv-2*m,'native undirected distance')
 pairs+=1

def polynomial_literal(n,h):
 pa,_,_=tree(n,h);xs=set().union(*(set(path(u,n)) for u in pa))
 return Counter(len(path(x,n))-1 for x in xs)
def polynomial_recurrence(n,h):
 if h==0:return Counter({0:1})
 children=kids(n);out=Counter({i:1 for i in range(max(a for _,a in children)+1)})
 for u,a in children:
  out.update({i+a+1:c for i,c in polynomial_recurrence(u,h-1).items()})
 return out
for n in (5,13,23,53):
 for h in range(5):check(polynomial_literal(n,h)==polynomial_recurrence(n,h),'source-dependent port polynomial')
check(polynomial_literal(5,1)==Counter({0:1,1:1,2:1,3:1,4:2,5:1,6:1}),'literal n5 h1 polynomial')

def first(a,b):
 pa,pb=path(a,5),path(b,5);lookup={v:i for i,v in enumerate(pb)}
 i=next(i for i,v in enumerate(pa) if v in lookup)
 return pa[i],i,lookup[pa[i]]
check(first(13,53)==(40,1,3),'smallest sibling hostile')
check(first(853,53)==(160,5,1),'equal phase unequal fixed-port response')
check(3*13+1==2**3*5 and 3*853+1==2**9*5 and 3%6==9%6,'legal equal exponent phase')

def det(A):
 A=[[F(v) for v in row] for row in A];result=F(1)
 for i in range(len(A)):
  pivot=next((k for k in range(i,len(A)) if A[k][i]),None)
  if pivot is None:return F(0)
  if pivot!=i:A[i],A[pivot]=A[pivot],A[i];result=-result
  v=A[i][i];result*=v
  for j in range(i+1,len(A)):
   factor=A[j][i]/v
   for k in range(i,len(A)):A[j][k]-=factor*A[i][k]
 return result
gram_tests=0
for m in range(1,5):
 for aa in combinations(range(1,11),m):
  expected=aa[0]
  for i in range(1,m):expected*=aa[i]-aa[i-1]
  check(det([[min(a,b) for b in aa] for a in aa])==expected,'min-kernel Gram determinant')
  gram_tests+=1
for weights,expected in (((1,0,0,2,0,0,1,0,0),2),((1,0,0,0,0,0,1,0,0),0)):
 sums=[sum(weights[:a]) for a in (3,5,9)]
 check(det([[min(a,b) for b in sums] for a in sums])==expected,'weighted positive or collapsed block')
check(det([[0]])==0,'depth-zero root path has zero Gram rank')
leaves=[u for u in parent if depth[u]==6]
leaf_edges={u:set(zip(paths[u],paths[u][1:])) for u in leaves}
for u in leaves:
 firstedge=(u,C(u))
 check(all((firstedge in leaf_edges[v])==(u==v) for v in leaves),'private first-edge identity minor')
pin=json.loads((pathlib.Path(__file__).parent/'SOURCE_PINS.json').read_text());data=(pathlib.Path(__file__).parent/pin['local_copy']).read_bytes()
check(len(data)==pin['bytes'] and hashlib.sha256(data).hexdigest()==pin['sha256'],'exact Git source copy pin')
print(json.dumps({'status':'PASS','checks':N,'odd_vertices':len(parent),'native_vertices':len(vertices),'native_edges':len(edges),'duplicate_vertex_occurrences':duplicates,'all_odd_source_pairs':pairs,'source_depth_polynomial_comparisons':20,'Gram_determinants':gram_tests,'positive_depth_leaf_Gram_rank':len(leaves),'depth_zero_leaf_Gram_rank':0,'hostiles':{'13_53':[40,1,3],'853_53':[160,5,1]},'scope':'J1-J4 reconstruction identities on supplied finite rooted odd trees; selected finite family is independent and not a full SF5-packet replay'},indent=2,sort_keys=True))
