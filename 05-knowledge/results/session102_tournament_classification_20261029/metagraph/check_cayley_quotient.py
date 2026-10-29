#!/usr/bin/env python3
"""Fresh fundamental-triangle quotient encoder and binary rank checks."""
import itertools as it,json,math
from collections import Counter
checks=0
def ck(t):
 global checks;checks+=1
 if not t:raise RuntimeError(checks)
def rank(rows):
 piv={}
 for x in rows:
  while x:
   i=x.bit_length()-1
   if i in piv:x^=piv[i]
   else:piv[i]=x;break
 return len(piv)
def setup(n):
 edges=list(it.combinations(range(n),2));idx={e:k for k,e in enumerate(edges)}
 def edge(u,v):return 1<<idx[tuple(sorted((u,v))) ]
 def path(p):
  x=0
  for u,v in zip(p,p[1:]):x^=edge(u,v)
  return x
 triangles=[edge(0,i)^edge(0,j)^edge(i,j) for i,j in it.combinations(range(1,n),2)]
 def quo(x):return sum(((x&t).bit_count()%2)<<i for i,t in enumerate(triangles))
 cuts=[sum(edge(v,u) for u in range(n) if u!=v) for v in range(n)]
 return edges,edge,path,quo,cuts
rows=[]
for n in range(3,11):
 edges,edge,path,quo,cuts=setup(n);m=len(edges);beta=(n-1)*(n-2)//2;J=(1<<m)-1
 ck(rank(cuts)==n-1);ck(all(quo(c)==0 for c in cuts));ck(rank(quo(1<<i) for i in range(m))==beta)
 if n<=8:
  words=[p for p in it.permutations(range(n)) if p[0]<p[-1]];universe='all undirected Hamilton paths'
 else:
  words=set()
  for a,b,c in it.permutations(range(n),3):
   rest=tuple(v for v in range(n) if v not in (a,b,c));words.add((a,b,c)+rest);words.add((b,a,c)+rest)
  universe='explicit all ordered-three-vertex wedge family (not full path census)'
 raw=[J^path(p) for p in words];gens=[quo(x) for x in raw];base=gens[0]
 r=rank(gens);dr=rank(g^base for g in gens);expected=beta-int(n%4==1)
 ck(r==expected);ck(dr==beta-(n%2));ck(rank(raw)==m-int(beta%2==0));ck(rank(raw+cuts)-(n-1)==r)
 ck(all(g!=0 for g in gens))
 rec={'n':n,'beta':beta,'path_universe':universe,'paths_tested':len(raw),'generator_rank':r,'difference_rank':dr,'components':1<<(beta-r)}
 if n<=8:
  counts=Counter(gens);ck(len(set(counts.values()))==1);rec.update(distinct_generators=len(counts),aggregate_weight=2*next(iter(counts.values())))
 rows.append(rec)
# Exhaustive gauge and aggregate incidence at n=3,4, with direct orientation bits.
aggregates=[]
for n in [3,4]:
 edges,edge,path,quo,cuts=setup(n);m=len(edges);J=(1<<m)-1;words=list(it.permutations(range(n)));C=Counter()
 def legal(x,p):return all(bool(x&edge(u,v))==(u>v) for u,v in zip(p,p[1:])) # x marks backward arcs
 def orient(x,p):
  d={p[0]:0}
  for u,v in zip(p,p[1:]):d[v]=d[u]^bool(x&edge(u,v))^(u>v)
  y=x
  for v in range(n):
   if d[v]:y^=cuts[v]
  return y
 for x in range(1<<m):
  for p in words:
   y=orient(x,p);ck(legal(y,p));ck(quo(x)==quo(y));ck(orient(y,p)==y)
   if legal(x,p):C[quo(x),quo(x^(J^path(p)))]+=1
 G=Counter(quo(J^path(p)) for p in words)
 for q in range(1<<((n-1)*(n-2)//2)):
  ck(sum(v for (a,b),v in C.items() if a==q)==math.factorial(n))
  for g,w in G.items():ck(C[q,q^g]==w)
 aggregates.append({'n':n,'ordered_native_incidence_count':sum(C.values()),'nonzero_aggregate_entries':len(C)})
# A genuine support-edge lift failure from a specified representative.
n=4;edges,edge,path,quo,cuts=setup(n);J=(1<<len(edges))-1
P=(0,1,2,3);Q=(0,1,3,2);gP=quo(J^path(P));gQ=quo(J^path(Q))
ck(gP!=gQ);ck(not legal(0,Q));ck(legal(orient(0,Q),Q))
lift={'n':4,'source_backward_mask':0,'only_source_Hamilton_path':P,'quotient_path_label':Q,'native_reachable_quotient':gP,'requested_quotient':gQ,'Q_gauged_representative_backward_mask':orient(0,Q),'reason':'transitive source has only P; quotient Q edge requires a different switching representative'}
# Full labelled native component census; this is not an unlabelled metagraph census.
native_census=[]
for n in range(3,7):
 edges,edge,path,quo,cuts=setup(n);m=len(edges);J=(1<<m)-1;size=1<<m;parent=list(range(size));deg=[0]*size
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 def unite(x,y):
  x,y=find(x),find(y)
  if x!=y:parent[y]=x
 for x in range(size):
  outs=[0]*n
  for u,v in edges:
   if x&edge(u,v):outs[v]|=1<<u
   else:outs[u]|=1<<v
  def visit(v,left,pm):
   if not left:
    y=x^J^pm;deg[x]+=1;unite(x,y)
    ck((x.bit_count()^y.bit_count())%2==((n-1)*(n-2)//2)%2)
    return
   choices=left&outs[v]
   while choices:
    bit=choices&-choices;choices-=bit;w=bit.bit_length()-1
    visit(w,left^bit,pm^edge(v,w))
  for v in range(n):visit(v,((1<<n)-1)^(1<<v),0)
 comps={}
 for x in range(size):comps.setdefault(find(x),[]).append(x)
 hist=Counter(map(len,comps.values()))
 if n==3:
  ck(sorted(map(len,comps.values()))==[4,4]);ck(all(sorted(deg[x] for x in c)==[1,1,1,3] for c in comps.values()))
 if n in [5,6]:ck(all(len({x.bit_count()%2 for x in c})==1 for c in comps.values()))
 native_census.append({'n':n,'labelled_tournaments':size,'native_edges':sum(deg)//2,'components':len(comps),'component_size_multiplicities':dict(sorted(hist.items())),'switching_quotient_components':2 if n%4==1 else 1,'method':'directed-path backtracking and union-find','component_parity_supports':sorted([sorted({x.bit_count()%2 for x in c}) for c in comps.values()])})
print(json.dumps({'status':'PASS','checks':checks,'rank_census':rows,'native_component_census':native_census,'aggregate_controls':aggregates,'lift_failure':lift},indent=2))
