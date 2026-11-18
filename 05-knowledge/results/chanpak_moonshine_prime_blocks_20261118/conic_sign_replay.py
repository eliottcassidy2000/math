from itertools import product,combinations,permutations

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def project(v,p):
 v=tuple(x%p for x in v)
 for a in v:
  if a:
   z=pow(a,-1,p);return tuple(x*z%p for x in v)
 raise ValueError('zero vector')
def Q(v,p):x,y,z=v;return (x*x+y*y-z*z)%p
def B(v,w,p):return (v[0]*w[0]+v[1]*w[1]-v[2]*w[2])%p
def chi(x,p):
 x%=p
 if not x:return 0
 return 1 if pow(x,(p-1)//2,p)==1 else -1
def graph(p,typ):
 pts=sorted({project(v,p) for v in product(range(p),repeat=3) if any(v)})
 vs=[v for v in pts if chi(Q(v,p),p)==typ]
 g=[{j for j,w in enumerate(vs) if i!=j and B(v,w,p)==0} for i,v in enumerate(vs)]
 return vs,g
for p in [5,7,11,19,23]:
 for t,u in product(range(p),repeat=2):
  v=(1-t*t,2*t,1+t*t);w=(1-u*u,2*u,1+u*u)
  need(Q(v,p)==0,'null cone')
  need(B(v,w,p)==(-2*(t-u)**2)%p,'polarization')
  if t!=u:need(chi(B(v,w,p),p)==chi(-2,p),'sign erased')
_,g=graph(5,-1);_,e=graph(5,1)
need(len(g)==10 and all(len(a)==3 for a in g),'Petersen size degree')
for i,j in combinations(range(10),2):need(len(g[i]&g[j])==(0 if j in g[i] else 1),'Petersen srg')
# Independent explicit Petersen model: pairs of five labels, edges disjoint.
pairs=list(combinations(range(5),2));kg=[{j for j,b in enumerate(pairs) if not set(a)&set(b)} for a in pairs]
# Backtracking isomorphism, not uniqueness of SRG assumed.
mapping={};used=set()
def iso():
 if len(mapping)==10:return dict(mapping)
 v=max((v for v in range(10) if v not in mapping),key=lambda v:sum(w in mapping for w in g[v]))
 for w in range(10):
  if w in used:continue
  if all((u in g[v])==(z in kg[w]) for u,z in mapping.items()):
   mapping[v]=w;used.add(w);r=iso()
   if r:return r
   del mapping[v];used.remove(w)
 return None
need(iso()!=None,'explicit Petersen isomorphism')
seen=set();comps=[]
for i in range(len(e)):
 if i in seen:continue
 c={i};todo=[i]
 while todo:
  v=todo.pop()
  for w in e[v]-c:c.add(w);todo.append(w)
 seen|=c;comps.append(c)
need(len(comps)==5 and all(len(c)==3 and all(len(e[v])==2 for v in c) for c in comps),'5K3')
def match(rem):
 if not rem:yield [];return
 v=min(rem)
 for w in g[v]&rem:
  for tail in match(rem-{v,w}):yield [(v,w)]+tail
ms=list(match(set(range(10))))
need(len(ms)==6,'six perfect matchings')
for m in ms:
 edges=set()
 owner={v:i for i,p in enumerate(m) for v in p}
 for v in range(10):
  for w in g[v]:
   if owner[v]!=owner[w]:edges.add(tuple(sorted((owner[v],owner[w]))))
 need(len(edges)==10,'K5 contraction')
for drop in range(10):
 h={v:set(g[v])-{drop} for v in range(10) if v!=drop}
 while any(len(a)==2 for a in h.values()):
  v=next(v for v,a in h.items() if len(a)==2);a,b=h.pop(v);h[a].remove(v);h[b].remove(v);h[a].add(b);h[b].add(a)
 need(len(h)==6 and all(len(a)==3 for a in h.values()),'K33 subdivision')
 color={next(iter(h)):0};todo=list(color)
 while todo:
  v=todo.pop()
  for w in h[v]:
   if w not in color:color[w]=1-color[v];todo.append(w)
   need(color[w]!=color[v],'bipartite')
 need(len(color)==6 and sum(color.values())==3,'K33 parts')
# Explicit A5 involution on the Petersen pair labels.
swap={0:1,1:0,2:3,3:2,4:4};perm=[pairs.index(tuple(sorted((swap[a],swap[b])))) for a,b in pairs]
need(any(i!=j for i,j in enumerate(perm)) and all(perm[perm[i]]==i for i in range(10)),'nontrivial involution')
for i in range(10):need({perm[j] for j in kg[i]}==kg[perm[i]],'A5 Petersen action')
print('PASS conic polarization at5,7,11,19,23; explicit Petersen isomorphism;5K3;6 K5 contractions;10 K33 subdivisions;A5 involution')
print('Petersen projective points',graph(5,-1)[0])
print('Petersen to 2-subsets mapping',iso())
