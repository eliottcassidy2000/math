"""Constructive reduction with independent native-relation reconstruction."""
import itertools as it, collections, functools, json, random, hashlib

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def connected(n,edges):
 seen={0}
 while True:
  more=seen|{v for u,v in edges if u in seen}|{u for u,v in edges if v in seen}
  if more==seen:return len(seen)==n
  seen=more

def reduce_graph(n,edges):
 incident=[[] for _ in range(n)];copies={};owner=[]
 for e,(u,v) in enumerate(edges):
  for x in (u,v):copies[x,e]=len(owner);incident[x].append(e);owner.append(x)
 constraints={}
 def add(u,v,E):
  key=tuple(sorted((u,v)));need(u!=v and key not in constraints,'simple constraint graph');constraints[key]=set(E)
 for e,(u,v) in enumerate(edges):add(copies[u,e],copies[v,e],{1,2})
 for u in range(n):
  for e,f in zip(incident[u],incident[u][1:]):add(copies[u,e],copies[u,f],{0})
 deg=collections.Counter(x for uv in constraints for x in uv)
 need(max(deg.values())<=3,'copy degree')
 odd=sorted(x for x in range(len(owner)) if deg[x]%2);need(len(odd)%2==0,'odd parity')
 for u,v in zip(odd[::2],odd[1::2]):
  z=len(owner);owner.append(None);add(u,z,{0,1,2});add(v,z,{0,1,2})
 root=0;s=len(owner);owner.append(None);t=len(owner);owner.append(None)
 add(root,s,{0,1,2});add(root,t,{0,1,2})
 cells=len(owner);adj=[set() for _ in owner]
 for u,v in constraints:adj[u].add(v);adj[v].add(u)
 need(max(map(len,adj))<=6 and [i for i,a in enumerate(adj) if len(a)%2]==[s,t],'Euler degrees')
 unused=[set(a) for a in adj];stack=[s];trail=[]
 while stack:
  u=stack[-1]
  if unused[u]:
   v=min(unused[u]);unused[u].remove(v);unused[v].remove(u);stack.append(v)
  else:trail.append(stack.pop())
 trail.reverse();need(trail[0]==s and trail[-1]==t and len(trail)==len(constraints)+1,'Euler trail')
 visits=collections.Counter(trail);need(max(visits.values())<=3,'three visits')
 used=collections.Counter();template=[]
 for u in trail:
  length=4-visits[u] if used[u]==0 else 1
  for _ in range(length):template.append((u,used[u]));used[u]+=1
 need(len(template)==3*cells and all(v==3 for v in used.values()),'spanning template')
 N=3*cells;a=[[0]*N for _ in range(N)]
 for u in range(N):
  for v in range(u+1,N):
   bit=int((v-u)%3==1) if u//3==v//3 else 1
   a[u][v]=bit;a[v][u]=1-bit
 boundary_seen=set()
 for (i,u),(j,v) in zip(template,template[1:]):
  if i==j:need((v-u)%3==1,'padding arc');continue
  key=tuple(sorted((i,j)));need(key not in boundary_seen,'unique block boundary');boundary_seen.add(key)
  E=constraints[key];D={(d+v-u)%3 for d in E}
  for x,y in it.product(range(3),repeat=2):
   bit=int((y-x)%3 in D);a[3*i+x][3*j+y]=bit;a[3*j+y][3*i+x]=1-bit
 need(boundary_seen==set(constraints),'all constraints represented')
 return owner,constraints,template,a,len(odd)//2,root

def native_relations(cells,template,a):
 relations={}
 for (i,u),(j,v) in zip(template,template[1:]):
  R={(x,y) for x,y in it.product(range(3),repeat=2) if a[3*i+(u+x)%3][3*j+(v+y)%3]}
  if i==j:need(all((x,x) in R for x in range(3)),'internal all rotations');continue
  if i>j:i,j=j,i;R={(y,x) for x,y in R}
  relations[i,j]=R
 return relations

def count_csp(cells,relations,pin=False):
 # Generic arc-consistent binary finite-domain solver, not original-graph coloring.
 rel={e:R for e,R in relations.items() if len(R)<9}
 active={x for e in rel for x in e};free=set(range(cells))-active
 factor=3**len(free);start=[7]*cells
 for i in free:start[i]=1
 if pin:start[0]=1
 @functools.lru_cache(None)
 def solve(state):
  d=list(state);change=True
  while change:
   change=False
   for (i,j),R in rel.items():
    di=sum(1<<x for x in range(3) if d[i]>>x&1 and any(d[j]>>y&1 and (x,y) in R for y in range(3)))
    dj=sum(1<<y for y in range(3) if d[j]>>y&1 and any(d[i]>>x&1 and (x,y) in R for x in range(3)))
    if not di or not dj:return 0
    if (di,dj)!=(d[i],d[j]):d[i],d[j]=di,dj;change=True
  choices=[i for i in active if d[i].bit_count()>1]
  if not choices:return 1
  i=min(choices,key=lambda i:(d[i].bit_count(),i));answer=0
  for x in range(3):
   if d[i]>>x&1:
    dd=d[:];dd[i]=1<<x;answer+=solve(tuple(dd))
  return answer
 return factor*solve(tuple(start))

cases=[]
for n in range(2,5):
 E=list(it.combinations(range(n),2))
 for mask in range(1<<len(E)):
  es=[e for k,e in enumerate(E) if mask>>k&1]
  if connected(n,es):cases.append((n,es))
rng=random.Random(430003)
for _ in range(12):
 while True:
  es=[e for e in it.combinations(range(5),2) if rng.randrange(2)]
  if connected(5,es):cases.append((5,es));break
cases.append((5,list(it.combinations(range(5),2))))
reports=[]
for n,es in cases:
 colors=[c for c in it.product(range(3),repeat=n) if all(c[u]!=c[v] for u,v in es)]
 owner,constraints,template,a,k,root=reduce_graph(n,es);cells=len(owner)
 relations=native_relations(cells,template,a)
 for (i,j),R in relations.items():need(R=={(x,y) for x,y in it.product(range(3),repeat=2) if (y-x)%3 in constraints[i,j]},'decoded phase constraints')
 actual=count_csp(cells,relations);pinned=count_csp(cells,relations,True)
 need(actual==len(colors)*3**(k+2),'counting lift');need(actual==3*pinned,'root owner pin')
 raw=None
 if cells<=8:
  raw=sum(all(a[3*i+(u+phase[i])%3][3*j+(v+phase[j])%3] for (i,u),(j,v) in zip(template,template[1:])) for phase in it.product(range(3),repeat=cells))
  need(raw==actual,'exhaustive literal phase check')
 witness=None
 if colors:
  phase=[colors[0][u] if u is not None else 0 for u in owner]
  path=[3*i+(u+phase[i])%3 for i,u in template]
  need(len(set(path))==3*cells and all(a[u][v] for u,v in zip(path,path[1:])),'literal witness')
  witness=path
 # Host still has a Hamilton path even when the prescribed family is empty.
 ordinary=[]
 for v in range(3*cells):
  pos=next((i for i,u in enumerate(ordinary) if a[v][u]),len(ordinary));ordinary.insert(pos,v)
 need(all(a[u][v] for u,v in zip(ordinary,ordinary[1:])),'ordinary tournament Hamilton path')
 # Regular doubling with a fixed genuine B suffix and actual owner pins.
 N=len(a);c=2*N;reg=[[0]*(2*N+1) for _ in range(2*N+1)]
 for u in range(N):
  for v in range(N):
   reg[u][v]=reg[N+u][N+v]=a[u][v]
   bit=a[v][u] if u!=v else 0
   reg[u][N+v]=bit;reg[N+v][u]=1-bit
  reg[u][c]=1;reg[c][N+u]=1
 need(all(sum(row)==N for row in reg),'regular doubling degree')
 partition=[[3*i+j for j in range(3)] for i in range(cells)]+[[N+3*i+j for j in range(3)] for i in range(cells)]+[[c]]
 for C in partition:
  for D in partition:need(len({sum(reg[u][v] for v in D) for u in C})==1,'equitable regular partition')
 regular_template=[3*i+u for i,u in template]+[c]+[N+u for u in ordinary]
 pin_positions=[next(pos for pos,v in enumerate(regular_template) if v==N+3*i) for i in range(cells)]
 need(len(set(pin_positions))==cells,'one genuine indexed owner per B cell')
 for trial in range(12):
  phase=[(trial+i*i+trial*i)%3 for i in range(cells)]
  pa=[3*i+(u+phase[i])%3 for i,u in template]
  path=pa+[c]+[N+u for u in ordinary]
  need(all(reg[u][v] for u,v in zip(path,path[1:]))==all(a[u][v] for u,v in zip(pa,pa[1:])),'regular template equivalence')
 reports.append({'original_vertices':n,'original_edges':es,'proper_colorings':len(colors),'cells':cells,'host_vertices':3*cells,'fresh_parity_vertices':k,'family_paths':actual,'root_pinned_paths':pinned,'raw_exhaustive_count':raw,'witness':witness,'ordinary_path_sha256':hashlib.sha256(bytes(ordinary)).hexdigest(),'regular_host_vertices':2*N+1,'regular_score':N,'pinned_B_cells':cells,'regular_equitable_cells':2*cells+1})
print(json.dumps({'status':'PASS','cases':len(reports),'scope':'actual constructed tournaments and prescribed shifted templates; ordinary Hamilton paths retained even for empty families','reports':reports},sort_keys=True,indent=2))
