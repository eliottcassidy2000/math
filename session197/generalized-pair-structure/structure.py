"""A local graph criterion and an explicit N/M/Q compiler for Gamma."""
from core import inventories,cycles,require,adjacency_map,glue
from pair_families import construct

def induced(T,vs):return [sum(1<<j for j,v in enumerate(vs)if T[u]>>v&1)for u in vs]
def converse(T):return [sum(1<<j for j in range(len(T))if T[j]>>i&1)for i in range(len(T))]
def components(T,vs):
 A=induced(T,vs);n=len(vs);reach=[x|1<<i for i,x in enumerate(A)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 parts=[];left=set(range(n))
 while left:
  i=min(left);p={j for j in left if reach[i]>>j&1 and reach[j]>>i&1};left-=p;parts.append(sorted(p))
 original=tuple(parts)
 parts.sort(key=lambda p:-sum(bool(A[p[0]]>>q[0]&1)for q in original if p!=q))
 return [[vs[i]for i in p]for p in parts]

def root_records(T):
 I=inventories(T);out=[]
 for r,s in enumerate(I['s']):
  if s!=2:continue
  blocks=components(T,[v for v in range(len(T))if v!=r]);U=blocks[0]
  out.append(dict(root=r,initial=U,initial_cycles=cycles(induced(T,U)),tail=blocks[1:]))
 return out

def compile_root(T,r):
 """Return a literal source/family/map, or None if r fails the local criterion."""
 I=inventories(T)
 if cycles(T)!=2 or I['s'][r]!=2:return None
 blocks=components(T,[v for v in range(len(T))if v!=r]);U=blocks[0]
 if cycles(induced(T,U))!=1:return None
 require(all(len(p)==1 for p in blocks[1:]),'start-two tail must be singleton SCCs')
 W=[p[0]for p in blocks[1:]];v=1+sum((1<<i)for i,x in enumerate(W)if T[x]>>r&1)
 E=[x for x in U if T[r]>>x&1];degree={x:sum(bool(T[x]>>y&1)for y in U)for x in U}
 if len(E)==2:
  adjacent=[(a,b)for a in E for b in E if a!=b and degree[a]==1 and T[a]>>b&1]
  if adjacent:
   a,b=adjacent[0];V=[x for x in U if x!=a];kind='N';mapping=V+[r,a]+W
   require(v>=2,'empty adjacent pair gives C1, not C2')
  else:
   choices=[(a,c,b)for a in E for c in E if a!=c for b in U if b not in E and degree[a]==degree[b]==1 and T[a]>>b&1 and T[b]>>c&1]
   require(len(choices)==1,'nonadjacent three-block')
   a,c,b=choices[0];V=[x for x in U if x not in(a,b)];kind='M';mapping=V+[b,r,a]+W;b=c
 elif len(E)==1:
  a=E[0]
  choices=[(y,b)for y in U if y!=a and degree[y]==1 for b in U if b not in(a,y)and T[y]>>b&1 and degree[a]==2 and T[a]>>y&1 and T[a]>>b&1]
  require(len(choices)==1,'first clone of sole-out-neighbor vertex')
  y,b=choices[0];V=[x for x in U if x not in(a,y)];kind='Q';mapping=V+[a,y,r]+W
 else:raise RuntimeError('start-two entry sum')
 source=induced(T,V);root=V.index(b)
 require(len(source)>=2 and inventories(source)['s'][root]==1,'extracted physical source')
 A,q=construct(kind,source,root,v)
 require(mapping[q]==r and adjacency_map(A,T,mapping),'exact full graph/root mapping')
 return dict(family=kind,source=source,source_root=root,parameter=v,vertex_map=mapping,root=r)

def compile_graph(T):
 if cycles(T)!=2:return None
 for r,s in enumerate(inventories(T)['s']):
  if s==2:
   out=compile_root(T,r)
   if out is not None:return out
 return None

def redundant_P(V,b,v):
 A,r=construct('P',V,b,v);J,a=glue(V,{b},[0]);B,q=construct('Q',J,a,v)
 require(A==B and r==q,'P equals Q after J, without relabeling')
 return True
