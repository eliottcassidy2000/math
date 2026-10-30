"""Literal integer minors and complete modular fibers; standard library only."""
from itertools import product, combinations
from math import gcd
from collections import Counter
import json

def require(ok,msg):
 if not ok: raise RuntimeError(msg)
def det(a):
 n=len(a)
 if not n:return 1
 a=[r[:] for r in a];sign=1;prev=1
 for k in range(n-1):
  if a[k][k]==0:
   p=next((i for i in range(k+1,n) if a[i][k]),None)
   if p is None:return 0
   a[k],a[p]=a[p],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):a[i][j]=(a[i][j]*pivot-a[i][k]*a[k][j])//prev
  for i in range(k+1,n):a[i][k]=0
  prev=pivot
 return sign*a[-1][-1]
def divisors(a,n):
 e=len(a);ds=[1]
 for k in range(1,min(e,n)+1):
  g=0
  for rs in combinations(range(e),k):
   for cs in combinations(range(n),k):g=gcd(g,abs(det([[a[i][j] for j in cs] for i in rs])))
  ds.append(g)
 return ds

def cycle_matrix(signs):
 n=len(signs);a=[[int(i==j) for j in range(n)] for i in range(n)]
 for i,s in enumerate(signs):a[i][(i+1)%n]-=s
 return a
cycle_fibers=cycle_matrices=0
for n in range(1,9):
 for signs in product([-1,1],repeat=n):
  eta=1
  for s in signs:eta*=s
  a=cycle_matrix(signs)
  require(det(a)==1-eta,('cycle determinant',signs))
  # Verify inverse formula as an integer matrix identity A sum Q^j=(1-eta)I.
  for col in range(n):
   v=[int(i==col) for i in range(n)];total=[0]*n
   for _ in range(n):
    total=[x+y for x,y in zip(total,v)];v=[signs[i]*v[(i+1)%n] for i in range(n)]
   require([sum(a[i][j]*total[j] for j in range(n)) for i in range(n)]==[(1-eta)*int(i==col) for i in range(n)],'cycle inverse identity')
  if n<=5:
   ds=divisors(a,n);require(ds==[1]*n+[abs(1-eta)],('cycle smith',signs,ds))
  cycle_matrices+=1
  if n>4:continue
  for m in range(2,6):
   counts=Counter(tuple(sum(a[i][j]*x[j] for j in range(n))%m for i in range(n)) for x in product(range(m),repeat=n))
   g=gcd(m,1-eta)
   lam=[1]
   for s in signs[:-1]:lam.append(lam[-1]*s)
   for b in product(range(m),repeat=n):
    expected=g if sum(x*y for x,y in zip(lam,b))%g==0 else 0
    require(counts[b]==expected,('cycle fiber',signs,m,b,counts[b],expected));cycle_fibers+=1

def graph_matrix(n,edges):
 a=[]
 for v,w,s in edges:
  row=[0]*n;row[v]+=1;row[w]-=s;a.append(row)
 return a

def forest(n,edges):
 adj=[[] for _ in range(n)]
 for k,(v,w,s) in enumerate(edges):
  adj[v].append((w,s,-s,k));adj[w].append((v,s,1,k))
 comp=[-1]*n;alpha=[0]*n;paths=[None]*n;components=[]
 for root in range(n):
  if comp[root]>=0:continue
  ci=len(components);verts=[];comp[root]=ci;alpha[root]=1;paths[root]=[0]*len(edges);stack=[root]
  while stack:
   v=stack.pop();verts.append(v)
   for w,s,c,k in adj[v]:
    if comp[w]<0:
     comp[w]=ci;alpha[w]=s*alpha[v];paths[w]=[s*z for z in paths[v]];paths[w][k]+=c;stack.append(w)
  components.append(verts)
 rows=[[] for _ in components]
 for k,(v,w,s) in enumerate(edges):
  coef=alpha[v]-s*alpha[w]
  rhs=[-paths[v][j]+s*paths[w][j]+int(j==k) for j in range(len(edges))]
  rows[comp[v]].append((coef,rhs))
 return components,rows

def test_graph(n,edges,moduli):
 global graph_fibers,graph_matrices
 a=graph_matrix(n,edges);components,rows=forest(n,edges)
 o=sum(any(coef for coef,rhs in block) for block in rows);b=len(components)-o
 units=n-len(components);rank=units+o
 ds=divisors(a,n)
 expected=[1]
 for k in range(1,min(len(edges),n)+1):expected.append(0 if k>rank else 2**max(0,k-units))
 require(ds==expected,('graph smith',n,edges,ds,expected));graph_matrices+=1
 for m in moduli:
  counts=Counter(tuple(sum(row[j]*x[j] for j in range(n))%m for row in a) for x in product(range(m),repeat=n))
  expected_count=m**b*gcd(m,2)**o
  for rhs in product(range(m),repeat=len(edges)):
   ok=True
   for block in rows:
    value=None
    for coef,vec in block:
     v=sum(x*y for x,y in zip(vec,rhs))%m
     if coef==0:
      if v:ok=False;break
     else:
      norm=(v if coef==2 else -v)%m
      if norm%gcd(m,2):ok=False;break
      if value is not None and value!=norm:ok=False;break
      value=norm
    if not ok:break
   require(counts[rhs]==(expected_count if ok else 0),('graph fibers',n,edges,m,rhs,counts[rhs],ok));graph_fibers+=1
  require(len(counts)*expected_count==m**n,'kernel image count')
  require(m**len(edges)//len(counts)==m**(len(edges)-n+b)*gcd(m,2)**o,'cokernel count')

graph_fibers=graph_matrices=0
for n in range(1,5):
 pairs=list(combinations(range(n),2))
 for pattern in product([0,-1,1],repeat=len(pairs)):
  edges=[(v,w,s) for (v,w),s in zip(pairs,pattern) if s]
  test_graph(n,edges,[2,3])
# Loops, parallel-edge conflicts and isolated coordinates are genuine boundaries.
for s in [-1,1]:test_graph(1,[(0,0,s)],[2,3,4,5])
for s,t in product([-1,1],repeat=2):test_graph(2,[(0,1,s),(0,1,t)],[2,3,4,5])

# Multiplication orbit formula, all unit actions including nonunit residue strata.
def phi(n):return sum(gcd(x,n)==1 for x in range(n)) if n>1 else 1
def order(u,n):
 if n==1:return 1
 x=u%n;h=1
 while x!=1:x=x*u%n;h+=1
 return h
arith=0
for n in range(1,81):
 for u in range(n):
  if gcd(u,n)!=1:continue
  seen=set();lens=[]
  for x in range(n):
   if x in seen:continue
   y=x;L=0
   while y not in seen:seen.add(y);L+=1;y=u*y%n
   lens.append(L)
  expected=sum(phi(d)//order(u,d) for d in range(1,n+1) if n%d==0)
  require(len(lens)==expected,'arithmetic cycle count')
  require(all(L%2 for L in lens)==(order(u,n)%2==1),'arithmetic odd-order criterion')
  arith+=1
print(json.dumps({'status':'PASS','cycle_matrices':cycle_matrices,'cycle_modular_rhs':cycle_fibers,'signed_graph_matrices':graph_matrices,'graph_modular_rhs':graph_fibers,'arithmetic_unit_actions':arith,'scope':'literal minors, complete modular fibers in stated finite universes; general proofs separate'},indent=2,sort_keys=True))
