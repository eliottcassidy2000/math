#!/usr/bin/env python3
"""Independent exact Fourier/forest/square controls; no imported producer code."""
from itertools import combinations, permutations
from collections import Counter
from math import factorial, comb
from pathlib import Path
import hashlib,json
checks=0
def require(b,msg):
 global checks
 checks+=1
 if not b: raise RuntimeError(msg)
def walsh(values):
 a=list(values); step=1
 while step<len(a):
  for base in range(0,len(a),2*step):
   for i in range(base,base+step):
    x,y=a[i],a[i+step];a[i],a[i+step]=x+y,x-y
  step*=2
 return a
def linear_forest(n,edges,mask):
 degree=[0]*n;parent=list(range(n));nontrivial=set()
 def find(x):
  while parent[x]!=x: x=parent[x]
  return x
 for i,(u,v) in enumerate(edges):
  if mask>>i&1:
   degree[u]+=1;degree[v]+=1
   if degree[u]>2 or degree[v]>2: return None
   a,b=find(u),find(v)
   if a==b: return None
   parent[a]=b;nontrivial.update((u,v))
 return len({find(v) for v in nontrivial})
def forest_sum(n,edges,Z):
 total=0;sub=Z
 while True:
  c=linear_forest(n,edges,sub)
  if c is not None:
   r=sub.bit_count(); total+=((-2)**r)*(2**c)*factorial(n-r)
  if sub==0:break
  sub=(sub-1)&Z
 return total
receipt={'status':'PASS','orders':[]}
for n in range(1,7):
 E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
 F=tuple(combinations(range(1,n),2));full=(1<<len(E))-1
 def bit(u,v):return 1<<idx[tuple(sorted((u,v)))]
 def project(x):return sum(1<<k for k,(i,j) in enumerate(F) if (bool(x&bit(i,j))^bool(x&bit(0,i))^bool(x&bit(0,j))))
 def eulerian(h):
  Z=0
  for k,(i,j) in enumerate(F):
   if h>>k&1: Z ^=bit(i,j)^bit(0,i)^bit(0,j)
  return Z
 path_words=[]; kernel=[0]*(1<<len(F))
 for p in permutations(range(n)):
  path=sum(bit(u,v) for u,v in zip(p,p[1:]))
  path_words.append(path)
  kernel[project(full^path)]+=1
 eigen=walsh(kernel)
 require(sum(kernel)==factorial(n),'aggregate row sum')
 require(sum(v*v for v in eigen)==len(kernel)*sum(w*w for w in kernel),'Fourier second moment')
 for h,lam in enumerate(eigen):
  Z=eulerian(h)
  require(all(sum(bool(Z&bit(v,w)) for w in range(n) if w!=v)%2==0 for v in range(n)),'character not Eulerian')
  direct=sum((-1)**((Z& (full^path)).bit_count()) for path in path_words)
  signed=sum((-1)**((Z&path).bit_count()) for path in path_words)
  require(lam==direct==(-1)**Z.bit_count()*signed,'signed path Fourier formula')
  if n<=5: require(signed==forest_sum(n,E,Z),'linear-forest expansion')
 generators=[g for g,w in enumerate(kernel) if w]
 if n<=2:
  require(kernel==[factorial(n)],'singleton loop boundary')
 else:
  require(kernel[0]==0 and sum(eigen)==0,'loopless first spectral moment')
  require(len(set(kernel[g] for g in generators))==1,'uniform generator multiplicity')
  w=kernel[generators[0]]
  require(sum(x*x for x in eigen)==len(kernel)*factorial(n)*w,'uniform second moment')
  require(eigen.count(factorial(n))==(2 if n%4==1 else 1),'top multiplicity')
  require(eigen.count(-factorial(n))==(1 if n%4==3 else 0),'bottom multiplicity')
  if n==3:
   require(generators==[1] and w==6,'K2 support bridge with six parallel incidences')
  else:
   require(len(generators)>=2,'insufficient generators for square')
   for g in generators:
    h=next(h for h in generators if h!=g)
    cycle=(0,g,g^h,h)
    require(len(set(cycle))==4,'square vertices collide')
    require(all(kernel[cycle[i]^cycle[(i+1)%4]]>0 for i in range(4)),'square missing edge')
  # A transitive tournament has exactly one directed Hamiltonian word.
  directed_paths=[p for p in permutations(range(n)) if all(u<v for u,v in zip(p,p[1:]))]
  require(len(directed_paths)==1,'transitive source not a leaf')
  unique_path=sum(bit(u,v) for u,v in zip(directed_paths[0],directed_paths[0][1:]))
  require((full^unique_path)!=0,'original leaf incident edge is a loop')
 receipt['orders'].append({'n':n,'quotient_vertices':len(kernel),'generator_count':len(generators),
  'positive_weight':kernel[generators[0]],'spectrum':dict(sorted(Counter(eigen).items()))})
# Independent signed Hamilton-path subset DP for the triangle eigenvalue.
def signed_path_dp(n,Z):
 E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
 sign=[[0 if i==j else (-1 if Z>>idx[tuple(sorted((i,j)))]&1 else 1) for j in range(n)] for i in range(n)]
 dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not(mask>>v&1):continue
   if mask==(1<<v):continue
   prior=mask^(1<<v)
   dp[mask][v]=sum(dp[prior][u]*sign[u][v] for u in range(n) if prior>>u&1)
 return sum(dp[-1])
triangle_rows=[]
for n in range(3,12):
 E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)}
 Z=sum(1<<idx[e] for e in combinations(range(3),2))
 signed=signed_path_dp(n,Z)
 formula=-factorial(n-2)*(n-4)*(n-9)
 require(-signed==formula,'triangle eigenvalue subset-DP control')
 require(signed==factorial(n)-12*factorial(n-1)+24*factorial(n-2),'triangle forest polynomial')
 triangle_rows.append({'n':n,'triangle_eigenvalue':formula})
# Character count, not a giant order-nine quotient census.
n=9;E=tuple(combinations(range(n),2));idx={e:i for i,e in enumerate(E)};full=(1<<len(E))-1
triangles={sum(1<<idx[e] for e in combinations(vs,2)) for vs in combinations(range(n),3)}
complements={full^Z for Z in triangles}
require(len(triangles)==84 and len(triangles|complements)==168,'order-nine zero-character multiplicity')
require(all(sum(bool(Z>>idx[tuple(sorted((v,u)))]&1) for u in range(n) if u!=v)%2==0 for v in range(n) for Z in triangles|complements),'order-nine zero characters not Eulerian')
receipt['triangle_modes']=triangle_rows
receipt['n9_nullity_lower_bound_full']=168
receipt['n9_nullity_lower_bound_each_component']=84
receipt['checks']=checks
receipt['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print(json.dumps(receipt,indent=2,sort_keys=True))
