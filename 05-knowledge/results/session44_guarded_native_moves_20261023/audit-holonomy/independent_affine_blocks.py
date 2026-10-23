"""Independent exact affine-Fourier interface audit; no producer imports."""
from fractions import Fraction as F
from itertools import product
import json
P=list(product(range(3),repeat=2));I=((1,0),(0,1))
# Reversible mixed-affine chain: symmetric cross identity, shear/inverse at
# vertex0, involutive affine reflection at vertex1, and holding probability1/3.
M=((1,1),(0,1));Mi=((1,2),(0,1));N=((2,0),(0,1))
E={0:[(0,I,(0,0),F(1,3)),(1,I,(0,0),F(1,3)),(0,M,(0,0),F(1,6)),(0,Mi,(0,0),F(1,6))],
   1:[(1,I,(0,0),F(1,3)),(0,I,(0,0),F(1,3)),(1,N,(1,0),F(1,3))]}
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
def conj(x):return (x[0]-x[1],-x[1])
def scale(w,x):return (w*x[0],w*x[1])
def chi(k,a):return [(1,0),(0,1),(-1,-1)][sum(x*y for x,y in zip(k,a))%3]
def act(A,c,a):return tuple((sum(A[i][j]*a[j] for j in range(2))+c[i])%3 for i in range(2))
def dual(A,k):return tuple(sum(A[i][j]*k[i] for i in range(2))%3 for j in range(2))
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
checks=0
adj={k:{dual(A,k) for v in E for _,A,c,w in E[v]} for k in P};orbits=[];todo=set(P)
while todo:
 s=min(todo);orbit={s};stack=[s]
 while stack:
  for t in adj[stack.pop()]:
   if t not in orbit:orbit.add(t);stack.append(t)
 todo-=orbit;orbits.append(orbit)
which={k:i for i,o in enumerate(orbits) for k in o};coeff={}
for u,v,k,l in product(range(2),range(2),P,P):
 literal=(0,0);formula=(0,0)
 for a in P:
  image=(0,0)
  for dest,A,c,w in E[u]:
   if dest==v:image=add(image,scale(w,chi(k,act(A,c,a))))
  literal=add(literal,scale(F(1,9),mul(conj(chi(l,a)),image)))
 for dest,A,c,w in E[u]:
  if dest==v and dual(A,k)==l:formula=add(formula,scale(w,chi(k,c)))
 need(literal==formula,'direct full-state Fourier coefficient versus affine mode transport')
 if which[k]!=which[l]:need(literal==(0,0),'separate frequency orbit blocks')
 coeff[(u,l,v,k)]=literal
for u,l,v,k in coeff:
 need(coeff[(u,l,v,k)]==conj(coeff[(v,k,u,l)]),'Hermitian block coefficients in uniform stationary gauge')
# This demonstrates an interface failure, not a call to a complete-cover compiler.
vertices=[(v,a) for v in range(2) for a in range(3)]
def comp(edges):
 todo=set(vertices);count=0;adj={s:set() for s in vertices}
 for s,t in edges:adj[s].add(t);adj[t].add(s)
 while todo:
  seen={next(iter(todo))};stack=list(seen);count+=1
  while stack:
   for t in adj[stack.pop()]:
    if t not in seen:seen.add(t);stack.append(t)
  todo-=seen
 return count
need(comp([((0,0),(1,0))])==5,'partial transition must not be silently completed')
need(comp([((0,a),(1,a)) for a in range(3)])==3,'full matching has a different component partition')
print(json.dumps({'status':'PASS','checks':checks,'frequency_orbit_sizes':sorted(map(len,orbits)),'coefficient_ring':'Q[w]/(w^2+w+1)','producer_imports':False,'partial_components':5,'completed_components':3},indent=2,sort_keys=True))
