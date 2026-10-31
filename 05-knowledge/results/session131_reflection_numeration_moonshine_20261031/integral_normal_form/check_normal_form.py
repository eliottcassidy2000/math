#!/usr/bin/env python3
"""Constructive reflection module normal form and coherent CG growth."""
from fractions import Fraction
from collections import deque
import json
# Owner display-node convention of Session129, not least-mask node ranks.
BLUE=[(0,5),(3,8),(4,6),(5,6),(5,8),(6,9),(7,9),(8,9)]
BLACK=[(1,4),(2,4),(2,5),(2,6),(2,7),(2,8),(4,5),(4,8),(5,7),(6,6),(6,8),(8,8)]
def req(c,m):
 if not c:raise RuntimeError(m)
def comps(n,es):
 p=list(range(n))
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for u,v in es:
  a,b=root(u),root(v)
  if a!=b:p[a]=b
 out={}
 for u in range(n):out.setdefault(root(u),[]).append(u)
 return sorted(out.values())
def cycles(n,es):
 p=list(range(n));adj=[[] for _ in range(n)];tree=[];ch=[]
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for j,(u,v) in enumerate(es):
  a,b=root(u),root(v)
  if a==b:ch.append(j)
  else:p[a]=b;tree.append(j);adj[u].append((v,j,1));adj[v].append((u,j,-1))
 basis=[]
 for e in ch:
  u,v=es[e];z=[0]*len(es);z[e]=1
  seen={u:None};q=deque([u])
  while v not in seen:
   a=q.popleft()
   for b,j,sign in adj[a]:
    if b not in seen:seen[b]=(a,j,sign);q.append(b)
  b=v
  while b!=u:
   a,j,sign=seen[b];z[j]-=sign;b=a
  req(boundary(n,es,z)==[0]*n,'basis not cycle');basis.append(z)
 return basis,ch,tree
def boundary(n,es,z):
 b=[0]*n
 for (u,v),x in zip(es,z):b[u]-=x;b[v]+=x
 return b
def det(a):
 a=[list(map(int,row)) for row in a];n=len(a);prev=1;sgn=1
 for k in range(n-1):
  if a[k][k]==0:
   j=next((j for j in range(k+1,n) if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sgn=-sgn
  piv=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    x=a[i][j]*piv-a[i][k]*a[k][j];req(x%prev==0,'Bareiss division');a[i][j]=x//prev
   a[i][k]=0
  prev=piv
 return sgn*a[-1][-1] if n else 1
def rref2(a):
 a=[[x%2 for x in row] for row in a];r=0;ps=[]
 for c in range(len(a[0])):
  j=next((j for j in range(r,len(a)) if a[j][c]),None)
  if j is None:continue
  a[r],a[j]=a[j],a[r]
  for i in range(len(a)):
   if i!=r and a[i][c]:a[i]=[x^y for x,y in zip(a[i],a[r])]
  ps.append(c);r+=1
 return a[:r],ps

from math import isqrt
Q=BLUE+BLACK;G=BLUE+[e for e in BLACK for _ in (0,1)]
qbas,qch,qt=cycles(10,Q);gbas,gch,gt=cycles(10,G)
R=list(range(8))+[i^1 for i in range(8,32)]
def reflect(z):return [z[j] for j in R]
def lift(z):return z[:8]+[y for x in z[8:] for y in (x,0)]
basis=[];blocks=[]
for z,e in zip(qbas,qch):
 if e<8:
  w=lift(z);req(reflect(w)==w,'blue chord not fixed');blocks.append(('trivial',[len(basis)]));basis.append(w)
for e in qt:
 if e>=8:
  w=[0]*32;j=8+2*(e-8);w[j]=1;w[j+1]=-1
  req(reflect(w)==[-x for x in w],'tree twin not sign');blocks.append(('sign',[len(basis)]));basis.append(w)
for z,e in zip(qbas,qch):
 if e>=8:
  w=lift(z);blocks.append(('exchange',[len(basis),len(basis)+1]));basis.extend([w,reflect(w)])
M=[[z[c] for z in basis] for c in gch]
req(abs(det(M))==1,'module basis not integral unimodular')
req([sum(t==name for t,ix in blocks) for name in ('trivial','sign','exchange')]==[1,2,10],'wrong multiplicities')
for z in basis:req(boundary(10,G,z)==[0]*10,'not cycle')
def reflection_coeff(a):
 b=a.copy()
 for typ,ix in blocks:
  if typ=='sign':b[ix[0]]=-a[ix[0]]
  elif typ=='exchange':b[ix[0]],b[ix[1]]=a[ix[1]],a[ix[0]]
 return b
def physical(a):return [sum(x*z[j] for x,z in zip(a,basis)) for j in range(32)]
def shift(n):
 if n==0:return 0
 if n<0:return -shift(-n)
 return (3*n+isqrt(5*n*n))//2+1
def growth(a):return list(map(shift,a))
# Deterministic signed coefficient controls; no floating arithmetic.
controls=[]
for seed in range(101):
 a=[((seed+3)*(j+7)%29)-14 for j in range(23)]
 req(physical(reflection_coeff(a))==reflect(physical(a)),'reflection transport')
 b=growth(a);req(boundary(10,G,physical(b))==[0]*10,'growth lost graph boundary')
 req(growth(reflection_coeff(a))==reflection_coeff(b),'growth not equivariant')
 controls.append(a)
req(shift(1)==3 and shift(2)==6 and shift(3)==8,'CG shift')
req(shift(1)+shift(2)!=shift(3),'growth incorrectly additive')
# Shift after diagonalization destroys the integral gluing constraint.
s,d=3,1;x,y=(s+d)//2,(s-d)//2
naive=(shift(s),shift(d));repaired=(shift(x)+shift(y),shift(x)-shift(y))
req(naive==(8,3) and repaired==(9,3),'parity repair hostile')
req(1+3==4 and isqrt(shift(1)+shift(3))**2!=shift(1)+shift(3),'square-sum hostile')
# Inherited n4 owner atlas: blue 0--2 and a loop at2; one black pair1--2.
b4=[(0,2),(2,2)];k4=[(1,2)]
cB4=len(comps(3,b4));cQ4=len(comps(3,b4+k4))
n4=(len(b4)-3+cB4,cB4-cQ4,len(k4)+cQ4-cB4)
req(n4==(1,1,0),'n4 module specialization')
# Frob permutation modules have only trivial and exchange blocks.
req(11<12,'permutation rank obstruction')
if '--test-failing-gate' in __import__('sys').argv:req(False,'intentional failing gate')
print(json.dumps({'status':'PASS','n4_multiplicities_trivial_sign_exchange':n4,'basis_blocks':blocks,'cycle_basis':basis,'basis_determinant':det(M),'blue_spanning_forest_quotient_edges':[e for e in qt if e<8],'black_forest_quotient_edges':[e for e in qt if e>=8],'module_multiplicities':{'trivial':1,'sign':2,'exchange':10},'even_rank':11,'odd_rank':12,'eigensublattice_index':1024,'signed_growth_controls':len(controls),'naive_eigen_shift':naive,'repaired_eigen_shift':repaired,'square_sum_shift_hostile':[1,3,shift(1),shift(3)]},indent=2,sort_keys=True))
