#!/usr/bin/env python3
"""Exact polynomial and literal-subspace controls; no external packages."""
from itertools import combinations,product
from functools import lru_cache
from math import comb,factorial
from collections import Counter
import json
CHECKS=0

def check(ok,ctx):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(ctx)

def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 while len(c)>1 and c[-1]==0:c.pop()
 return tuple(c)

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 while len(c)>1 and c[-1]==0:c.pop()
 return tuple(c)

def shift(a,k):return (0,)*k+tuple(a)
def evalp(a,q):return sum(x*q**i for i,x in enumerate(a))
def qpowerplus(k,s=1):return (s,)+(0,)*(k-1)+(1,)
@lru_cache(None)
def gauss(n,k):
 if k<0 or k>n:return (0,)
 if k in (0,n):return (1,)
 return add(gauss(n-1,k),shift(gauss(n-1,k-1),n-k))

def iso(n,k):
 p=gauss(n,k)
 for j in range(n-k+1,n+1):p=mul(p,qpowerplus(j))
 return p

def lag(n):return iso(n,n)
def sp(n):
 p=(1,)
 for j in range(1,n+1):p=mul(p,qpowerplus(2*j,-1))
 return shift(p,n*n)

def taylor_at_one(p):
 return tuple(sum(p[i]*comb(i,j) for i in range(j,len(p))) for j in range(len(p)))

def rank(rows,p):
 a=[list(r) for r in rows];h=0
 if not a:return 0
 for j in range(len(a[0])):
  i=next((i for i in range(h,len(a)) if a[i][j]%p),None)
  if i is None:continue
  a[h],a[i]=a[i],a[h];v=pow(a[h][j],-1,p)
  a[h]=[(x*v)%p for x in a[h]]
  for i in range(len(a)):
   if i!=h:
    v=a[i][j]%p;a[i]=[(x-v*y)%p for x,y in zip(a[i],a[h])]
  h+=1
  if h==len(a):break
 return h

def subspaces(N,k,p):
 for pivots in combinations(range(N),k):
  slots=[(i,j) for i,pi in enumerate(pivots) for j in range(pi+1,N) if j not in pivots]
  for values in product(range(p),repeat=len(slots)):
   a=[[int(j==pi) for j in range(N)] for pi in pivots]
   for (i,j),v in zip(slots,values):a[i][j]=v
   yield a

def symp(x,y,p):return sum(x[i]*y[i+1]-x[i+1]*y[i] for i in range(0,len(x),2))%p

def isotropic(a,p):return all(symp(x,y,p)==0 for x in a for y in a)

for m in range(7):
 for n in range(7):
  total=(0,)
  for r in range(min(m,n)+1):
   sector=mul(mul(iso(m,m-r),iso(n,n-r)),sp(r));total=add(total,sector)
   t=taylor_at_one(sector)
   check(all(x==0 for x in t[:r]),('vanishing order',m,n,r))
   expected=2**(m+n-r)*factorial(r)*comb(m,r)*comb(n,r)
   check(t[r]==expected,('first coefficient',m,n,r,t[r],expected))
   check(all(x>=0 for x in t),('positive t expansion',m,n,r))
  check(total==lag(m+n),('composition polynomial',m,n))
  check(evalp(lag(m),1)==2**m,('thin Lagrangians',m))
  flags=(1,)
  for i in range(1,m+1):flags=mul(flags,(1,)*(2*i))
  check(evalp(flags,1)==2**m*factorial(m),('Weyl count',m))

universes=[]
for p,n in [(2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(5,2)]:
 counts=[];visited=0
 for k in range(n+1):
  got=0
  for a in subspaces(2*n,k,p):
   visited+=1
   if isotropic(a,p):got+=1
  check(got==evalp(iso(n,k),p),('literal isotropic',p,n,k,got))
  counts.append(got)
 universes.append({'prime':p,'half_dimension':n,'subspaces_examined':visited,'isotropic_counts':counts})

coupled=[]
for p,m,n in [(2,1,1),(3,1,1),(5,1,1),(2,1,2),(3,1,2)]:
 k=m+n;counts=Counter()
 for a in subspaces(2*k,k,p):
  if not isotropic(a,p):continue
  dimA=k-rank([row[2*m:] for row in a],p)
  dimB=k-rank([row[:2*m] for row in a],p)
  r=m-dimA
  check(r==n-dimB,('balanced quotient',p,m,n,a))
  counts[r]+=1
 for r in range(min(m,n)+1):
  expected=evalp(iso(m,m-r),p)*evalp(iso(n,n-r),p)*evalp(sp(r),p)
  check(counts[r]==expected,('literal sector',p,m,n,r,counts[r],expected))
 coupled.append({'prime':p,'halves':[m,n],'rank_counts':dict(counts)})

# Both pairings take values 0 or 1/2 modulo Z; matrix entries are doubled.
H=((0,1),(1,0));I=((1,0),(0,1));vectors=list(product(range(2),repeat=2))
def pair(x,y,M):return sum(x[i]*M[i][j]*y[j] for i in range(2) for j in range(2))%2
for M in (H,I):
 check(all(pair(x,y,M)==pair(y,x,M) for x in vectors for y in vectors),'antisymmetric mod Q/Z')
 check([x for x in vectors if all(pair(x,y,M)==0 for y in vectors)]==[(0,0)],'nondegenerate')
check(sum(pair(x,x,H)==0 for x in vectors)==4,'hyperbolic self-isotropy')
check(sum(pair(x,x,I)==0 for x in vectors)==2,'nonalternating self-isotropy')
check([x for x in vectors if all(pair(y,y,I)==pair(y,x,I) for y in vectors)]==[(1,1)],'characteristic vector')
check(pair((1,1),(1,1),I)==0,'square but nonalternating')
print(json.dumps({'status':'PASS','checks':CHECKS,'literal_subspaces':universes,'coupled_sectors':coupled,'same_order_pairing_hostile':{'group_order':4,'alternating_isotropic_elements':4,'nonalternating_isotropic_elements':2}},sort_keys=True))
