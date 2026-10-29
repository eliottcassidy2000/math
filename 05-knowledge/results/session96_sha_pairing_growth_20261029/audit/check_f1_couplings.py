#!/usr/bin/env python3
"""Independent polynomial and literal symplectic controls; no producer imports."""
from itertools import product
from math import comb,factorial
from functools import lru_cache
import json
N=0

def ck(v,c):
 global N
 N+=1
 if not v:raise RuntimeError(c)

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def shift(a,k):return [0]*k+a
@lru_cache(None)
def qb(n,k):
 if k<0 or k>n:return (0,)
 if k==0 or k==n:return (1,)
 return tuple(add(list(qb(n-1,k)),shift(list(qb(n-1,k-1)),n-k)))
def L(n):
 a=[1]
 for i in range(1,n+1):a=mul(a,[1]+[0]*(i-1)+[1])
 return a
def isotropic(n,k):
 a=list(qb(n,k))
 for i in range(n-k+1,n+1):a=mul(a,[1]+[0]*(i-1)+[1])
 return a
def sp(r):
 a=shift([1],r*r)
 for i in range(1,r+1):a=mul(a,[-1]+[0]*(2*i-1)+[1])
 return a
def at_one(a):return [sum(a[i]*comb(i,j) for i in range(j,len(a))) for j in range(len(a))]
for m,n in product(range(5),repeat=2):
 total=[0]
 for r in range(min(m,n)+1):
  c=mul(mul(isotropic(m,m-r),isotropic(n,n-r)),sp(r));total=add(total,c);t=at_one(c)
  ck(all(x>=0 for x in t),('positive coefficients',m,n,r))
  ck(all(x==0 for x in t[:r]) and t[r]==2**(m+n-r)*factorial(r)*comb(m,r)*comb(n,r),('exact order',m,n,r))
 ck(total==L(m+n),('composition polynomial',m,n))

def rref(rows,p):
 A=[list(row) for row in rows];rank=0
 for j in range(4):
  i=next((i for i in range(rank,len(A)) if A[i][j]%p),None)
  if i is None:continue
  A[i],A[rank]=A[rank],A[i];a=pow(A[rank][j],-1,p);A[rank]=[a*x%p for x in A[rank]]
  for i in range(len(A)):
   if i!=rank:
    a=A[i][j];A[i]=[(x-a*y)%p for x,y in zip(A[i],A[rank])]
  rank+=1
 return tuple(tuple(row) for row in A[:rank])
def omega(v,w,p):return (v[0]*w[1]-v[1]*w[0]+v[2]*w[3]-v[3]*w[2])%p
rows=[]
for p in (2,3):
 vectors=list(product(range(p),repeat=4));planes=set()
 for v in vectors[1:]:
  for w in vectors[1:]:
   B=rref([v,w],p)
   if len(B)==2 and omega(v,w,p)==0:planes.add(B)
 split=mixed=0
 for B in planes:
  span={tuple((a*B[0][i]+b*B[1][i])%p for i in range(4)) for a,b in product(range(p),repeat=2)}
  lu=sum(x[2:]==(0,0) for x in span);lw=sum(x[:2]==(0,0) for x in span)
  ck((lu,lw) in ((p,p),(1,1)),('coupling classification',p,B))
  if lu==p:split+=1
  else:mixed+=1
 ck(split==(p+1)**2 and mixed==p*(p*p-1),'literal sectors')
 rows.append({'field_order':p,'lagrangians':len(planes),'separable':split,'mixed':mixed})
print(json.dumps({'status':'PASS','checks':N,'literal_four_dimensional_spaces':rows},sort_keys=True))
