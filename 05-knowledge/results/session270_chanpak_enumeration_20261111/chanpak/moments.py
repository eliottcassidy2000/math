#!/usr/bin/env python3
"""Exact all-color moment transform; this counts presentations, not classes."""
import json,math
from fractions import Fraction
from itertools import product
from spine import value,graph,count

def covers(core):
 k=len(core);N=1<<k
 H=[0]*N;H[0]=1
 for mask in range(1,N):
  ids=[i for i in range(k) if mask>>i&1]
  rows=[sum(((core[i]>>j)&1)<<p for p,j in enumerate(ids)) for i in ids]
  H[mask]=count(rows)
 old=[0]*N;old[0]=1;ans=[]
 for r in range(1,k+1):
  new=[0]*N
  for mask in range(1,N):
   s=mask
   while s:
    new[mask]+=H[s]*old[mask^s];s=(s-1)&mask
  ans.append(new[-1]);old=new
 return ans

def onto(m,t):
 return sum((-1)**(t-j)*math.comb(t,j)*j**m for j in range(t+1))

def coefficient(k,m,r):
 return sum(Fraction(math.comb(2,b)*onto(m,r-1+b),2**(2*(r-1)+b)) for b in range(3))*2**(k*m)

def mass(core,m):
 q=covers(core);ans=sum(coefficient(len(core),m,r)*x for r,x in enumerate(q,1))
 if ans.denominator!=1:raise RuntimeError('nonintegral mass')
 return int(ans)

def recover(k,A):
 q=[]
 for m,a in enumerate(A):
  x=(a-sum(coefficient(k,m,r)*v for r,v in enumerate(q,1)))/coefficient(k,m,m+1)
  if x.denominator!=1:raise RuntimeError('nonintegral recovery')
  q.append(int(x))
 return q

def tests():
 cores=[(0,),(2,0),(2,4,1),(14,12,8,0),(2,4,9,3)]
 data=[]
 for C in cores:
  k=len(C);checks=[]
  for m in range(4 if k<4 else 2):
   direct=sum(value(w,C) for w in product(range(1<<k),repeat=m));a=mass(C,m)
   if direct!=a:raise RuntimeError((C,m,direct,a))
   checks.append([m,a])
  A=[mass(C,m) for m in range(k)]
  if recover(k,A)!=covers(C):raise RuntimeError('tomography')
  data.append({'core':C,'covers':covers(C),'direct_mass_checks':checks,'recovered_covers':recover(k,A)})
 A=[mass((2,0),m) for m in range(20)]
 for m in range(4,20):
  if A[m]!=24*A[m-1]-176*A[m-2]+384*A[m-3]:raise RuntimeError('recurrence')
 return {'core_controls':data,'two_core_mass':A,'recurrence':[24,-176,384],'recurrence_starts_at':4}
if __name__=='__main__':print(json.dumps(tests(),indent=2))
