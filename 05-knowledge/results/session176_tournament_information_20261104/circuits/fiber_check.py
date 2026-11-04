#!/usr/bin/env python3
"""Literal independent cycle-orientation choices; exact endpoint tables."""
import json,itertools
from explore import kernel

def require(x,m):
 if not x:raise RuntimeError(m)

def build(k,bits):
 a=[]
 for i in range(k):
  c=[4,1,2] if bits>>i&1 else [2,4,1]
  later=((1<<(3*k))-1)^((1<<(3*i+3))-1)
  a.extend((r<<(3*i))|later for r in c)
 return a

def signature(a):
 H,K=kernel(a);n=len(a);w={(u,v):(K[u][v]+K[v][u])//2 for u,v in itertools.combinations(range(n),2)}
 d=tuple(sum(w[min(u,v),max(u,v)] for v in range(n)if a[u]>>v&1)for u in range(n))
 return H,tuple(w.values()),d,K

out=[]
for k in range(1,5):
 old=None
 for bits in range(1<<k):
  a=build(k,bits);s=signature(a)
  if old is None:old=s[:3]
  require(s[:3]==old,'response collision')
  H,_,_,K=s;require(H==3**k,'H')
  for u,v in itertools.combinations(range(3*k),2):
   i,j=u//3,v//3
   if i==j:require(K[u][v]+K[v][u]==2*3**(k-1),'internal capacity')
   else:
    t=j-i-1
    require(K[u][v]==(3**(k-2) if t==0 else 0),'forward cover')
    require(K[v][u]==3**k*4**t,'reverse cover')
 out.append({'blocks':k,'vertices':3*k,'distinct_tournaments_with_identical_full_response':1<<k,'H':3**k})
print(json.dumps({'status':'PASS','fibers':out},indent=2,sort_keys=True))
