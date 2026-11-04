"""Subset DP and attachment helpers inherited from S185 closure/verify.py.
Source: eliottcassidy2000/math at 8880585ffd98072e93aa39d97476c7a1d4c9129d.
Only the helper prefix is used; no source discovery code or startup is run.
"""
"""Literal strong-suffix escape certificate; Python standard library only."""
import json,itertools
from pathlib import Path

def need(x,msg):
 if not x:raise RuntimeError(msg)
def tables(r):
 n=len(r);N=1<<n;h=[0]*N;h[0]=1;e=[[0]*n for _ in range(N)]
 for v in range(n):e[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:e[S][v]=sum(e[T][u] for u in range(n) if r[u]>>v&1)
  h[S]=sum(e[S])
 return h,e
def direct(r,W):
 n=len(r);h,e=tables(r);rev=[sum(1<<j for j in range(n) if r[j]>>i&1) for i in range(n)];_,s=tables(rev);N=(1<<n)-1
 ep=[sum(x[i] for i in W) for x in e];sp=[sum(x[i] for i in W) for x in s]
 return(h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(ep[x]*h[N^x] for x in range(N+1)),sum(h[x]*sp[N^x] for x in range(N+1)),sum(ep[x]*sp[N^x] for x in range(N+1)),ep[N],sp[N])
def strong(r):
 N=(1<<len(r))-1
 for v in range(len(r)):
  S=1<<v
  while True:
   T=S
   for i in range(len(r)):
    if S>>i&1:T|=r[i]
   if T==S:break
   S=T
  if S!=N:return False
 return True
def tournament(r):
 n=len(r)
 return all(not(r[i]>>i&1) and r[i]<(1<<n) for i in range(n)) and all(bool(r[i]>>j&1)!=bool(r[j]>>i&1) for i in range(n) for j in range(i))
def attach(r,W,t):
 n=len(r);root=t+n
 rows=[sum(1<<j for j in range(i+1,t+n)) for i in range(t)]
 rows += [(x<<t)|((1<<root) if i in W else 0) for i,x in enumerate(r)]
 rows += [((1<<t)-1)|sum(1<<(t+i) for i in range(n) if i not in W)]
 return rows


