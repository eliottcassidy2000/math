#!/usr/bin/env python3
"""Exact labelled-tournament decision; deliberately exponential. No floating arithmetic."""
import sys
from functools import lru_cache

def f(n):
 q,r=divmod(n,3)
 return (3*5**(q-1),5**q,9*5**(q-1))[r]

def bound(h):
 n=0
 while 3**(n+1)<=h**3:n+=1
 return n

def paths(adj,n):
 dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not (s>>v)&1: continue
   a=dp[s][v]
   for u in range(n):
    if not (s>>u)&1 and (adj[v]>>u)&1:dp[s|1<<u][u]+=a
 return sum(dp[-1])

def strong(adj,n):
 def reach(start):
  seen=1<<start
  while True:
   nex=seen
   for i in range(n):
    if seen>>i&1:nex|=adj[i]
   if nex==seen:return seen==(1<<n)-1
   seen=nex
 return all(reach(i) for i in range(n))

def find_strong(h):
 n=3
 while f(n)<=h:
  pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
  for mask in range(1<<len(pairs)):
   a=[0]*n
   for k,(i,j) in enumerate(pairs):
    if (mask>>k)&1:a[i]|=1<<j
    else:a[j]|=1<<i
   if strong(a,n) and paths(a,n)==h:return (n,mask)
  n+=1
 return None

def divisors(h):return [d for d in range(3,h+1,2) if h%d==0]

def solve(h):
 if h==1:return []
 if h<1 or h%2==0:return None
 # An exact but expensive enumeration; only counts dividing h can occur.
 witnesses={d:find_strong(d) for d in divisors(h)}
 print('strong witnesses by divisor:',witnesses)
 @lru_cache(None)
 def factor(x):
  if x==1:return ()
  for d,w in witnesses.items():
   if w is not None and x%d==0:
    rest=factor(x//d)
    if rest is not None:return ((d,w),)+rest
  return None
 return factor(h)

if __name__=='__main__':
 for h in (map(int,sys.argv[1:]) if len(sys.argv)>1 else (13,21)):
  print('target',h,'global bound',bound(h),'result',solve(h))
