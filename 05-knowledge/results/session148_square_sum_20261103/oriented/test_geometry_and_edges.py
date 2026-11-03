#!/usr/bin/env python3
from collections import deque
from functools import lru_cache
from math import isqrt
from orientation import graph,paths

def require(b,m):
 if not b:raise ValueError(m)

def square(x):return isqrt(x)**2==x

def distance(adj,start,end):
 Q=deque([(start,0)]);seen={start}
 while Q:
  v,d=Q.popleft()
  if v==end:return d
  for w in adj[v]-seen:seen.add(w);Q.append((w,d+1))
 return None

for n in range(10,41):
 _,adj=graph(n);want=4 if n<15 else 3
 require(distance(adj,8,10)==want,'distance n=%d'%n)
require([8+1,1+15,15+10]==[9,16,25],'three labels')
for t in range(2,101):
 x,z=t*t-1,t*t+1; a,b=8*t*t+1,8*t*t-1
 require(len({x,a,b,z})==4 and min(x,a,b,z)>0,'bad vertices')
 require([x+a,a+b,b+z]==[9*t*t,16*t*t,9*t*t],'bad path')
 require(not square(x+z),'adjacent endpoints')
print('Verified 8,10 prefix distances n=10..40 and t=2..100 witness family.')
for n in (23,25):
 E,P=paths(n)
 def dp(exclude=None):
  adj={v:set() for v in range(1,n+1)}
  for a,b in E:
   if (a,b)!=exclude:adj[a].add(b);adj[b].add(a)
  @lru_cache(None)
  def f(mask,x):
   old=mask^(1<<(x-1))
   if not old:return 1
   return sum(f(old,y) for y in adj[x] if old>>(y-1)&1)
  full=(1<<n)-1
  return sum(f(full,x) for x in adj)
 total=dp()
 for e in E:
  enumerated=sum(e in {tuple(sorted(ab)) for ab in zip(p,p[1:])} for p in P)
  require(total-dp(e)==2*enumerated,'incidence mismatch %d %s'%(n,e))
 print('n=%d independently verified every edge incidence, %d directed paths.'%(n,total),flush=True)
print('All distance and independent edge-incidence checks passed.')
