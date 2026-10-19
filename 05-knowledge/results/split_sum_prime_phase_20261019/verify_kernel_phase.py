#!/usr/bin/env python3
"""Finite exact controls; universal arguments are in KERNEL_PHASE_RECOVERY.md."""
from itertools import combinations_with_replacement
from collections import defaultdict
import json
count=0

def require(ok,msg):
 if not ok: raise RuntimeError(msg)

def kernel(ts,h):
 d=[0]*h; d[0]=1
 for t in ts:d=[d[x]+d[(x-t)%h] for x in range(h)]
 return tuple(d)

for h in range(1,17):
 for k in range(0,6):
  for ts in combinations_with_replacement(range(h),k):
   d=kernel(ts,h);T=sum(ts)%h
   per=[a for a in range(h) if all(d[x]==d[(x+a)%h] for x in range(h))]
   centers=[c for c in range(h) if all(d[x]==d[(c-x)%h] for x in range(h))]
   require(set(centers)=={(T+a)%h for a in per},('centers',h,ts))
   require((1<<k)%len(per)==0,('period order',h,ts))
   if h%2==1:require(per==[0] and centers==[T],('odd recovery',h,ts))
   if h%2==0:
    sylow=h&-h
    require((len({c%2 for c in centers})==1)==(len(per)<sylow),('parity',h,ts))
   if h>1 and h&(h-1)==0:
    uniform=len(set(d))==1
    if not uniform:require(all(c%2==T%2 for c in centers),('recover',h,ts))
    elif d[0]==1:require(T%2==1,('unit multiplicity',h,ts))
   count+=1
for s in range(1,10):
 h=1<<s
 for k in range(s+1,s+7):
  base=[1<<v for v in range(s)]
  a=base+[0]*(k-s);b=base+[1]+[0]*(k-s-1)
  require(kernel(a,h)==kernel(b,h)==(1<<(k-s),)*h,('ambiguity',s,k))
  require(sum(a)%2!=sum(b)%2,('opposite',s,k))

def prime_kernel(ps,q):
 ds=[1];n=1
 for p in ps:ds +=[d*p for d in ds];n*=p
 d=tuple(sum(min(a%q,q-a%q)==i for a in ds) for i in range(1,(q+1)//2))
 z=sum((a+n//a)%q==0 for a in ds)
 return d,pow((-n)%q,(q-1)//2,q),z
for q,a,b in [(5,[11],[19]),(17,[2,43],[2,59])]:
 ka,pa,za=prime_kernel(a,q);kb,pb,zb=prime_kernel(b,q)
 require(ka==kb and pa==pb==1 and za==0 and zb==2,('root hostile',q))
print(json.dumps({'status':'PASS','cyclic_lists':count,'uniform_ambiguity_pairs':54,'root_hostiles':2},sort_keys=True,indent=2))
