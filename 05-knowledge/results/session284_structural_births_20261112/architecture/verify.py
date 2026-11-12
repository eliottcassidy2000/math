#!/usr/bin/env python3
"""Independent small controls for rooted/unrooted grammar (not large census)."""
import itertools,json,math
from pathlib import Path

def rows(n,c):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   a,b=(i,j) if c>>k&1 else (j,i);r[a]|=1<<b;k+=1
 return r

def code(r,p):
 c=k=0
 for i in range(len(p)):
  for j in range(i+1,len(p)):
   c|=((r[p[i]]>>p[j])&1)<<k;k+=1
 return c

def canon(r):return min(code(r,p) for p in itertools.permutations(range(len(r))))
def autos(r):
 old=code(r,tuple(range(len(r))))
 return [p for p in itertools.permutations(range(len(r))) if code(r,p)==old]
def cycles(p):
 seen=set();n=0
 for x in range(len(p)):
  if x not in seen:
   n+=1
   while x not in seen:seen.add(x);x=p[x]
 return n
def add(r,z):
 n=len(r);s=list(r)+[0]
 for i in range(n):
  if z>>i&1:s[n]|=1<<i
  else:s[i]|=1<<n
 return s
def check(cond,msg):
 if not cond:raise RuntimeError(msg)
G={0};out=[]
for n in range(1,6):
 G={canon(add(rows(n-1,c),z)) for c in G for z in range(1<<(n-1))}
 check(len(G)==[1,1,2,4,12][n-1],'class recurrence')
 out.append({'n':n,'classes':len(G)})
 if n<=4:
  for c in G:
   r=rows(n,c);au=autos(r)
   orbit={min(sum(((z>>p[i])&1)<<i for i in range(n)) for p in au) for z in range(1<<n)}
   E=sum(2**cycles(p) for p in au)//len(au)
   check(E==len(orbit),'Burnside')
   children={canon(add(r,z)) for z in range(1<<n)}
   check(len(children)<=E<=(n+1)*len(children),'forget root bound')
def prime(r):
 n=len(r)
 for mask in range(1,1<<n):
  if mask.bit_count() in (1,n):continue
  inside=[i for i in range(n) if mask>>i&1]
  if all(len({(r[v]>>u)&1 for u in inside})==1 for v in range(n) if not(mask>>v&1)):return False
 return True
w=[sum(1<<j for j in range(i+1,4)) for i in range(4)]+[0]
for i in range(4):
 if i%2==0:w[4]|=1<<i
 else:w[i]|=1<<4
check(prime(w) and len(autos(w))==1,'W5 prime rigidity')
pc=[add(w,z) for z in range(32) if prime(add(w,z))]
check(len(pc)==20,'prime-parent attachment classification')
pct=len({canon(r) for r in pc})
check(pct>=math.ceil(20/6),'prime child lower bound')
result={'prime_W5_children':pct,'prime_W5_attachments':len(pc),'status' :'PASS','scope':'all classes through5; orbit/forgetting controls for all parents through4','levels':out}
print(json.dumps(result,indent=2))
