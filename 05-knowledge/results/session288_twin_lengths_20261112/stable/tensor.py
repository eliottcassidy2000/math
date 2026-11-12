"""Exact compressed skeleton compiler for acyclic unselected ports."""
from collections import defaultdict
from itertools import product
from math import factorial,comb

def kernel(rows,O,ports):
 n=len(rows);t=len(ports);fixed=[v for v in range(n) if v not in ports];k=len(fixed);fidx={v:i for i,v in enumerate(fixed)};pidx={v:i for i,v in enumerate(ports)}
 # Reject cyclic port sets, for which the finite run bound is unavailable.
 if sorted(sum(bool(rows[v]>>w&1) for w in ports) for v in ports)!=list(range(t)):raise ValueError('ports must be acyclic')
 D=k+2;out=defaultdict(int);states={}
 for v in range(n):
  c=[0]*t;mask=0
  if v in pidx:c[pidx[v]]=1
  else:mask=1<<fidx[v]
  states[(mask,tuple(c),v,0)]=1
 for length in range(1,k+t*D+1):
  nxt=defaultdict(int)
  for (mask,c,last,br),number in states.items():
   if mask==(1<<k)-1 and all(c):out[c]+=number*(1 if not br else -1)
   for v in range(n):
    if v==last:continue
    nm=mask;cc=list(c)
    if v in fidx:
     bit=1<<fidx[v]
     if mask&bit:continue
     nm|=bit
    else:
     i=pidx[v]
     if cc[i]>=D:continue
     cc[i]+=1
    key=(nm,tuple(cc),v,br)
    if rows[last]>>v&1:nxt[key]+=number
    if not br and O>>last&1 and O>>v&1:nxt[(nm,tuple(cc),v,1)]+=number
  states=nxt
  if not states:break
 return {c:v for c,v in out.items() if v}

def expo(K):
 ans=defaultdict(int)
 for counts,c in K.items():
  for bases in product(*(range(1,r+1) for r in counts)):
   a=c
   for r,b in zip(counts,bases):a*=(-1)**(r-b)*comb(r,b)
   ans[bases]+=a
 return {b:c for b,c in ans.items() if c}

def anchored(E):
 # b^m = b + (b^m-b), b>=2; 1^m=1.
 A=defaultdict(int)
 for bases,c in E.items():
  choices=[(0,b) if b>1 else (0,) for b in bases]
  for labels in product(*choices):
   val=c
   for b,label in zip(bases,labels):
    if not label:val*=b
   A[labels]+=val
 return {a:c for a,c in A.items() if c}

def stable_faces(A,t):
 out=[]
 for mask in range(1<<t):
  if all(not c or not any(a) or any(x and not(mask>>i&1) for i,x in enumerate(a)) for a,c in A.items()):out.append(mask)
 return out

def evaluate(E,lengths):return sum(c*__import__('math').prod(b**m for b,m in zip(bs,lengths)) for bs,c in E.items())
