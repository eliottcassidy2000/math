#!/usr/bin/env python3
import itertools,json,hashlib
from collections import defaultdict
from functools import lru_cache

def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def arcs(n,mask):
 a=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def paths(a,S):
 return [p for p in itertools.permutations(S) if all(a[u]>>v&1 for u,v in zip(p,p[1:]))]
def merge(a,P,Q):
 i=j=0;out=[]
 while i<len(P) and j<len(Q):
  if a[P[i]]>>Q[j]&1:out.append(P[i]);i+=1
  else:out.append(Q[j]);j+=1
 return tuple(out)+P[i:]+Q[j:]
def overlap_merge(a,P,Q,S):
 order=tuple(x for x in P if x in S)
 check(order==tuple(x for x in Q if x in S),'incompatible input')
 out=();i=j=0
 for x in order:
  u=P.index(x);v=Q.index(x)
  out+=merge(a,P[i:u],Q[j:v])+(x,);i=u+1;j=v+1
 return out+merge(a,P[i:],Q[j:])

tournaments=overlap_cases=pair_checks=0
for n in range(1,6):
 for mask in range(1<<(n*(n-1)//2)):
  a=arcs(n,mask);V=set(range(n))
  hp={m:paths(a,tuple(v for v in range(n) if m>>v&1)) for m in range(1<<n)}
  full=set(hp[(1<<n)-1]);tournaments+=1
  for u in range(n):
   for v in range(n):
    if not(a[u]>>v&1):continue
    f=sum(p.index(u)<p.index(v) for p in full)
    if f==1:check(len(full)==1,'forward-one rigidity')
    if f==2:check(len(full)==3,'forward-two rigidity')
  # Each vertex is X-only, Y-only, or shared. Both private parts are nonempty.
  for marks in itertools.product(range(3),repeat=n):
   if 0 not in marks or 1 not in marks:continue
   X={v for v in V if marks[v]!=1};Y={v for v in V if marks[v]!=0};S=X&Y
   mx=sum(1<<v for v in X);my=sum(1<<v for v in Y)
   image=set();expected=0
   for P in hp[mx]:
    for Q in hp[my]:
     if tuple(v for v in P if v in S)!=tuple(v for v in Q if v in S):continue
     R=overlap_merge(a,P,Q,S)
     check(R in full,'merged output is not a path')
     check(tuple(v for v in R if v in X)==P,'X restriction failed')
     check(tuple(v for v in R if v in Y)==Q,'Y restriction failed')
     check(R not in image,'injection failed')
     image.add(R);expected+=1;pair_checks+=1
   check(len(image)==expected and expected<=len(full),'overlap lower bound failed')
   if len(S)<=1:check(expected==len(hp[mx])*len(hp[my]),'product count failed')
   overlap_cases+=1
 print('completed order',n,flush=True)
sharp=[10,28,25,16,1]
check(len(paths(sharp,range(5)))==9,'sharp graph count')
check(len(paths(sharp,(0,1,2)))==len(paths(sharp,(0,3,4)))==3,'sharp components')
hostile=arcs(6,91); pp=paths(hostile,range(6));f=sum(p.index(0)<p.index(1) for p in pp);r=len(pp)-f
check(hostile[0]>>1&1,'hostile edge')
check((f,r)==(7,8),'hostile order profile')
result={'status':'PASS','universe':'all labeled tournaments through order five; all ordered covering pairs with both private parts nonempty; every compatible Hamilton-path pair','tournaments':tournaments,'overlap_cases':overlap_cases,'compatible_pair_checks':pair_checks,'sharp_rows':sharp,'sharp_H':9,'majority_hostile_rows':hostile,'majority_hostile_forward':f,'majority_hostile_reverse':r}
print(json.dumps(result,sort_keys=True,indent=2))
