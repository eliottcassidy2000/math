#!/usr/bin/env python3
"""Independent pair-mask fixed-path DP and half-skeleton enumeration."""
import itertools as it,json,math
from common import *
def half_skeleton(r):
 q=r//2;counts={}
 def walk(last,used,k):
  endpoint=(last=='L' or isinstance(last,int) and last<q) if r%2 else (last=='A' or isinstance(last,int) and last<q)
  if len(used)==q and k and endpoint:counts[k]=counts.get(k,0)+1
  for x in ['A','L',*range(r)]:
   if isinstance(x,int) and (x==r//2 and r%2 or min(x,r-1-x) in used):continue
   if x in ('A','L') and k>=r+1:continue
   ok=last is None or last=='A' and x=='L' or last=='L' and isinstance(x,int) or isinstance(last,int) and x=='A' or isinstance(last,int) and isinstance(x,int) and last<x
   if ok:walk(x,used|{min(x,r-1-x)} if isinstance(x,int) else used,k+(x in ('A','L')))
 walk(None,set(),0)
 return counts
def fixed(a,r):
 g=wrap(a,r);n=len(a);N=len(g)
 sigma=list(range(n,2*n))+list(range(n))+list(reversed(range(2*n,N)))
 pairs=[(i,sigma[i]) for i in range(N) if i<sigma[i]];q=len(pairs)
 pair_id={v:k for k,pair in enumerate(pairs) for v in pair};dp={}
 for v,k in pair_id.items():dp[1<<k,v]=1
 for mask in range(1,1<<q):
  for v in pair_id:
   val=dp.get((mask,v),0)
   if not val:continue
   for w,k in pair_id.items():
    if not(mask>>k&1) and g[v]>>w&1:dp[mask|(1<<k),w]=dp.get((mask|(1<<k),w),0)+val
 centers=[i for i in range(N) if sigma[i]==i]
 require(len(centers)==N%2,'center count')
 count=sum(dp.get(((1<<q)-1,v),0) for v in pair_id if (g[v]>>(centers[0] if centers else sigma[v]))&1)
 if N<=8:
  brute=0
  for path in it.permutations(range(N)):
   if all(g[path[i]]>>path[i+1]&1 for i in range(N-1)) and tuple(sigma[v] for v in reversed(path))==path:brute+=1
  require(count==brute,'literal fixed-path count')
 coeff=half_skeleton(r);p=cover(a)
 require(count==sum(c*math.factorial(k)*(p[k] if k<len(p) else 0) for k,c in coeff.items()),'independent half-skeleton expansion')
 if r==1:require(count==evalF(p,2)-evalF(p,1),'r1 fixed formula')
 if r==2:require(count==2*evalF(p,3)-3*evalF(p,2)+3*evalF(p,1),'r2 fixed formula')
 h=paths(g)[-1];require(count%2==h%2==1,'inherited odd fixed-path parity')
 return {'module':list(a),'r':r,'order':N,'H':h,'fixed_paths':count}
print(json.dumps({'half_skeletons':{r:half_skeleton(r) for r in range(1,5)},'tests':[fixed(a,r) for a in [(0,),(2,0),(2,4,1),*B] for r in range(1,5)]},sort_keys=True,indent=2))
