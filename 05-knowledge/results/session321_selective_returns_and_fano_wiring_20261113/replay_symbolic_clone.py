from functools import lru_cache
from replay_clone_profiles import C,P,X,Y,U,V,add
def compress_cut(nq,n,cut):
 b=cut&1
 if (cut&((1<<n)-1)) != (b*((1<<n)-1)):
  raise ValueError('nonuniform source')
 return ((cut>>n)<<1)|b
def layouts(R):
 n=len(R);N=(1<<n)-2-(1<<(n-1))
 @lru_cache(None)
 def f(S,last,i,j):
  from collections import Counter
  d=Counter()
  if S==N and i and j:d[i,j]+=1
  for v in range(n):
   if last!=-1 and not R[last]>>v&1:continue
   if v==0:
    if i<n:d.update(f(S,0,i+1,j))
   elif v==n-1:
    if j<n:d.update(f(S,n-1,i,j+1))
   elif not S>>v&1:d.update(f(S|(1<<v),v,i,j))
  return d
 return f(0,-1,0,0)
for r,n,k in ((C,9,7167),(P,10,6144),(P,10,2048),(P,10,5120)):
 t=compress_cut(len(r),n,k);R=add(r,t);L=layouts(R)
 want={7167:(29,15,24,12,4,2),6144:(53,27,48,24,12,6),
       2048:(43,15,40,15,10,4),5120:(75,61,12,72,66,14,18,18,4)}[k]
 if tuple(v for _,v in sorted(L.items()))!=want:raise ValueError((k,R,L))
 print('marked',k,R,sorted(L.items()))
