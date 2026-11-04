from math import isqrt
from collections import Counter
import time

def sq(n):return {v:{w for w in range(1,n+1) if v!=w and isqrt(v+w)**2==v+w} for v in range(1,n+1)}

def direct(G):
 V=sorted(G);N=len(V);out=Counter();witness={}
 if not N:return out,witness
 if N==1:return Counter({(V[0],V[0]):1}),{(V[0],V[0]):V}
 ix={v:i for i,v in enumerate(V)};A=[sum(1<<ix[w] for w in G[v]) for v in V]
 if any(not a for a in A):return out,witness
 leaves=[i for i,a in enumerate(A) if a.bit_count()==1]
 if len(leaves)>2:return out,witness
 starts=leaves[:1] if leaves else range(N)
 allbits=(1<<N)-1
 def dfs(p,left):
  if not left:
   if leaves or p[0]<p[-1]:
    pair=tuple(sorted((V[p[0]],V[p[-1]])));out[pair]+=1
    if pair not in witness:witness[pair]=[V[i] for i in p]
   return
  a=A[p[-1]]&left
  while a:
   b=a&-a;a-=b;w=b.bit_length()-1
   dfs(p+[w],left^b)
 for s in starts:dfs([s],allbits^(1<<s))
 return out,witness
if __name__=='__main__':
 for n in range(15,34):
  t=time.time();C,W=direct(sq(n));print(n,sum(C.values()),len(C),time.time()-t,flush=True)
