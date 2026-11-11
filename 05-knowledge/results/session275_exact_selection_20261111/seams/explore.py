from itertools import permutations
from functools import lru_cache

def H(a):
 n=len(a)
 @lru_cache(None)
 def f(mask,v):
  if mask==1<<v:return 1
  return sum(f(mask^(1<<v),u) for u in range(n) if mask>>u&1 and a[u]>>v&1)
 return sum(f((1<<n)-1,v) for v in range(n))
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def flip(a,i,j):
 b=a[:];b[i]^=1<<j;b[j]^=1<<i;return b
if __name__=='__main__':
 for n in range(3,12):print(n,H(flip(tt(n),0,n-1)))
def fans(a,r,m,n):
 # m new source-type (a->s->old other); n sink-type (old other->b->a)
 b=a[:]
 for j in range(m):
  k=len(b); b.append(((1<<k)-1)^(1<<r)); b[r]|=1<<k
 for j in range(n):
  k=len(b); b.append(1<<r)
  for u in range(k):
   if u!=r:b[u]|=1<<k
 return b
if __name__=='__main__':
 for a in [tt(2),flip(tt(3),0,2),flip(tt(4),0,3)]:
  for r in range(len(a)):
   arr=[[H(fans(a,r,m,n)) for n in range(1,4)] for m in range(1,4)]
   print(a,r,arr)
