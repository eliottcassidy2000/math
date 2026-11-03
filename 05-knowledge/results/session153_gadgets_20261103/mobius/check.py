from itertools import permutations,combinations
from functools import lru_cache
from math import factorial

def calc(n,mask):
 a=[[0]*(n+1) for _ in range(n+1)]
 for i in range(1,n+1): a[0][i]=a[i][0]=1
 for bit,(i,j) in enumerate(combinations(range(1,n+1),2)):
  if mask>>bit&1:a[i][j]=1
  else:a[j][i]=1
 @lru_cache(None)
 def per(s):
  vs=[j for j in range(n+1) if s>>j&1]
  return sum(all(a[i][j] for i,j in zip(vs,ps)) for ps in permutations(vs))
 @lru_cache(None)
 def polys(s):
  if not s:return (1,)
  first=s&-s; rem=s^first; sub=rem; out=[0]*(n+2)
  while True:
   block=sub|first; v=per(block)
   if v:
    for k,c in enumerate(polys(s^block)[:n+1]):out[k+1]+=v*c
   if not sub:break
   sub=(sub-1)&rem
  return tuple(out)
 p=polys((1<<(n+1))-1)
 terms=[((-1)**(k-1))*factorial(k-1)*p[k] for k in range(1,n+2)]
 h=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(1,n+1)))
 if sum(terms)!=h: raise RuntimeError((n,mask,terms,h))
 return h,terms
for n in range(1,6):
 hist={}
 for mask in range(1<<(n*(n-1)//2)):
  h,_=calc(n,mask);hist[h]=hist.get(h,0)+1
 print('n',n,'tournaments',sum(hist.values()),'spectrum',sorted(hist.items()))
n=7; mask=0
for bit,(i,j) in enumerate(combinations(range(1,n+1),2)):
 if (j-i)%n in (1,2,3): mask|=1<<bit
print('cyclic 7',calc(n,mask))
