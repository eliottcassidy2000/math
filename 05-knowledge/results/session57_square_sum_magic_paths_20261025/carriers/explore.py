from itertools import combinations
from math import isqrt
Q=lambda p:{x*x%p for x in range(p)}
def arr(c,a,b):return [a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a]
for p in [29,59]:
 q=Q(p)
 for c in sorted(q):
  found=False
  for a in sorted(q):
   for b in sorted(q):
    m=arr(c,a,b)
    if len(set(m))==9 and set(m)<=q:
     print('nowrap',p,m, 'total',3*c);found=True;break
   if found:break
  if found:break
# false triangle in fixed-total hypergraph shadow
for T in range(1,100):
 q=[x*x for x in range(isqrt(T)+1)]
 h=[t for t in combinations(q,3) if sum(t)==T]
 es=set(e for t in h for e in combinations(t,2))
 for t in combinations(q,3):
  if all(e in es for e in combinations(t,2)) and sum(t)!=T:
   print('false triangle',T,t,h);raise SystemExit
