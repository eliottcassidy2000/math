from itertools import combinations,permutations
from collections import Counter
for n,targets in [(4,[(5,9)]),(5,[(5,9),(19,23)]),(6,[(19,23)])]:
 P=list(permutations(range(n))); edge=list(combinations(range(n),2)); ind={e:i for i,e in enumerate(edge)}
 d=[(sum(1<<ind[(a,b)] for a,b in zip(p,p[1:]) if a<b),sum(1<<ind[(b,a)] for a,b in zip(p,p[1:]) if b<a)) for p in P]
 k=1<<len(edge); Ms=[tuple(i for i,(a,b) in enumerate(d) if a&m==a and b&m==0) for m in range(k)]
 for a,b in targets:
  C=Counter();ex={}
  for m,M in enumerate(Ms):
   if len(M)!=a:continue
   for i,(x,y) in enumerate(edge):
    q=m^(1<<i); N=Ms[q]
    if len(N)!=b:continue
    o=len(set(M)&set(N)); C[o]+=1;ex[o]=(m,q,x,y)
  assert sorted(C.items())=={(4,5,9):[],(5,5,9):[(3,240),(4,240)],(5,19,23):[],(6,19,23):[(11,1440),(14,1440),(17,1440)]}[n,a,b]
  print(n,a,b,sorted(C.items()),ex)
