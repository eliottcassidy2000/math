#!/usr/bin/env python3
"""Independent Jones degree exclusions, via dimension-wise binary search.
No group census; theorem applicability still depends on certified primitivity
and the cited CFSG-dependent Jones Theorem1.2(k=0).
"""
from math import isqrt,prod
import json
checks=0
def check(p,s):
 global checks
 checks+=1
 if not p:raise RuntimeError(s)
def fibonacci(n):
 x,y=0,1
 for _ in range(n):x,y=y,x+y
 return x
def repunit(q,d):return sum(q**i for i in range(d))
rows=[]
for n,factors in [(25,[11,101,151]),(55,[11,11,199,331,39161])]:
 N=fibonacci(n-1)+fibonacci(n+1)
 check(prod(factors)==N,'complete factorization')
 for p in factors:check(p>=2 and all(p%d for d in range(2,isqrt(p)+1)),'factor primality')
 check(len(factors)>1 and N not in (11,23),'non-affine and nonsporadic degree')
 q=N-1;check(q%2==0 and q.bit_count()>1,'dimension2 not prime power')
 excluded=[];d=3
 while repunit(2,d)<=N:
  lo,hi=2,isqrt(N)+1
  while lo<hi:
   mid=(lo+hi)//2
   if repunit(mid,d)<N:lo=mid+1
   else:hi=mid
  check(repunit(lo,d)!=N,'integer repunit degree excluded')
  check(repunit(lo,d)>N and (lo==2 or repunit(lo-1,d)<N),'binary search boundary')
  excluded.append(d);d+=1
 check(repunit(2,d)>N,'all higher dimensions excluded')
 # Reflection quotient is path0..m with a loop forbidding its last vertex.
 m=(n-1)//2;zero,one=1,1
 for i in range(1,m+1):zero,one=zero+one,(zero if i<m else 0)
 fixed=zero+one
 check(fixed==fibonacci(m+2),'reflection fixed word DP')
 check((N-fixed)%4==0 and N%2==1,'both generators even')
 rows.append({'period':n,'N':N,'reflection_fixed':fixed,'transpositions':(N-fixed)//2,'projective_dimensions_checked':excluded,'conclusion':'A_N, given sparse primitivity and Jones k=0 classification'})
print(json.dumps({'checks':checks,'rows':rows},sort_keys=True,indent=2))
