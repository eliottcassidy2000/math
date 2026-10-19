#!/usr/bin/env python3
import json,math

def need(b,m):
 if not b:raise RuntimeError(m)
def sigma(n):
 ans=0
 for d in range(1,math.isqrt(n)+1):
  if n%d==0:ans+=d+(n//d if d*d!=n else 0)
 return ans
def prime(n):
 return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
ps=[p for p in range(5,2001) if prime(p)]
count=0;twins=0
for C in (6,28,496):
 need(sigma(C)==2*C,('perfect',C))
 for p in ps:
  if C%p==0:continue
  need(sigma(C*p)-C*p==C*(p+2),('edge',C,p));count+=1
for p in ps:
 n=6*p;q=p+2
 need(sigma(n)-n==6*q,('six',p))
 if prime(q):
  twins+=1;need(not prime(q+2),('short lane',p))
  need(sigma(6*q)-6*q==6*(q+2),('second edge',p))
need(sigma(18)-18!=6*(3+2),'coprimality hostile')
need(sigma(150)-150!=6*(25+2),'composite hostile')
need(sigma(54)-54==6*(9+2) and not prime(9),'accidental agreement is not primality')
print(json.dumps({'status':'PASS','guarded_edges':count,'twin_edges_lower_prime_5_to_2000':twins,'native_guard_hostiles':3},sort_keys=True,indent=2))
