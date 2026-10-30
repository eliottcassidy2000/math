from math import isqrt,gcd
from collections import defaultdict
import json

def sigma(n):
    return sum(d+n//d if d*d!=n else d for d in range(1,isqrt(n)+1) if n%d==0)
def kernel(n):
    n//=n&-n
    k=1;p=3
    while p*p<=n:
        a=0
        while n%p==0:n//=p;a^=1
        if a:k*=p
        p+=2
    return k*n
def parity(n):return (sigma(n)-n)%2
checks=0
for n in range(1,301):
    for d in range(1,101):
        predicted=int(kernel(n)==kernel(d))^((n*d)%2)
        if predicted!=parity(n*d):raise RuntimeError((n,d))
        checks+=1
primes=[3,5,7,11,13,17,19]
rows=[]
for p in primes:
    n=2**p; successor=2**p-1;k=kernel(successor)
    if sigma(n)-n!=successor or k==1:raise RuntimeError(p)
    rows.append((p,n,successor,k))
for i,(_,n,s,k) in enumerate(rows):
    for j,(_,nn,ss,kk) in enumerate(rows):
        if i==j:continue
        if gcd(s,ss)!=1 or kernel(s*k)!=1 or kernel(ss*k)==1:raise RuntimeError((i,j))
        checks+=1
print(json.dumps(dict(checks=checks,rows=rows),indent=2))
