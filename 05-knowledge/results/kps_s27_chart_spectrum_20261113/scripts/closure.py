# permanent omissions of the chart class below a Busch bound B from a miss list computed through order n0
# usage: python closure.py missfile B   (B = f(n0+1): every strong member with H < B has order <= n0)
import sys
miss=set(map(int,open(sys.argv[1]).read().split())); B=int(sys.argv[2])
att=[v for v in range(1,B,2) if v not in miss]
S=set(att)
changed=True
while changed:
    changed=False
    for a in list(S):
        for b in list(S):
            if a>1 and b>1 and a*b<B and a*b not in S: S.add(a*b); changed=True
perm=[v for v in range(1,B,2) if v not in S]
from sympy import isprime, factorint
print('permanent omissions below',B,':',perm)
print('  primes:',[v for v in perm if isprime(v)])
print('  7p:',[v for v in perm if v%7==0 and isprime(v//7)])
print('  other:',[(v,factorint(v)) for v in perm if not isprime(v) and not (v%7==0 and isprime(v//7))])
pr=[p for p in range(3,B,2) if isprime(p)]
print('  primes below B:',len(pr),' omitted:',sum(1 for p in pr if p not in S))
