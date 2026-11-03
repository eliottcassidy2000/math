#!/usr/bin/env python3
from math import gcd, isqrt
from collections import Counter

def require(flag, message):
    if not flag: raise RuntimeError(message)

def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def count(p,q):
    n=p*q
    holes={t for t in range(n) if t%p==0 or t%q==0}
    units=set(range(n))-holes
    edges={(t,(t+1)%n) for t in units if (t+1)%n in units}
    hole_edges={(t,(t+1)%n) for t in holes if (t+1)%n in holes}
    starts=[t for t in units if (t-1)%n not in units]
    lengths=[]
    for t in starts:
        k=0
        while t in units:
            k+=1; t=(t+1)%n
        lengths.append(k)
    require(len(units)==(p-1)*(q-1),'vertex count')
    require(len(edges)==(p-2)*(q-2),'edge count')
    require(len(holes)==p+q-1 and len(hole_edges)==2,'holes')
    require(len(lengths)==p+q-3 and sum(lengths)==len(units),'components')
    small,large=sorted((p,q))
    expected={k:2 for k in range(1,small-1)}
    expected[small-1]=large-small+1
    require(dict(sorted(Counter(lengths).items()))==expected,'full path-length distribution')
    return dict(sorted(Counter(lengths).items()))

for p,q in [(3,5),(3,7),(5,29)]:
    print('p,q',p,q,'component lengths',count(p,q))
P=[n for n in range(3,50,2) if prime(n)]
for i,p in enumerate(P):
    for q in P[i+1:]: count(p,q)
print('PASS:',len(P)*(len(P)-1)//2,'distinct odd-prime pairs below 50; explicit checks survive -O')
