from math import gcd,isqrt
from collections import Counter

def require(v,msg):
    if not v: raise ValueError(msg)

def prime(x):
    return x >= 2 and all(x%d for d in range(2,isqrt(x)+1))

primes=[p for p in range(3,100,2) if prime(p)]
count=0; examples={}
for i,p in enumerate(primes):
  for q in primes[i+1:]:
    m=2*p*q; odds=list(range(1,m,2)); V={r for r in odds if gcd(r,m)==1}; H=set(odds)-V
    EE=[(r,(r+2)%m) for r in odds if r in V and (r+2)%m in V]
    HH=[(r,(r+2)%m) for r in odds if r in H and (r+2)%m in H]
    require(len(H)==p+q-1,('holes',p,q))
    require(len(HH)==2 and len(set(HH))==2,('hole adjacencies',p,q,HH))
    require(len(V)==(p-1)*(q-1) and len(EE)==(p-2)*(q-2),('counts',p,q))
    U={frozenset(e) for e in EE}
    require(len(U)==len(EE) and all(len(e)==2 for e in U),('directed edge duplication',p,q))
    adj={r:set() for r in V}
    for u,v in EE: adj[u].add(v);adj[v].add(u)
    unseen=V.copy(); components=[]
    while unseen:
        stack=[unseen.pop()]; C=[]
        while stack:
            u=stack.pop(); C.append(u)
            for v in adj[u]&unseen: unseen.remove(v);stack.append(v)
        components.append(C)
    require(len(components)==p+q-3 and max(map(len,components))<=min(p-1,q-1),('paths',p,q))
    histogram=dict(Counter(map(len,components)))
    expected={length:2 for length in range(1,p-1)}
    expected[p-1]=q-p+1
    require(histogram==expected,('full distribution',p,q,histogram,expected))
    require(all(sum(len(adj[v]) for v in C)==2*(len(C)-1) and all(len(adj[v])<=2 for v in C) for C in components),('not paths',p,q))
    # For odd units, the r=m-1 -> 1 edge must occur exactly once and in the +2 direction.
    require(EE.count((m-1,1))==1 and (1,m-1) not in EE,('wrap',p,q))
    # Resolve hole locations by CRT independently of the production t labels.
    require(all((u%p==0 and v%q==0) or (u%q==0 and v%p==0) for u,v in HH),('CRT correspondence',p,q))
    if (p,q) in [(3,5),(3,7),(5,29)]: examples[p,q]=dict(sorted(Counter(map(len,components)).items()))
    count+=1
require(examples=={(3,5):{1:2,2:3},(3,7):{1:2,2:5},(5,29):{1:2,2:2,3:2,4:25}},examples)
print('PASS',count,'distinct unordered odd-prime pairs below 100; counts, simple directed edges, CRT hole edges, path components, maximum lengths and wrap verified')
print('examples',examples)
