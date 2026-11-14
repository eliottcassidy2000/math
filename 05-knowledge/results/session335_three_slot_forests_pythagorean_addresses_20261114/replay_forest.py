from forest import *
from counter35 import H
from itertools import combinations_with_replacement

def literal(w):
    n=len(w); r=[0]*(n+2)
    for i in range(n):
        r[i]=sum(1<<j for j in range(i+1,n))
        for h in (0,1):
            if bit(w[i],h):r[n+h]|=1<<i
            else:r[i]|=1<<(n+h)
    r[n]|=1<<(n+1)
    return r

def orders(w):
    r=literal(w)
    return [p for p in permutations(range(len(r))) if all(r[a]>>b&1 for a,b in zip(p,p[1:]))]

types=tuple(product(range(4),repeat=2))
seven=[]
for t in range(4):
    for q in combinations_with_replacement(types,t):
        if count(q)==7:seven.append(q)
assert seven==[((1,0),(3,0)),((3,0),(3,1))]
witnesses=[]
for w in product(range(4),repeat=4):
    for a,b in (((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))):
        t=tuple(sorted([(w[a[0]],w[a[1]]),(w[b[0]],w[b[1]])]))
        if t in seven:
            p=orders(w)
            raw=0
            for order in p:
                runs=[[]]
                for v in order:
                    if v<4:runs[-1].append(v)
                    else:runs.append([])
                raw+=sorted(tuple(x) for x in runs if x)==sorted((a,b))
            assert raw==7
            assert H(w)==len(p)
            witnesses.append((len(p),w,a,b))
assert min(witnesses)==(23,(1,0,3,0),(0,1),(2,3))
from collections import Counter
dist=Counter(sum(z>0 for z in fibre_load(t)) for t in product(types,repeat=3))
assert sorted(dist.items())==[(1,16),(2,48),(4,288),(5,384),(6,1680),(7,576),(8,816),(9,192),(10,96)]
print('Local seven-gluing endpoint types',seven)
print('Exact all-order minimum for seven-gluing witness',min(witnesses))
print('16-chart support distribution',sorted(dist.items()))
