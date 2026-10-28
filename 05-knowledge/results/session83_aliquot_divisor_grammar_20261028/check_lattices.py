#!/usr/bin/env python3
"""Exact finite controls for LATTICE_COMPANIONS.md; Python standard library."""
from collections import Counter
from itertools import product
from math import gcd, isqrt, prod
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                a += 1
                n //= p
            out.append((p, a))
        p += 1 if p == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def sigma(n):
    return prod(sum(p**i for i in range(a+1)) for p, a in factors(n))


def divisors(n):
    low, high = [], []
    for d in range(1, isqrt(n)+1):
        if n % d == 0:
            low.append(d)
            if d*d != n:
                high.append(n//d)
    return low + list(reversed(high))


def aliquot(n):
    return sigma(n)-n


def c(n):
    return n*(n-1)


def star(a, b):
    x, u = a
    y, v = b
    return x*y, x*v+y*u+u*v


box_pairs = 0
for n in range(1, 2001):
    fac = factors(n)
    vectors = list(product(*(range(a+1) for _, a in fac)))
    evaluate = lambda e: prod(p**j for (p, _), j in zip(fac, e))
    weighted = {e: evaluate(e) for e in vectors}
    ds = divisors(n)
    check(sorted(weighted.values()) == ds, ('box bijection', n))
    check(sum(ds) == sigma(n), ('box sum', n))
    check(sum(ds[:-1]) == aliquot(n), ('proper divisor sum', n))
    if n >= 2:
        top = tuple(a for _, a in fac)
        support = set(vectors) - {top}
        unit = lambda i: tuple(int(j == i) for j in range(len(fac)))
        active = [i for i in range(len(fac)) if unit(i) in support]
        edges = []
        for i in active:
            for j in active:
                if i < j:
                    pair = tuple(int(k in (i,j)) for k in range(len(fac)))
                    if pair not in support:
                        edges.append((i,j))
        independent = set()
        for bits in product((0,1), repeat=len(active)):
            v = [0]*len(fac)
            for i,bit in zip(active,bits):
                v[i]=bit
            if all(not(v[i] and v[j]) for i,j in edges):
                independent.add(tuple(v))
        check((support == independent) == (sum(a for _,a in fac)<=2),
              ('active singleton graph boundary',n))
        if len(fac)==1 and fac[0][1]==1:
            check(active==[] and support=={(0,)}, ('prime ghost boundary',n))
    for e in vectors:
        for f in vectors:
            meet = tuple(min(x, y) for x, y in zip(e, f))
            join = tuple(max(x, y) for x, y in zip(e, f))
            d, q = weighted[e], weighted[f]
            check(weighted[meet] == gcd(d,q), ('meet', n,e,f))
            check(weighted[join] == d*q//gcd(d,q), ('join', n,e,f))
            box_pairs += 1

for a in range(1, 81):
    for b in range(1, 81):
        actual = Counter(d*e for d in divisors(a) for e in divisors(b))
        af, bf = dict(factors(a)), dict(factors(b))
        primes = sorted(set(af)|set(bf))
        expected = {}
        for exps in product(*(range(af.get(p,0)+bf.get(p,0)+1) for p in primes)):
            d = prod(p**k for p,k in zip(primes,exps))
            multiplicity = prod(max(0,min(af.get(p,0),k)-max(0,k-bf.get(p,0))+1)
                                for p,k in zip(primes,exps))
            expected[d] = multiplicity
        check(actual == expected, ('multiplicity', a,b))
        check(set(actual) == set(divisors(a*b)), ('surjection',a,b))
        check((len(actual) == sum(actual.values())) == (gcd(a,b)==1), ('injectivity',a,b))
        check((sigma(a)*sigma(b) == sigma(a*b)) == (gcd(a,b)==1), ('sigma guard',a,b))
        if gcd(a,b)==1:
            check(star((a,aliquot(a)),(b,aliquot(b))) == (a*b,aliquot(a*b)), ('star',a,b))

for a,b,d in product(range(1,7), repeat=3):
    x,y,z = (a,a+1),(b,b+2),(d,d+3)
    check(star(star(x,y),z) == star(x,star(y,z)), ('associativity',x,y,z))

for x in range(2,101):
    prime = factors(x)==[(x,1)]
    for k in range(1,5):
        g = sum(x**j for j in range(k+1))
        check(sigma(x**k)>=g, ('grouping inequality',x,k))
        check((sigma(x**k)==g)==prime, ('grouping equality',x,k))

check(divisors(27)[:-1] == [1,3,9], 'cube shape')
check(divisors(35)[:-1] == [1,5,7], 'semiprime shape')
check(aliquot(27)==aliquot(35)==13, 'fiber collision')
check((aliquot(54),aliquot(70))==(66,74), 'fiber operation failure')
check(aliquot(64)==63 and 1+4+16==21, 'composite base hostile')
proper = divisors(1001)[:-1]
columns = [0,0,0]
for d in proper:
    j=0
    while d:
        columns[j] += d%7
        d//=7
        j+=1
check(columns==[14,19,4], 'base 7 columns')
check(aliquot(1001)==343 and aliquot(343)==57 and sigma(343)==400, 'anchor chain')
for x in range(-20,21):
    check(14+19*x+4*x*x-x**3 == (7-x)*(x+1)*(x+2), ('carry polynomial',x))

N=10000
sieve=[0]*(N+1)
for d in range(1,N//2+1):
    for n in range(2*d,N+1,d):
        sieve[n]+=d
edge_images=[]
sigma_commuting=[]
aliquot_commuting=[]
for n in range(2,N+1):
    t,u=sieve[n],sieve[n-1]
    sc=aliquot(c(n))
    check(sc == (n-1)*t+n*u+t*u, ('companion law',n))
    check(sc-c(t)==t*(n-t)+(n+t)*u, ('companion defect',n))
    check(sigma(c(n))==(n+t)*(n-1+u), ('sigma product',n))
    check((sigma(c(n))==c(sigma(n))) == (t==u), ('sigma commutation',n))
    total=t+sc
    check(total==(n+t)*(1+sigma(n-1))-n*n, ('edge total',n))
    if n<=1000:
        check(sc==sum(divisors(c(n))[:-1]), ('independent companion divisors',n))
    if n>=3 and t>0 and sc>0 and t!=sc and isqrt(total)**2==total:
        edge_images.append([n,c(n),t,sc,isqrt(total)])
    if t==u:
        sigma_commuting.append(n)
    if sc==c(t):
        aliquot_commuting.append(n)
        check(t>n, ('abundance exclusion',n))
        check(factors(n-1)!=[(n-1,1)], ('prime predecessor exclusion',n))
check(edge_images==[[264,69432,456,120648,348]], 'companion image census')
check(aliquot(3)+aliquot(6)==7, 'minimal edge hostile')
check(aliquot(1001)+aliquot(1001000)==2144303, 'anchor edge failure')
check(sigma_commuting==[3], 'finite sigma commuting census')
check(aliquot_commuting==[], 'finite aliquot commuting census')

print(json.dumps({
    'status':'PASS',
    'box_sources':2000,
    'box_meet_join_pairs':box_pairs,
    'divisor_product_pairs':6400,
    'grouped_base_cases':396,
    'source_collision':{'sources':[27,35],'aliquot':13,'doubled_aliquot':[66,74]},
    'anchor':{'proper_divisors':proper,'base7_column_sums':columns,'aliquot':343,'next':57},
    'companion_range':[2,N],
    'companion_direct_divisor_range':[2,1000],
    'square_edge_image_survivors':edge_images,
    'sigma_commuting_finite':sigma_commuting,
    'aliquot_commuting_finite':aliquot_commuting,
    'scope':'Finite lists are complete only within the displayed intervals; universal statements use the proofs.'
}, indent=2, sort_keys=True))
