#!/usr/bin/env python3
"""Exact independent checks for ARITHMETIC.md; no third-party dependencies."""
from fractions import Fraction
from functools import lru_cache
from itertools import product, permutations
from math import gcd, isqrt, prod
import json

checks = 0

def check(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

@lru_cache(None)
def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

def R(p, a):
    return (p**(a+1)-1)//(p-1)

@lru_cache(None)
def sigma(n):
    return prod(R(p, a) for p, a in factor(n).items())

def square(n):
    return isqrt(n)**2 == n

@lru_cache(None)
def Q(n):
    return frozenset(p for p, a in factor(sigma(n)).items() if a % 2)

def D(m,n):
    return Q(m*n)^Q(m)^Q(n)

for n in range(1,2001):
    check(sigma(n) == sum(d for d in range(1,n+1) if n%d==0), f'divisor sum {n}')

for m,n in product(range(1,101), repeat=2):
    fm,fn=factor(m),factor(n)
    ratio=prod((Fraction(R(p,fm[p]+fn[p]), R(p,fm[p])*R(p,fn[p])) for p in fm.keys() & fn.keys()), start=Fraction(1))
    check(Fraction(sigma(m*n),sigma(m)*sigma(n))==ratio, f'port defect {m,n}')
    if gcd(m,n)==1:
        check(Q(m*n)==Q(m)^Q(n), f'coprime map {m,n}')
        s=lambda x:sigma(x)-x
        check(s(m*n)==s(m)*s(n)+m*s(n)+n*s(m), f'aliquot identity {m,n}')

for l,m,n in product(range(1,21), repeat=3):
    check(D(l,m)^D(l*m,n)==D(m,n)^D(l,m*n), f'cocycle {l,m,n}')

primes=[p for p in range(2,201) if factor(p)=={p:1}]
hits=[(p,a,isqrt(R(p,a))) for p in primes for a in range(1,17) if square(R(p,a))]
check(hits==[(3,1,2),(3,4,11),(7,3,20)], 'bounded zero columns')

# Test local updates by exact rational identities; independently test cocycle by telescoping.
for p in (3,5,7,19):
    for a,b,c in product(range(7),repeat=3):
        check(Fraction(R(p,b),R(p,a))*Fraction(R(p,c),R(p,b))==Fraction(R(p,c),R(p,a)), f'update telescope {p,a,b,c}')
    for a,b,u in product(range(7),range(7),range(1,30)):
        if gcd(p,u)==1:
            check(Fraction(sigma(p**b*u),sigma(p**a*u))==Fraction(R(p,b),R(p,a)), f'update {p,a,b,u}')

rows=[]
for bits in product(range(2),repeat=3):
    i,j,k=bits
    n=3**i*5**j*19**k
    check(sigma(n)==4**i*6**j*20**k, f'cube value {bits}')
    check(Q(n)==(frozenset((2,3)) if j else frozenset())^(frozenset((5,)) if k else frozenset()), f'cube class {bits}')
    check(square(sigma(n))==(j==k==0), f'cube kernel {bits}')
    rows.append({'bits':bits,'n':n,'sigma':sigma(n),'class':sorted(Q(n))})

states=list(product((0,1,4),(0,3)))
for u in (1,5,11,19,25,55,95):
    values=[3**a*7**b*u for a,b in states]
    check(len(set(values))==6, 'six distinct states')
    check(all(Q(n)==Q(u) for n in values), f'six-state invariant {u}')
    actions=[]
    for ap,bp in product(permutations((0,1,4)),permutations((0,3))):
        amap=dict(zip((0,1,4),ap)); bmap=dict(zip((0,3),bp))
        action=tuple(3**amap[a]*7**bmap[b]*u for a,b in states)
        check(set(action)==set(values), f'action bijective {u}')
        actions.append(action)
    check(len(set(actions))==12,'faithful S3 x C2')

for n,value in [(9,13),(25,31),(361,381),(125,156),(6859,7240),(243,364),(115,144),(1501,1600)]:
    check(sigma(n)==value,f'certificate {n}')
check(Q(9)!=Q(3)^Q(3), 'multiplicativity hostile')
check(not square(sigma(243)), 'shared 3 port hostile')
check(square(sigma(3)-3) and not square(sigma(81)-81), 'aliquot hostile')
check(sigma(1001)-1001==343 and sigma(343)-343==57, 'aliquot segment')
check(sorted(3**a*7**b for a,b in states)==[1,3,81,343,1029,27783], 'neutral products')
check(sigma(5*19*29)==3600, 'three-way nonzero charge cancellation')
check(Q(9)==frozenset((13,)) and Q(243)==frozenset((7,13)), 'shared-port differing output classes')
check(pow(13,2,5)==4 and pow(91,2,5)==1, 'quadratic character at 5')
check(Q(3)==Q(81)==frozenset(), 'same global squareclass before raw observer')
check(pow(sigma(3),5,11)==1 and pow(sigma(81),5,11)==0, 'raw Legendre observer at 11 distinguishes unit from zero')
print(json.dumps({'status':'PASS','explicit_checks':checks,'prime_power_scan':{'prime_bound':200,'exponent_bound':16,'square_hits':hits},'requested_cube':rows,'six_state_group_order':12,'classification_scope':'Infinite statement uses cited square-case theorem; scan is finite only.'},indent=2))
