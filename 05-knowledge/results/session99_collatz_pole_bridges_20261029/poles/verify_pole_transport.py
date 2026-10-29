#!/usr/bin/env python3
"""Finite exact controls for POLE_TRANSPORT.md; only Python standard library."""
from fractions import Fraction as Q
from math import comb, factorial, gcd, prod
from collections import Counter
import json


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def pow2(a):
    return Q(2**a) if a >= 0 else Q(1, 2**(-a))


def cores(n, least=3):
    if n == 0:
        yield ()
    for p in range(least, n+1, 2):
        for tail in cores(n-p, p):
            yield (p,)+tail


def data(mu):
    m = sum(mu)
    c = len(mu)
    E = sum((p-1)//2 for p in mu)
    E += sum(gcd(mu[i],mu[j]) for i in range(c) for j in range(i+1,c))
    L = comb(m,2)-E
    z = prod(p**t*factorial(t) for p,t in Counter(mu).items())
    return m, (m-c)//2, Q(factorial(m),z)*pow2(-L)


def coefficient(mu,n):
    m,k,C = data(mu)
    return C*comb(n,m)*pow2(-2*k*(n-m)) if n>=m else Q(0)


checks = 0
core_count = 0
for mass in range(19):
    for mu in cores(mass):
        core_count += 1
        m,k,C = data(mu)
        for p in (3,5,7,9):
            hp=sum(gcd(p,a) for a in mu)
            factor=pow2(-m+hp-(p-1)**2//2)/ (p*(mu.count(p)+1))
            child=tuple(sorted(mu+(p,)))
            for n in range(36):
                r=n-p
                rhs=(factor*prod(range(r+1,r+p+1))*coefficient(mu,r)
                     *pow2(-(p-1)*r)) if r>=0 else Q(0)
                require(coefficient(child,n)==rhs, ('insertion',mu,p,n))
                checks+=1

for k in range(11):
    for n in range(36):
        r=n-3
        rhs=(Q(1,12)*prod(range(r+1,r+4))*coefficient((3,)*k,r)
             *pow2(-2*r)) if r>=0 else Q(0)
        require((k+1)*coefficient((3,)*(k+1),n)==rhs,('tower',k,n))
        checks+=1



for ell in (3,5,7,9):
    prefactor=pow2(-(ell-1)**2//2)/ell
    for k in range(6):
        m,pole,C=data((ell,)*k)
        expected=Q(factorial(ell*k),ell**k*factorial(k))*pow2(-(ell-1)*(ell*k*k-k)//2)
        require(C==expected,('general C',ell,k))
        amplitude=C*pow2((ell-1)*ell*k*k)
        expected_A=Q(factorial(ell*k),ell**k*factorial(k))*pow2((ell-1)*(ell*k*k+k)//2)
        require(amplitude==expected_A,('general amplitude',ell,k))
        for n in range(56):
            r=n-ell
            rhs=(prefactor*prod(range(r+1,r+ell+1))*coefficient((ell,)*k,r)
                 *pow2(-(ell-1)*r)) if r>=0 else Q(0)
            require((k+1)*coefficient((ell,)*(k+1),n)==rhs,('general tower',ell,k,n))
            checks+=1


def F(n):
    if n<0 or n%3:
        return Q(0)
    k=n//3
    return pow2(-3*k*k+k)/ (3**k*factorial(k))

for n in range(61):
    rhs=Q(1,4)*F(n-2)*pow2(-2*(n-2)) if n>=2 else Q(0)
    require((n+1)*F(n+1)==rhs,('pantograph',n))
    checks+=1

u,v=(3,7),(5,5)
require(data(u)[:2]==data(v)[:2], 'equal pole labels')
hu,hv=(sum(gcd(3,a) for a in mu) for mu in (u,v))
require((hu,hv)==(4,2),'distinct insertion responses')
require(pow2(-sum(u)+hu)/pow2(-sum(v)+hv)==4,'hostile ratio')
print(json.dumps({'status':'PASS','arithmetic':'exact fractions',
                  'core_mass_max':18,'core_count':core_count,
                  'inserted_lengths':[3,5,7,9],'Taylor_degree_max':35,
                  'top_core_k_max':10,'general_lengths':[3,5,7,9],
                  'general_k_max':5,'general_Taylor_degree_max':55,'pantograph_degree_max':60,
                  'coefficient_equalities':checks,
                  'hostile_same_labels':{'cores':[u,v],'k':4,'m':10,
                                         'h3':[hu,hv],'response_ratio':4}},indent=2))
