#!/usr/bin/env python3
"""Independent, standard-library-only finite checks. Not proofs of infinitude."""
import json
from math import isqrt

def factors(n):
    out=[]
    d=2
    while d*d<=n:
        while n%d==0:
            out.append(d); n//=d
        d+=1
    if n>1: out.append(n)
    return out

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def digits(n,b):
    s=0
    while n: s+=n%b; n//=b
    return s

def smith(n,b):
    fs=factors(n)
    return len(fs)>1 and digits(n,b)==sum(digits(p,b) for p in fs)

N=80000
phis=list(range(N+1))
for p in range(2,N+1):
    if phis[p]==p:
        for n in range(p,N+1,p): phis[n]-=phis[n]//p
qs=[q for q in range(5,101) if prime(q)]
for q in qs:
    actual=[n for n in range(1,N+1) if phis[n]==2*q]
    expected=[2*q+1,4*q+2] if prime(2*q+1) else []
    assert actual==expected,(q,actual,expected)
# q=3 is genuinely outside the exact-fiber classification.
assert [n for n in range(1,30) if phis[n]==6]==[7,9,14,18]
pairs=[]
for q in range(3,32,2):
    M=2**q-1; W=(2**q+1)//3
    if prime(M) and prime(W):
        pairs.append(q)
        if q>5: assert not prime(2*q+1)
assert pairs[:6]==[3,5,7,13,17,19],pairs
from fractions import Fraction
for n in range(1,101):
    for t in range(1,51):
        rhs=Fraction(t*phis[n])
        for p in set(factors(t)):
            if n%p: rhs*=Fraction(p-1,p)
        assert rhs==phis[t*n]
assert phis[3]==phis[4] and phis[6]!=phis[8]
print("PASS: complete inverse fibers at 2q for prime q=5..97; all preimages included by n<=2m^2<=80000.")
print("PASS: simultaneous-prime q<=31:", pairs)
print("PASS: support-aware scaling n=1..100,t=1..50; sharp q=3,5 and scaling hostile controls.")
