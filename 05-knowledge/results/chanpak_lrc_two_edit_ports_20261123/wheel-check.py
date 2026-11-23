#!/usr/bin/env python3
"""Finite local and CRT replay of the six-read physical wheel."""
from math import isqrt,prod,gcd

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def forms(x,d=1):
    return (40*x+d,52*x+3*d,52*x+5*d,312*x+67*d,840*x+19*d,840*x+23*d)
def crt(rs):
    m=prod(p for p,a in rs)
    return sum(a*(m//p)*pow(m//p,-1,p) for p,a in rs)%m,m
def measured(x,d=1):
    f=list(forms(x,d));assert f[2]%3==0;f[2]//=3;return f

assert [q for q in range(5) if any(f%5==0 for f in forms(q))]==[0,1,4]
assert all(x%3 for x in measured(1))
for p in [p for p in range(7,501) if prime(p)]:
    rootsets=[{q for q in range(p) if forms(q)[i]%p==0} for i in range(6)]
    assert all(len(R)<=1 for R in rootsets)
    assert len(set.union(*rootsets))<=6<p
    three=set.union(*(rootsets[i] for i in (0,4,5)))
    assert len(three)==(1 if p==7 else 3)
for y in (5,7,17,43,97):
    ls=[p for p in range(5,y+1) if prime(p)]
    Q,M=crt([(9,1)]+[(p,next(q for q in range(p) if all(f%p for f in forms(q)))) for p in ls])
    for j in range(5):
        x=Q+j*M
        assert all(gcd(m,2*M)==1 for m in measured(x))
print("Exact local root counts through 499, and six simultaneous CRT wheels through 97")
