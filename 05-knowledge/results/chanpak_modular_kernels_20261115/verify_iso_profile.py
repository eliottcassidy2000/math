#!/usr/bin/env python3
from verify_hyperbolic import rows,covers,equiv
from fractions import Fraction

pairs=[(8,217,12,27280129),(8,223,12,10760449),
       (9,695,14,10342203650),(9,707,14,31344104449)]
for n,h,k,code in pairs:
    a=rows(n,code)
    assert covers(a)[:2]==(h,2**k)
    for arcs in (a,equiv(a)):
        seen={0}
        for _ in range(n):
            seen|=set.union(*(arcs[i] for i in seen))
        assert len(seen)==n

u=Fraction(707,695)
v=Fraction(223,217)
assert 1<u<v<u*u and min(u,v/u,u*u/v)>Fraction(201,200)
for r in range(2,15):
    q=2**r
    u=707*pow(695,-1,q)%q
    v=223*pow(217,-1,q)%q
    m=2**(r-2)-1
    scale=217*pow(695,m,q)%q
    seq=[]
    a=1
    for j in range(m+1):
        seq.extend(((-scale*a)%q,(-scale*a*v)%q))
        a=a*u%q
    assert len(set(seq))==len(seq)==q//2
print('4 literal strong profiles, real separation inequalities, and complete odd residues through 2^14 verified')
