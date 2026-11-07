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

N=20000
phis=list(range(N+1))
for p in range(2,N+1):
    if phis[p]==p:
        for n in range(p,N+1,p): phis[n]-=phis[n]//p
qs=[q for q in range(5,501) if prime(q)]
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
for p in range(5,501):
    if prime(p):
        assert smith(p*p,p-1) and smith(2*p,p-1)
assert not smith(9,2) and not smith(6,2)
repunit_cases=0
for b in range(2,41):
    for k in range(2,8):
        R=(b**k-1)//(b-1)
        if prime(R):
            for c in range(2,b):
                repunit_cases+=1
                assert smith(c*R,b)==(c==2 and k==2),(b,k,c)
for b in range(2,1001):
    assert not smith(b*b-1,b),b
print(json.dumps({'status':'PASS','inverse_totient':{'prime_q_range':[5,500],'x_range':[1,N],'tested_q':len(qs)},'simultaneous_prime_q_up_to31':pairs,'smith_prime_range':[5,500],'prime_repunit_cases':repunit_cases,'repunit_base_range':[2,40],'repunit_length_range':[2,7],'b_squared_minus_one_bases':[2,1000]},indent=2))
# Multi-digit multiplier outside the one-digit obstruction.
assert factors(1540)==[2,2,5,7,11]
assert sum(digits(p,10) for p in factors(1540))==18
for k in range(4,201):
    R=(10**k-1)//9
    assert digits(1540*R,10)==k+18,k
# Dyadic repair criterion checked against an independent totient sieve.
for q in range(3,101,2):
    if not prime(q): continue
    for a in range(1,7):
        m=(2**a)*q
        predicted=any(prime((2**j)*q+1) for j in range(1,a+1)) or any(q-1==2**j for j in range(1,a+1))
        # Positive constructions always fit this sieve: n<=4*m.
        actual=any(phis[n]==m for n in range(1,min(N,4*m)+1))
        assert predicted==actual,(q,a)
print('Additional checks PASS: 1540 identity k=4..200; dyadic criterion odd prime q<=100, a=1..6.')
def weight(n,b): return sum(digits(p,b) for p in factors(n))
def raw_carry(x,y,b):
    def ds(n):
        out=[]
        while n: out.append(n%b); n//=b
        return out or [0]
    dx,dy=ds(x),ds(y)
    raw=[0]*(len(dx)+len(dy))
    for i,u in enumerate(dx):
        for j,v in enumerate(dy): raw[i+j]+=u*v
    carry=mass=0
    for value in raw:
        carry=(value+carry)//b
        mass+=carry
    while carry:
        carry//=b; mass+=carry
    return mass

def grouping(n,b,h):
    out=i=0
    while n:
        out+=(n%b)*(b**(i%h)-1)//(b-1)
        n//=b; i+=1
    return out
for b in range(2,13):
    for x in range(1,31):
        for y in range(1,31):
            C=raw_carry(x,y,b)
            assert digits(x*y,b)==digits(x,b)*digits(y,b)-(b-1)*C
    for n in range(1,301):
        for h in range(1,5):
            assert digits(n,b**h)==digits(n,b)+(b-1)*grouping(n,b,h)
            assert digits(n,b**h)-weight(n,b**h)==digits(n,b)-weight(n,b)+(b-1)*(grouping(n,b,h)-sum(grouping(p,b,h) for p in factors(n)))
for k in range(1,9):
    for c in range(1,2**k+1): assert digits(c*(2**k-1),2)==k
assert smith(238,2)
for n in range(4,301):
    if prime(n): continue
    for b in range(2,n+3):
        if smith(n,b): assert n==4 or ((n-sum(factors(n)))%(b-1)==0 and b<=n)
print('Additional checks PASS: raw-convolution carries; grouped-base transport; binary obstruction; finite all-base sieve.')
