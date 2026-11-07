#!/usr/bin/env python3
"""Exact bounded tests. Standard library; integer trial-division primality only."""
from functools import lru_cache
from math import isqrt
import json
from pathlib import Path

def prime(n):
    if n<2:return False
    if n%2==0:return n==2
    return all(n%d for d in range(3,isqrt(n)+1,2))

def divs(m):
    ds=[]
    for d in range(1,isqrt(m)+1):
        if m%d==0:
            ds.append(d)
            if d*d!=m:ds.append(m//d)
    return sorted(ds)

def inverse_count(m):
    ps=[d+1 for d in divs(m) if prime(d+1)]
    @lru_cache(None)
    def rec(i,r):
        if i==len(ps):return int(r==1)
        p=ps[i]; out=rec(i+1,r); c=p-1
        while c<=r:
            if r%c==0:out+=rec(i+1,r//c)
            c*=p
        return out
    return rec(0,m)

# Pure powers require no large primality test: h not a power of two gives
# algebraic factorization, h=32 gives the displayed F5 factorization.
fermat=[3,5,17,257,65537]
assert all(prime(p) for p in fermat)
assert prime(641) and prime(6700417)
assert 641*6700417==(1<<32)+1
assert 6700416==(1<<7)*3*17449
assert prime(17449)
assert pow(2,32,641)==640 and pow(2,32,6700417)==6700416
assert (641-1)*(6700417-1)==(1<<14)*15*17449

def coefficient(a,weights):
    # Enumerate allowed odd support; each remainder specifies one 2-part,
    # except remainder zero, which allows 2-part 1 or 2.
    sums=[0]
    for w in weights:sums += [s+w for s in sums]
    return sum(2 if s==a else 1 for s in sums if s<=a)

pure=[coefficient(a,[1,2,4,8,16]) for a in range(64)]
assert pure==[a+2 if a<32 else 32 for a in range(64)]
assert all(inverse_count(1<<a)==pure[a] for a in range(19))

def predicted(q,a):
    weights=[1,2,4,8,16]
    out=sum(coefficient(a-j,weights) for j in range(1,a+1) if prime((1<<j)*q+1))
    if q in fermat:
        h=(q-1).bit_length()-1
        if a>=h:out+=coefficient(a-h,[w for w in weights if w!=h])
    return out

rows=[]
for q in [3,5,7,13,17,19,31]:
    vals=[]
    for a in range(1,19):
        actual=inverse_count((1<<a)*q)
        assert actual==predicted(q,a),(q,a,actual,predicted(q,a))
        vals.append(actual)
    if q not in fermat:
        cs=[0]+vals
        ss=[0]+[int(prime((1<<j)*q+1)) for j in range(1,19)]
        for a in range(1,19):
            assert cs[a]-2*cs[a-1]+(cs[a-2] if a>=2 else 0)==2*ss[a]-ss[a-1]
    rows.append({'q':q,'supplier_heights_through_18':[j for j in range(1,19) if prime((1<<j)*q+1)],'counts_at_a_1_through_18':vals})
# Independent complete sieve controls on small targets, using phi(n)^2 >= n/2.
limit=2*100**2
phi=list(range(limit+1))
for p in range(2,limit+1):
    if phi[p]==p:
        for n in range(p,limit+1,p):phi[n]-=phi[n]//p
freq=[0]*101
for x in range(1,limit+1):
    if phi[x]<=100:freq[phi[x]]+=1
assert all(freq[m]==inverse_count(m) for m in range(1,101))
assert [m for m in range(2,77,2) if not freq[m]]==[14,26,34,38,50,62,68,74,76]
# A007374 minima certified only while the minimum lies within this complete window.
mins={str(n):next(m for m in range(1,101) if freq[m]==n) for n in range(2,12)}
assert list(mins.values())==[1,2,4,8,12,32,36,40,24,48]
receipt={'status':'PASS','pure_power_counts_a0_to63':pure,'prime_q_rows':rows,'a007374_minima_multiplicities_2_to11':mins,'inverse_count_checks':19+126+100,'sieve_source_bound':limit,'F5_factors':[641,6700417],'phi_F5':(641-1)*(6700417-1),'limits':'Pure powers through a63 proved by Fermat factorization; independent inverse recursion for a<=18. Prime-q test a<=18, seven q values. No unbounded primality assertion.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
