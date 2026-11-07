#!/usr/bin/env python3
"""Exact, dependency-free finite controls. Writes small JSON to stdout."""
import json, math

def s(n,b):
    out=0
    while n: out+=n%b; n//=b
    return out

def fac(n):
    out=[];p=2
    while p*p<=n:
        while n%p==0:out.append(p);n//=p
        p+=1
    if n>1:out.append(n)
    return out

def A(n,b):return sum(s(p,b) for p in fac(n))
def prime(n):return n>=2 and fac(n)==[n]
def smith(n,b):return n>=4 and not prime(n) and s(n,b)==A(n,b)
def R(b,k):return (b**k-1)//(b-1)
def G(n,b,h):
    ans=0;i=0
    while n:
        ans+=(n%b)*(b**(i%h)-1)//(b-1);n//=b;i+=1
    return ans

checks={}
for n in range(4,1001):
    if prime(n):continue
    F=sum(fac(n))
    for b in range(2,n+3):
        if smith(n,b):assert n==4 or (n-F)%(b-1)==0
checks['base_sieve']={'composite_n':[4,1000],'bases':'2 through n+2'}
assert [b for b in range(2,20) if smith(4,b)]==list(range(5,20))
primes=[p for p in range(5,300) if prime(p)]
for p in primes:assert smith(2*p,p-1) and smith(p*p,p-1)
checks['moving_base_prime_limit']=299
for b in range(2,31):
    for k in range(2,8):
        r=R(b,k)
        if prime(r):
            for c in range(2,b):assert smith(c*r,b)==(c==k==2)
checks['one_digit_classification']={'bases':[2,30],'lengths':[2,7]}
for b in range(2,501):assert not smith(b*b-1,b)
checks['maximal_two_digit_nonsmith']={'bases':[2,500]}
for b in range(2,12):
    for t in range(1,60):
        m=1
        while R(b,m)<t:m+=1
        c=(b-1)*t+1
        for k in range(m,m+5):assert s(c*R(b,k),b)==k+s(t,b)+s(R(b,m)-t,b)-m
checks['separated_boundary']={'bases':[2,11],'t':[1,59],'k':'m through m+4'}
assert fac(1540)==[2,2,5,7,11] and A(1540,10)==18
for k in range(4,301):assert s(1540*R(10,k),10)==k+18
checks['decimal_seed_digit_identity_lengths']=[4,300]
for k in range(2,12):
    for c in range(1,2**k+1):assert s(c*(2**k-1),2)==k
checks['binary_small_multiplier']={'lengths':[2,11],'c':'1 through 2^k'}
for b in range(2,10):
    for h in range(2,5):
        for n in range(2,400):
            assert s(n,b**h)==s(n,b)+(b-1)*G(n,b,h)
            assert s(n,b**h)-A(n,b**h)==s(n,b)-A(n,b)+(b-1)*(G(n,b,h)-sum(G(p,b,h) for p in fac(n)))
checks['grouping']={'bases':[2,9],'h':[2,4],'n':[2,399]}
for b in range(2,21):
    for x in range(1,60):
        for y in range(1,60):
            raw=s(x,b)*s(y,b)-s(x*y,b)
            assert raw>=0 and raw%(b-1)==0
checks['carry_defect']={'bases':[2,20],'x_y':[1,59]}
for b in range(2,15):
    for k in range(2,12):
        rr=R(b,k)
        assert math.gcd(rr,b-1)==math.gcd(k,b-1)
        if prime(rr):assert prime(k) and (b-1)%k!=0
        if prime(k):
            for q in set(fac(rr)):
                j=1
                while pow(b,j,q)!=1:j+=1
                assert (q==k and b%q==1) or (j==k and (q-1)%k==0)
checks['repunit_order']={'bases':[2,14],'lengths':[2,11]}
examples={str(n):{'factors':fac(n),'smith_bases':[b for b in range(2,n+1) if smith(n,b)]} for n in [14,26,34,38,22,63,238]}
print(json.dumps({'status':'PASS','checks':checks,'examples':examples,'seed_hostile_k6':{'n':1540*R(10,6),'digit_sum':s(1540*R(10,6),10),'factor_weight':A(1540*R(10,6),10)},'scope':'Finite exact controls support, but do not replace, the elementary proofs.'},indent=2))
