#!/usr/bin/env python3
from math import gcd, prod
from fractions import Fraction
from sympy import isprime, factorint, divisors, totient

def inv(m):
    ps=[d+1 for d in divisors(m) if isprime(d+1)]
    out=[]
    def rec(i,n,v):
        if i==len(ps):
            if v==m: out.append(n)
            return
        p=ps[i]
        rec(i+1,n,v)
        a=p; f=p-1
        while m%(v*f)==0:
            rec(i+1,n*a,v*f)
            a*=p; f*=p
    rec(0,1,1)
    return sorted(out)

M=200
N=2*M*M
phi=list(range(N+1))
for p in range(2,N+1):
    if phi[p]==p:
        for n in range(p,N+1,p): phi[n]-=phi[n]//p
fib={m:[] for m in range(1,M+1)}
for n in range(1,N+1):
    if phi[n]<=M: fib[phi[n]].append(n)
for m in range(1,M+1): assert inv(m)==fib[m],(m,inv(m),fib[m])
print('Complete inverse fibers m=1..200: support enumeration agrees with rigorous-bound sieve')
for q in range(5,101):
    if isprime(q):
        assert inv(2*q)==([2*q+1,4*q+2] if isprime(2*q+1) else [])
print('Prime-index fiber classification checked q=5..97')
rows=[]
for q in range(3,32,2):
    a=2**q-1; b=(2**q+1)//3
    if isprime(a) and isprime(b):
        assert q<=5 or not inv(2*q)
        rows.append((q,2*q,a,b,inv(2*q)))
print('Double-prime rows (q,k,Mq,Wq,inverse_phi(k)):',rows)
for n in range(1,101):
    for t in range(1,51):
        z=Fraction(t*int(totient(n)))
        for p in factorint(t):
            if n%p:z*=Fraction(p-1,p)
        assert z==totient(t*n)
print('Support-aware scaling 5000 exact tests: PASS')
assert inv(14)==[] and 29 in inv(28)
assert inv(34)==inv(68)==[]
assert totient(3)==totient(4) and totient(6)!=totient(8)
assert inv(74)==[] and not isprime(2**37-1)
print('Hostile controls: PASS; no global sequence equivalence asserted')
