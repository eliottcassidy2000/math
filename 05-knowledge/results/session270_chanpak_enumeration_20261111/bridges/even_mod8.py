#!/usr/bin/env python3
"""Exact even-order A000568 mod-8 query from a supplied factorization."""
from math import gcd, factorial, prod
from fractions import Fraction
from collections import Counter
import json

def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=1
    return out

def phi_from(f):
    return prod((p-1)*p**(e-1) for p,e in f.items())

def query(f):
    n=prod(p**e for p,e in f.items())
    assert n>=4 and n%2==0
    terms=[(1,1)]
    for p in f:terms += [(d*p,-mu) for d,mu in terms[:]]
    A=sum(mu*d*((n-1)//(2*d))*(((n-1)//(2*d))+1)//2 for d,mu in terms)
    c=phi_from(f)//2
    correction=0
    if n%3==0 and n>=12:
        g=f.copy();g[3]-=1
        if not g[3]:del g[3]
        correction=4*(phi_from(g)//2)
    return (n*A-c+correction+6*(n==4)+2*(n==6))%8

def parts(n,low=1):
    if not n:yield ();return
    for a in range(low,n+1,2):
        for rest in parts(n-a,a):yield (a,)+rest

def direct(n):
    ans=Fraction()
    for p in parts(n):
        e=sum(x//2 for x in p)+sum(gcd(x,y) for i,x in enumerate(p) for y in p[:i])
        z=prod(p)*prod(factorial(m) for m in Counter(p).values())
        ans+=Fraction(1<<e,z)
    assert ans.denominator==1 and ans.numerator%(1<<(n//2))==0
    return (ans.numerator>>(n//2))%8

if __name__=='__main__':
    rows=[]
    for n in range(4,41,2):
        got=query(factor(n));expected=direct(n)
        assert got==expected,(n,got,expected)
        rows.append({'n':n,'normalized_mod8':got})
    huge=[{2:1000},{2:1,3:1000},{2:2,5:1000},{2:1,1000000007:2}]
    print(json.dumps({'status':'FINITE-EXACT','direct_orders':rows,'factorization_queries':[{'factors':f,'n_bits':prod(p**e for p,e in f.items()).bit_length(),'normalized_mod8':query(f)} for f in huge]},indent=2))
