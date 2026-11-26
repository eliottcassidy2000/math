"""Exact height census checks for primitive twist bases (excluding j=0,1728)."""
from collections import Counter
from math import isqrt
from functools import lru_cache

def cubert(n):
    a,b=0,1
    while b**3<=n:b*=2
    while b-a>1:
        c=(a+b)//2
        if c**3<=n:a=c
        else:b=c
    return a

def valuation(a,p):
    c=0
    while a%p==0:
        a//=p;c+=1
    return c

def primes_to(N):
    a=[]
    for d in range(2,N+1):
        if all(d%p for p in a if p*p<=d):a.append(d)
    return a

def mobius_table(N):
    M=[1]*(N+1)
    for p in primes_to(N):
        for k in range(p,N+1,p):M[k]*=-1
        for k in range(p*p,N+1,p*p):M[k]=0
    return M

@lru_cache(None)
def base(X,power):
    a,b=cubert(X),isqrt(X)
    D=cubert(isqrt(X)) if power==1 else isqrt(cubert(isqrt(X)))
    M=mobius_table(D)
    total=2*sum(M[k]*(a//k**(2*power))*(b//k**(3*power)) for k in range(1,D+1))
    # Singular short forms with B>0 are (-3*t^2,2*t^3); the minimum excludes
    # nonsquarefree t when power=2 and every t>1 when power=1.
    N=cubert(isqrt(X//27))
    return total-(int(N>0) if power==1 else sum(bool(x) for x in mobius_table(N)[1:]))

def partition(X):
    D=cubert(isqrt(X))
    M=mobius_table(D)
    return {d:2*base(X//d**6,1) for d in range(1,D+1) if M[d]}

def brute(X):
    A,B=cubert(X),isqrt(X)
    P=primes_to(cubert(B))
    C=Counter()
    for a in range(-A,A+1):
        if a==0:continue
        for b in range(1,B+1):
            if 4*a**3+27*b*b==0:continue
            d=1;ok=True
            for p in P:
                if p*p<=abs(a) and p**3<=b:
                    k=min(valuation(a,p)//2,valuation(b,p)//3)
                    if k>1:ok=False;break
                    if k:d*=p
            if ok:C[d]+=2
    return dict(C)

for X in (32,4096,24576,2**18):
    P=partition(X)
    assert {d:v for d,v in P.items() if v}==brute(X)
    assert sum(P.values())==2*base(X,2)
for X in (10**6,10**12,10**18):
    D=partition(X)
    assert sum(D.values())==2*base(X,2)
    print(X, sum(D.values()), D[1], D[1]/sum(D.values()))
