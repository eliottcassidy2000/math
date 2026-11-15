#!/usr/bin/env python3
"""Literal controls for the digit-reversal row limit."""
from math import log, sqrt
from collections import Counter
from bisect import bisect_right
from random import Random

def p_values(N):
    P=[0]*N
    for n in range(1,N):
        u=n//2
        P[n]=P[u]+(0 if n&1 else u+(u&1))
    return P

def rev(n,L):
    r=0
    for _ in range(L):
        r=2*r+(n&1); n>>=1
    return r

def cdf(t,c,k):
    if t<=0:return 0
    if t>=c+k:return 1
    lo=max(c,t)
    return max(0,t-c)/k+t*log((c+k)/lo)/k

def main():
    rng=Random(109); tests=0
    for L in (8,9,12,13,14):
        N=1<<L; P=p_values(N); R=[rev(n,L) for n in range(N)]
        assert all(abs(N*P[n]-n*(N-R[n]))<=N*(L+1) for n in range(N))
        for s in range(1,L//2+1):
            C=Counter((n>>(L-s),r>>(L-s)) for n,r in enumerate(R))
            assert len(C)==1<<(2*s) and set(C.values())=={1<<(L-2*s)}
        for c,k,phi,b in ((189,11,11,11),(18,3,2,9)):
            v=sorted(c*n+k*P[n]+phi*(n&1)+b for n in range(1,N))
            assert all(x&1 for x in v)
            for j in range(200):
                y=rng.randrange((c+k+2)*N)
                count=bisect_right(v,y)
                mass=N*cdf(y/N,c,k)
                assert abs(mass-count)<=8*sqrt(N)+4*(L+1)
                tests+=1
    print('exact digit-prefix certificates at five lengths;',tests,'row distribution checks')
if __name__=='__main__':main()
