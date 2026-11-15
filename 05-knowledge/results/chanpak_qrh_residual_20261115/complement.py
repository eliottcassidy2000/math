#!/usr/bin/env python3
"""Preserve the physical three-window while improving uniform bracketing."""
from math import isqrt
from random import Random

K = (135,137,139)

def p(n):
    r = 0
    while n:
        b = n & 1
        n >>= 1
        if not b:
            r += n + (n & 1)
    return r

def baseline(x,y):
    return 432*x*y + 129*x + 258*y

def checksum(x,y):
    return 45*(p(x)+p(y)) + 15*(x&1) + 30*(y&1)

def value(x,y):
    v=baseline(x,y)+checksum(x,y)
    return tuple(v+z for z in K)

def first_across(lo,hi,side):
    """Return adjacent i,j with side(i) false and side(j) true."""
    assert not side(lo) and side(hi)
    while hi-lo>1:
        m=(lo+hi)//2
        if side(m):
            hi=m
        else:
            lo=m
    return lo,hi

def pair(N):
    assert N >= 10_000
    g=1
    while 432*(2*g)**6 <= N:
        g*=2
    a=isqrt(N//108)
    s0=a+(108*a*a<N)
    S=s0+(g-1-s0)%g
    def xy(k):
        x=(S+1)//2+k
        return x,S-x
    def V(k):
        return value(*xy(k))
    def lower(k):
        return baseline(*xy(k)) + 45*S + 45 + K[-1] <= N
    assert V(0)[0]>N
    cap=S-(S+1)//2-1
    assert not lower(0) and lower(cap)
    _,end=first_across(0,cap,lower)
    i,_=first_across(0,end,lambda k:V(k)[0]<=N)
    _,j=first_across(0,end,lambda k:V(k)[-1]<=N)
    upper=(xy(i),V(i))
    lower=(xy(j),V(j))
    assert lower[1][-1] <= N < upper[1][0]
    assert all((x+y)%g==g-1 and x*y>0 for (x,y),_ in (lower,upper))
    return g,S,lower,upper

def transfer(g,b):
    sig=0
    s=g.bit_length()-1
    for j in range(s):
        sig|=(1-((b>>(s-1-j))&1))<<j
    return sig,(0 if not s else 1-((b>>(s-1))&1))

def controls():
    count=0
    for s in range(8):
        g=1<<s
        for b in range(g):
            r=g-1-b
            A,z=transfer(g,b)
            B,w=transfer(g,r)
            assert A+B==g-1
            for u in range(40):
                assert p(g*u+b)==p(u)+A*u+p(b)+z*(u&1)
                count+=1
    print('transfer cases',count)
    R=Random(491065)
    tested=0
    for N in range(10_000,50_000):
        pair(N);tested+=1
    for m in range(14,513):
        for N in [2**m-1,2**m+1, R.randrange(2**m,2**(m+1))]:
            pair(N);tested+=1
    print('brackets including 512-bit targets',tested)
    for N in (35_099,1_000_001,10**12+1,10**36+1):
        print(N,pair(N))
if __name__=='__main__':
    import sys
    if len(sys.argv)>1:
        for q in sys.argv[1:]:
            print(q,pair(int(q)))
    else:
        controls()
