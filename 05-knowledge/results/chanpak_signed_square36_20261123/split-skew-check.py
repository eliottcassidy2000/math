#!/usr/bin/env python3
"""Finite replay of the centered branches, split reachability and cycle criterion."""
from random import Random
from collections import deque

def v2(x):
    assert x
    return (x&-x).bit_length()-1
def address(word,finish=1):
    z=finish;M=4
    for A in reversed(word):
        M<<=A
        z=((z<<A)-1)*pow(3,-1,M)%M
    return z,M
def transfer(x,word):
    for A in word:
        assert x%2 and v2(3*x+1)==A
        x=(3*x+1)>>A
    return x
def qmap(x,A,q):return (3*x+1)*pow(2,-A,q)%q
def word_q(x,y,q):
    maybe={x:()}; Q=deque([x])
    while Q:
        z=Q.popleft()
        if z==y:return maybe[z]
        for A in (2,3):
            zz=qmap(z,A,q)
            if zz not in maybe:
                Q.append(zz);maybe[zz]=maybe[z]+(A,)
    raise AssertionError("target not reached")

rng=Random(29);count=0
for q0 in (5,7,11,13,25,35):
    for k in (0,1,2,3):
        q=q0*3**k
        for _ in range(10):
            a=rng.randrange(q); b=rng.randrange(q)
            while k and b%3==0:b=rng.randrange(q)
            tail=[2]*(k-1) if k else []
            if k:
                A=next(A for A in range(2,2+2*3**(k-1)) if b*pow(2,A,3**k)%3**k==4%3**k)
                tail.append(A)
            aim=b%q0
            for A in reversed(tail):aim=(pow(2,A,q0)*aim-1)*pow(3,-1,q0)%q0
            w=list(word_q(a%q0,aim,q0))+tail
            x,M=address(w)
            x+=M*((a-x)*pow(M,-1,q)%q)
            for t in (-2,2):
                z=x+t*M*q; y=transfer(z,w)
                assert y%q==b and y%4==1 and y*z>0
                count+=1
assert count==480
print(count,'signed literal split-only CRT itineraries')
verified=[]
for A in range(2,30):
    for k in range(70):
        d=2**(A+k)-3**(k+1)
        if (2**A-2)%d==0:
            assert not(A%2==k%2==1)
            x=3**k*(2**A-2)//d-1
            assert transfer(x,[A]+[1]*k)==x
            verified.append((A,k,x))
assert verified==[(2,0,1),(2,1,-7)]
print('1820 one-split signatures tested; only 1 and -7/-5 in that bounded scan')
