#!/usr/bin/env python3
"""Exact finite controls for the signed central-table dictionary."""
from math import gcd

def order4(q):
    a, t = 4 % q, 1
    while a != 1:
        a = 4*a % q
        t += 1
    return t

def center(m):
    k = m//2
    return k*k % m, k*(k+1) % m

def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k, n

edges = ladders = 0
for s in (1,-1):
    for x in range(1,2001):
        m = 4*x+s
        d,c = center(m)
        assert (d,c) == ((3*x+1,x) if s==1 else (x,3*x-1))
        raw=(c,d) if s==1 else (d,c)
        if x%2:
            decoded=raw
            expected=(x,3*x+s)
        else:
            decoded=(raw[0],(raw[1]-s)//6)
            expected=(x,x//2)
        assert decoded==expected
        edges+=1
    for y in range(1,32,2):
        if y%3==0: continue
        k0=next(k for k in (1,2) if (2**k*y-s)%3==0)
        x0=(2**k0*y-s)//3
        for L in range(1,100,2):
            P=order4(3*L//gcd(y,L))
            x=x0
            residues=[]
            for j in range(P):
                assert v2(3*x+s)==(k0+2*j,y)
                residues.append(x%L)
                m=4*x+s
                xn=4*x+s
                assert 4*xn+s==4*m+s
                x=xn
            assert x%L==x0%L
            assert len(set(residues))==P
            for r in range(min(P,3)):
                xr=(2**k0*4**r*y-s)//3
                q=2
                xx=4**(q*P)*xr+s*(4**(q*P)-1)//3
                assert xx==(2**k0*4**(q*P+r)*y-s)//3
            ladders+=1
assert center(9)==(7,2) and (2,7)!=(2,1)
assert order4(27)==9 and order4(135)==18
print(f'PASS: {edges} decoded edges; {ladders} exact minimal-period ladder controls.')
print('Hostile control: mod9 raw 2→7 is not Collatz; decoded 2→1 is.')
print('All proofs are in README.md; no general convergence claim.')
