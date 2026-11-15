#!/usr/bin/env python3
"""Exact digit/lattice identities and target-aware physical carrier inversion."""
from collections import defaultdict
from math import isqrt
from random import Random

def checksum(n):
    q = 0
    while n:
        b, n = n & 1, n >> 1
        if not b:
            q += n + (n & 1)
    return q

def rev(n, bits):
    r = 0
    for _ in range(bits):
        r, n = (r << 1) + (n & 1), n >> 1
    return r

def row(n, c, k, phi):
    return c*n + k*checksum(n) + phi*(n & 1)

def F(x, y):
    return 100*x*y + 89*x + 11*(checksum(x)+checksum(y)+1+(x & 1))

def ceildiv(n, q):
    return -(-n//q)

def hosts(T):
    """Yield every positive (x,y) in the physical carrier F with F(x,y)=T."""
    k = 11
    for a in range(1, isqrt(T//100)+1):
        pa = checksum(a)
        for x_fixed in (True, False):
            if x_fixed:
                c, phi = 100*a, 0
                Q = T - 89*a - k*(pa+1+(a & 1))
            else:
                c, phi = 100*a+89, 11
                Q = T - k*(pa+1)
            lo = max(a + (not x_fixed), ceildiv(Q-phi, c+k))
            hi = Q//c
            if lo > hi:
                continue
            bbits = max(1, (isqrt(T)//a).bit_length()-1)
            K = 1 << bbits
            residual_max = (c+k)*(K-1)+phi+k
            for w in range(lo//K, hi//K+1):
                L = (c*K+k*(K-1))*w + k*checksum(w)
                if w:
                    e0 = max(0, ceildiv(L-Q, k*w))
                    e1 = min(K-1, (L+residual_max-Q)//(k*w))
                    ns = (K*w+rev(e,bbits) for e in range(e0,e1+1))
                else:
                    ns = range(lo,min(hi,K-1)+1)
                for n in ns:
                    if lo <= n <= hi and row(n,c,k,phi)==Q:
                        yield (a,n) if x_fixed else (n,a)

def check_lattice(cases=100000):
    rng = Random(731)
    for _ in range(cases):
        s, j = rng.randint(1,18), rng.randint(1,16)
        y = rng.randrange(1 << (s-1))
        x = rng.randrange(1 << (s-1), 1 << s)
        x = (x & ~1) | (1-(y & 1))
        if x < (1 << (s-1)): # s=1: only possible x=1, y=0
            raise AssertionError((s,x,y))
        w = rng.randrange(1<<18)
        t = rng.randrange(1<<j)
        c,k,ph = rng.randint(1,301),rng.randint(1,31),rng.randint(-20,20)
        u = (w << (s+j))+(x << j)+t
        v = (w << (s+j))+(y << j)+t
        a,d,r = x-y,rev(x,s)-rev(y,s),rev(t,j)
        eps = (x & 1)-(y & 1)
        right = ((c+k)*(1<<j)-k)*a+k*(checksum(x)-checksum(y))
        left = k*(d*w+a*r+(w & 1)-eps*(1-(r & 1)))
        if row(u,c,k,ph)-row(v,c,k,ph) != right-left:
            raise AssertionError((s,j,x,y,w,t))
    print(f"Passed {cases} signed digit-lattice checks")

def check_decoder():
    upper = 130000
    all_values = defaultdict(set)
    for x in range(1,upper//100+1):
        for y in range(1,upper//(100*x)+1):
            T = F(x,y)
            if T <= upper:
                all_values[T].add((x,y))
    test_points = set(all_values)
    rng = Random(239)
    test_points.update(rng.randrange(1,upper+1) for _ in range(250))
    for T in test_points:
        out=list(hosts(T))
        assert len(out)==len(set(out))
        assert set(out)==all_values[T], (T, out, all_values[T])
    for _ in range(15):
        x,y = rng.randint(1,1000),rng.randint(1,1000)
        assert (x,y) in set(hosts(F(x,y)))
    print(f"Passed complete target replay at {len(test_points)} checked target values plus 15 larger witnesses")


def check_orientation_and_rays():
    rows = [(18,3,2),(189,11,11),(2,1,0)]
    for c,k,phi in rows:
        fibers=defaultdict(list)
        for n in range(1<<15):
            v=row(n,c,k,phi)
            for m in fibers[v]:
                assert rev(n,15)>rev(m,15), (n,m,c,k,phi)
            fibers[v].append(n)
    for j in range(1,11):
        B=1<<j
        for r in range(0,B,2):
            t=rev(r,j)
            w=7*B-r
            assert row((2*w+1)*B+t,18,3,2)==row(2*w*B+t,18,3,2)
            if (200*B-8-11*r)%7==0:
                w=(200*B-8-11*r)//7
                assert w%2==0
                assert row((16*w+13)*B+t,189,11,11)==row((16*w+2)*B+t,189,11,11)
    assert row(2,1,1,3)==row(1,1,1,3)
    print("Physical collision orientation and persistent primitive rays passed")

if __name__=='__main__':
    check_lattice()
    check_decoder()
    check_orientation_and_rays()
