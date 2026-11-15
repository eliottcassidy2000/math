#!/usr/bin/env python3
"""Finite checks of the binary split, derivative weights and affine collisions."""
from collections import Counter
from functools import lru_cache
from itertools import product
from random import Random

@lru_cache(None)
def p(n):
    if not n:
        return 0
    w, bit = divmod(n, 2)
    return p(w) + (0 if bit else w + (w & 1))

def rev(n, width):
    out = 0
    for _ in range(width):
        n, bit = divmod(n, 2)
        out = 2*out + bit
    assert n == 0
    return out

def split(n, n1, v):
    u0,b0 = divmod(n,v)
    u1,b1 = divmod(n1,v)
    assert abs(u1-u0) <= 1
    a = v
    if u0 != u1:
        u = min(u0,u1)
        j = (u ^ (u+1)).bit_length()-1
        a *= 1 << (j+1)
    t,b = divmod(n,a)
    t1,b1 = divmod(n1,a)
    assert t == t1
    return a,t,b,b1

def checksum_coeff(a,b,b1):
    width = a.bit_length()-1
    d = rev(b1,width)-rev(b,width)
    boundary = int(2*b1<a)-int(2*b<a) if a>1 else 0
    return -d, p(b1)-p(b), boundary

def f(x,y):
    return 100*x*y+89*x+11*(p(x)+p(y)+1+(x&1))

def test_signed():
    for r in range(9):
        signatures = {}
        for s in product((-1,0,1),repeat=r):
            h = sum(c<<i for i,c in enumerate(s))
            d = sum(c<<(r-1-i) for i,c in enumerate(s))
            assert (h,d) not in signatures
            signatures[h,d] = 2**s.count(0)
        counts = Counter((u-v,rev(u,r)-rev(v,r))
                         for u in range(1<<r) for v in range(1<<r))
        assert counts == signatures
        assert sum(w*w for w in counts.values()) == 6**r

def test_affine():
    rg = Random(20261115)
    for _ in range(80000):
        xx,yy = rg.randrange(10),rg.randrange(10)
        vx,vy = 1<<rg.randrange(10),1<<rg.randrange(10)
        x,y = rg.randrange(1,(1<<xx)+1),rg.randrange(1,(1<<yy)+1)
        x1 = x+rg.randint(-min(vx,x-1),vx)
        y1 = y+rg.randint(-min(vy,y-1),vy)
        a,u,b,b1=split(x,x1,vx)
        c,v,d,d1=split(y,y1,vy)
        sx,hx,px=checksum_coeff(a,b,b1)
        sy,hy,py=checksum_coeff(c,d,d1)
        assert p(x1)-p(x) == sx*u+hx+px*(u&1)
        assert p(y1)-p(y) == sy*v+hy+py*(v&1)
        alpha = 100*a*(d1-d)+11*sx
        beta = 100*c*(b1-b)+11*sy
        const = (100*(b1*d1-b*d)+89*(x1-x)
                 +11*(hx+hy+px*(u&1)+py*(v&1)+(x1&1)-(x&1)))
        assert f(x1,y1)-f(x,y)==alpha*u+beta*v+const
        ar,br=alpha-89*a,beta+89*c
        ct=const+89*(d1-b1)+11*((y1&1)-(x1&1))
        assert f(y1,x1)-f(x,y)==ar*u+br*v+ct
        assert ar*br!=0
        if (x,y)!=(x1,y1):
            assert alpha*beta != 0
    # Exhaust all small adjacent-quotient cases, including through a carry.
    for a in (1,2,4,8,16):
        for x in range(1,256):
            for x1 in range(max(1,x-a),x+a+1):
                b,u,v,w = split(x,x1,a)
                assert x == b*u+v and x1 == b*u+w

def test_orbit():
    for x in range(1,80):
        for y in range(1,80):
            assert f(2*x,y)==f(x,2*y+1)
            assert f(8*x,y)==f(x,8*y+7)

if __name__ == "__main__":
    test_signed()
    test_affine()
    test_orbit()
    print("PASS: binary derivative weights, affine coefficients and fundamental transport")
