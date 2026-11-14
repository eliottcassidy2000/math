"""Literal path/colour checks and exact weighted-grammar certificate."""
from functools import cache
from itertools import combinations, product
from fractions import Fraction as Q
from math import isqrt


def transitive(n):
    return [((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n)]


def parallel(rows, h):
    n = len(rows)
    out = transitive(n + 2*h)
    later = ((1 << (n + 2*h)) - (1 << (n+h)))
    for j, row in enumerate(rows):
        out[h+j] = (row << h) | later
    for i in range(h):
        a, b = i, h+n+i
        out[a] ^= 1 << b
        out[b] ^= 1 << a
    return out


def paths(rows):
    @cache
    def end(S,v):
        R = S ^ (1 << v)
        if not R:
            return 1
        return sum(end(R,u) for u,r in enumerate(rows)
                   if (R >> u)&1 and (r >> v)&1)
    @cache
    def H(S):
        return sum(end(S,v) for v in range(len(rows)) if (S >> v)&1) if S else 1
    return H


def whole(rows):
    return paths(rows)((1 << len(rows)) - 1)


def colours(rows,m):
    H = paths(rows)
    w = [0]*m
    total = 0
    for c in product(range(m),repeat=len(rows)):
        w[:] = [0]*m
        for i,j in enumerate(c):
            w[j] |= 1 << i
        weight = 1
        for mask in w:
            weight *= H(mask)
        total += weight
    return total


def explicit(n, mask):
    rows = transitive(n)
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if mask >> k & 1:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return rows


def grammar():
    def A(z): return z+z**3+z**4/(1-2*z)
    def U(z,u): return A(z)+z**2/(1-z**2)*u/(1-u)
    s,t = 119,125
    z3,z4,z5,z6 = (Q(9,25),Q(7,25),Q(11,50),Q(19,100))
    for k,z in zip((3,4,5,6),(z3,z4,z5,z6)):
        assert z**t*k**s>1
    u4,u5 = Q(2,5),Q(7,25)
    assert U(z4,u4)<u4 and U(z5,u5)<u5
    q3 = A(z3)+z3*(u4/(1-u4))+z6/(1-z3**2)*(u5/(1-u5))
    assert q3<1
    print('Grammar replay: F3 atom weight bound',q3,'=',float(q3))


def cross():
    known=[5,13,33,85,225,613,1713,4885]
    for n,N in enumerate(known):
        assert whole(parallel(transitive(n),2))==N==2*3**n+4*2**n-1
    checked=0
    for n in range(5):
        for mask in range(1<<(n*(n-1)//2)):
            B=explicit(n,mask)
            F1,F2,F3=tuple(colours(B,m) for m in (1,2,3))
            P=parallel(B,2)
            H=whole(P)
            assert H==2*F3+4*F2-F1
            if n: assert H%3==(F1%3 if n%2 else 0)
            if n<=3:
                for m in (1,2):
                    assert colours(P,m)>=2*m*colours(B,m+2)
            checked+=1
    print('Cross replay:',checked,'literal reservoirs through order four; interpolation seeds and colour inequality checked')
    B0=[18,8,3,5,14]
    assert tuple(colours(B0,m) for m in (1,2,3))==(13,104,471)
    B=[r|96 for r in B0]+[64,0]
    assert whole(parallel(B,2))==10129==7*1447
    assert all(1447%d for d in range(2,isqrt(1447)+1))
    for r in range(24):
        N=942*3**r+416*2**r-13
        assert (N%7==0)==(r%6 in (2,3))
        assert (N%3!=0)==(r%2==0)
        assert 3**r*942<=N+13<3**(r+1)*942
    print('7-prime witness and exact ternary padding decoder checked')


if __name__=='__main__':
    grammar()
    cross()
