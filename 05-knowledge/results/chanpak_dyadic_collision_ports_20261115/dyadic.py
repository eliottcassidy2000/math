#!/usr/bin/env python3
"""Exact controls for the dyadic-rank certificate and a physical test box."""
from itertools import product
from collections import Counter
from math import exp,pi

def p(n):
    if n==0:return 0
    u,r=divmod(n,2)
    return p(u)+(0 if r else u+(u%2))

def coefficient(d):
    return 1 if d==1 else 1<<(d-2)

def q_bits(n):
    s=[(n>>i)&1 for i in range(n.bit_length())]
    return sum((1-s[i])*s[j]*coefficient(j-i) for j in range(len(s)) for i in range(j))

def pairs(n,d):
    w=d+1
    return [(i,i+w) for i in range(n-w) if i%(2*w)<w]

def rank_check():
    tested=0
    for d in range(1,7):
        M=1<<d
        for n in range(1,33):
            rows=pairs(n,d)
            A=[x for x,y in rows];B=[y for x,y in rows]
            assert len(set(A+B))==2*len(A)
            assert len(A)*2>=max(0,n-d-1)
            matrix=[]
            for i,a in enumerate(A):
                row=[]
                for j,b in enumerate(B):
                    c=coefficient(abs(a-b))%M
                    if i==j: assert c==M//2
                    if i<j: assert c==0
                    row.append(c)
                matrix.append(row)
            if len(rows)<=9:
                # Explicit orthogonality: for every difference of binary rows,
                # some column contributes -1, canceling the full inner product.
                for delta in product((-1,0,1),repeat=len(rows)):
                    if any(delta):
                        assert any(sum(delta[i]*matrix[i][j] for i in range(len(rows)))%M==M//2 for j in range(len(rows)))
            tested+=1
    return tested

def modulus_diagnostic():
    # On a whole aligned interval with a fixed external high prefix, and
    # x fixed modulo 2^(m+1), compare actual integer histograms.
    T=1<<15
    out=[]
    for m in range(1,6):
        Q=1<<m;M=2*Q
        stats=[Counter() for _ in range(M)]
        for x in range(T,2*T):
            stats[x%M][p(x)//2%Q]+=1
        lo=min(row[j] for row in stats for j in range(Q))
        hi=max(row[j] for row in stats for j in range(Q))
        out.append((m,lo,hi))
    return out

def carrier(x,y):
    return 100*x*y+89*x+11*(p(x)+p(y)+1+(x&1))

def physical_box():
    # Verify the reduction to the two checksum histograms independently
    # of the direct modular output counts.
    r=7;s=8;m=3;Q=1<<m;M=2*Q
    X=range(3<<r,4<<r);Y=range(5<<s,6<<s)
    def hist(I):return Counter((x%M,p(x)//2%Q) for x in I)
    hx,hy=hist(X),hist(Y)
    C=Counter()
    for (x,u),a in hx.items():
        for (y,v),b in hy.items():
            C[(100*x*y+89*x+11*(2*u+2*v+1+(x&1)))%M]+=a*b
    E=Counter(carrier(x,y)%M for x in X for y in Y)
    assert C==E and len(E)==Q and all(k%2 for k in E)
    return sorted(E.items())

if __name__=='__main__':
    assert all(2*q_bits(n)==p(n) for n in range(1<<13))
    print(rank_check(),'dyadic rank templates')
    print('one-tape histogram extrema',modulus_diagnostic())
    print('carrier box',physical_box())
