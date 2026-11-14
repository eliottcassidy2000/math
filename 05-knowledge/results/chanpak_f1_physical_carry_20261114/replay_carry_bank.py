"""Literal DP: all labelled punctures through order 4, then sample large ones."""
from itertools import combinations
from random import Random

def counts(rows):
    n=len(rows)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,c in enumerate(d[s]):
            if not c:
                continue
            x=rows[i]&~s
            while x:
                b=x&-x
                x^=b
                d[s|b][b.bit_length()-1]+=c
    return [1]+[sum(z) for z in d[1:]]

def cases(n, masks):
    P=list(combinations(range(n),2))
    for m in masks:
        T=[0]*n
        for j,(a,b) in enumerate(P):
            if m>>j&1:
                a,b=b,a
            T[a]|=1<<b
        yield T

def extension(T,U):
    n=len(T)
    return [r | (((U>>i&1)^1)<<n) for i,r in enumerate(T)]+[U]

def run():
    expected={2:{(1,3),(3,5)},4:{(1,5),(5,9)}}
    seen={2:set(),4:set()}
    checked=0
    rng=Random(20261114)
    for n in range(1,9):
        rows=range(1<<(n*(n-1)//2)) if n<=4 else (
          rng.randrange(1<<(n*(n-1)//2)) for _ in range(12))
        for T in cases(n,rows):
            ht=counts(T)
            full=(1<<n)-1
            masks=range(1<<n) if n<=5 else (
                rng.randrange(1<<n) for _ in range(10))
            for U in masks:
                W=full^U
                n1=counts(extension(T,U))[-1]
                n0=ht[-1]
                d=n1-n0
                k=sum((T[i]&W).bit_count() for i in range(n) if U>>i&1)
                base=ht[U]*ht[W]
                assert d>=0 and d%2==0
                assert (d==0)==(k==0)
                if d:
                    lg=d.bit_length()-1
                    assert 2*k<=d and k*base<=d
                    assert n0-base<=d*lg
                    assert n0<=d*(1+lg)
                    assert n0*k<=d*(k+1)
                    if d>=6:
                        assert n0<=2*d-3
                if d in expected:
                    assert (n0,n1) in expected[d]
                    seen[d].add((n0,n1))
                checked+=1
    for r in range(2,8):
        W=[sum(1<<j for j in range(i+1,r+1)) | (0 if i==0 else 1<<(r+1))
           for i in range(r+1)]+[1]
        b=len(W)-1
        # a is the added old cut tail, and it only beats b.
        T=extension(W,1<<b)
        G=extension(T,1<<len(W))
        h=counts(T)[-1]
        d=counts(G)[-1]-h
        assert (h,d)==(2*(2+2**r)-3,2+2**r)
    assert seen==expected
    print(f"PASS: {checked} literal extensions; least-gap pairs {seen}")

if __name__=='__main__':
    run()
