#!/usr/bin/env python3
from pair_structural import C7, ends
from math import factorial
from functools import cache
@cache
def S(n,k):
    if n==k==0:return 1
    if k==0 or n==0:return 0
    return k*S(n-1,k)+S(n-1,k-1)
def skew(n):
    T,U,V,p=C7
    c=[i for i in range(7) for _ in range(n if i==4 else 1)]
    R=tuple(sum(1<<b for b,j in enumerate(c) if (j!=i and T[i]>>j&1) or (j==i and a<b))
            for a,i in enumerate(c))
    A=sum(1<<i for i,j in enumerate(c) if U>>j&1)
    B=sum(1<<i for i,j in enumerate(c) if V>>j&1);N=len(R);F=(1<<N)-1
    e=ends(R);s=ends(tuple(F^r^(1<<i) for i,r in enumerate(R)))
    def D(U,V):
        return sum(sum(e[m][i] for i in range(N) if U>>i&1)*
                   sum(s[F^m][i] for i in range(N) if V>>i&1) for m in range(1,F))
    return D(B,A)-D(A,B)
if __name__=='__main__':
    d=[]
    for n in range(1,9):
        x=skew(n)
        q,r=divmod(x-sum(d[j-1]*factorial(j)*S(n,j) for j in range(1,n)),factorial(n))
        assert r==0
        d.append(q)
        print(n,x)
    assert d==[2,17,12,0,0,0,0,0]
    print('PASS: complete fixed-vertex slot coefficients',d)
