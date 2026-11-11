#!/usr/bin/env python3
"""Inherited Session264 fusion carrier helpers; copied with attribution.
The replay counts literal paths independently of matrix multiplication.
"""
from itertools import combinations, product
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json

@lru_cache(None)
def hp(T):
    n=len(T); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    H=[1]+[0]*(size-1)
    for S in range(1,size):
        H[S]=sum(dp[S])
        for v in range(n):
            if not dp[S][v]:continue
            rem=T[v]&~S
            while rem:
                bit=rem&-rem; rem-=bit
                dp[S|bit][bit.bit_length()-1]+=dp[S][v]
    return tuple(H)

def op(T):return tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T)))

def rooted(T):
    """Root is 0; return weighted counts by actual path vertex length."""
    n=len(T);H=hp(T);full=(1<<n)-1;dp=[[0]*n for _ in range(1<<n)];dp[1][0]=1
    P=[0]*(n+1)
    for S in range(1,1<<n):
        P[S.bit_count()]+=sum(dp[S])*H[full^S]
        for v in range(n):
            if not dp[S][v]:continue
            rem=T[v]&~S
            while rem:
                bit=rem&-rem;rem-=bit;dp[S|bit][bit.bit_length()-1]+=dp[S][v]
    return P

def profile(T):
    Hs=hp(T);full=(1<<len(T))-1;R=full^1;U=rooted(T);V=rooted(op(T));h=Hs[R]
    Q=k=0;S=R
    while True:
        val=Hs[S]*Hs[R^S];Q+=val;k+=(-1)**S.bit_count()*val
        if not S:break
        S=(S-1)&R
    d=sum((-1)**(j+1)*x for j,x in enumerate(U));e=sum((-1)**(j+1)*x for j,x in enumerate(V))
    return dict(H=Hs[full],h=h,Q=Q,k=k,U=sum(U),V=sum(V),d=d,e=e)

def matrix(p):
    h=p['h']
    return [[h,p['U']-h,p['d']-h,p['H']-h],[0,p['Q'],0,p['V']-h],[0,0,p['k'],p['e']-h],[0,0,0,h]]

def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def from_matrix(M):
    h=M[0][0]
    return dict(H=M[0][3]+h,h=h,Q=M[1][1],k=M[2][2],U=M[0][1]+h,V=M[1][3]+h,d=M[0][2]+h,e=M[2][3]+h)

def fuse(A,B):
    """Identify roots 0, orient all A-root vertices toward all B-root vertices."""
    n=len(A);m=len(B);out=[0]*(n+m-1)
    for T,inds in [(A,list(range(n))),(B,[0]+list(range(n,n+m-1)))]:
        for i in range(len(T)):
            for j in range(len(T)):
                if T[i]>>j&1:out[inds[i]]|=1<<inds[j]
    for i in range(1,n):
        for j in range(n,n+m-1):out[i]|=1<<j
    return tuple(out)

def tournaments(n):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        out=[0]*n
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:out[a]|=1<<b
            else:out[b]|=1<<a
        yield tuple(out)

