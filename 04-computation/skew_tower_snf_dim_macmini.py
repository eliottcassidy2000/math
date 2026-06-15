#!/usr/bin/env python3
"""
Independent confirmation that dim C(H_{2^k}) = n/2 via:
  (i) direct GF(2) rank of B=(J-H)/2  (already done elsewhere, repeated here), and
  (ii) the integer Smith Normal Form of M = H (= I+S, the tournament gauge matrix),
       whose mod-2 reduction controls the GF(2) rank of (J-H)/2.

Note (J-H)/2 = (J - I - S)/2... but more directly: over F_2, B=(J-H)/2 has
B_{ij} = 1 iff H_{ij}=-1. The GF(2) rank of B is the code dimension. We also
compute SNF(H) over Z to show the '2-rank' structure (number of odd elementary
divisors) and that exactly n/2 of the elementary divisors are odd (=1), which
forces rank_{F2}(H mod 2) = n/2; and B relates to H by B = (J - H)/2, an integer
matrix. We report the SNF and the F_2 ranks of both H mod 2 and B.

Pure-Python integer SNF (fraction-free, with full pivoting). No project reuse.
"""
from fractions import Fraction

def transpose(M):
    n=len(M); m=len(M[0]); return [[M[i][j] for i in range(n)] for j in range(m)]
def double(H):
    m=len(H); HT=transpose(H)
    top=[H[i]+H[i] for i in range(m)]
    bot=[[-HT[i][j] for j in range(m)]+HT[i][:] for i in range(m)]
    return top+bot
def skew_tower(k):
    H=[[1]]
    for _ in range(k): H=double(H)
    return H

def f2_rank(rows_bitmasks):
    pivots={}
    for v in rows_bitmasks:
        x=v
        for pb in sorted(pivots,reverse=True):
            if (x>>pb)&1: x^=pivots[pb]
        if x: pivots[x.bit_length()-1]=x
    return len(pivots)

def mat_to_f2_bitmasks(M):
    n=len(M); out=[]
    for i in range(n):
        v=0
        for j in range(n):
            if M[i][j]%2!=0: v|=(1<<j)
        out.append(v)
    return out

def B_matrix(H):
    n=len(H)
    return [[(1 - H[i][j])//2 for j in range(n)] for i in range(n)]  # (J-H)/2

def smith_normal_form(Mat):
    """Integer SNF, returns list of diagonal elementary divisors. Pure Python,
    standard algorithm with row/col ops and pivot-on-smallest-|entry|. OK for
    n up to ~64 given the matrix is well-conditioned (entries +-1)."""
    import copy
    A=[row[:] for row in Mat]
    n=len(A); m=len(A[0]) if n else 0
    divisors=[]
    def find_pivot(t):
        best=None
        for i in range(t,n):
            for j in range(t,m):
                if A[i][j]!=0:
                    if best is None or abs(A[i][j])<abs(A[best[0]][best[1]]):
                        best=(i,j)
        return best
    t=0
    import math
    while t<min(n,m):
        p=find_pivot(t)
        if p is None: break
        pi,pj=p
        # move pivot to (t,t)
        A[t],A[pi]=A[pi],A[t]
        for r in range(n): A[r][t],A[r][pj]=A[r][pj],A[r][t]
        # reduce until column t and row t are cleared
        done=False
        while not done:
            done=True
            # clear column t
            for i in range(t+1,n):
                if A[i][t]!=0:
                    q=A[i][t]//A[t][t]
                    for j in range(m): A[i][j]-=q*A[t][j]
                    if A[i][t]!=0:
                        A[t],A[i]=A[i],A[t]; done=False
            # clear row t
            for j in range(t+1,m):
                if A[t][j]!=0:
                    q=A[t][j]//A[t][t]
                    for i in range(n): A[i][j]-=q*A[i][t]
                    if A[t][j]!=0:
                        for r in range(n): A[r][t],A[r][j]=A[r][j],A[r][t]
                        done=False
        # ensure pivot divides everything in the lower-right block
        d=A[t][t]
        bad=False
        for i in range(t+1,n):
            for j in range(t+1,m):
                if A[i][j]%d!=0:
                    # add row i to row t
                    for jj in range(m): A[t][jj]+=A[i][jj]
                    bad=True; break
            if bad: break
        if bad:
            continue
        divisors.append(abs(A[t][t]))
        t+=1
    return divisors

if __name__=="__main__":
    import sys
    from collections import Counter
    # SNF over Z only for small orders (8,16); F2-rank for ALL orders.
    for k in (3,4):
        n=1<<k
        H=skew_tower(k)
        B=B_matrix(H)
        rB=f2_rank(mat_to_f2_bitmasks(B))
        rHmod2=f2_rank(mat_to_f2_bitmasks(H))
        snf=smith_normal_form(H)
        odd_divs=sum(1 for d in snf if d%2==1)
        sig=Counter(snf)
        print(f"order {n}: F2rank(B=(J-H)/2)={rB} (n/2={n//2})  F2rank(H mod2)={rHmod2}", flush=True)
        print(f"          SNF(H) multiset = {dict(sorted(sig.items()))}  #odd divisors={odd_divs}", flush=True)
    # orders 32, 64: F2 rank only (SNF over Z too slow in pure python)
    for k in (5,6):
        n=1<<k
        H=skew_tower(k)
        B=B_matrix(H)
        rB=f2_rank(mat_to_f2_bitmasks(B))
        print(f"order {n}: F2rank(B=(J-H)/2)={rB} (n/2={n//2})   [SNF over Z skipped for speed]", flush=True)
