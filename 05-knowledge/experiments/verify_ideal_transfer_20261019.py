#!/usr/bin/env python3
"""Tiny fixed-instance controls for the ideal-transfer lemmas; not a census."""
from fractions import Fraction as F
from math import comb
from itertools import combinations
import json


def layers(p):
    n=len(p);ans=[[] for _ in range(n+1)]
    for s in range(1<<n):
        if all(not(s>>i&1) or p[i]&s==p[i] for i in range(n)):
            ans[s.bit_count()].append(s)
    return ans


def kernel(p,k,length,sortkey=None):
    lay=layers(p);rows=lay[k];cols=lay[k+length]
    if sortkey:
        rows=sorted(rows,key=sortkey);cols=sorted(cols,key=sortkey)
    K=[]
    for I in rows:
        v={I:1}
        for _ in range(length):
            w={}
            for s,z in v.items():
                for i in range(len(p)):
                    if not(s>>i&1) and p[i]&s==p[i]:
                        t=s|1<<i;w[t]=w.get(t,0)+z
            v=w
        K.append([v.get(J,0) for J in cols])
    return rows,cols,K


def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matrix_rank(M):
    A=[[F(x) for x in row] for row in M];rank=0
    for j in range(len(A[0])):
        pivot=next((i for i in range(rank,len(A)) if A[i][j]),None)
        if pivot is None:continue
        A[rank],A[pivot]=A[pivot],A[rank];z=A[rank][j];A[rank]=[v/z for v in A[rank]]
        for i in range(len(A)):
            if i!=rank:
                z=A[i][j];A[i]=[v-z*w for v,w in zip(A[i],A[rank])]
        rank+=1
    return rank


def incomparability_edges(p):
    return [(i,j) for i,j in combinations(range(len(p)),2) if not(p[i]>>j&1 or p[j]>>i&1)]


def connected(n,edges):
    reached={0}
    while True:
        new=reached|{y for x,y in edges if x in reached}|{x for x,y in edges if y in reached}
        if new==reached:return len(reached)==n
        reached=new


def two_chains(m):
    return tuple([(1<<i)-1 for i in range(m)]+[((1<<i)-1)<<m for i in range(m)])


def lag_two(m):
    return tuple([((1<<i)-1)|(((1<<max(0,i-1))-1)<<m) for i in range(m)]
                 +[(((1<<i)-1)<<m)|((1<<max(0,i-1))-1) for i in range(m)])


def run():
    sharp=[]
    for D in range(2,8):
        p=two_chains(D);r=D//2;rows,cols,K=kernel(p,r,D)
        C=comb(D,r)
        assert min(min(row) for row in K)>0
        cross=max(F(K[i][j]*K[k][l],K[i][l]*K[k][j])
                  for i in range(len(rows)) for k in range(len(rows))
                  for j in range(len(cols)) for l in range(len(cols)))
        assert cross==C*C
        U=0
        for I in rows:U|=I
        assert U.bit_count()-r==r
        assert all(I&~J==0 for I in rows for J in cols)
        # Delay sharpness in a (D+1)-antichain.
        a=(0,)*(D+1);I=1;J=((1<<(D+1))-1)^1
        rr,cc,H=kernel(a,1,D-1)
        assert H[rr.index(I)][cc.index(J)]==0
        sharp.append({'D':D,'central_binomial':C,'max_cross_ratio':str(cross),
                      'contraction_bound':str(F(C-1,C+1))})
    m=6;p=lag_two(m);edges=incomparability_edges(p)
    deg=[sum(i==a or i==b for a,b in edges) for i in range(2*m)]
    assert max(deg)==3 and connected(2*m,edges)
    assert all(connected(2*m,[e for e in edges if e!=deleted]) for deleted in edges)
    key=lambda I:(I&((1<<m)-1)).bit_count()
    _,_,T=kernel(p,4,2,key)
    expected=[[1,1,0],[1,2,1],[0,1,1]]
    assert T==expected and matrix_rank(T)==2
    _,_,T2=kernel(p,4,4,key)
    assert T2==matmul(T,T)==[[2,3,1],[3,6,3],[1,3,2]]
    assert min(min(row) for row in T2)>0 and matrix_rank(T2)==2
    _,_,K3=kernel(p,4,3,key)
    assert min(min(row) for row in K3)>0
    print(json.dumps({'status':'fixed controls passed; all-size conclusions use the proofs',
                      'sharp_contraction_examples':sharp,
                      'singular_bridgeless_example':{'rail_length':m,'vertices':2*m,
                        'max_incomparability_degree':max(deg),'incomparability_bridges':0,
                        'T':T,'T_squared':T2,'rank_T':2,'rank_T_squared':2}},indent=2))

if __name__=='__main__':run()
