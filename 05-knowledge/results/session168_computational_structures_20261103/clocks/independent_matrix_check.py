#!/usr/bin/env python3
"""Independent F3 matrix verification; imports neither the field nor main checker."""
import json


def need(v,msg):
    if not v:raise RuntimeError(msg)


def eye(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def multiply(A,B):
    n=len(A)
    out=[[0]*n for _ in range(n)]
    for i in range(n):
        for k,a in enumerate(A[i]):
            if a:
                for j,b in enumerate(B[k]):
                    out[i][j]+=a*b
        out[i]=[x%3 for x in out[i]]
    return out


def mpow(A,e):
    out=eye(len(A))
    while e:
        if e&1:out=multiply(out,A)
        A=multiply(A,A)
        e//=2
    return out


def rank(A):
    A=[r[:] for r in A]
    n=len(A);r=0
    for j in range(n):
        pivot=next((i for i in range(r,n) if A[i][j]%3),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r]
        inv=A[r][j]%3
        A[r]=[(x*inv)%3 for x in A[r]]
        for i in range(n):
            if i!=r and A[i][j]%3:
                a=A[i][j]%3
                A[i]=[(x-a*y)%3 for x,y in zip(A[i],A[r])]
        r+=1
    return r


def vector(A,v):
    return tuple(sum(a*b for a,b in zip(row,v))%3 for row in A)


def check(k):
    # Coordinate order s,u,d. Literal C(s,u,d)=(u,s+u+d,s).
    C=((0,1,0),(1,1,1),(1,0,0))
    n=3*k
    A=[[0]*n for _ in range(n)]
    for i in range(k):
        for r in range(3):
            for c in range(3):A[3*i+r][3*i+c]+=C[r][c]
            A[3*i+r][3*((i+1)%k)+r]+=1
    A=[[x%3 for x in row] for row in A]
    I=eye(n);period=26*k
    seed=tuple(int(j==2) for j in range(n))
    need(mpow(A,period)==I,'matrix return')
    for p in (2,13)+( (3,) if k>1 else () ):
        B=mpow(A,period//p)
        need(B!=I and vector(B,seed)!=seed,'exact matrix and delta orders')
    fixed=[]
    q=1
    while q<=k:
        B=mpow(A,26*q)
        D=[[(B[i][j]-I[i][j])%3 for j in range(n)] for i in range(n)]
        nullity=n-rank(D)
        need(nullity==3*q,'fixed-space dimension')
        fixed.append({'period_dividing':26*q,'nullity':nullity,'fixed_states':3**nullity})
        q*=3
    need(rank([[(A[i][j]-I[i][j])%3 for j in range(n)] for i in range(n)])==n,'zero unique fixed point')
    return {'sites':k,'exact_order':period,'fixed_spaces':fixed}


if __name__=='__main__':
    print(json.dumps([check(k) for k in (1,3,9,27)],sort_keys=True,indent=2))
