#!/usr/bin/env python3
"""Exact nested-contact classification on all small outside tournaments."""
from itertools import combinations, permutations, product
from collections import Counter
from random import Random

def need(ok,message):
    if not ok: raise ArithmeticError(message)

def paths(A):
    n=len(A)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n): D[1<<i][i]=1
    for B in range(1,1<<n):
        for i,val in enumerate(D[B]):
            if val:
                for j in range(n):
                    if not (B>>j&1) and A[i][j]:
                        D[B|1<<j][j]+=val
    return sum(D[-1]) if n else 1

def make(r,mask,modes):
    n=r+2
    A=[[0]*n for _ in range(n)]
    def arc(i,j): A[i][j]=1
    arc(0,1)
    for k,(i,j) in enumerate(combinations(range(2,n),2)):
        if mask>>k&1: arc(i,j)
        else: arc(j,i)
    for i,g in enumerate(modes,2):
        if g==2: arc(0,i);arc(1,i)
        elif g==1: arc(0,i);arc(i,1)
        else: arc(i,0);arc(i,1)
    return A

def nested():
    hit=Counter()
    tested=0
    for r in range(5):
        for mask in range(1<<(r*(r-1)//2)):
            for modes in product(range(3),repeat=r):
                A=make(r,mask,modes)
                lower=paths(A)
                A[0][1],A[1][0]=0,1
                diff=(paths(A)-lower)//2
                D=[i for i,g in enumerate(modes,2) if g==1]
                need(diff>=len(D),"nested positivity")
                if len(D)==1:
                    t=D[0]
                    rem=[i for i in range(2,r+2) if i!=t]
                    order=sorted(rem,key=lambda i:-sum(A[i][j] for j in rem))
                    transitive=all(A[i][j] for i,j in combinations(order,2))
                    unit=transitive and all(A[t][i]<=A[1][j] and A[1][i]<=A[t][j]
                                            for i,j in combinations(order,2))
                else: unit=False
                need((diff==1)==unit,"nested tape classification")
                if unit: hit[lower]+=1
                tested+=1
    need(set(hit)=={1,3,9,13},"four lower counts")
    print(f"{tested} nested pairs; exact unit lower values {sorted(hit)}")

def cycles(A,vertices):
    for k in range(3,len(vertices)+1,2):
        for S in combinations(sorted(vertices),k):
            for tail in permutations(S[1:]):
                p=(S[0],*tail)
                if all(A[u][v] for u,v in zip(p,p[1:]+p[:1])):
                    yield p

def multiply(a,b):
    c=Counter()
    for S,x in a.items():
        for T,y in b.items():
            if not S&T:
                c[S|T]+=x*y
    return {S:x for S,x in c.items() if x}

def product_poly(C,activity):
    z={0:1}
    for p in C:
        k=sum(1<<i for i in p)
        z=multiply(z,{0:1,k:activity})
    return z

def overlaps():
    rng=Random(887)
    for _ in range(200):
        V=set(range(6)); A=[[0]*6 for _ in V]
        for i,j in combinations(sorted(V),2):
            if rng.randrange(2): A[i][j]=1
            else: A[j][i]=1
        left={0,1,2,3};right={1,2,3,4,5};mid=left&right
        allc=list(cycles(A,V))
        mixed=[c for c in allc if not set(c)<=left and not set(c)<=right]
        expected=product_poly(allc,2)
        got={0:1}
        for C,v in ((cycles(A,left),2),(cycles(A,right),2),
                    (cycles(A,mid),-2),(mixed,2)):
            got=multiply(got,product_poly(C,v))
        need(got==expected,"squarefree shared-cycle correction")
        need(sum(got.values())==paths(A),"packing/HP check")
    print("200 cyclic-overlap polynomial replays against ordinary Hamilton counts")
if __name__=="__main__":
    nested()
    overlaps()
