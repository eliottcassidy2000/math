#!/usr/bin/env python3
"""Exact bounds and small examples for the edge-rooted satellite grid."""
from fractions import Fraction as F
from math import factorial, prod
from functools import cache

w=tuple(F(x,1000) for x in (1000,1000,1000,1000,1000,1000,1000,1071,1149,1233,1324,1421,1525,1637,1758,1887,2025,2174,2334,2505,2689,2886,3098,3326,3570,3832,4113,4415,4740,5088,5461,5862,6293,6755,7251,7784,8356,8970,9629,10336,11095,11910,12786,13725,14732,15816,16978,18225,19565,21003,22546,24203,25984,27891,29939,32142))
wa=F(33844,1000)
e0=F(27182818,10**7)
e1=F(271828183,10**8)
partial=sum(F(1,factorial(i)) for i in range(14))
assert e0<partial and partial+F(1,13*factorial(13))<e1

# The exact layer bounds also used in the native affine-threshold note.
f=1
for d in range(1,3138):
    f*=d
    assert f*e1.numerator**d<((d+5)*e1.denominator)**d
    if d<=59:
        assert f*(500*e1.numerator)**d<((510*d+1087)*e1.denominator)**d

def lower(d,k,bits):
    b=k+1-4*(d-1)-3
    bad=k*(d-1)+1
    num=factorial(k)*(d-1)**((d-1)*k)*b**bad
    den=d**((d-2)*k)*(b+1)**bad
    assert num>den<<bits

lower(28,1568,19663)
lower(3,18,53)

def ceil(x):
    return -((-x.numerator)//x.denominator)

def upper_grid(mate):
    m=len(w)
    n=m*m
    assert w==tuple(sorted(w)) and all(x>=1 for x in (*w,wa))
    scale=10**6
    apex=m*sum(wa/t for t in w)+5*wa
    f=ceil(apex/e0*scale)
    for i,x in enumerate(w):
        edge_backward=[x/y for y in w[:i]]
        edge_forward=[x/y for y in w[i+1:]]
        r=sum(edge_backward)+(m-1)*sum(edge_forward)+x/wa
        penalty_at_high=sum(max(y-1,0) for y in edge_backward+[x/wa])
        for a in range(m):
            total=r+(m-1-a)
            high=penalty_at_high
            if mate:
                if a==0:
                    total-=1 # lost next, within row
                    if i>0: # gained previous, across row boundary
                        v=x/w[i-1];total+=v;high+=v-1
                if a==m-1:
                    total+=1 # gained previous, within row
                    if i<m-1:total-=x/w[i+1] # lost next
            factor=total+F(1413,500)+F(1087,500)*x+F(1,50)*high
            f*=ceil(factor/e0*scale)
    assert f<(scale**(n+1)<<25627)

for mate in (False,True):
    upper_grid(mate)

def count(A):
    n=len(A)
    @cache
    def go(i,todo):
        if not todo:return 1
        return sum(go(j,todo^(1<<j)) for j in range(n) if todo>>j&1 and A[i][j])
    return sum(go(i,((1<<n)-1)^(1<<i)) for i in range(n))

def grid(V,edges,root):
    def edge(i,j):
        return frozenset((i,j)) in edges
    u,v=root
    P=[x for x in V if x!=v and edge(x,u)]
    Q=[x for x in V if x!=u and edge(x,v)]
    m=len(P)
    rows=[]
    col={}
    for i,p in enumerate(P):
        row=[]
        for j,q in enumerate(Q):
            sites=[x for x in V if edge(x,p) and edge(x,q)]
            assert len(sites)==1
            x=sites[0];row.append(x);col[x]=j
        if rows:
            if edge(rows[-1][-1],row[0]):
                row[0],row[1]=row[1],row[0]
            assert not edge(rows[-1][-1],row[0])
        rows.append(row)
    X=[x for row in rows for x in row]
    assert len(set(X))==m*m
    A=[[edge(x,y) for y in X] for x in X]
    for i in range(m*m):
        for j in range(m*m):
            s=sum(a and b for a,b in zip(A[i],A[j]))+A[i][j]
            expected=m*(i==j)+1-(i//m==j//m)-(col[X[i]]==col[X[j]])
            assert s==expected
    def orient(mate):
        return [[i!=j and (i<j) != (A[i][j] or (mate and abs(i-j)==1))
                 for j in range(m*m)] for i in range(m*m)]
    if m<=2:
        return [count(orient(mate)) for mate in (False,True)]
    # A simple exact Bregman certificate for the concrete HS controls.
    @cache
    def bf(d):
        scale=10000
        a,b=1,d*scale+1
        while b-a>1:
            c=(b+a)//2
            if c**d>=factorial(d)*scale**d:b=c
            else:a=c
        return b
    tops=[]
    for mate in (False,True):
        T=orient(mate)
        factors=[bf(sum(row)+1) for row in T]+[bf(m*m)]
        x=prod(factors); den=10000**len(factors)
        top=0
        while x>=den<<top:top+=1
        tops.append(top)
    return tops

def control(m):
    if m==1:
        V=list(range(5))
        E=[(i,(i+1)%5) for i in V]
        root=(0,1)
    elif m==2:
        V=list(range(10))
        E=[(i,(i+1)%5) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)]+[(i,5+i) for i in range(5)]
        root=(0,1)
    elif m==6:
        V=[(side,x,y) for side in (0,1) for x in range(5) for y in range(5)]
        E=[((side,x,y),(side,x,(y+side+1)%5)) for side,x,y in V]
        E += [((0,i,x),(1,j,(x+i*j)%5)) for i in range(5) for x in range(5) for j in range(5)]
        root=((0,0,0),(0,0,1))
    return grid(V,set(map(frozenset,E)),root)
assert control(1)==[1,1]
assert control(2)==[5,5]
print("grid counts: C5 1,1; Petersen 5,5; HS upper exponents",control(6))
print("degree 57: both greater than 2^19663 and less than 2^25627")

