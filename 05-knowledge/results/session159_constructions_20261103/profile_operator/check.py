#!/usr/bin/env python3
"""Finite controls for the triangle operator and its mod-9 certificate."""
from itertools import product
from math import comb

def check(ok,msg):
    if not ok: raise RuntimeError(msg)

def tournament(n,mask):
    arcs=[[False]*n for _ in range(n)]
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    for b,(i,j) in enumerate(edges):
        if mask>>b&1: arcs[j][i]=True
        else: arcs[i][j]=True
    return arcs

def G(a):
    n=len(a); b=[[False]*(n+2) for _ in range(n+2)]
    for i in range(n):
        for j in range(n): b[i][j]=a[i][j]
        b[i][n]=True; b[n+1][i]=True
    b[n][n+1]=True
    return b

def h_subsets(a):
    n=len(a); N=1<<n; dp=[[0]*n for _ in range(N)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,N):
        for j in range(n):
            if mask>>j&1:
                rem=mask^(1<<j)
                if rem: dp[mask][j]=sum(dp[rem][i] for i in range(n) if a[i][j])
    return [1]+[sum(row) for row in dp[1:]]

def F(a,r):
    h=h_subsets(a); n=len(a); total=0
    for assignment in product(range(r),repeat=n):
        slots=[0]*r
        for i,j in enumerate(assignment):slots[j]|=1<<i
        v=1
        for s in slots:v*=h[s]
        total+=v
    return total

def transform(p):
    q=[0]*(len(p)+2)
    for j,v in enumerate(p):
        if j: q[j-1]+=j*(j-1)*v
        q[j]+=j*(j+2)*v
        q[j+1]+=(2*j+1)*v
        q[j+2]+=v
    return q

def profile(a):
    n=len(a); vals=[F(a,r) for r in range(1,n+1)]
    # F_r=sum_j p_j (r)_j, recover Newton coefficients.
    vals=[0]+vals
    from math import factorial
    p=[sum((-1)**(j-i)*comb(j,i)*vals[i] for i in range(j+1))//factorial(j) for j in range(n+1)]
    return p

cases=0
for n in range(1,5):
 for mask in range(1<<(n*(n-1)//2)):
    a=tournament(n,mask); b=G(a); p=profile(a);q=transform(p)
    for r in range(1,4):
        lhs=F(b,r)
        rhs=r*F(a,r+1)+r*r*F(a,r)-(r*F(a,r-1) if r>1 else 0)
        calc=sum(q[k]*__import__('math').prod(range(r-k+1,r+1)) for k in range(min(r,len(q)-1)+1))
        check(lhs==rhs==calc,(n,mask,r,lhs,rhs,calc))
    cases+=1
print('all triangle identities verified:', cases, 'input tournaments; r=1..3')
p=[0,1]; expected=[1,3,15,123,1479,23955,496479]
for m,h in enumerate(expected):
    check(p[1]==h,('iterated count',m,p[1],h));p=transform(p)
print('iterated G Hamilton counts:',expected)
# A closed 8-dimensional state over Z/9Z, independent of p_9 because 8*9=0.
def step(v):
    a=[0]+list(v)+[0]
    return tuple(((k+1)*k*a[k+1]+k*(k+2)*a[k]+(2*k-1)*a[k-1]+(a[k-2] if k>=2 else 0))%9 for k in range(1,9))
v=(1,0,0,0,0,0,0,0);states=[]
for i in range(12): states.append(v);v=step(v)
check(states[5]==states[11], 'mod 9 cycle mismatch')
check(all(v[0] in (3,6) for v in states[1:11]), 'valuation mismatch')
print('mod-9 8-coordinate states, steps 0..11:')
for i,v in enumerate(states):print(i, *v)
print('state 5 = state 11; v_3(h)=1 for every iteration m>=1')

def F_weighted(a,r,marked,t):
    n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,N):
        for j in range(n):
            if mask>>j&1:
                rem=mask^(1<<j)
                if rem:dp[mask][j]=sum(dp[rem][i]*(t if (i,j) in marked else 1) for i in range(n) if a[i][j])
    h=[1]+[sum(row) for row in dp[1:]];total=0
    for assn in product(range(r),repeat=n):
        slots=[0]*r
        for i,j in enumerate(assn):slots[j]|=1<<i
        v=1
        for s in slots:v*=h[s]
        total+=v
    return total
for mask in range(8):
    a=tournament(3,mask); b=G(a)
    arcs=[(i,j) for i in range(3) for j in range(3) if a[i][j]]
    for marking in range(8):
        marked={e for i,e in enumerate(arcs) if marking>>i&1}
        for t in (2,3):
            for r in (1,2,3):
                lhs=F_weighted(b,r,marked,t)
                rhs=r*F_weighted(a,r+1,marked,t)+r*r*F_weighted(a,r,marked,t)
                if r>1:rhs-=r*F_weighted(a,r-1,marked,t)
                check(lhs==rhs,('weighted',mask,marking,t,r,lhs,rhs))
print('weighted literal controls: all order-3 orientations and arc markings, t=2,3, r=1,2,3')
