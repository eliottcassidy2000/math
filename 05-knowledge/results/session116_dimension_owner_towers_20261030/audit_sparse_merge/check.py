#!/usr/bin/env python3
"""Independent C1--C5 audit. Explicit Baez convention, no imported checker."""
from functools import cache
from collections import Counter
from fractions import Fraction
from math import gcd
import json, hashlib, pathlib

@cache
def sign(d,i,j):
    if d == 0: return 1
    h=1<<(d-1); u=i%h; v=j%h; c=1 if u==0 else -1
    if i<h:
        return sign(d-1,u,v)*(1 if j<h else c)
    return sign(d-1,v,u)*(1 if j<h else -c)

@cache
def predicted(d,i,j):
    if i==0 or j==0: return 0
    h=1<<(d-1);u=i%h;v=j%h
    if (i<h)==(j<h): return 2*predicted(d-1,u,v)
    if u==0 or v==0 or u==v: return 0
    return h-4-2*predicted(d-1,u,v)

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def rankmod(A,p):
    A=[[v%p for v in row] for row in A]; r=0
    for c in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][c]),None)
        if pivot is None: continue
        A[r],A[pivot]=A[pivot],A[r]
        z=pow(A[r][c],-1,p); A[r]=[(v*z)%p for v in A[r]]
        for i in range(r+1,len(A)):
            z=A[i][c]
            A[i]=[(v-z*w)%p for v,w in zip(A[i],A[r])]
        r+=1
    return r

rows=[];prior=None
for d in range(1,8):
    n=1<<d; h=n//2; S=[[sign(d,i,j) for j in range(n)] for i in range(n)]
    census=Counter(); rectangles=0
    for i in range(n):
        for j in range(i+1,n):
            t=i^j;k=0;rc=0
            for v in range(n):
                w=v^t
                if v>w: continue
                eps=S[j][v]*S[i][w]*S[j][w]*S[i][v]
                k+=eps==1
                coefficient=2*(S[i][v]*S[j][w]+S[i][w]*S[j][v])
                require(coefficient in (-4,0,4),'rectangle magnitude')
                rc+=coefficient!=0
            require(rc==k,'rectangle/block identification')
            require(k==predicted(d,i,j),f'recurrence {d,i,j,k}')
            census[k]+=1; rectangles+=rc
    require(sum(census.values())==n*(n-1)//2,'pair total')
    if d>=3:
        require(census[0]==(6*d-5)*(1<<(d-2))+2,'zero census')
        require(set(census)==set(range(0,h-3,4)),'support')
        require(rectangles==(h-4)*(h-1)*(h-2),'first moment')
        rec=Counter({0:3*h-2})
        for k,c in prior.items():
            rec[2*k]+=2*c
            count=c-(h-1 if k==0 else 0)
            if count:rec[h-4-2*k]+=2*count
        require(census==rec,'polynomial recurrence')
    rows.append({'dimension':n,'pairs':n*(n-1)//2,'census':dict(sorted(census.items())),
                 'violated_rectangles':rectangles,'sign_norm_variance':str(Fraction(16*rectangles,n**4))})
    prior=census

# Direct modular elimination, not block formula, for the named odd-prime support.
ranks=[]
for d in (4,5,6):
    n=1<<d
    for p in (3,5,7,11):
        M=[[0]*n for _ in range(n)]
        for v in range(n):
            M[1^v][v]+=sign(d,1,v)
            M[10^v][v]+=(p-1)*sign(d,10,v)
        r=rankmod(M,p)
        require(r==3*n//4,'odd-prime rank')
        ranks.append([n,p,r])
# Exhaustive individual 2x2-block kernels for composite moduli, zero and signed coefficients.
block_tests=0
for a,b in ((0,2),(3,0),(1,1),(1,-1),(2,4),(3,5),(-2,3)):
    g=gcd(a,b)
    for s in (-1,1):
        for u in (-1,1):
            for m in range(1,13):
                actual=sum((a*x+b*s*y)%m==0 and (b*u*x+a*y)%m==0 for x in range(m) for y in range(m))
                expected=gcd(m,g)*gcd(m,(a*a-b*b*s*u)//g)
                require(actual==expected,'Smith kernel formula'); block_tests+=1
print(json.dumps({'status':'PASS','scope':'Pinned incoming C1-C5 only; Baez multiplication explicitly implemented',
 'source_commit':'7f70dcd6012439b7cf75fb43f26850f3807126b9','dimensions':rows,
 'modular_ranks_N_p_rank':ranks,'exhaustive_block_kernel_tests':block_tests,
 'checker_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},indent=2))
