#!/usr/bin/env python3
"""Independent finite audit: no imports of production checkers; optimization-safe."""
from fractions import Fraction as F
from itertools import product, combinations
from math import factorial, comb, gcd
import json

def need(p,m):
    if not p: raise RuntimeError(m)
def comp(n,j):
    if j==1: yield (n,); return
    for a in range(n+1):
        for t in comp(n-a,j-1): yield (a,)+t

def legal(w):
    q=1
    for i,a in enumerate(w):
        q+=2 if a%3==0 else -1
        if q<=0: return i==len(w)-1 and q==0
    return False

def carry(w,b=3):
    C=D=0
    for a in w: C=b*C+2**D; D+=a
    return C,D

def native(n,b):
    x=b*n+1; a=(x&-x).bit_length()-1
    return x>>a,a

cylinder_count=0
for b in (3,5,7):
    for J in range(1,5):
        for w in product(range(1,5),repeat=J):
            C,D=carry(w,b); mod=2**(D+1)
            n=((2**D-C)*pow(b**J,-1,mod))%mod
            need(n%2==1,'odd source')
            for offset in (0,1,3):
                source=n+offset*mod; x=source
                for a in w:
                    x,actual=native(x,b); need(a==actual,'cylinder exactness')
                need((x<<D)==b**J*source+C,'affine endpoint')
            cylinder_count+=1

# Direct ordered words versus independent convolution of rational factors.
N=14; coefficient_count=0
for k in range(3):
    J=3*k+1
    sigs=[]
    for pos in combinations(range(J),k):
        w=tuple(3 if i in pos else 1 for i in range(J))
        if legal(w): sigs.append(tuple(a==3 for a in w))
    need(len(sigs)==comb(J,k)//J,'shape count')
    direct=[[0]*(N+1) for _ in range(4)]
    for D in range(J,N+1):
        for t in comp(D-J,J):
            w=tuple(a+1 for a in t)
            if legal(w):
                C,_=carry(w)
                for r in range(4): direct[r][D]+=C**r
    for r in range(4):
        ans=[0]*(N+1)
        for sig in sigs:
            for m in comp(r,J):
                scale=factorial(r)//__import__('math').prod(factorial(a) for a in m)*3**sum((J-i-1)*m[i] for i in range(J))
                arr=[1]+[0]*N
                for j,branch in enumerate(sig):
                    h=sum(m[j+1:]); out=[0]*(N+1)
                    for d in range(1,N+1):
                        if (d%3==0)==branch:
                            for t in range(N+1-d): out[t+d]+=arr[t]*2**(h*d)
                    arr=out
                for d in range(N+1): ans[d]+=scale*arr[d]
        need(ans==direct[r],f'coefficient agreement {k,r}')
        coefficient_count+=N+1

# Finite bank reconstruction and uniqueness of its highest pole.
def parts(n,lo=1):
    if n==0: yield (); return
    for d in range(lo,n+1):
        for t in parts(n-d,d): yield (d,)+t
banks=[]
for k in range(11):
    cores=[tuple(2*d+1 for d in ds) for ds in parts(k)]
    top=[mu for mu in cores if sum(mu)==3*k]
    need(top==[(3,)*k],'unique maximal pole order')
    A=F(factorial(3*k),3**k*factorial(k))*2**(3*k*k+k)
    if k:
        need(A==prev*(3*k-2)*(3*k-1)*2**(6*k-2),'leading amplitude ratio')
    prev=A; banks.append({'k':k,'cores':len(cores),'order':3*k+1})

# A native positive 5x+1 cycle, not merely a formal valuation word.
x=13; cycle=[]; digits=[]
for _ in range(3):
    cycle.append(x); x,a=native(x,5); digits.append(a)
need(x==13 and cycle==[13,33,83] and digits==[1,1,5],'hostile 5x+1 cycle')
x=1; unit=[]
for _ in range(2): unit.append(x); x,_=native(x,5)
need(x==1 and not set(unit)&set(cycle),'disjoint native cycles')
print(json.dumps({'status':'PASS','source_cylinders':cylinder_count,'coefficient_comparisons':coefficient_count,'banks':banks,'hostile_5x_plus_1':{'cycle':cycle,'digits':digits,'unit_cycle':unit}},indent=2))
