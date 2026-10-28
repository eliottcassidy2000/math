#!/usr/bin/env python3
"""Fixed-excess prime-depth diagonals, using only P through order9.
Independent ordinary-series implementation. Exact fractions, no third-party deps.
"""
from fractions import Fraction as F
import json
from math import gcd, factorial
from collections import Counter
H=15; R=6; N=2*H+1+R
P={
(1,1,1):F(1,3),(3,):F(2,3),
(1,)*5:F(11,5),(5,):F(4,5),
(1,)*6:F(43,3),(3,3):F(2,3),
(1,)*7:F(4057,21),(1,3,3):F(8,3),(7,):F(8,7),
(1,)*8:F(4000),(1,1,3,3):F(8),
(1,)*9:F(1229279,9),(1,1,1,3,3):F(48),(3,3,3):F(172,9),(9,):F(4,3)}
def mul(a,b):
    out=[F(0)]*(N+1)
    for i,x in enumerate(a):
        if not x: continue
        for j,y in enumerate(b[:N+1-i]):
            if y: out[i+j]+=x*y
    return out
def compose(a,deriv=False):
    out=[F(0)]*(N+1); cache={}
    for part,v in P.items():
        if deriv:
            if 1 not in part: continue
            v*=part.count(1); part=list(part); part.remove(1)
        term=[v]+[F(0)]*N
        for d in part:
            if d not in cache: cache[d]=[a[i//d] if i%d==0 else F(0) for i in range(N+1)]
            term=mul(term,cache[d])
        out=[x+y for x,y in zip(out,term)]
    return out
def step(a):
    s=compose(a);s[1]+=1; t=[F(1)]+[F(0)]*N
    for n in range(1,N+1): t[n]=sum(s[k]*t[n-k] for k in range(1,n+1))
    t[0]=0;return t
levels=[[F(0)]+[F(1)]*N]
for h in range(1,H+1): levels.append(step(levels[-1]))
a=levels[-1]; B=a.copy();B[0]=1
K=mul(mul(B,B),compose(a,deriv=True)); L=K[2:2+R+1]
# Independent pointed-Burnside identity K=1-B^2/Dpoint.
def oddparts(n,lo=1):
    if n==0: yield ()
    for k in range(lo,n+1,2):
        for p in oddparts(n-k,k): yield (k,)+p
def weight(p):
    e=sum((i-1)//2 for i in p)+sum(gcd(i,j) for k,i in enumerate(p) for j in p[k+1:])
    z=1
    for i,m in Counter(p).items(): z*=i**m*factorial(m)
    return F(2**e,z)
burn=[F(0)]*(N+1); point=burn.copy()
for n in range(R+3):
    for p in oddparts(n):
        w=weight(p);burn[n]+=w;point[n]+=w*2**len(p)
square=mul(burn,burn); ratio=[F(0)]*(N+1)
for n in range(R+3): ratio[n]=square[n]-sum(point[k]*ratio[n-k] for k in range(1,n+1))
for n in range(R+3):
    independent=(1 if n==0 else 0)-ratio[n]
    if independent!=K[n]: raise RuntimeError(('Pointed kernel mismatch',n,independent,K[n]))

D={h:[levels[h][2*h+1+r]-levels[h-1][2*h+1+r] for r in range(R+1)] for h in range(1,H+1)}
# Nilpotency filtration directly checks the single-Jordan-block carrier.
q=[F(0)]+L[1:]+[F(0)]*(N-R)
power=[F(1)]+[F(0)]*N
for j in range(R+2):
    if j<=R and (any(power[k] for k in range(j)) or power[j]!=4**j):
        raise RuntimeError(('Jordan basis leading diagonal',j))
    if j==R+1 and any(power[k] for k in range(R+1)):
        raise RuntimeError('Nilpotency bound')
    power=mul(power,q)

rows={}
for r in range(R+1):
    h0=max(1,(r+2)//2)
    for h in range(h0+1,H+1):
        predicted=sum(L[j]*D[h-1][r-j] for j in range(r+1))
        if predicted!=D[h][r]: raise RuntimeError(('Transfer mismatch',r,h,predicted,D[h][r]))
    vals=[D[h][r] for h in range(h0,H+1)]; differences=[];work=vals
    for degree in range(r+2):
        differences.append(work[0]);work=[b-a for a,b in zip(work,work[1:])]
    if differences[-1]!=0 or differences[-2]!=4**r: raise RuntimeError(('Polynomial degree',r,differences))
    if any(x.denominator!=1 for x in vals): raise RuntimeError('Integrality')
    rows[r]={'valid_from_h':h0,'binomial_basis_in_h_minus_h0':list(map(int,differences[:-1])),'values_through_h15':list(map(int,vals))}
print(json.dumps({'checks':'PASS: exact transfer recurrence; independent pointed-Burnside kernel; single Jordan block; polynomial degrees0..6; leading coefficient4^r/r!; h through15','L_coefficients_0_to_6':list(map(str,L)),'diagonals':rows},indent=2))
