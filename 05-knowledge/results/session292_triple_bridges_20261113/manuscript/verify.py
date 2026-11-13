#!/usr/bin/env python3
"""Exact controls for the elementary reconstruction lemmas; no RH verification."""
from fractions import Fraction as F
from math import prod, lcm
import json
from pathlib import Path

def check(p, message):
    if not p:
        raise RuntimeError(message)

B=(1,2,3,4,6,8,9,12)
w=[F(1,prod(b-c for c in B if c!=b)) for b in B]
L=lcm(*(x.denominator for x in w))
W=[int(L*x) for x in w]
check(W==[-216,1188,-2464,2079,-924,594,-264,7], 'first witness')
values=[sum(x*b**m for x,b in zip(W,B)) for m in range(10)]
check(values[:7]==[0]*7 and values[7]==L, 'sharp delay')
w2=[x/b**2 for x,b in zip(w,B)]
L2=lcm(*(x.denominator for x in w2)); W2=[int(L2*x) for x in w2]
check(W2==[-186624,256608,-236544,112266,-22176,8019,-2816,42], 'shifted witness')
check(all(sum(x*b**m for x,b in zip(W2,B))==0 for m in range(2,9)), 'physical-length delay')
check(sum(x*b**9 for x,b in zip(W2,B))==L2, 'physical-length separation')
# Build monic characteristic polynomial in ascending powers.
poly=[1]
for b in B:
    new=[0]*(len(poly)+1)
    for j,c in enumerate(poly):
        new[j]-=b*c; new[j+1]+=c
    poly=new
check(all(sum(c*b**j for j,c in enumerate(poly))==0 for b in B), 'annihilator')
for coeff in (W,W2,list(range(8)),[1]*8):
    a=[sum(c*b**m for c,b in zip(coeff,B)) for m in range(25)]
    check(all(sum(poly[j]*a[m+j] for j in range(9))==0 for m in range(17)), 'recurrence')
# Rational Gaussian-complex pair operations certify quartet roots exactly.
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,c): return (a[0]*c,a[1]*c)
controls=0
for R in (1,2,10,100):
    delta=F(1,4)
    for sign in (-1,1):
        for imsign in (-1,1):
            z=(sign*delta,F(imsign*R)); z2=mul(z,z); z4=mul(z2,z2)
            numerator=add(add(((F(R*R)+delta*delta)**2,F(0)),scale(z2,2*(F(R*R)-delta*delta))),z4)
            check(numerator==(0,0),'quartet root'); controls+=1
# Unbounded-delay systems: positive series with arbitrary common prefix.
for length in range(1,30):
    p=[1]* (length+2); q=[1+(m>=length) for m in range(length+2)]
    check(p[:length]==q[:length] and p[length]!=q[length],'unbounded delay')
check(prod(F(p,p-1) for p in (2,3,5))==F(15,4),'235 support')
check(prod(F(p,p-1) for p in (2,3,7))==F(7,2),'237 support')
# Sparse wheel-surviving observations also form a complete eight-mode basis.
samples=[2,4,5,6,8,10,12,14]
small_primes=[p for p in range(2,100) if all(p%d for d in range(2,int(p**0.5)+1))]
check(len(small_primes)==25,'prime list')
check(all(all((99*pow(2,m,p)+23)%p for p in small_primes) for m in samples),'sample survival')
mat=[[F(b**m) for b in B] for m in samples]
determinant=F(1)
for j in range(8):
    pivot=next((r for r in range(j,8) if mat[r][j]),None)
    check(pivot is not None,'sparse observation rank')
    if pivot!=j: mat[j],mat[pivot]=mat[pivot],mat[j]; determinant=-determinant
    v=mat[j][j]; determinant*=v
    for r in range(j+1,8):
        ratio=mat[r][j]/v
        for c in range(j,8):mat[r][c]-=ratio*mat[j][c]
check(determinant!=0 and determinant.denominator==1,'integral determinant')
result={'status':'PASS' ,'bases':B,'weights':W,'moments_0_to_9':values,'shifted_weights':W2,
        'shifted_first_difference':L2,'annihilator_ascending':poly,'quartet_root_controls':controls,'sparse_surviving_lengths':samples,'sparse_determinant':int(determinant),
        'scope':'Exact elementary controls; no physical realization or RH conclusion.'}
print(json.dumps(result,indent=2))
