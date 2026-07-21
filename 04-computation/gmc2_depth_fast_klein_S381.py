#!/usr/bin/env python3
"""klein-S381 -- FAST: detection depth on {-1,0,1} at d=0 and d=1 only, m<=8, exact."""
import sympy as sp
from math import factorial
from sympy import groebner
Z,Zb=sp.symbols('Z Zb'); ZZb=Z*Zb
def Ew(poly):
    poly=sp.expand(poly)
    if poly==0: return sp.Integer(0)
    t=sp.Integer(0)
    for (a,b),c in sp.Poly(poly,Z,Zb).as_dict().items():
        if a==b: t+=c*factorial(a)
    return t
def buildP(d):
    A=[sp.symbols(f'A{k}') for k in range(d+1)]; B=[sp.symbols(f'B{k}') for k in range(d+1)]; C=[sp.symbols(f'C{k}') for k in range(d+1)]
    P=sum(A[k]*ZZb**k*Z for k in range(d+1))+sum(B[k]*ZZb**k for k in range(d+1))+sum(C[k]*ZZb**k*Zb for k in range(d+1))
    return P,A,B,C
def depth(d,mmax=8):
    P,A,B,C=buildP(d); allc=A+B+C
    moms=[sp.nsimplify(sp.expand(Ew(sp.expand(P**m)))) for m in range(1,mmax+1)]
    targets=[a*c for a in A for c in C]+list(B)
    for D in range(1,mmax+1):
        gens=[g for g in moms[:D] if g!=0]
        if not gens: continue
        G=groebner(gens,*allc,order='grevlex')
        if all(any(sp.simplify(G.reduce(t**k)[1])==0 for k in range(1,7)) for t in targets):
            return D,moms
    return None,moms
for d in (0,1):
    D,moms=depth(d)
    print(f"d={d}: #coeffs={3*(d+1)}  detection depth = {D}")
    for i in range(min(4,len(moms))): print(f"    E[P^{i+1}] = {moms[i]}")
    print()
