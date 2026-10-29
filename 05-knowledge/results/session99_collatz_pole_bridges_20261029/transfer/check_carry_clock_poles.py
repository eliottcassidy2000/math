#!/usr/bin/env python3
"""Exact symbolic and native checks; no assert-dependent verification."""
import itertools, json, math
import sympy as s
z=s.symbols('z')
def need(p,msg):
    if not p: raise RuntimeError(msg)
def compositions(n,k):
    if k==1:
        yield (n,); return
    for x in range(n+1):
        for tail in compositions(n-x,k-1): yield (x,)+tail
def shapes(k):
    J=3*k+1
    for indices in itertools.combinations(range(J),k):
        sig=tuple(i in indices for i in range(J)); q=1
        for i,b in enumerate(sig):
            q+=2 if b else -1
            if q<=0 and i<J-1: break
        else:
            if q==0: yield sig
def V(b,x): return x**3/(1-x**3) if b else (x+x*x)/(1-x**3)
def formula(k,r):
    J=3*k+1; terms=[]
    for sig in shapes(k):
        for m in compositions(r,J):
            coefficient=math.factorial(r)*3**sum((J-i-1)*mi for i,mi in enumerate(m))
            coefficient//=math.prod(math.factorial(mi) for mi in m)
            terms.append(coefficient*s.prod(V(b,2**sum(m[j+1:])*z) for j,b in enumerate(sig)))
    return s.cancel(sum(terms))
def direct(k,r,D):
    J=3*k+1; total=0
    for aa in compositions(D-J,J):
        w=tuple(a+1 for a in aa); q=1; C=clock=0
        for i,a in enumerate(w):
            q+=2 if a%3==0 else -1
            if q<=0 and i<J-1: break
            C=3*C+2**clock; clock+=a
        else:
            if q==0: total+=C**r
    return total
rows=[]
for k,r in [(0,0),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]:
    J=3*k+1; ck=math.comb(J,k)//J; H=formula(k,r)
    need(len(list(shapes(k)))==ck,'shape count')
    series=s.series(H,z,0,J+8).removeO().expand()
    for D in range(J,J+8): need(series.coeff(z,D)==direct(k,r,D),f'coefficient {k,r,D}')
    lead=s.limit((1-z**3)**J*H,z,1)
    need(lead==3**(r*(J-1))*ck*2**(2*k+1),'unshifted leading pole')
    if k:
        for e in range(1,r+1):
            lead=s.limit((1-2**(3*e)*z**3)**(J-1)*H,z,s.Rational(1,2**e))
            expected=ck*math.comb(r,e)*3**((J-1)*(r-e))*2**(2*k)*V(False,s.Rational(1,2**e))
            need(lead==expected and lead>0,'shifted leading pole')
    else: need(s.cancel(H-V(False,z))==0,'one leaf')
    rows.append({'k':k,'r':r,'J':J,'coefficients_checked':8})
hostiles=[]
for w in [(3,1,1,2),(3,1,2,1)]:
    C=D=0
    for a in w: C=3*C+2**D; D+=a
    J=len(w); modulus=2**(D+1); source=((2**D-C)*pow(3**J,-1,modulus))%modulus
    n=source
    for a in w:
        v=3*n+1; actual=(v&-v).bit_length()-1
        need(actual==a,'source guard'); n=v//2**a
    need(2**D*n==3**J*source+C,'affine endpoint')
    hostiles.append({'word':w,'carry':C,'source':source,'modulus':modulus,'endpoint':n})
print(json.dumps({'status':'PASS','rows':rows,'hostiles':hostiles},indent=2))
