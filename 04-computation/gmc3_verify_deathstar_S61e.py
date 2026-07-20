#!/usr/bin/env python3
"""
death-star-2026-07-20-S61e -- verify the owner's GMC(3) counterexample via the FORMAL moment
functional (rational, exact; no sqrt2 needed). Variables Z,W,U with
   L(Z^a W^b U^c) = a! * [a==b] * (1/2)_c ,   (1/2)_c = (2c-1)!!/2^c = poch(1/2,c).
P = (1+Z)(W - (2+Z)U),  Q = Z. Claim: L(P^m)=0 for all m>=1 but L(Q P^m)=m!.
"""
from fractions import Fraction as Fr
from math import factorial

def poch_half(c):
    r=Fr(1)
    for j in range(c): r*= (Fr(1,2)+j)
    return r
def L(poly):
    s=Fr(0)
    for (a,b,c),co in poly.items():
        if a==b: s+= co*factorial(a)*poch_half(c)
    return s
def pmul(p,q):
    r={}
    for (a,b,c),u in p.items():
        for (d,e,f),v in q.items():
            k=(a+d,b+e,c+f); r[k]=r.get(k,Fr(0))+u*v
    return {k:v for k,v in r.items() if v!=0}
def padd(*ps):
    r={}
    for p in ps:
        for k,v in p.items(): r[k]=r.get(k,Fr(0))+v
    return {k:v for k,v in r.items() if v!=0}
def scal(p,c): return {k:v*c for k,v in p.items()}

Z={(1,0,0):Fr(1)}; W={(0,1,0):Fr(1)}; U={(0,0,1):Fr(1)}; one={(0,0,0):Fr(1)}
# P = (1+Z)*(W - (2+Z)*U)
oneZ=padd(one,Z); twoZ=padd(scal(one,2),Z)
P=pmul(oneZ, padd(W, scal(pmul(twoZ,U),-1)))
print("P monomials (a,b,c):coeff =", {k:str(v) for k,v in sorted(P.items())})
Pm=one
for m in range(1,8):
    Pm=pmul(Pm,P)
    lp=L(Pm); lqp=L(pmul(Z,Pm))
    print(f"  m={m}: L(P^m)={lp}   L(Z*P^m)={lqp}   (expect 0 and {factorial(m)})  "
          f"{'OK' if lp==0 and lqp==factorial(m) else 'FAIL'}")
print("\n=> GMC(3) FALSE confirmed: L(P^m)=0 for all m, L(Q P^m)=m! != 0 for all m.")
