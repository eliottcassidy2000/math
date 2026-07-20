#!/usr/bin/env python3
"""
death-star-2026-07-20-S61e -- (A) rigidity of the GMC(3) mechanism; (B) search for a GMC(2)
counterexample. Formal functional L(Z^a W^b U^c) = a![a==b]*mu_c on Q[Z,W,U].
"""
from fractions import Fraction as Fr
from math import factorial

def poch_half(c):
    r=Fr(1)
    for j in range(c): r*=(Fr(1,2)+j)
    return r
def make_L(mu):   # mu: c -> Fraction (moments of U)
    def L(poly):
        s=Fr(0)
        for (a,b,c),co in poly.items():
            if a==b: s+= co*factorial(a)*mu(c)
        return s
    return L
def pmul(p,q):
    r={}
    for k1,u in p.items():
        for k2,v in q.items():
            k=tuple(x+y for x,y in zip(k1,k2)); r[k]=r.get(k,Fr(0))+u*v
    return {k:v for k,v in r.items() if v!=0}
def padd(*ps):
    r={}
    for p in ps:
        for k,v in p.items(): r[k]=r.get(k,Fr(0))+v
    return {k:v for k,v in r.items() if v!=0}
def scal(p,c): return {k:v*c for k,v in p.items() if v*c!=0}
Z={(1,0,0):Fr(1)}; W={(0,1,0):Fr(1)}; U={(0,0,1):Fr(1)}; one={(0,0,0):Fr(1)}

def test_family(c_const, mu, label, M=6):
    # P = (1+Z)(W - (c+Z)U)
    P=pmul(padd(one,Z), padd(W, scal(pmul(padd(scal(one,c_const),Z),U),-1)))
    L=make_L(mu)
    Pm=one; res=[]
    for m in range(1,M+1):
        Pm=pmul(Pm,P); res.append(L(Pm))
    allzero=all(r==0 for r in res)
    print(f"  {label}: E[P^m] m=1..{M} = {[str(r) for r in res]}  {'ALL ZERO' if allzero else 'nonzero'}")
    return allzero

print("=== (A) rigidity: only c=2 AND mu_c=(1/2)_c gives the counterexample ===")
test_family(2, poch_half,          "c=2, mu=(1/2)_c  [the counterexample]")
test_family(3, poch_half,          "c=3, mu=(1/2)_c  [wrong constant]")
test_family(2, lambda c: Fr(factorial(c)), "c=2, mu=c!       [U~exponential/ZW, wrong EGF]")
test_family(2, lambda c: Fr(1),    "c=2, mu=1        [U~point mass]")
test_family(1, poch_half,          "c=1, mu=(1/2)_c")
test_family(4, poch_half,          "c=4, mu=(1/2)_c")

print("\n=== (B) GMC(2): search for P in C[Z,W] (no U) with E[P^m]=0 all m, E[Q P^m]!=0 large m ===")
L2=make_L(lambda c: Fr(1) if c==0 else Fr(0))  # ignore U (c must be 0); effectively C[Z,W]
def EP(P,M):   # returns [E[P^m] for m=1..M]
    Pm=one; out=[]
    for m in range(1,M+1): Pm=pmul(Pm,P); out.append(L2(Pm))
    return out
def EQP(Q,P,M):
    Pm=one; out=[]
    for m in range(1,M+1):
        Pm=pmul(Pm,P); out.append(L2(pmul(Q,Pm)))
    return out
import itertools
# ansatz families to test (Z,W polynomials), c_00=0 forced
ansatze = {
  "(1+Z)(W-(2+Z)W) = -(1+Z)^2 W": pmul(padd(one,Z), padd(W, scal(pmul(padd(scal(one,2),Z),W),-1))),
  "(1+Z)(W-(2+Z)Z)":             pmul(padd(one,Z), padd(W, scal(pmul(padd(scal(one,2),Z),Z),-1))),
  "W - Z*W (=(1-Z)W)":           padd(W, scal(pmul(Z,W),-1)),
  "W - W^2":                      padd(W, scal(pmul(W,W),-1)),
  "Z*W - Z":                      padd(pmul(Z,W), scal(Z,-1)),
  "W + Z*W - Z^2*W":              padd(W, pmul(Z,W), scal(pmul(pmul(Z,Z),W),-1)),
}
found=False
for name,P in ansatze.items():
    ep=EP(P,8)
    if all(r==0 for r in ep):
        eqpW=EQP(W,P,8); eqpZ=EQP(Z,P,8)
        big = any(eqpW[m]!=0 for m in range(4,8)) or any(eqpZ[m]!=0 for m in range(4,8))
        print(f"  E[P^m]=0 for '{name}': True; E[W P^m] m=1..8={[str(x) for x in eqpW]}; large-m nonzero: {big}")
        if big: found=True
    else:
        pass
# random search over small integer-coeff P in Z,W up to degree 3
import random
rng=random.Random(0)
mons=[(a,b,0) for a in range(4) for b in range(4) if 1<=a+b<=3]
tries=0; hits=0
for _ in range(4000):
    P={m:Fr(rng.randint(-2,2)) for m in mons if rng.random()<0.5}
    P={k:v for k,v in P.items() if v!=0}
    if not P: continue
    tries+=1
    ep=EP(P,7)
    if all(r==0 for r in ep):
        hits+=1
        eqpW=EQP(W,P,9)
        if any(eqpW[m]!=0 for m in range(5,9)):
            print(f"  *** GMC(2) CANDIDATE: P={ {k:str(v) for k,v in P.items()} }, E[W P^m] large={ [str(eqpW[m]) for m in range(5,9)] }")
            found=True
print(f"  random search: {tries} tried, {hits} had E[P^m]=0 for m<=7, GMC(2)-counterexample found: {found}")
if not found:
    print("  => NO GMC(2) counterexample found (degree<=3, integer coeffs). Evidence GMC(2) is TRUE.")
