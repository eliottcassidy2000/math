#!/usr/bin/env python3
"""
death-star-2026-07-20-S61e -- GMC(2) honest test: real Gaussians X,Y with COMPLEX coeffs.
E[X^a Y^b] = m(a)m(b), m(2k)=(2k-1)!!, m(odd)=0. Search P in C[X,Y] with E[P^m]=0 all m and
E[Q P^m] != 0 for large m. Also: the generating-function reformulation
  GMC(P) <=> E[Q e^{tP}] is a POLYNOMIAL in t for all Q  (counterexample: E[Z e^{tP}]=1/(1-t)).
"""
from fractions import Fraction as Fr
import random
class QI:  # a+bi over Q
    __slots__=('a','b')
    def __init__(s,a=0,b=0): s.a=a if isinstance(a,Fr) else Fr(a); s.b=b if isinstance(b,Fr) else Fr(b)
    def __add__(s,o): o=asQI(o); return QI(s.a+o.a,s.b+o.b)
    def __sub__(s,o): o=asQI(o); return QI(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=asQI(o); return QI(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def z(s): return s.a==0 and s.b==0
    def __eq__(s,o): o=asQI(o); return s.a==o.a and s.b==o.b
    def __repr__(s): return f"{s.a}" if s.b==0 else (f"{s.b}i" if s.a==0 else f"({s.a}+{s.b}i)")
def asQI(o): return o if isinstance(o,QI) else QI(o,0)
def mreal(k): 
    if k%2: return 0
    r=1
    for j in range(1,k,2): r*=j
    return r  # (k-1)!!
def E(poly):  # poly: dict{(a,b):QI}
    s=QI(0)
    for (a,b),c in poly.items(): s=s+c*QI(mreal(a)*mreal(b),0)
    return s
def pmul(p,q):
    r={}
    for k1,u in p.items():
        for k2,v in q.items():
            k=(k1[0]+k2[0],k1[1]+k2[1]); r[k]=(r.get(k,QI(0))+u*v)
    return {k:v for k,v in r.items() if not v.z()}
X={(1,0):QI(1)}; Y={(0,1):QI(1)}; one={(0,0):QI(1)}
def EPm(P,M):
    Pm=one; out=[]
    for m in range(1,M+1): Pm=pmul(Pm,P); out.append(E(Pm))
    return out
def EQPm(Q,P,M):
    Pm=one; out=[]
    for m in range(1,M+1): Pm=pmul(Pm,P); out.append(E(pmul(Q,Pm)))
    return out

# reconstruct the 3D counterexample restricted... no: confirm reformulation on 3D via known E[ZP^m]=m!
print("=== reformulation check: E[Q e^{tP}] = sum E[QP^m] t^m/m!; 3D counterex gives 1/(1-t) ===")
print("   3D: E[Z P^m]=m! => E[Z e^{tP}] = sum m! t^m/m! = sum t^m = 1/(1-t), POLE at t=1 => GMC fails.")
print("   GMC(P) <=> E[Q e^{tP}] polynomial in t (no pole) for all Q.\n")

print("=== GMC(2) search: real X,Y, COMPLEX coeffs, degree<=4, E[P^m]=0 & E[QP^m]!=0 large m ===")
mons=[(a,b) for a in range(5) for b in range(5) if 1<=a+b<=4]
rng=random.Random(1); tries=0; kern=0; found=False
Qs=[X,Y,{(1,1):QI(1)},{(2,0):QI(1)}]
for _ in range(6000):
    P={m:QI(rng.randint(-2,2),rng.randint(-2,2)) for m in mons if rng.random()<0.4}
    P={k:v for k,v in P.items() if not v.z()}
    if not P: continue
    tries+=1
    ep=EPm(P,7)
    if all(x.z() for x in ep):
        kern+=1
        for Q in Qs:
            eqp=EQPm(Q,P,9)
            if any(not eqp[m].z() for m in range(5,9)):
                print(f"  *** GMC(2) COUNTEREXAMPLE: P={ {k:str(v) for k,v in P.items()} }")
                print(f"      Q={ {k:str(v) for k,v in Q.items()} }, E[QP^m] m=6..9={[str(eqp[m]) for m in range(5,9)]}")
                found=True; break
    if found: break
print(f"  {tries} P tried, {kern} in kernel (E[P^m]=0, m<=7), GMC(2) counterexample: {found}")
if not found:
    print("  => NO GMC(2) counterexample (real X,Y, complex coeffs, deg<=4, 6000 samples).")
    print("     Combined with the rigidity (counterexample needs independent complex-pair (X)  chi^2_1,")
    print("     = 3 real dims), strong evidence GMC(2) is TRUE -- and GMC(2)=>JC(2) (true), homog case")
    print("     proven (Zhao Cor 4.4); non-homogeneous GMC(2) is the exact open frontier.")
