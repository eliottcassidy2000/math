# New unconditional test: span-6 charge family {+3,+1,-1,-3}, P = a Z^3 + b Z + c Zbar + d Zbar^3.
# E[Z^A Zbar^B] = A! [A=B].  E[P^m] = sum over configs (i3+,i1,i1-,i3-) with total charge 0.
# Nullcone: E[P^m]=0 all m.  GMC(2)/NC2 here: forces (a,b)=0 OR (c,d)=0 (one-sided).
from fractions import Fraction as Fr
from itertools import product as iproduct
from math import factorial
# term charges and (A-contribution): aZ^3:(+3,Z-power3), bZ:(+1,1), cZbar:(-1, Zbar1), dZbar^3:(-3,Zbar3)
# monomial in coeffs (a,b,c,d) -> Fraction coeff of E[P^m]
def moments(M):
    res=[]
    # each term: (charge, zpow, zbarpow, varindex)
    terms=[(+3,3,0,0),(+1,1,0,1),(-1,0,1,2),(-3,0,3,3)]
    for m in range(1,M+1):
        poly={}  # (i0,i1,i2,i3) exponent tuple -> coeff
        for combo in iproduct(range(m+1),repeat=4):
            if sum(combo)!=m: continue
            A=sum(combo[k]*terms[k][1] for k in range(4))   # total Z power
            B=sum(combo[k]*terms[k][2] for k in range(4))   # total Zbar power
            if A!=B: continue
            mult=factorial(m)
            for k in range(4): mult//=factorial(combo[k])
            coeff=Fr(mult*factorial(A))
            key=combo
            poly[key]=poly.get(key,Fr(0))+coeff
        res.append(poly)
    return res
def show(poly):
    names=['a','b','c','d']
    parts=[]
    for key,co in sorted(poly.items()):
        mon="".join(f"{names[k]}^{key[k]}" if key[k]>1 else (names[k] if key[k]==1 else "") for k in range(4))
        parts.append(f"{co}{mon}")
    return " + ".join(parts) if parts else "0"
Ms=moments(8)
for m in (2,4,6,8):
    print(f"E[P^{m}] = {show(Ms[m-1])}")
# numerical emptiness probe: search for a TWO-SIDED (a or b !=0) AND (c or d !=0) nullcone solution
import random
def val(poly,v):
    s=0.0
    for key,co in poly.items():
        t=float(co)
        for k in range(4): t*=v[k]**key[k]
        s+=t
    return s
rng=random.Random(0); found=False
for _ in range(300000):
    v=[rng.uniform(-2,2) for _ in range(4)]
    if abs(v[0])+abs(v[1])<1e-3 or abs(v[2])+abs(v[3])<1e-3: continue
    if all(abs(val(Ms[m-1],v))<1e-6 for m in (2,4,6,8)):
        found=True; print("TWO-SIDED nullcone candidate:",v); break
print("two-sided nullcone solution found in 300k search:",found,"(False => consistent with GMC(2))")
