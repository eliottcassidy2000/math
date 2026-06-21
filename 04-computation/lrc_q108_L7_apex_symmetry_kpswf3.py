#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Symmetries of D_{p,q}: transpose p<->q, and the true reflection that DOES preserve |D|."""
from fractions import Fraction as Fr
from math import gcd
P=7
def sector(yf): return int(P*yf)
def mu_full(p,q):
    bp={Fr(0),Fr(1)}
    for f in (p,q):
        for t in range(0,P*f): bp.add(Fr(t,P*f))
    vs=sorted(bp); cell={}
    for a,b in zip(vs,vs[1:]):
        mid=(a+b)/2; k=(sector((q*mid)%1),sector((p*mid)%1))
        cell[k]=cell.get(k,Fr(0))+(b-a)
    return cell
def D_pq(p,q):
    cell=mu_full(p,q); inv=Fr(1,49)
    return sum(abs(cell.get((i,j),Fr(0))-inv) for i in range(P) for j in range(P))

# (A) transpose symmetry D_{p,q} = D_{q,p}
tfail=0
for q in range(1,30):
    for p in range(1,30):
        if gcd(p,q)!=1: continue
        if D_pq(p,q)!=D_pq(q,p): tfail+=1
print("Transpose symmetry D_{p,q}=D_{q,p}:", "HOLDS" if tfail==0 else f"{tfail} FAIL")

# (B) does D depend only on (p mod 7, q mod 7)? test: same residues -> same D?
from collections import defaultdict
byres=defaultdict(set)
for q in range(1,40):
    for p in range(1,40):
        if gcd(p,q)!=1: continue
        byres[(p%7,q%7)].add(D_pq(p,q))
nconst=sum(1 for v in byres.values() if len(v)==1)
print(f"D depends only on (p mod7, q mod7)? {nconst}/{len(byres)} residue classes have a single D value")
print("  => NO; the integer size of p,q matters (D ~ 1/q). residue sets the cell PATTERN's support.")

# (C) the reflection that DOES fix |D|: which residue map preserves the FULL D-spectrum?
# test g: (p,q)->(g*p mod7-lift, q) impossible directly; instead test whether the SET of
# D-values for slope s equals that for slope s' under s'=QR-multiple m*s.
def slope_spectrum(target_s):
    out=set()
    for q in range(1,30):
        if q%7==0: continue
        si=pow(q%7,-1,7)
        for p in range(1,60):
            if gcd(p,q)!=1 or p%7==0: continue
            if (p*si)%7==target_s: out.add(D_pq(p,q))
    return out
spec={s:slope_spectrum(s) for s in range(1,7)}
print("\nWhich m in (Z/7)* satisfy spec[s]=spec[m*s mod7] for ALL s? (a multiplicative symmetry of |D|)")
for m in range(1,7):
    ok=all(spec[s]==spec[(m*s)%7] for s in range(1,7))
    print(f"  m={m} ({'QR' if m in {1,2,4} else 'NQR'}): {'SYMMETRY' if ok else 'no'}")
