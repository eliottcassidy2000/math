#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast numeric corroboration of the full 6^6 coset Tr scan (rounds to nearest int).
Tr(D7(c)) is a rational integer; numeric D7 sum over dilation orbit -> round.
Confirms: sum Tr = 0, max|Tr| = 9240, per-#QR-class sum = 0.
"""
import itertools, math, cmath
from collections import defaultdict

Z=cmath.exp(2j*math.pi/7)
QR={1,2,4}
Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:(-1)**len(T) for T in Tlist}
SIG={(T,m):sum(Z**((-m*t)%7) for t in T) for T in Tlist for m in range(1,7)}
PREF={m:(1-Z**((-m)%7)) for m in range(1,7)}
_c={}
def D7(c):
    v=_c.get(c)
    if v is not None: return v
    pref=1+0j
    for cj in c: pref*=PREF[cj]
    acc=0j
    for T in Tlist:
        p=1+0j
        for cj in c: p*=SIG[(T,cj)]
        acc+=SGN[T]*p
    v=pref*acc; _c[c]=v; return v

def Tr(c):
    s=0j
    for a in range(1,7):
        s+=D7(tuple((a*cj)%7 for cj in c))
    return s

sumTr=0; maxabs=0; byqr=defaultdict(int); maxc=None
imag_max=0.0
for c in itertools.product(range(1,7),repeat=6):
    t=Tr(c)
    imag_max=max(imag_max, abs(t.imag))
    ti=round(t.real)
    assert abs(t.real-ti)<1e-6, f"non-integer Tr at {c}: {t}"
    sumTr+=ti
    if abs(ti)>maxabs: maxabs=abs(ti); maxc=c
    byqr[sum(1 for x in c if x in QR)]+=ti
print(f"max imag(Tr) over all cosets = {imag_max:.2e} (should ~0, confirms rational)")
print(f"sum Tr over all 6^6 cosets = {sumTr}  (claim 0)")
print(f"max|Tr| = {maxabs} at c={maxc}  (claim 9240 = {2**3*3*5*7*11})")
print("per #QR-class sum Tr:")
for k in sorted(byqr): print(f"  #QR={k}: {byqr[k]}")
print(f"all per-class zero: {all(v==0 for v in byqr.values())}")
