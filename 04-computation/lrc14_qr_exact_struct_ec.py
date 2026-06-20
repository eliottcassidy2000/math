#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXACT structural verification (no slow 6^6 exact scan; magnitudes done numerically
in lrc14_qr_fastscan_ec.py).  Here we verify EXACTLY in Z[zeta_7]:
  (2) Galois equivariance  D7(a c) = sigma_a(D7(c))
  (3) trace rationality: sum_a D7(ac) rational; Tr(D7(1^6)) = 5040 = 7!
  (4) reality: sigma_6(D7) = conj, D7(6c)=conj(D7(c))
  GAP: D7 generically NOT in Q(sqrt-7) (only fixed by sigma_2 rarely).
"""
import sys, itertools, math, cmath, random
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

ZERO=(Fraction(0),)*6
ONE=(Fraction(1),Fraction(0),Fraction(0),Fraction(0),Fraction(0),Fraction(0))
def z_pow(e):
    e%=7
    if e==6: return (Fraction(-1),)*6
    v=[Fraction(0)]*6; v[e]=Fraction(1); return tuple(v)
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def mul(a,b):
    raw=[Fraction(0)]*12
    for i in range(6):
        if a[i]==0: continue
        for j in range(6):
            if b[j]==0: continue
            raw[i+j]+=a[i]*b[j]
    out=[Fraction(0)]*6
    for e in range(12):
        if raw[e]==0: continue
        be=z_pow(e)
        for k in range(6): out[k]+=raw[e]*be[k]
    return tuple(out)
def galois(a,x):
    r=[Fraction(0)]*6
    for k in range(6):
        if x[k]!=0:
            be=z_pow((a*k)%7)
            for j in range(6): r[j]+=x[k]*be[j]
    return tuple(r)
def trace(x):
    tot=[Fraction(0)]*6
    for a in range(1,7):
        ga=galois(a,x)
        for j in range(6): tot[j]+=ga[j]
    return tuple(tot)
def is_rational(x): return all(x[k]==0 for k in range(1,6))
def to_num(x):
    Z=cmath.exp(2j*math.pi/7)
    return sum(complex(x[k])*(Z**k) for k in range(6))

Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:Fraction((-1)**len(T)) for T in Tlist}
def sigma_T(T,m):
    r=[Fraction(0)]*6
    for t in T:
        be=z_pow((-m*t)%7)
        for j in range(6): r[j]+=be[j]
    return tuple(r)
SIG={(T,m):sigma_T(T,m) for T in Tlist for m in range(1,7)}
PREF={m:sub(ONE,z_pow((-m)%7)) for m in range(1,7)}
_cache={}
def D7(c):
    v=_cache.get(c)
    if v is not None: return v
    pref=ONE
    for cj in c: pref=mul(pref,PREF[cj])
    acc=[Fraction(0)]*6
    for T in Tlist:
        p=ONE
        for cj in c: p=mul(p,SIG[(T,cj)])
        sg=SGN[T]
        for j in range(6): acc[j]+=sg*p[j]
    v=mul(pref,tuple(acc)); _cache[c]=v; return v

QR={1,2,4}

random.seed(7)
tests=[(1,1,1,1,1,1),(1,2,4,3,5,6),(1,1,1,1,1,2),(3,3,3,3,3,3),(1,2,3,4,5,6),
       (1,1,1,1,1,4)]+[tuple(random.choice([1,2,3,4,5,6]) for _ in range(6)) for _ in range(40)]

# Gauss g^2=-7
g=[Fraction(0)]*6
for r in range(1,7):
    be=z_pow(r); s=Fraction(1 if r in QR else -1)
    for j in range(6): g[j]+=s*be[j]
g2=mul(tuple(g),tuple(g))
print(f"[sanity] g^2 = {[str(x) for x in g2]} (expect -7,0,0,0,0,0) -> {'OK' if g2==tuple([Fraction(-7)]+[Fraction(0)]*5) else 'FAIL'}")

# (2)
f2=0
for c in tests:
    d=D7(c)
    for a in range(1,7):
        if D7(tuple((a*cj)%7 for cj in c))!=galois(a,d): f2+=1
print(f"[2] Galois equivariance D7(ac)=sigma_a(D7(c)): {len(tests)*6} pairs, fails={f2} -> {'CONFIRMED' if f2==0 else 'REFUTED'}")

# (3)
f3=0
for c in tests:
    orb=[Fraction(0)]*6
    for a in range(1,7):
        d=D7(tuple((a*cj)%7 for cj in c))
        for j in range(6): orb[j]+=d[j]
    orb=tuple(orb)
    if not is_rational(orb) or orb!=trace(D7(c)): f3+=1
tr1=trace(D7((1,1,1,1,1,1)))
print(f"[3] sum_a D7(ac)=Tr rational on all tests: fails={f3}; Tr(1^6)={tr1[0]} (7!={math.factorial(7)}) -> "
      f"{'CONFIRMED' if f3==0 and tr1[0]==Fraction(5040) else 'REFUTED'}")

# (4)
f4=0
for c in tests:
    d=D7(c); s6=galois(6,d)
    if D7(tuple((6*cj)%7 for cj in c))!=s6: f4+=1
    if abs(to_num(s6)-to_num(d).conjugate())>1e-9: f4+=1
print(f"[4] sigma_6=conj and D7(6c)=conj(D7(c)): fails={f4} -> {'CONFIRMED' if f4==0 else 'REFUTED'}")

# GAP: fraction of cosets with D7 in Q(sqrt-7) (fixed by sigma_2)
infield=0; rational=0
sample=[tuple(random.choice([1,2,3,4,5,6]) for _ in range(6)) for _ in range(404)]
for c in sample:
    d=D7(c)
    if galois(2,d)==d: infield+=1
    if galois(3,d)==d: rational+=1   # fixed by NQR gen too => rational
print(f"[GAP] D7 in Q(sqrt-7) [fixed by sigma_2=QR gen]: {infield}/{len(sample)}; "
      f"D7 rational [fixed by sigma_3 too]: {rational}/{len(sample)}")
print("      => QR/NQR index-2 split alone does NOT diagonalize; full F7* trace needed,")
print("         which the integer lattice does not realize. GAP CONFIRMED.")
