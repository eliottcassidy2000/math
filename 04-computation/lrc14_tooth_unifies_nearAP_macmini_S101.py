#!/usr/bin/env python3
"""mac-mini-S101c: does THM-751's tooth-narrowing (at the core's tight point, ANY q -- not just q<=13)
close the near-AP-with-far families the SHADOW tile (k<=13) MISSED (S99)? The near-AP core's tight
point is q=14 (the AP 1/14), and lcm-carriers like 182=14*13 are ALIGNED there (14|182). If THM-751's
bound mu0*wm/(wm+pmax) closes them, it UNIFIES the shadow + near-AP tiles under one mechanism."""
from fractions import Fraction as F
from math import gcd
import numpy as np
def M_t0(S,dens=4_000_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def best_rat(t,qmax=60):
    b=None
    for q in range(2,qmax+1):
        a=round(t*q)
        if 0<a<q and gcd(a,q)==1:
            d=abs(t-a/q)
            if b is None or d<b[0]: b=(d,a,q)
    return b
# the S99 escapee family (near-AP {1..13}\{6} + 182m) that the k<=13 shadow MISSED
tests=[("{1..13}\\{6}+182m", [v for v in range(1,14) if v!=6], 182),
       ("{1..13}\\{6,7}+ 14m? core+169", [v for v in range(1,14) if v not in (6,7)], 169),
       ("{1..12}+182m (control, shadow-closed)", list(range(1,13)), 182)]
one14=F(1,14)
for name,core,w in tests:
    mu0,t0=M_t0(core); br=best_rat(t0); q=br[2] if br else 0
    aligned = (w % q==0) if q else False
    pmax=max(core)
    print(f"\n{name}: core M={mu0:.5f}, t0~{t0:.5f}~{br[1]}/{br[2]} (q={q}); w={w}; ALIGNED(q|w)? {aligned}; pmax={pmax}")
    print(f"   THM-751 bound mu0*wm/(wm+pmax) vs actual M(core u {{{w}m}}):")
    murat=F(round(mu0*100000),100000)  # rational approx of mu0 for display; use exact core M below if q clean
    for m in range(1,5):
        wm=w*m; Mm,_=M_t0(core+[wm],dens=2_000_000)
        bound=mu0*wm/(wm+pmax)  # float bound
        print(f"     m={m}: bound~{bound:.5f}  actual M={Mm:.5f}  bound<=M? {bound<=Mm+1e-5}  M>=1/14? {Mm>=float(one14)-1e-6}")
print()
print("=> if aligned (q|w) and the bound holds for the near-AP q=14 cores, THM-751 (tooth-narrowing at")
print("ANY q) closes the near-AP-with-far families the k<=13 SHADOW MISSED (S99) -- UNIFYING the shadow")
print("(q<=13) and near-AP (q=14) tiles as ONE mechanism: an aligned tooth at the core's tight point.")
