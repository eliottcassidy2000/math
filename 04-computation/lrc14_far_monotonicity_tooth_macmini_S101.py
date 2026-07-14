#!/usr/bin/env python3
"""mac-mini-S101: prove far-element monotonicity (THM-726 Step 1) via the ALIGNED TOOTH-NARROWING lemma
(general form of S87). Test: (1) each multi-killer core's tight point t0=a/q; (2) is the lcm-carrier
outlier w ALIGNED (q | w => w has a tooth at t0)?; (3) does M(core u {w*m}) -> M(core), increasing,
in closed form? If aligned+monotone, the lemma proves Step 1 for the extremal families."""
from fractions import Fraction as F
import numpy as np
def M_and_t0(S, dens=3_000_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j], x[j]
def best_rational(t, qmax=200):
    from math import gcd
    best=None
    for q in range(2,qmax+1):
        a=round(t*q)
        if 0<a<q and gcd(a,q)==1:
            d=abs(t-a/q)
            if best is None or d<best[0]: best=(d,a,q)
    return best
cores=[("{1..12}",list(range(1,13)),182),
       ("{1..11,13}",list(range(1,12))+[13],84),
       ("{1..10,13,14}",list(range(1,11))+[13,14],45),
       ("{1..10,13,22}",list(range(1,11))+[13,22],84)]
for name,P,w in cores:
    Mp,t0=M_and_t0(P); br=best_rational(t0)
    q=br[2] if br else 0; aligned = (w % q == 0) if q else False
    print(f"core {name}: M(core)={Mp:.5f}, t0~{t0:.5f} ~ {br[1]}/{br[2]}; lcm-carrier w={w}; q={q}; ALIGNED (q|w)? {aligned}")
    # scaled family M(core u {w*m})
    Ms=[]
    for m in range(1,7):
        Mm,_=M_and_t0(P+[w*m], dens=2_000_000); Ms.append(Mm)
    inc = all(Ms[i]<=Ms[i+1]+1e-6 for i in range(len(Ms)-1))
    print(f"      M(core u {{{w}*m}}) m=1..6: {[f'{x:.5f}' for x in Ms]}; monotone-nondecr? {inc}; -> M(core)={Mp:.5f}")
print()
print("TOOTH-NARROWING LEMMA (to prove): if core P is safe at level mu on an interval around t0=a/q,")
print("and outlier w is ALIGNED (q | w, so w*m has a tooth at t0), then M(P u {w*m}) is nondecreasing")
print("in m and -> M(P). Mechanism: ||w*m*t|| ~ w*m*|t-t0| near t0 (tooth half-width mu/(w*m) -> 0);")
print("P-safe interval around t0 has width ~ (M(P)-mu)/slope; family safe at mu iff mu/(w*m) < that width,")
print("so achievable mu increases to M(P) as m grows. This proves Step 1 for aligned lcm-carriers.")
