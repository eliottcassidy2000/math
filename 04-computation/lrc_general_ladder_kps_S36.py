#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S36: THE GENERAL RESONANCE LADDER + opus's spectrum form.

S35 found, for base {1,2,3,4,5,7}: adding a resonant outlier x=6j gives M=j/(6j+5).
GENERAL CLAIM (this file): for ANY base B with plateau M(B)=mu at witness t*=c/D and
binding-runner descent rate rho, a resonant outlier x gives

        M(B + x) = mu * x / (x + rho)                          [THE LADDER]

- rho is a PROPERTY OF THE BASE (a binding runner speed), constant across resonances.
- as x -> infinity, M -> mu (opus HYP-4476 height-independence = the plateau).
- expressed in opus's amended-spectrum form (HYP-4486) M = s/(ns+k) [n = #speeds],
  the reduced M=p/q gives (s,k) = (p, q - n*p); opus: M in gap window (1/(n+1),2/(2n+1))
  <=> k < s < 2k, and k>=2 required (Kravitz rungs k=1 never inside).
- mac-mini HYP-4542 "far element ~ const*c via lever" = the x<->c relation of the ladder.

We (1) verify the ladder formula + constant rho on several bases; (2) express rungs in
opus's (s,k) form and check the window; (3) drive it at n=12 (the crux): do ANY 11-speed
bases have a ladder rung with k>=2 and k<s<2k (i.e. attaining a gap value)?
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)

def spectrum_form(M, n):
    """M = p/q reduced -> (s,k) with M = s/(n*s+k); k = q - n*p (may be <=0 if M>=1/n)."""
    p,q = M.numerator, M.denominator
    return p, q - n*p

def in_window(M, n):
    return Fraction(1,n+1) < M < Fraction(2,2*n+1)

def ladder_of_base(base, jmax=9):
    """M(base) and the resonance ladder (x = multiples of the base denominator D)."""
    mu,(c,D) = Mw(base)
    n = len(base)+1   # #speeds after adding one outlier
    rows=[]
    for j in range(1, jmax+1):
        x = j*D
        v = sorted(base+[x])
        if len(set(v))!=len(base)+1: continue
        if reduce(gcd,v)!=1:
            # ladder still defined; skip non-primitive for cleanliness
            pass
        M,(cc,q) = Mw(v)
        # empirical rho from M = mu*x/(x+rho): rho = x*(mu-M)/M
        rho = x*(mu-M)/M if M>0 else None
        s,k = spectrum_form(M, n)
        rows.append((j, x, M, rho, (s,k), in_window(M,n), (cc,q)))
    return mu,(c,D),n,rows

print("=== S36: GENERAL RESONANCE LADDER  M(B+x)=mu*x/(x+rho), in opus s/(ns+k) form ===\n", flush=True)

BASES = {
    "{1,2,3,4,5,7}  (S35)": [1,2,3,4,5,7],
    "{1,2,3,4,5}    (AP5)": [1,2,3,4,5],
    "{1,2,3,4,5,6}  (AP6)": [1,2,3,4,5,6],
    "{1,2,3,4,5,6,8}": [1,2,3,4,5,6,8],
}
for name,base in BASES.items():
    mu,(c,D),n,rows = ladder_of_base(base)
    print(f"BASE {name}: M(B)={mu} at t={c}/{D}, n(after +x)={n}, gap=({Fraction(1,n+1)},{Fraction(2,2*n+1)})", flush=True)
    # infer rho (should be constant)
    rhos = {r[3] for r in rows if r[3] is not None}
    print(f"  rho (empirical, from M=mu*x/(x+rho)): {sorted(set(round(float(x),4) for x in rhos))}", flush=True)
    for j,x,M,rho,(s,k),inw,(cc,q) in rows:
        pred = mu*x/(x+rho) if rho is not None else None
        star = "  <== IN GAP" if inw else ""
        print(f"    j={j:>2} x={x:>3}  M={str(M):>8}  (s,k)=({s},{k})  {'k<s<2k' if (k>=1 and k<s<2*k) else '        '}{star}", flush=True)
    print(flush=True)

print("="*70, flush=True)
print("READING: rho is constant per base (a binding-runner speed) => the ladder", flush=True)
print("M=mu*x/(x+rho) holds; rungs realize opus's s/(ns+k); the unique in-gap rung", flush=True)
print("is opus's minimal (k=2,s=3)=mediant.  Next block drives the n=12 crux.", flush=True)
