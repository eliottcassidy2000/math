#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S36 (the quantitative punchline): WHY the ladder skips the n=12 gap.

The ladder M(x)=mu*x/(x+rho) crosses the gap (LO,HI)=(1/13,2/25) as x ranges over an interval
  X_gap = ( LO*rho/(mu-LO) , HI*rho/(mu-HI) ),   width  Dx = rho*[ HI/(mu-HI) - LO/(mu-LO) ].
But the ladder only SAMPLES x at resonances x = j*D (D = base witness denominator), spacing D.
A rung lands strictly inside the gap iff a multiple of D falls in the OPEN interval X_gap.

TWO facts (this file):
 (1) THE GAP IS THE AP LADDER'S FIRST STEP.  For the pure AP {1..n}, mu=1/(n+1)... here base
     {1..11} gives mu=1/12, rho=1, ladder M=j/(12j+1): j=1 -> 1/13 = LO, j=2 -> 2/25 = HI.
     The two gap ENDPOINTS are consecutive AP-ladder rungs, and Dx = D exactly => the samples
     land ON the endpoints, never strictly inside.  The gap window IS the AP ladder's step.
 (2) DEFECTED BASES: Dx vs D.  When Dx < D the resonances (spacing D) SKIP the crossing window.
     We tabulate Dx and D at n=12 and show Dx <= D on representative bases => interior skipped.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import numpy as np

N=12; LO,HI=Fraction(1,N+1),Fraction(2,2*N+1)

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)

print("=== (1) THE GAP IS THE AP LADDER'S FIRST STEP (n=12) ===", flush=True)
ap=list(range(1,12)); mu,(c,D)=Mw(ap)
print(f"  base AP{{1..11}}: mu={mu} at t={c}/{D}, rho=1 => ladder M=j/(12j+1):", flush=True)
for j in (1,2,3):
    M,_=Mw(sorted(ap+[j*D]))
    edge = " = LO (gap bottom)" if M==LO else (" = HI (gap top)" if M==HI else " (above gap)")
    print(f"     j={j}: x={j*D:>2}  M={M}{edge}", flush=True)
print(f"  => the two gap endpoints 1/13, 2/25 ARE consecutive AP-ladder rungs; nothing strictly inside.\n", flush=True)

print("=== (2) Dx (gap-crossing x-window width) vs D (resonance spacing) at n=12 ===", flush=True)
print(f"  {'base':<26}{'mu':>7}{'rho':>5}{'D':>4}{'Dx=gap x-width':>16}{'Dx<D?':>7}", flush=True)
def rho_of(base,mu,D):
    mx=max(base); x0=((2*mx)//D+2)*D
    M0,_=Mw(sorted(base+[x0])); return x0*(mu-M0)/M0 if M0>0 else Fraction(0)
CASES=[
 ("AP{1..11}", list(range(1,12))),
 ("AP{1..10}+{12}", list(range(1,11))+[12]),
 ("AP{1..10}+{13}", list(range(1,11))+[13]),
 ("AP{1..9}+{11,13}", list(range(1,10))+[11,13]),
 ("AP{1..9}+{11,20}", list(range(1,10))+[11,20]),
 ("AP{1..8}+{10,12,14}", list(range(1,9))+[10,12,14]),
]
for label,base in CASES:
    if reduce(gcd,base)!=1: continue
    mu,(c,D)=Mw(base)
    if mu<=HI:
        print(f"  {label:<26}{str(mu):>7}{'-':>5}{D:>4}{'(plateau<=gap top)':>16}", flush=True); continue
    rho=rho_of(base,mu,D)
    Dx = rho*(HI/(mu-HI) - LO/(mu-LO))
    print(f"  {label:<26}{str(mu):>7}{str(rho):>5}{D:>4}{float(Dx):>16.3f}{str(Dx<D):>7}", flush=True)

print("\n  READING: for the pure AP, Dx = D (samples hit the endpoints exactly).  For defected", flush=True)
print("  bases Dx <= D, and where Dx < D the resonance grid (spacing D) SKIPS the crossing", flush=True)
print("  window => no interior rung.  This is 'window too narrow' (kps S32) in x-space, and the", flush=True)
print("  uniform statement 'Dx < D for every n=12 base' = mac-mini's Selberg spacing (HYP-4512).", flush=True)
