#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S35: THE PLATEAU + RESONANCE LADDER (corrects S35's first guess).

Slice {1,2,3,4,5,7,x} (base B={1,2,3,4,5,7}, M(B)=1/6 at t=1/6):
  * PLATEAU: generic outlier x (x !=0 mod 6) => M=1/6 (opus HYP-4476 height-independence, concrete).
  * RESONANCE: x = 6j has dist 0 at t=1/6 => KILLS the base witness => M drops onto a ladder.
  * CLOSED FORM (j>=3, x=6j additively isolated): M = j/(6j+5), witness t=(j+1)/(6j+5).
      LOWER BOUND is closed-form provable (residue table below); equality is computational.
      Denominators 6j+5 = 23,29,35,41,... form an AP (step 6) rising to the plateau 1/6.
  * UNIQUE GAP RUNG: j/(6j+5) in (1/8, 2/15) <=> 2j>5 and 3j<10 <=> j=3 => x=18, M=3/23 (mediant, at wall 3k+2=23).

The gap (window) catches EXACTLY ONE rung of the resonance ladder -- structure x width made explicit.
"""
from fractions import Fraction
import numpy as np

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)

base=[1,2,3,4,5,7]
gap=(Fraction(1,8),Fraction(2,15))
print(f"base {base}: M={Mw(base)[0]} at t=1/6  (generic denominator 6)\n", flush=True)

print("PLATEAU test (x NOT = 0 mod 6): M should be 1/6", flush=True)
plateau_ok=True
for x in range(13,50):
    if x%6==0: continue
    M,_=Mw(sorted(base+[x]))
    if M!=Fraction(1,6): plateau_ok=False; print(f"  x={x}: M={M} (off plateau)", flush=True)
print(f"  plateau holds for all non-resonant x in 13..49: {plateau_ok}\n", flush=True)

print("RESONANCE LADDER (x=6j): closed form M=j/(6j+5) at t=(j+1)/(6j+5)", flush=True)
print(f"  {'j':>2} {'x=6j':>5} {'M(computed)':>12} {'j/(6j+5)':>10} {'match':>6} {'t':>9} {'in gap?':>8}", flush=True)
for j in range(3,12):
    x=6*j; M,(c,q)=Mw(sorted(base+[x]))
    cf=Fraction(j,6*j+5); ingap = gap[0]<M<gap[1]
    print(f"  {j:>2} {x:>5} {str(M):>12} {str(cf):>10} {str(M==cf):>6} {c}/{q:<7} {'YES' if ingap else '':>8}", flush=True)

print("\nCLOSED-FORM LOWER BOUND at t=(j+1)/(6j+5), q=6j+5, c=j+1 -- residue-distances:", flush=True)
print("  runner 1 -> j+1 ;  2 -> 2j+2 ;  3 -> 3j+2 ;  4 -> 2j+1 ;  5 -> j ;  7 -> j+2 ;  6j -> j", flush=True)
print("  => min = j (runners 5 and 6j bind), so M >= j/(6j+5).  [equality computational above]", flush=True)
for j in range(3,8):
    q=6*j+5; c=j+1
    ds={r: min((r*c)%q, q-(r*c)%q) for r in base+[6*j]}
    assert min(ds.values())==j, (j,ds)
print("  verified min-dist = j for j=3..7.", flush=True)

print("\nUNIQUENESS: j/(6j+5) in (1/8,2/15)  <=>  8j>6j+5 (j>=3) AND 15j<12j+10 (j<=3)  <=>  j=3.", flush=True)
print("  => the ONLY gap member in this slice is x=6*3=18, M=3/23 (mediant, denom 23=3k+2 wall).", flush=True)
print("  The gap window catches exactly ONE rung of the resonance ladder.", flush=True)
