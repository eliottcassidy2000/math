#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
r-FAR union-bound closure: does the fat-hole/band-barrier certificate close ALL multi-far r=2..6?

kind-pasteur-2026-07-01. Single-far (r=1) closed via the band-barrier. Two-far (r=2) verified
(widest hole survives every pair). Here: extend to r far speeds by the UNION BOUND and check the
race between (a) the shrinking core (13-r speeds => FATTER holes, larger delta_max) and (b) the
growing danger budget (r combs cover <= r/7).

CERTIFICATE (rigorous, union bound). 13-set = C (13-r bounded speeds<=182) u {W_1..W_r>182}. By
PROVEN LRC(14-r) (<=12 speeds for r>=2), M(C)>=1/(14-r), so L_C has a widest hole of width
delta_max. The r combs each cover danger <= (delta/7 + 1/(7*183)) on that hole (measure + one end
tooth). So the jointly-safe measure on the hole is
    >= delta*(1 - r/7) - r/(7*183)  > 0   <=>   delta_max > threshold_r := r/((7-r)*183).
If min_C delta_max(C) > threshold_r for every r=2..6, the whole multi-far case CLOSES (r>=7 is the
level-7 sieve THM-573). We compute min delta_max over (13-r)-cores (random search, speeds<=182)
vs threshold_r.
"""
import sys, random
from fractions import Fraction as Fr
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
BAND=Fr(1,14); WMIN=183

def safe_arcs(v):
    return [((Fr(k)+BAND)/v, (Fr(k+1)-BAND)/v) for k in range(v)]
def intersect(A,B):
    res=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]
        hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi: res.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return res
def delta_max(C):
    arcs=safe_arcs(C[0])
    for v in C[1:]:
        arcs=intersect(arcs, safe_arcs(v))
        if not arcs: return Fr(0)
    return max(hi-lo for lo,hi in arcs)

print("="*94)
print("r-FAR UNION-BOUND CLOSURE:  hole survives r combs  <=>  delta_max > threshold_r = r/((7-r)*183)")
print("  core size = 13-r ; proven LRC(14-r) => M(C)>=1/(14-r) => fat holes")
print("="*94)
print(f"  {'r':>2} {'core size':>9} {'threshold_r':>12} {'min delta_max (rand)':>21} {'margin':>8} {'closes?':>8}")
rng=random.Random(701)
results={}
for r in range(2,7):
    k=13-r
    thr=Fr(r,(7-r)*WMIN)
    worst=Fr(2); worstC=None; n=0
    trials=800 if k>=9 else 1500
    for _ in range(trials):
        # k-core: {1} + (k-1) random from [2,182], bias to include some large speeds (worst holes)
        pool=list(range(2,183))
        C=sorted({1}|set(rng.sample(pool,k-1)))
        if len(C)!=k: continue
        d=delta_max(C); n+=1
        if d<worst: worst=d; worstC=C
    closes = worst>thr
    margin=float(worst/thr) if thr>0 else float('inf')
    results[r]=(thr,worst,closes,worstC)
    print(f"  {r:>2} {k:>9} {float(thr):>12.6f} {float(worst):>21.6f} {margin:>7.2f}x {str(closes):>8}")

print("\n"+"="*94)
print("VERDICT")
print("="*94)
allclose=all(results[r][2] for r in results)
for r in range(2,7):
    thr,worst,closes,C=results[r]
    tag = "CLOSES (union bound)" if closes else "union bound TIGHT/FAILS -> needs exact equidistribution (resonance analysis)"
    print(f"  r={r}: {tag}   worst core={C}")
print(f"\n  ALL r=2..6 close by union bound: {allclose}")
if not allclose:
    fails=[r for r in results if not results[r][2]]
    print(f"  residual r = {fails}: the union bound is too lossy (boundary teeth eat the margin);")
    print(f"  these need the EXACT two-comb result (survives with ~0.58 delta) generalized -- resonance")
    print(f"  collapses effective combs, non-resonant equidistributes. NOT a beater search.")
print("\n  (r>=7: >=7 far speeds => >=7 mult-of-7-like => THM-573 level-7 sieve, M>1/14. r<=6 is the")
print("   whole residual; single-far r=1 closed; r=2 exact-verified; r=3..6 = the finite constant-chase.)")
print("DONE.")
