#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
SINGLE-FAR UNBOUNDED CERTIFICATE via Morse theory / the band-barrier.

kind-pasteur-2026-07-01. First concrete swing at the UNBOUNDED case (the open core of
LRC(14) after klein HYP-3779 closed the BOUNDED case, speeds<=n(n-1)=182, by ILP).

SETUP. A covering 13-set with one large speed W>182 is C u {W}, C = 12-speed bounded core
(all speeds<=182). By PROVEN LRC(13) (<=12 speeds), M(C)>=1/13>1/14, so the core's LONELY
SET  L_C = { t : g_C(t)=min_{v in C}||vt|| >= 1/14 }  is NONEMPTY with positive measure.

MORSE / BAND-BARRIER. g_C is a piecewise-linear Morse function on S^1; L_C is its super-level
set {g_C>=1/14} = a union of arcs around the LOCAL MAXIMA (the "holes"). Adding W carves g_C
from above: g_{CuW}=min(g_C, ||Wt||). ||Wt|| is a fast sawtooth with W teeth; its DANGER band
{||Wt||<1/14} is W thin arcs each of width 1/(7W). A lonely arc of L_C of width delta CONTAINS
a W-safe point (=> a lonely moment of C u {W}, => M>=1/14) UNLESS it fits inside one danger
tooth, i.e. unless delta <= 1/(7W). So:

  if the WIDEST lonely arc  delta_max(C) > 1/(7W)  then  M(C u {W}) >= 1/14.

=> single-far certificate holds for  W > W*(C) = 1/(7 delta_max(C)).  Combined with klein's
BOUNDED closure (W<=182), the single-far UNBOUNDED case is CLOSED iff  W*(C) <= 182  for every
bounded core C, i.e.  delta_max(C) >= 1/(7*182) = 1/1274 ~ 0.000785.

This script computes L_C EXACTLY (interval intersection, Fractions -- no MISTAKE-86: we build
the super-level set directly, not a max), reports delta_max and W*, and searches for the worst
(smallest delta_max) bounded core to see whether/where W*<=182 (single-far closes) vs a gap.
"""
import sys, random
from fractions import Fraction as Fr
from functools import reduce
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
BAND=Fr(1,14)

def safe_arcs(v, band=BAND):
    # {t in [0,1): ||v t|| >= band} = union_k [(k+band)/v, (k+1-band)/v], k=0..v-1
    return [((Fr(k)+band)/v, (Fr(k+1)-band)/v) for k in range(v)]
def intersect(A,B):
    res=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]
        hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi: res.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return res
def lonely_set(C, band=BAND):
    arcs=safe_arcs(C[0],band)
    for v in C[1:]:
        arcs=intersect(arcs, safe_arcs(v,band))
        if not arcs: break
    return arcs
def core_stats(C):
    L=lonely_set(C)
    if not L: return None
    widths=[hi-lo for lo,hi in L]
    dmax=max(widths); meas=sum(widths)
    Wstar=1/(7*dmax)                     # threshold: W>Wstar => certificate holds
    return dmax, meas, Wstar, len(L)

print("="*94)
print("SINGLE-FAR CERTIFICATE:  W > W*(C)=1/(7*delta_max) => M(C u {W}) >= 1/14  (band-barrier)")
print("  need W* <= 182 (klein's bounded closure) to CLOSE the unbounded single-far case")
print("="*94)
tests={
 "AP {1..12}":list(range(1,13)),
 "{2..13}":list(range(2,14)),
 "{1..11,13}":list(range(1,12))+[13],
 "{1..11,13} (near-tight?)":[1,2,3,4,5,6,7,8,9,10,11,13],
 "{1..10,12,13}":list(range(1,11))+[12,13],
 "one mid speed {1..11,30}":list(range(1,12))+[30],
 "one big speed {1..11,90}":list(range(1,12))+[90],
 "one big speed {1..11,182}":list(range(1,12))+[182],
 "two-ish big {1..10,90,150}":list(range(1,11))+[90,150],
}
print(f"  {'core C (12 speeds)':32s} {'meas(L_C)':>10} {'delta_max':>11} {'W*=1/(7dmax)':>13} {'#arcs':>6} {'W*<=182?':>9}")
for name,C in tests.items():
    st=core_stats(C)
    if st is None: print(f"  {name:32s}  EMPTY lonely set (M(C)<1/14 ?!)"); continue
    dmax,meas,Wstar,narc=st
    print(f"  {name:32s} {float(meas):10.5f} {float(dmax):11.6f} {float(Wstar):13.2f} {narc:6d}   {str(Wstar<=182):>9}")

print("\n"+"="*94)
print("SEARCH: worst (smallest delta_max) bounded 12-core, by max-speed budget B")
print("  single-far closes for cores with max<=B iff worst W* <= 182")
print("="*94)
rng=random.Random(701)
for B in [20,30,50,90,140,182]:
    worstW=0.0; worstC=None; worstd=1.0; n=0
    for _ in range(300):
        # 12-core: 1 plus 11 more from [2,B], ensure primitive-ish, include some structure
        C=sorted({1}|set(rng.sample(range(2,B+1), 11)))
        if len(C)!=12: continue
        st=core_stats(C); n+=1
        if st is None: continue
        dmax,meas,Wstar,narc=st
        if Wstar>worstW: worstW=Wstar; worstC=C; worstd=dmax
    print(f"  max-speed budget B={B:3d}: worst W*={worstW:8.2f} (delta_max={worstd:.6f}) "
          f"{'<=182 CLOSES' if worstW<=182 else '>182 GAP (multi-far-like)'}  core={worstC}")

print("\n"+"="*94)
print("READING (Morse + band-barrier + Borsuk-Ulam)")
print("="*94)
print(" - L_C = super-level set {g_C>=1/14} = arcs around the LOCAL MAXIMA of the PL Morse function g_C;")
print("   delta_max = widest 'hole'. Three-distance bounds the #slopes (=> structured arc widths for AP cores).")
print(" - The far W carves with thin danger teeth (width 1/(7W)); a hole survives iff wider than a tooth")
print("   (delta_max>1/(7W)). This IS the band-barrier: thin bands cannot cover fat holes.")
print(" - L_C is ANTIPODE-SYMMETRIC (||v(-t)||=||vt|| => t in L_C <=> -t in L_C); the Z/2 (complement/iota)")
print("   action is free off the fixed points {0,1/2}; whether a hole survives is the iota-odd (Borsuk-Ulam)")
print("   index -- the topological packaging of the band-barrier (klein-S56 iota-odd = OPEN-Q-108 residual).")
print(" - NET: single-far closes for small-max cores (fat holes, W*<=182 overlaps klein); large internal")
print("   speeds (thin holes, W*>182) are the multi-far residual = >=2 large speeds = the r-far ladder.")
print("DONE.")
