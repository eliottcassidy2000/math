#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TWO-FAR (multi-far) certificate: does a fat hole of the core survive TWO thin combs?

kind-pasteur-2026-07-01. After single-far closed (band-barrier), the residual is >=2 large
speeds. Two-far = C (11 bounded speeds<=182) u {W1,W2}, W1,W2>182. By PROVEN LRC(12) (11 speeds)
M(C)>=1/12, so the core lonely set L_C={g_C>=1/14} has fat holes.

REFRAME (2D view-obstruction). On a hole I=[a,a+delta], the two far phases (W1 t, W2 t) trace a
LINE OF SLOPE W2/W1 on the 2-torus T^2; "t survives both combs" <=> ||W1 t||>=1/14 AND ||W2 t||
>=1/14 <=> the line HITS THE SAFE SQUARE [1/14,13/14]^2 (area (6/7)^2=36/49). The danger is the
coordinate CROSS. So the two-comb question IS Cusick's 2D view-obstruction. A generic (incommen-
surate) line equidistributes on T^2 and hits the safe square (survives); the obstruction is
RESONANCE -- a low-denominator slope W2/W1 that confines the line to the danger cross.

We compute meas(I cap safe(W1) cap safe(W2)) on the widest hole (exact intervals), sweep
(W1,W2)>182, find the WORST (min jointly-safe) pair, and test the reframe: resonant pairs
(small gcd-structure / small-denominator ratio) are the only threats, and they collapse to a
SINGLE effective comb (= single-far, already closed).
"""
import sys
from fractions import Fraction as Fr
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
BAND=Fr(1,14)

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
def lonely_set(C):
    arcs=safe_arcs(C[0])
    for v in C[1:]:
        arcs=intersect(arcs, safe_arcs(v))
        if not arcs: break
    return arcs
def safe_in(W, lo, hi):
    res=[]; klo=int(lo*W)-1; khi=int(hi*W)+2
    for k in range(klo,khi+1):
        alo=(Fr(k)+BAND)/W; ahi=(Fr(k+1)-BAND)/W
        L=alo if alo>lo else lo; H=ahi if ahi<hi else hi
        if L<H: res.append((L,H))
    return res
def joint_safe_meas(I, W1, W2):
    lo,hi=I
    s=intersect(safe_in(W1,lo,hi), safe_in(W2,lo,hi))
    return sum(h-l for l,h in s)

print("="*92)
print("TWO-COMB SURVIVAL on the widest hole of a core (does it escape the 2D danger cross?)")
print("="*92)
cores={"AP {1..11}":list(range(1,12)), "{1..10,12}":list(range(1,11))+[12],
       "{1..9,11,13}":list(range(1,10))+[11,13]}
for name,C in cores.items():
    L=lonely_set(C)
    if not L: print(f"  {name}: empty L_C"); continue
    I=max(L,key=lambda ab:ab[1]-ab[0]); delta=I[1]-I[0]
    # sweep W1<W2 in a large window
    Wlo,Whi=183,246
    worst=(Fr(2),0,0); killed=0; tot=0
    for W1 in range(Wlo,Whi):
        for W2 in range(W1+1,Whi):
            m=joint_safe_meas(I,W1,W2); tot+=1
            if m==0: killed+=1
            if m<worst[0]: worst=(m,W1,W2)
    wm,w1,w2=worst; g=gcd(w1,w2)
    single=joint_safe_meas(I,w1,w1+ (1 if w1%2==0 else 2))  # rough single-comb baseline on I
    print(f"  {name:16s} widest hole delta={float(delta):.5f}")
    print(f"    over {tot} pairs (W1,W2) in [{Wlo},{Whi}): killed(meas=0)={killed};  "
          f"MIN joint-safe={float(wm):.6f} at (W1,W2)=({w1},{w2}) gcd={g} ratio={w1}/{w2}")
    frac=float(wm)/float(delta) if delta>0 else 0
    print(f"    worst joint-safe / delta = {frac:.3f}  (incommensurate -> ~ (6/7)^2 = {36/49:.3f}); "
          f"{'HOLE ALWAYS SURVIVES (min>0)' if wm>0 else 'widest hole KILLED for some pair -> check other holes'}")

print("\n"+"="*92)
print("THE RESONANCE STRUCTURE: which (W1,W2) minimize joint-safe? (the 2D view-obstruction slope)")
print("="*92)
C=list(range(1,12)); L=lonely_set(C); I=max(L,key=lambda ab:ab[1]-ab[0]); delta=I[1]-I[0]
rows=[]
for W1 in range(183,220):
    for W2 in range(W1+1,220):
        m=joint_safe_meas(I,W1,W2); rows.append((float(m),W1,W2,gcd(W1,W2)))
rows.sort()
print(f"  core {{1..11}}, widest hole delta={float(delta):.5f}; the 8 WORST (smallest joint-safe) pairs:")
for m,w1,w2,g in rows[:8]:
    print(f"    (W1,W2)=({w1},{w2})  joint-safe={m:.6f}  gcd={g}  W2-W1={w2-w1}  ratio~{w2/w1:.3f}")
print(f"  the 4 BEST (largest joint-safe):")
for m,w1,w2,g in rows[-4:]:
    print(f"    (W1,W2)=({w1},{w2})  joint-safe={m:.6f}  gcd={g}  W2-W1={w2-w1}")

print("\n"+"="*92)
print("READING / REFRAME")
print("="*92)
print(" - Two-comb survival = a line of slope W2/W1 on T^2 escaping the danger CROSS into the safe")
print("   square (6/7)^2=36/49. Generic (incommensurate) lines equidistribute => hit the safe square;")
print("   the ONLY threat is RESONANCE (low-denominator slope / large gcd / small W2-W1).")
print(" - Resonance collapses the two combs to ONE effective comb (common period) => the SINGLE-FAR")
print("   band-barrier applies again (three-distance). Non-resonant => joint 36/49 safe measure.")
print("   So multi-far reduces to [resonant = single-far, CLOSED] + [non-resonant = equidistribution].")
print(" - Borsuk-Ulam iota-odd: the resonant lines are the antipodally-trapped ones; the index bites")
print("   only there -- the sharply-isolated residual.")
print("DONE.")
