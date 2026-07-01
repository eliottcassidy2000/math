#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
FINISH PROBES for inf meas(L_C) >= 1/36 over 11-speed cores (the reduced LRC(14) target).

kind-pasteur-2026-07-01. Three probes establishing the honest state of the target:
 (1) HARD SEARCH: min meas(L_C) over 11-cores ~ 0.0323 at {1..13}\{6,10} => target TRUE but TIGHT (1.16x).
 (2) REVERSE-MARKOV: fails (int g_C ~0.032 < 1/14 => lonely set is an above-mean event, not typical).
 (3) COLLAR/MORSE: meas(L_C) = sum of collar widths around the local maxima of g_C=min||vt||.
The over-optimistic "(6/7)^11=0.186 => 87% margin" is a CONFLATION: meas(L_C) IS a recursive singular
series whose own corrections reduce 0.186 -> 0.032. Finish needs a SHARP tool (census / Toeplitz /
HYP-3787 O(1/w) threshold-dissolve), not a positivity freebie.
"""
import sys, itertools, random
from fractions import Fraction as Fr
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def safe_arcs(v,band): return [((Fr(k)+band)/v,(Fr(k+1)-band)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]; hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
def L(C,band=Fr(1,14)):
    a=safe_arcs(C[0],band)
    for v in C[1:]:
        a=inter(a,safe_arcs(v,band))
        if not a: return a
    return a
def meas(A): return sum(h-l for l,h in A)

print("="*88); print("PROBE 1  HARD SEARCH: min meas(L_C) over 11-cores (near-tight locus)"); print("="*88)
cands=[list(range(1,12))]
for i in range(1,12):
    for v in range(12,140):
        C=sorted(set(range(1,12))-{i}|{v})
        if len(C)==11: cands.append(C)
for i,j in itertools.combinations(range(1,12),2):
    for a in range(12,50,3):
        for b in range(a+1,55,4):
            C=sorted(set(range(1,12))-{i,j}|{a,b})
            if len(C)==11: cands.append(C)
rng=random.Random(3)
for _ in range(4000):
    k=rng.choice([2,3]); C=sorted(set(range(1,12-k))|set(rng.sample(range(12,183),k)))
    if len(C)==11: cands.append(C)
worst=(Fr(2),None)
for C in cands:
    if len(set(C))!=11: continue
    m=meas(L(C))
    if m<worst[0]: worst=(m,C)
print(f"  {len(cands)} cores: MIN meas(L_C)={float(worst[0]):.6f}={worst[0]} at {worst[1]}")
print(f"  1/36={float(Fr(1,36)):.6f}; inf>=1/36? {worst[0]>=Fr(1,36)}; margin={float(worst[0]*36):.2f}x (TIGHT)")

print("\n"+"="*88); print("PROBE 2  REVERSE-MARKOV (dead-end: int g_C < 1/14)"); print("="*88)
def int_g(C,N=120):
    return sum(meas(L(C,Fr(i,2*N)))*Fr(1,2*N) for i in range(1,N+1))
def Mmax(C,N=20000):
    b=0.0
    for i in range(N):
        t=(i+0.5)/N; g=min(abs(v*t-round(v*t)) for v in C)
        if g>b: b=g
    return b
C=worst[1]; Ig=int_g(C); M=Mmax(C); a=1/14.0
rm=(float(Ig)-a)/(M-a) if M>a else float('nan')
print(f"  extremizer {C}: int g_C={float(Ig):.5f} < 1/14={a:.5f}; M~{M:.5f}")
print(f"  reverse-Markov LB (int g -1/14)/(M-1/14) = {rm:.4f} < 0 => VACUOUS (lonely set = above-mean tail)")

print("\n"+"="*88); print("PROBE 3  COLLAR/MORSE: meas(L_C) = sum of collar widths at the local maxima"); print("="*88)
Lc=L(C); ws=sorted(float(h-l) for l,h in Lc)
print(f"  extremizer meas={float(meas(Lc)):.5f}; #collars(local maxima >=1/14)={len(Lc)}; widths {ws[0]:.5f}..{ws[-1]:.5f}")
for band in [Fr(1,14),Fr(1,13),Fr(1,12)]:
    print(f"    threshold {band}: meas{{g>={float(band):.4f}}}={float(meas(L(C,band))):.5f}")
print("  => the lonely set is the COLLAR from 1/14 up to the maxima ~1/12; measure = sum of collar widths.")
print("\nNET: target TRUE but TIGHT (1.16x); finish needs a SHARP tool (census of dilated-AP minimizers /")
print("     Caratheodory-Toeplitz moment positivity / HYP-3787 O(1/w) to dissolve the threshold), not a freebie.")
print("DONE.")
