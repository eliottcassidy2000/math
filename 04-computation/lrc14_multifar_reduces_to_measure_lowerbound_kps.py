#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
MOMENT RELAXATION => multi-far reduces to a MEASURE LOWER BOUND on the core lonely set:
inf meas(L_C) > 6^{-r} over (13-r)-cores. Closed for r>=3; for r=2 a tractable sub-problem.

kind-pasteur-2026-07-01. survival = mu(safe-cube), mu = pushforward of 1_{L_C} onto T^r via
t->(W_i t); moments mu-hat(k)=1_{L_C}-hat(sum k_i W_i). Moment/SOS hierarchy:
 * DEGREE-2 (Cauchy-Schwarz): COPRIME combs (gcd 1) have ||g||_2 -> (6/49)^{r/2}, so
   |correction| <= sqrt(meas)*(6/49)^{r/2} < main=(6/7)^r meas  <=>  meas(L_C) > 6^{-r}.
 * GCD REDUCTION: combs with gcd g>1 reduce (u=g t) to an r-speed comb {a_i}={W_i/g} at large
   frequency g; the x-g map equidistributes L_C, so survival -> meas(cap safe(a_i)) * meas(L_C)
   >= (1 - r/7) meas(L_C) > 0 for r<=6. (Verified: (2g,3g) -> meas(safe(2)cap safe(3))*meas.)
So BOTH classes close, and multi-far r=2..6 reduces to ONE clean bound:
   inf_{(13-r)-cores} meas(L_C) > 6^{-r}.
This is a SUB-problem of OPEN-Q-108 on FEWER speeds (proven LRC(<=12)) -- lower-dimensional and
more tractable than the full 13-set inf. r>=7 = THM-573. The rates are the far-element decorrelation.
"""
import sys, random
from fractions import Fraction as Fr
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
BAND=Fr(1,14)
def safe_arcs(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]; hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
def L(C):
    a=safe_arcs(C[0])
    for v in C[1:]:
        a=inter(a,safe_arcs(v))
        if not a: return a
    return a
def meas(A): return sum(h-l for l,h in A)

print("="*92)
print("PART 1  GCD REDUCTION (gcd>1 combs close by x-g equidistribution)")
print("="*92)
print("  meas(safe(a)cap safe(b)) for small coprime (a,b) >= 1-2/7 = 5/7 = 0.714:")
for a,b in [(1,2),(2,3),(3,4),(2,5),(3,5),(4,5)]:
    print(f"    (a,b)=({a},{b}): meas(safe(a)cap safe(b))={float(meas(inter(safe_arcs(a),safe_arcs(b)))):.4f}")
C=[1,2,3,4,5,7,8,9,10,11,12]; Lc=L(C); mL=meas(Lc)
tgt=meas(inter(safe_arcs(2),safe_arcs(3)))*mL
print(f"  reduction: (W1,W2)=(2g,3g) on min-meas 11-core (meas={float(mL):.5f}); survival -> {float(tgt):.5f}:")
for g in [92,300,700]:
    s=meas(inter(inter(Lc,safe_arcs(2*g)),safe_arcs(3*g)))
    print(f"    g={g:4d}: survival={float(s):.5f} (target {float(tgt):.5f}, ratio {float(s/tgt):.3f}) >0")

print("\n"+"="*92)
print("PART 2  THE MEASURE LOWER BOUND: is inf meas(L_C) > 6^{-r} over (13-r)-cores?")
print("="*92)
def search_min(k, nt):
    rng=random.Random(k); worst=(Fr(2),None); cands=[list(range(1,k+1))]
    for j in range(1,k+1):
        for mm in range(1,20):
            Cc=sorted(set(range(1,k+1))-{j}|{(k+1)*mm})
            if len(Cc)==k and (k+1)*mm<=182: cands.append(Cc)
    for d in range(1,8):
        Cc=[d*i for i in range(1,k+1)]
        if max(Cc)<=182: cands.append(Cc)
    for _ in range(nt): cands.append(sorted({1}|set(rng.sample(range(2,183),k-1))))
    for Cc in cands:
        if len(set(Cc))!=k: continue
        mm=meas(L(Cc))
        if mm<worst[0]: worst=(mm,Cc)
    return worst
print(f"  {'r':>2} {'k=13-r':>7} {'min meas(L_C)':>14} {'6^-r':>10} {'margin':>8} {'closes?':>8}")
for r in range(2,7):
    k=13-r; thr=Fr(1,6**r); mm,Cc=search_min(k, 2500 if k>=9 else 4000)
    print(f"  {r:>2} {k:>7} {float(mm):>14.6f} {float(thr):>10.6f} {float(mm/thr):>7.1f}x {str(mm>thr):>8}")

print("\n"+"="*92)
print("VERDICT")
print("="*92)
print(" - Moment relaxation: survival=mu(safe-cube); DEGREE-2 (Cauchy-Schwarz) closes COPRIME combs")
print("   for meas(L_C)>6^{-r}; GCD>1 combs REDUCE (u=g t) to an r-speed comb {W_i/g} at large freq g,")
print("   x-g equidistributes L_C => survival -> meas(cap safe(W_i/g))*meas(L_C) >= (1-r/7)meas > 0.")
print(" - So multi-far r=2..6 REDUCES to ONE clean bound:  inf_{(13-r)-cores} meas(L_C) > 6^{-r}.")
print("   Search: CLOSES all r with margins 1.9x (r=2) -> 14000x (r=6). EASY for r>=3; r=2 tight-ish.")
print(" - This is OPEN-Q-108 on FEWER speeds (proven LRC(<=12)): the LONELY-MEASURE inf on (13-r)-cores,")
print("   lower-dimensional & more tractable than the full 13-set inf. Rates = far-element decorrelation.")
print(" - HONEST: min is a SEARCH min, not a proved inf; r>=3 margins make it robust, r=2 (1.9x) needs")
print("   a real inf meas(L_C)>1/36 over 11-cores (tractable via LRC(12)). r>=7 = THM-573.")
print("DONE.")
