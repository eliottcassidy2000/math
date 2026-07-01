#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ATTEMPT the facility-location PoA bound for inf meas >= 1/36, and find it REDUCES to the signed resonance =
the Gauss certificate. Plus the 21-Frobenius subgroup.

kind-pasteur-2026-07-01-S25. CORRECTION to S24: the r=2 residual is over 11-CORES (not 13 speeds).
 * meas(L_C) = (6/7)^k * (1 + R), R = the resonance correction. inf meas >= 1/36 <=> 1+R >= (1/36)/(6/7)^k.
 * ATTEMPT via 2nd moment / Gram: the congestion C(t)=#dangerous runners has mean k/7 > 1, so meas{C=0} is a
   BELOW-MEAN TAIL; second-moment/Paley-Zygmund/Markov give UPPER bounds on meas (WRONG direction) => the lower
   bound cannot come from an unsigned/absolute estimate.
 * The resonance sum R converges only with SIGNED cancellation (Sum|.| diverges, MISTAKE-078) => the PoA bound
   IS the signed resonance sum = the ι-odd Gauss certificate (S20) + tight-locus finiteness (THM-523). The
   game-theoretic and arithmetic routes are the SAME bound.
 * 21-Frobenius F_21 = C_7 x| C_3 = Aut(Paley_7) < PSL_2(F_7): the apex symmetry; 21 = 3*7 = compositum modulus.
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
BAND=Fr(1,14)
def safe(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def danger(v): return [((Fr(k)-BAND)/v,(Fr(k)+BAND)/v) for k in range(v)]  # arcs around k/v, width 2r
def inter(A,B):
    r=[];i=j=0
    A=sorted(A); B=sorted(B)
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def meas(A): return sum(h-l for l,h in A)
def L_of(S):
    a=safe(S[0])
    for v in S[1:]:
        a=inter(a,safe(v))
        if not a: return a
    return a
def pair_danger_meas(v,w):
    # meas{ ||vt||<r AND ||wt||<r } on [0,1): intersect danger arcs mod 1
    dv=[(a%1,b%1) if b<=1 else (a,b) for a,b in danger(v)]
    # simpler: sample-free exact via intersecting danger(v) and danger(w) as subsets of [0,1)
    def norm(arcs):
        out=[]
        for a,b in arcs:
            a=a-(a.numerator//a.denominator); # frac part of a
            if a<0: a+=1
            bb=a+(b-arcs[0][0] if False else (Fr(1,7)/v if False else b - ( ( (arcs.index((a,b)) if False else 0) )) ))
        return arcs
    Dv=danger(v); Dw=danger(w)
    # wrap into [0,1): each danger arc ((k-r)/v,(k+r)/v); for k=0 it's (-r/v, r/v) -> split
    def wrap(arcs):
        out=[]
        for a,b in arcs:
            if a<0: out.append((a+1,Fr(1))); out.append((Fr(0),b))
            elif b>1: out.append((a,Fr(1))); out.append((Fr(0),b-1))
            else: out.append((a,b))
        return sorted(out)
    return meas(inter(wrap(Dv),wrap(Dw)))

pent=[1,2,3,4,5,7,8,9,11,12,13]  # (Z/10)* extremizer 11-core
k=len(pent)
mP=meas(L_of(pent)); indep=(Fr(6,7))**k; R=mP/indep-1; target=Fr(1,36)/indep-1
print("="*96); print(" CORRECTED PoA (11-CORE, the r=2 residual): meas = (6/7)^k (1+R)"); print("="*96)
print(f"  pentagon 11-core {pent}: k={k}, meas={float(mP):.6f}, (6/7)^{k}={float(indep):.6f}")
print(f"  1+R = meas/indep = {float(mP/indep):.6f}  => R = {float(R):.6f}  (the resonance correction)")
print(f"  inf meas>=1/36 <=> 1+R >= (1/36)/(6/7)^{k} = {float(Fr(1,36)/indep):.6f} <=> R >= {float(target):.6f}")
print(f"  actual R={float(R):.4f} >= target {float(target):.4f}? {R>=target}  (clears by {float(R-target):.4f}, TIGHT)")
print(f"  PRICE OF ANARCHY = (6/7)^{k}/meas = {float(indep/mP):.3f}; target PoA <= (6/7)^{k}*36 = {float(indep*36):.3f}")
print(f"  [S24 CORRECTION: used k=13 (whole set) -> 4.18/4.85; the RESIDUAL is 11-cores -> {float(indep/mP):.2f}/{float(indep*36):.2f}]")

print("\n"+"="*96); print(" WHY 2nd-moment FAILS: meas{C=0} is a BELOW-MEAN TAIL (wrong direction)"); print("="*96)
EC=k*2*BAND
print(f"  congestion C(t)=#dangerous runners: E[C] = k*2r = {k}*(1/7) = {float(EC):.4f} > 1")
print(f"  => meas{{C=0}} is a BELOW-mean event; Markov/Paley-Zygmund/2nd-moment give UPPER bounds on meas (P(C>0)")
print(f"     >= E[C]^2/E[C^2]), i.e. LOWER bounds FAIL. The lower bound CANNOT be an unsigned/absolute estimate.")
# danger-Gram G[v,w]=I_{v,w}=meas{both dangerous}; trace = E[C]
import numpy as np
G=np.zeros((k,k))
for i in range(k):
    for j in range(k):
        G[i,j]=float(pair_danger_meas(pent[i],pent[j]))
tr=np.trace(G); eig=sorted(np.linalg.eigvalsh(G),reverse=True)
print(f"  danger-Gram trace = {tr:.4f} (=E[C]={float(EC):.4f}); top eigenvalues {[round(e,4) for e in eig[:4]]}")
print(f"  (mac-mini Bridge-2: the Gram 2nd moment; its spectrum bounds |R| in MAGNITUDE but NOT the signed value.)")

print("\n"+"="*96); print(" THE REDUCTION: the PoA bound = the SIGNED resonance = the Gauss certificate (S20)"); print("="*96)
print("  R = Sum_{nonzero resonances Sum m_v v=0} prod_v [ghat(m_v)/(6/7)].  Sum|.| DIVERGES (MISTAKE-078); R")
print("  converges only with SIGNED cancellation. So bounding R >= -0.849 needs the signed structure = the ι-odd")
print("  Lefschetz/Gauss certificate (S20, g_7=i sqrt7) + the tight-locus finiteness (THM-523 {AP,GW}). The")
print("  facility-location PoA bound and the arithmetic certificate are the SAME bound seen two ways.")
print("  => HONEST: not a discrepancy proof; the tight 1/36 route is (a) the FINITE near-tight census (THM-522")
print("     quantization bounds it) + (b) singular-series positivity (THM-501). PoA<=6.61 is the clean restatement.")

print("\n"+"="*96); print(" THE 21-FROBENIUS SUBGROUP + owner's faces"); print("="*96)
print("  F_21 = C_7 x| C_3 = Aut(Paley_7), order 21 = 3*7 = the compositum modulus (S23 sqrt21). F_21 < PSL_2(F_7)")
print("   (168 = 8*21): the apex Frobenius sits inside the LTC expander group. C_7 = rotation (the runners'")
print("   cyclic loop), C_3 = the QR multiplier (the tournament orientation). 21 = |F_21| = forbidden H (THM-079).")
print("  UNIFYING sqrt-p (odd/certificate): Gauss i sqrt p / Paley skew {0,±i sqrt p} / Ramanujan 2 sqrt p / Q(sqrt-p).")
print("  POCHHAMMER pi (even/measure): f(n)=(1/2)_{n-2}/(n-2)! ~ 1/sqrt(pi n) -- the (6/7)^k independent side.")
print("  So the PoA decomposition meas=(6/7)^k[1+R] IS the pi-side (independent) x the sqrt-p-side (signed R).")
print("DONE.")
