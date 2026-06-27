#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spec_tail_control_kpswf15.py  (kind-pasteur 2026-06-27, TOOL 3)

THE UNIFORM TAIL CONTROL: hybrid (exact LOW + L2 Cauchy-Schwarz HIGH).

This closes the wide-V residual that the 1/n total-variation tail left open.
The cross-correlation, on the resonance lattice n=14N:
   SPEC = SPEC_low + SPEC_high,
   SPEC_low  = sum_{1<=N<=M} 2 Re[chat(14N) conj(ghat(14N))]   (finite, EXACT),
   SPEC_high = sum_{N>M}    2 Re[chat(14N) conj(ghat(14N))].

UNIFORM HIGH-TAIL BOUND (L2 Cauchy-Schwarz, restricted to the resonance lattice):
   |SPEC_high| <= 2 * ( sum_{N>M} |chat(14N)|^2 )^{1/2} ( sum_{N>M} |ghat(14N)|^2 )^{1/2}
              =: 2 * Tc(M) * Tg(M),
where
   Tc(M)^2 = sum_{N>M} |chat(14N)|^2 = ResL2c - sum_{N<=M}|chat(14N)|^2,
   ResL2c  = sum_{N!=0} |chat(14N)|^2   (L2 mass of chat ON the resonance lattice 14Z),
and similarly Tg(M), ResL2g.   The factor 2 collects +-N.

KEY POINT (why this beats the TV bound):  ResL2c is the L2 mass of chat RESTRICTED to
14Z, NOT all of Z.  We have the EXACT Parseval ceiling
   sum_{all n!=0}|chat(n)|^2 = meas(Rsafe) - meas(Rsafe)^2 = var_R,
and ResL2c <= var_R (restriction can only shrink).  But the resonance restriction is the
crucial gain: most of chat's L2 mass is OFF 14Z (on low residues), so ResL2c << var_R.
We compute ResL2c exactly enough by summing |chat(14N)|^2 to large N (geometric tail).

THE FLOOR.  R' = 1 + SPEC/baseline >= 1 + (SPEC_low - |SPEC_high|)/baseline
                                     >= 1 + (SPEC_low - 2 Tc(M) Tg(M))/baseline.
We tabulate, over the multi-far family and increasing M, the resulting certified floor.
If for some fixed modest M the certified floor is uniformly > 0, the uniform R' >= c is
closed by ELEMENTARY harmonic analysis (finite low part + L2 tail), no BV/EH.
"""
import sys, math
from fractions import Fraction as F
from math import gcd, pi, sqrt
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    meas, intersect, safe_set, fourier_num_of_arcs,
)

def coeff(arcs, n):
    if n==0: return complex(float(meas(arcs)))
    return fourier_num_of_arcs(arcs,n)/(2j*pi*n)

def Qlon_t(Q):
    return safe_set([14*m for m in Q])

def resonance_L2(arcs, Nmax):
    """sum_{1<=N<=Nmax} |coeff(arcs,14N)|^2  (one-sided; double for +-)."""
    s=0.0
    parts=[]
    for N in range(1,Nmax+1):
        v=abs(coeff(arcs,14*N))**2
        s+=v
        parts.append((N,s))
    return s, parts

def full_L2_partial(arcs, nmax):
    """sum_{1<=n<=nmax} |hat f(n)|^2  (one-sided; full circle, ALL freqs).  Exact via numeric coeff."""
    s=0.0
    for n in range(1,nmax+1):
        s += abs(coeff(arcs,n))**2
    return s

def analyze(R,Q,Mlist=(5,10,20,40),Nmax=2000,label=""):
    R=sorted(set(int(x) for x in R if x!=0))
    Q=sorted(set(int(x) for x in Q if x!=0))
    Rsafe=safe_set(R); Qlon=Qlon_t(Q)
    mR=meas(Rsafe); mQ=meas(Qlon); mB=meas(intersect(Rsafe,Qlon))
    baseline=mR*mQ; bf=float(baseline)
    SPEC_exact=float(mB-baseline); Rprime=float(mB/baseline) if baseline>0 else 1.0
    var_R=float(mR-mR*mR); var_Q=float(mQ-mQ*mQ)

    # resonance-lattice L2 masses, numeric, to track how much of var lands on 14Z (diagnostic).
    L2c_num,_ = resonance_L2(Rsafe, Nmax)
    L2g_num,_ = resonance_L2(Qlon, Nmax)

    print("="*94)
    print(f"  ({label})  R={R}  Q={Q}")
    print(f"  baseline={bf:.6f}  SPEC={SPEC_exact:+.6f}  R'={Rprime:.5f}")
    print(f"  var_R/2 (one-sided all-freq L2)={var_R/2:.6f}  resonance-L2c(14Z, numeric)={L2c_num:.6f}  frac on 14Z={L2c_num/(var_R/2):.4f}")
    print(f"  var_Q/2 (one-sided all-freq L2)={var_Q/2:.6f}  resonance-L2g(14Z, numeric)={L2g_num:.6f}  frac on 14Z={L2g_num/(var_Q/2):.4f}")
    print("="*94)
    # RIGOROUS tail ceiling (full-circle Parseval): one-sided sum_{n>=1}|hat f(n)|^2 = var/2 EXACT.
    #   tail_resonance(M) = sum_{N>M}|chat(14N)|^2 <= sum_{n>14M}|chat(n)|^2 = var_R/2 - full_L2_partial(R,14M).
    best=(None,-9)
    for M in Mlist:
        SPEC_low=0.0
        for N in range(1,M+1):
            SPEC_low += 2.0*(coeff(Rsafe,14*N)*coeff(Qlon,14*N).conjugate()).real
        # rigorous one-sided tail ceilings
        Sc=full_L2_partial(Rsafe, 14*M)   # sum_{n<=14M}|chat|^2 (one-sided)
        Sg=full_L2_partial(Qlon, 14*M)
        tail_c2=max(var_R/2 - Sc, 0.0)    # >= sum_{N>M}|chat(14N)|^2
        tail_g2=max(var_Q/2 - Sg, 0.0)
        Tc=sqrt(tail_c2); Tg=sqrt(tail_g2)
        high_bound=2.0*Tc*Tg              # |SPEC_high| <= 2 (sum_{N>M}|chat|^2)^{1/2}(...)^{1/2}
        floor=1.0 + (SPEC_low - high_bound)/bf
        print(f"   M={M:3d}: SPEC_low={SPEC_low:+.6f}  RIGOROUS |SPEC_high|<=2*Tc*Tg={high_bound:.6f}  "
              f"=> certified R' >= {floor:+.5f}")
        if floor>best[1]: best=(M,floor)
    print(f"   best certified floor over M in {list(Mlist)}: R' >= {best[1]:+.5f} at M={best[0]}")
    return dict(R=R,Q=Q,baseline=baseline,Rprime=Rprime,bf=bf,label=label,
                best_floor=best[1],best_M=best[0],
                L2c=L2c_num,L2g=L2g_num,var_R=var_R,var_Q=var_Q)

def main():
    print("#"*94)
    print("# TOOL 3 : UNIFORM TAIL CONTROL (exact LOW + L2 Cauchy-Schwarz HIGH on resonance lattice)")
    print("#   R' >= 1 + (SPEC_low - 2 Tc(M) Tg(M))/baseline.  Is the certified floor uniformly > 0?")
    print("#"*94)
    cases = [
        ([1,2,3,4,5,6,7,8,9,10,11,13],[1],          "R={1..13}\\{12}, Q={1} r=1"),
        ([1,2,3,4,5,6,7,8,9,10,11,12],[1,2],        "R={1..12}, Q={1,2} r=2 (worst R'=0.70)"),
        ([1,2,3,4,5,6,7,8,9,10,11],  [1,2,3],      "R={1..11}, Q={1,2,3} r=3"),
        ([1,2,3,4,5,6,7,8,9,10],     [1,2,3,4],    "R={1..10}, Q={1,2,3,4} r=4"),
        ([1,2,3,4,5,6,7,8,9],        [1,2,3,4,5],  "R={1..9}, Q={1..5} r=5"),
        ([1,2,3,4,5,6,7,8],          [1,2,3,4,5,6],"R={1..8}, Q={1..6} r=6 max"),
        ([2,4,6,8,10,12],            [1,3,5],      "R gcd=2, Q={1,3,5} (worst SPEC=-0.029)"),
        ([3,6,9,12],                 [1,2,5],      "R gcd=3, Q={1,2,5}"),
        ([5,9,11,13],                [1,2,3],      "R coprime-ish, Q={1,2,3}"),
        # extra stress: small R-safe (deep tight perturbation), single far
        ([1,2,3,4,5,6,7,8,9,10,11,13],[2],          "R={1..13}\\{12}, Q={2} (far=28 only)"),
        ([1,2,3,4,5,6,7,8,9,10,12,13],[1],          "R={1..13}\\{11}, Q={1}"),
    ]
    res=[analyze(R,Q,Mlist=(5,10,20,40,80),Nmax=3000,label=lab) for R,Q,lab in cases]
    print("="*94)
    print("SUMMARY: best certified uniform floor (exact low + L2 tail)")
    print("="*94)
    print(f"{'label':<44}{'R-prime':>9}{'cert.floor':>12}{'at M':>6}{'L2c/varR':>10}{'L2g/varQ':>10}")
    for r in res:
        print(f"{r['label'][:44]:<44}{r['Rprime']:>9.4f}{r['best_floor']:>12.5f}{r['best_M']:>6}"
              f"{r['L2c']/r['var_R']:>10.4f}{r['L2g']/r['var_Q']:>10.4f}")
    minfloor=min(r['best_floor'] for r in res)
    print(f"\n  MIN certified floor over family = {minfloor:+.5f}")
    if minfloor>0:
        print(f"  ==> UNIFORM R' >= {minfloor:.4f} > 0 CERTIFIED by elementary (finite low + L2 tail).")
    else:
        print(f"  ==> certified floor still touches/below 0 at modest M; needs larger M or finer low part.")
    print("DONE.")

if __name__=="__main__":
    main()
