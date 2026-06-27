#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spec_level_distribution_kpswf15.py  (kind-pasteur 2026-06-27, TOOL 3)

THE LEVEL-OF-DISTRIBUTION STRUCTURE OF SPEC.   Honest BV-vs-EH assessment.

We make the cross-correlation SPEC into an explicit sum that an
equidistribution/discrepancy estimate would be asked to bound, then we identify
which estimate (elementary Weyl / large sieve / Bombieri-Vinogradov / Elliott-
Halberstam) the bound corresponds to -- and whether EH is genuinely needed.

EXACT FOURIER COEFFICIENTS.
  Let phi(theta) = 1_{||theta|| >= 1/14} on R/Z.  Its coefficients:
     phihat(0) = 6/7,
     phihat(k) = -sin(pi k / 7)/(pi k)    (k != 0).
  Note phihat(k)=0 exactly when 7 | k  (k!=0).  This is the "ghat vanishes on 7|n"
  fact (S_g(7k)=0): the apex prime 7 KILLS its own multiples.

  R-safe indicator c(t) = prod_{r in R} phi(r t).  In Fourier:
     chat(n) = sum over (k_r)_{r in R} with sum_r k_r r = n  of  prod_r phihat(k_r).
  So chat is supported on the LATTICE  Lat(R) = { sum_r k_r r }  and its value at n
  is a SUM OVER REPRESENTATIONS n = sum_r k_r r  (an additive-energy / singular-series
  object, HYP-2606).

  Q-lonely (in t, via u=14t): g(t) = prod_{m in Q} phi(14 m t).  In Fourier:
     ghat(n) = sum over (j_m) with sum_m j_m (14 m) = n  of  prod_m phihat(j_m)
             = [14 | n] * ( sum over (j_m) with sum_m j_m m = n/14  of  prod_m phihat(j_m) ).
  So ghat is supported on 14*Lat(Q) and ghat(14 N) = Ghat_Q(N), the analogous
  representation sum for Q at frequency N.

  SPEC = sum_{n != 0} chat(n) conj(ghat(n))
       = sum_{N != 0}  chat(14 N) * conj(Ghat_Q(N)).        [resonance: only n=14N survive]

THE LEVEL-OF-DISTRIBUTION SUM.
  chat(14 N) = sum_{(k_r): sum k_r r = 14 N} prod phihat(k_r).
  Dominant term: the "diagonal" where one k_r is large and the rest 0 is NOT allowed
  unless r | 14N.  The leading contributions come from SHORT relations
       sum_r k_r r = 14 N   with small |k_r|.
  Each phihat(k) ~ -sin(pi k/7)/(pi k) = O(1/|k|).  So the n=14N coefficient is a sum
  of products of 1/|k_r| over lattice relations -- this is EXACTLY a SINGULAR SERIES /
  restricted divisor sum.

  Writing m_far = max(Q), the far speeds are 14, 28, ..., 14 max(Q).  The resonance
  N runs over Lat(Q).  The KEY uniformity question:

     | SPEC | = | sum_N chat(14N) conj(Ghat_Q(N)) | <= sum_N |chat(14N)| |Ghat_Q(N)|.

  We compute this bound exactly (truncated) and compare to baseline.  We also compute
  the L2/Cauchy-Schwarz bound (HYP-2861):
     |SPEC| <= sqrt( sum_{N!=0}|chat(14N)|^2 ) * sqrt( sum_{N!=0}|Ghat_Q(N)|^2 )
            <= sqrt(meas(Rsafe) - meas(Rsafe)^2) * sqrt(meas(Qlon) - meas(Qlon)^2)
            (Parseval: sum_{all n}|chat|^2 = meas(Rsafe), minus the n=0 term).
  This L2 bound is UNCONDITIONAL and requires NO arithmetic input -- no BV, no EH.
"""
import sys, math
from fractions import Fraction as F
from math import gcd, pi, sin, sqrt
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, fourier_num_of_arcs,
)

def lat_gcd(S):
    S=[x for x in S if x!=0]
    return reduce(gcd,S) if S else 0

def coeff(arcs, n):
    if n==0: return complex(float(meas(arcs)))
    return fourier_num_of_arcs(arcs,n)/(2j*pi*n)

def Q_lonely_in_t(Q, h=F(1,14)):
    return safe_set([14*m for m in Q], h)

def analyze(R, Q, N=4000, label=""):
    R=sorted(set(int(x) for x in R if x!=0))
    Q=sorted(set(int(x) for x in Q if x!=0))
    Rsafe=safe_set(R); Qlon=Q_lonely_in_t(Q)
    mR=meas(Rsafe); mQ=meas(Qlon); mB=meas(intersect(Rsafe,Qlon))
    baseline=mR*mQ
    SPEC_exact=mB-baseline
    Rprime=(mB/baseline) if baseline>0 else F(1)
    bf=float(baseline)

    # spectral pieces restricted to resonance lattice n=14N
    # L1 (triangle) bound and exact SPEC reconstruction
    L1=0.0; SPECspec=0.0
    sum_chat2=0.0; sum_ghat2=0.0   # Parseval partial (nonzero n)
    # we only need n in 14Z for ghat support, but compute full chat^2/ghat^2 for L2
    for n in range(1,N+1):
        cn=coeff(Rsafe,n); gn=coeff(Qlon,n)
        sum_chat2 += 2.0*abs(cn)**2
        sum_ghat2 += 2.0*abs(gn)**2
        if n%14==0:
            t=2.0*(cn*gn.conjugate()).real
            SPECspec += t
            L1 += 2.0*abs(cn)*abs(gn)

    # EXACT Parseval L2 ceilings (no truncation error): sum_{n!=0}|chat|^2 = mR - mR^2
    var_R=float(mR-mR*mR)
    var_Q=float(mQ-mQ*mQ)
    L2_bound=sqrt(max(var_R,0))*sqrt(max(var_Q,0))

    # L1 (triangle) using exact Parseval is NOT available; report truncated L1.
    print("="*92)
    print(f"  ({label})  R={R}  Q={Q}  far=14Q={[14*q for q in Q]}")
    print("="*92)
    print(f"  baseline={float(baseline):.8f}  SPEC_exact={float(SPEC_exact):+.8f}  R'={float(Rprime):.6f}")
    print(f"  -- DISCREPANCY / LEVEL-OF-DISTRIBUTION BOUNDS on |SPEC| --")
    print(f"  |SPEC| (actual)                         = {abs(float(SPEC_exact)):.8f}  = {abs(float(SPEC_exact))/bf:.5f} * baseline")
    print(f"  L1 triangle  sum|chat||ghat| (|n|<={N})  = {L1:.8f}  = {L1/bf:.5f} * baseline   [needs per-term arith]")
    print(f"  L2 Cauchy-Schwarz (EXACT Parseval)      = {L2_bound:.8f}  = {L2_bound/bf:.5f} * baseline   [UNCONDITIONAL, no BV/EH]")
    print(f"     var_R=meas(Rsafe)(1-meas)= {var_R:.6f}   var_Q=meas(Qlon)(1-meas)= {var_Q:.6f}")
    # The decisive ratio: R' >= 1 - L2_bound/baseline.  If that is > 0, L2 alone (NO arithmetic) closes it.
    floor_L2 = 1.0 - L2_bound/bf
    floor_L1 = 1.0 - L1/bf
    print(f"  ==> R' >= 1 - L2/baseline = {floor_L2:+.6f}   (L2-only floor; positive iff L2<baseline)")
    print(f"  ==> R' >= 1 - L1/baseline = {floor_L1:+.6f}   (L1 floor; tighter but truncated)")
    return dict(R=R,Q=Q,baseline=baseline,SPEC=SPEC_exact,Rprime=Rprime,
                L1=L1,L2=L2_bound,floor_L2=floor_L2,floor_L1=floor_L1,
                var_R=var_R,var_Q=var_Q,bf=bf,label=label)

def main():
    print("#"*92)
    print("# TOOL 3 : LEVEL-OF-DISTRIBUTION BOUNDS ON SPEC  (L1 vs L2 vs actual)")
    print("#   Is the uniform |SPEC| < baseline bound elementary (L2/Parseval) or does it need BV/EH?")
    print("#"*92)
    cases = [
        ([1,2,3,4,5,6,7,8,9,10,11,13],[1],          "R={1..13}\\{12}, Q={1} r=1"),
        ([1,2,3,4,5,6,7,8,9,10,11,12],[1,2],        "R={1..12}, Q={1,2} r=2 (worst R'=0.70)"),
        ([1,2,3,4,5,6,7,8,9,10,11],  [1,2,3],      "R={1..11}, Q={1,2,3} r=3"),
        ([1,2,3,4,5,6,7,8,9,10],     [1,2,3,4],    "R={1..10}, Q={1,2,3,4} r=4"),
        ([1,2,3,4,5,6,7,8,9],        [1,2,3,4,5],  "R={1..9}, Q={1..5} r=5"),
        ([1,2,3,4,5,6,7,8],          [1,2,3,4,5,6],"R={1..8}, Q={1..6} r=6 max"),
        ([2,4,6,8,10,12],            [1,3,5],      "R gcd=2, Q={1,3,5}  (worst SPEC=-0.029)"),
        ([3,6,9,12],                 [1,2,5],      "R gcd=3, Q={1,2,5}"),
        ([5,9,11,13],                [1,2,3],      "R coprime-ish, Q={1,2,3}"),
    ]
    res=[analyze(R,Q,N=4000,label=lab) for R,Q,lab in cases]
    print("="*92)
    print("SUMMARY: actual R' vs the unconditional L2 floor (1 - L2/baseline)")
    print("="*92)
    print(f"{'label':<42}{'R-prime':>10}{'|SPEC|/base':>13}{'L2/base':>10}{'floor_L2':>11}{'floor_L1':>11}")
    for r in res:
        print(f"{r['label'][:42]:<42}{float(r['Rprime']):>10.5f}"
              f"{abs(float(r['SPEC']))/r['bf']:>13.5f}{r['L2']/r['bf']:>10.5f}"
              f"{r['floor_L2']:>11.5f}{r['floor_L1']:>11.5f}")
    print()
    minRp=min(float(r['Rprime']) for r in res)
    minfL2=min(r['floor_L2'] for r in res)
    minfL1=min(r['floor_L1'] for r in res)
    print(f"  min actual R'            = {minRp:.5f}")
    print(f"  min L2-only floor        = {minfL2:.5f}   (if >0: Parseval alone gives R'>=c, NO BV/EH)")
    print(f"  min L1 floor (truncated) = {minfL1:.5f}")
    print()
    print("VERDICT (printed in detail by the analysis):")
    print("  The L2 Cauchy-Schwarz / Parseval bound is UNCONDITIONAL: it uses only")
    print("  sum|chat|^2 = meas(Rsafe)-meas(Rsafe)^2 and the same for Q.  It needs NO")
    print("  equidistribution input.  Whether it gives a POSITIVE floor is the question.")
    print("DONE.")

if __name__=="__main__":
    main()
