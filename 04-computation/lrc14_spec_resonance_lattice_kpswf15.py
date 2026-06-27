#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spec_resonance_lattice_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 3 / EH-equidistribution)

THE RESONANCE-LATTICE DECOMPOSITION of SPEC for the r=2..6 MULTI-FAR family.

SETUP (post-substitution u=14t, HYP-2968 / HYP-3121 frame):
  S = R u 14Q,  R = 14-free small part,  14Q = r in {2..6} multiples of 14.
  R-safe   = { t : ||r t|| >= 1/14  for all r in R }        (= G_R in t-circle)
  Q-lonely = { u : ||m u|| >= 1/14  for all m in Q }, u=14t   (= G_Q in u-circle)
  loneliness = meas( R-safe  cap  Q-lonely ).

  chat(n) = Fourier coeffs of 1_{R-safe}  (t-circle).
  ghat(n) = Fourier coeffs of 1_{Q-lonely}(14 t) = 1_{Q-lonely} pulled back, freqs are 14*Z.
  R' = meas(both)/(meas(R-safe)*meas(Q-lonely)) = 1 + SPEC/baseline,
  SPEC = sum_{n!=0} chat(n) conj(ghat(n)).

KEY SUPPORT FACTS (this is the EH/equidistribution content):
  (A) chat(n) != 0 only for n in Lat(R) = additive group <r : r in R>.  Each factor
      1_{||r t|| in [1/14,13/14]} has ONLY frequencies that are multiples of r.  So the
      product's spectrum lives on  sum_r m_r r  -- the relation lattice of R (HYP-2606).
  (B) ghat(n) != 0 only for n in 14*Lat(Q) = 14*<m : m in Q>.  Pullback by u=14t multiplies
      all Q-lonely frequencies by 14.  So Q-lonely's t-spectrum lives on 14*Lat(Q).
  (C) RESONANCE CONDITION:  a frequency n contributes to SPEC iff
            n in Lat(R)  AND  n in 14*Lat(Q).
      Since 14 | n for the ghat side, ONLY MULTIPLES OF 14 resonate.  And n must also be
      expressible in the R-lattice.  This is the RESONANCE LATTICE  L = Lat(R) cap 14 Lat(Q).

  This script computes, for representative multi-far rows:
    - the exact SPEC (real-space, Fraction) and its spectral reconstruction
    - the per-frequency term(n), restricted to the resonance lattice n in 14 Z
    - the residue-0 trunk (here: n in 98 Z, since 14*7) vs the rest
    - the DISCREPANCY interpretation: term(n) = chat(n) * conj(ghat(n)) is a product of two
      Weyl/equidistribution discrepancies; |chat(n)|,|ghat(n)| <= V/(2 pi |n|).

OUTPUT is the resonance decomposition + the level-of-distribution sum structure that
BV/EH would be asked to bound.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from math import gcd, pi
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, danger_arcs, safe_set, fourier_num_of_arcs,
)

def lat_gcd(S):
    """gcd of S = generator of the additive lattice <S> in Z."""
    S = [x for x in S if x != 0]
    if not S: return 0
    return reduce(gcd, S)

def coeff(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)))
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def safe_set_h(P, h=F(1,14)):
    return safe_set(P, h)

def Q_lonely_in_t(Q, h=F(1,14)):
    """Q-lonely in u-circle pulled back to t-circle via u=14t:
       { t : ||14 m t|| >= 1/14  for all m in Q } = G_{14Q} in t.
       This is exactly safe_set(14*Q)."""
    return safe_set([14*m for m in Q], h)

def analyze_multifar(R, Q, N=6000, label=""):
    """R = 14-free small part (speeds), Q = multipliers so far speeds are 14*Q.
       Computes SPEC for S = R u 14Q via the t-circle factorization."""
    R = sorted(set(int(x) for x in R if x != 0))
    Q = sorted(set(int(x) for x in Q if x != 0))
    Rsafe = safe_set(R)                 # chat support: Lat(R)
    Qlon  = Q_lonely_in_t(Q)            # ghat support: 14 Lat(Q)
    both  = intersect(Rsafe, Qlon)
    mR = meas(Rsafe); mQ = meas(Qlon); mB = meas(both)
    baseline = mR * mQ
    SPEC_exact = mB - baseline
    Rprime = (mB / baseline) if baseline > 0 else F(1)

    latR = lat_gcd(R)                   # chat lives on latR * Z
    latQ = lat_gcd(Q)                   # Q-lonely lives on 14*latQ * Z
    res_step = 14 * latQ               # ghat nonzero only on multiples of this
    # resonance lattice generator: lcm(latR, res_step) ... but chat is on latR*Z and ghat on res_step*Z
    # both nonzero => n in (latR Z) cap (res_step Z) = lcm(latR,res_step) Z
    L = (latR * res_step) // gcd(latR, res_step) if latR and res_step else 0

    # Per-frequency terms restricted to multiples of 14 (ghat support).
    terms = []
    SPEC_spec = 0.0
    for n in range(1, N+1):
        cn = coeff(Rsafe, n)
        gn = coeff(Qlon, n)
        t = 2.0*(cn*gn.conjugate()).real
        SPEC_spec += t
        if abs(t) > 1e-12:
            terms.append((n, t, abs(cn), abs(gn)))
    # trunk: n in 98 Z (=14*7) -> the apex-prime trunk after the 14-pullback
    trunk = sum(t for (n,t,_,_) in terms if n % 98 == 0)
    nonzero_shell = sum(t for (n,t,_,_) in terms if n % 98 != 0)
    bf = float(baseline) if baseline != 0 else float('nan')

    print("="*94)
    print(f"  ({label})  R(14-free)={R}  Q(multipliers, far=14Q)={Q}")
    print(f"     full set S = {sorted(R + [14*q for q in Q])}")
    print("="*94)
    print(f"  meas(R-safe)            = {mR} = {float(mR):.8f}")
    print(f"  meas(Q-lonely)          = {mQ} = {float(mQ):.8f}")
    print(f"  baseline = product      = {baseline} = {float(baseline):.8f}")
    print(f"  meas(both)              = {mB} = {float(mB):.8f}")
    print(f"  SPEC (exact)            = {SPEC_exact} = {float(SPEC_exact):+.8f}")
    print(f"  R' = both/baseline      = {Rprime} = {float(Rprime):.8f}")
    print(f"  SPEC (spectral |n|<={N}) = {SPEC_spec:+.8f}  |err|={abs(SPEC_spec-float(SPEC_exact)):.2e}")
    print(f"  -- LATTICE STRUCTURE --")
    print(f"     gcd(R) = {latR}  =>  chat supported on {latR} Z")
    print(f"     gcd(Q) = {latQ}  =>  Q-lonely (in t) supported on {res_step} Z  (=14*gcd(Q))")
    print(f"     RESONANCE LATTICE L = lcm = {L} Z   (only these n give nonzero term)")
    print(f"     trunk (n in 98 Z, apex-7 after pullback)  = {trunk:+.8f}  ({trunk/bf:+.5f} * baseline)")
    print(f"     nonzero shell (n in 14Z \\ 98Z)            = {nonzero_shell:+.8f}  ({nonzero_shell/bf:+.5f} * baseline)")
    # top resonant frequencies
    top = sorted(terms, key=lambda r: -abs(r[1]))[:12]
    print(f"  -- TOP RESONANT FREQUENCIES (n, term, |chat|, |ghat|, n/14) --")
    for (n,t,acn,agn) in top:
        print(f"     n={n:5d} (={n//14 if n%14==0 else n/14:>5}*14 mod7={n//14%7 if n%14==0 else '-'})  term={t:+.7f}  |chat|={acn:.6f} |ghat|={agn:.6f}")
    return dict(R=R, Q=Q, mR=mR, mQ=mQ, baseline=baseline, SPEC=SPEC_exact, Rprime=Rprime,
                latR=latR, latQ=latQ, L=L, trunk=trunk, shell=nonzero_shell, terms=terms)

def main():
    print("#"*94)
    print("# TOOL 3 : RESONANCE-LATTICE DECOMPOSITION OF SPEC  (multi-far r=2..6)")
    print("#   chat on Lat(R), ghat on 14 Lat(Q)  =>  SPEC supported on resonance lattice")
    print("#   Question: does BV/EH bound SPEC uniformly over the multi-far family?")
    print("#"*94)

    # Multi-far rows: R = 14-free small part (drop one residue from {1..13}), Q = {2..r} far multipliers.
    # The "few-apex" branch is 1<=|Q|<=6.  We test r=2..6 (|Q| from 1 to 5) with various R.
    # GENUINE few-apex covering rows: R must be 14-FREE and leave R-safe > 0 (so it is NOT
    # already a covering set on its own); the far multiples 14Q then help cover.  The classic
    # tight-perturbation family drops ONE residue from {1..13} (R-safe gets positive measure)
    # and adds r=|Q| far speeds.  |Q| ranges 1..6 (the few-apex branch).
    cases = [
        # (R, Q, label)  -- R = {1..13} minus one residue, plus far multipliers Q
        ([1,2,3,4,5,6,7,8,9,10,11,13],[1],          "R={1..13}\\{12} (12->14), Q={1} (r=1 far)"),
        ([1,2,3,4,5,6,7,8,9,10,11,12],[1,2],        "R={1..12}, Q={1,2} (far 14,28; drop 13)"),
        ([1,2,3,4,5,6,7,8,9,10,11],  [1,2,3],      "R={1..11}, Q={1,2,3} (r=3 far)"),
        ([1,2,3,4,5,6,7,8,9,10],     [1,2,3,4],    "R={1..10}, Q={1,2,3,4} (r=4 far)"),
        ([1,2,3,4,5,6,7,8,9],        [1,2,3,4,5],  "R={1..9}, Q={1..5} (r=5 far)"),
        ([1,2,3,4,5,6,7,8],          [1,2,3,4,5,6],"R={1..8}, Q={1..6} (r=6 MULTI-FAR, max)"),
        # genuinely separated far multipliers (not consecutive) -- harder resonance:
        ([1,2,3,4,5,6,7,8,9],        [1,3,5],      "R={1..9}, Q={1,3,5} (spread far)"),
        ([1,2,3,4,5,6,7,8,9],        [2,3,5],      "R={1..9}, Q={2,3,5} (no q=1, gcd(Q)=1)"),
        # divisor-loaded / gcd>1 R to probe RESONANCE-LATTICE thinning:
        ([2,4,6,8,10,12],            [1,3,5],      "R even gcd=2 {2..12}, Q={1,3,5}"),
        ([3,6,9,12],                 [1,2,5],      "R gcd=3 {3,6,9,12}, Q={1,2,5}"),
        ([5,9,11,13],                [1,2,3],      "R coprime-ish {5,9,11,13}, Q={1,2,3}"),
    ]
    res = []
    for R,Q,lab in cases:
        d = analyze_multifar(R,Q,N=3000,label=lab)
        d['label'] = lab
        res.append(d)
        print()
    print("="*94)
    print("SUMMARY: R', SPEC, resonance-lattice generator L, trunk vs shell")
    print("="*94)
    print(f"{'label':<46}{'R-prime':>10}{'SPEC':>12}{'L':>8}{'trunk':>11}{'shell':>11}")
    for r in res:
        lab = (r['label'] if 'label' in r else '')[:46]
        print(f"{lab:<46}{float(r['Rprime']):>10.5f}{float(r['SPEC']):>12.6f}{r['L']:>8}{r['trunk']:>11.6f}{r['shell']:>11.6f}")
    Rs=[float(r['Rprime']) for r in res]
    print(f"\n  R' range over multi-far family: [{min(Rs):.5f}, {max(Rs):.5f}]")
    print("DONE.")

if __name__ == "__main__":
    # attach labels
    import types
    main()
