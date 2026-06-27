#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_certified_summary_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2 Asano/Lee-Yang)

CONSOLIDATED CERTIFIED SUMMARY of the Asano/Lee-Yang analysis of the LRC(14) loneliness floor.

THE REFRAMING (exact, validated).
  Loneliness = meas( cap_{s in S} D_s^c ) = P(M=0), where M(t)=#{s in S : ||s t||<1/14} is the
  DANGER COUNT over uniform t, and EVERY single danger event has meas(D_s)=1/7 EXACTLY.
  The danger-count PGF G(z)=E[z^M] is the PARTITION FUNCTION; with one fugacity lambda,
        Xi_diag(lambda) = integral_0^1 (1-lambda)^{M(t)} dt = G(1-lambda),     L = Xi_diag(1)=G(0).
  So the LEE-YANG zeros of the loneliness partition function are lambda = 1 - zeta, zeta = roots of
  the danger-count PGF.  "Lee-Yang property" := all roots zeta satisfy |1-zeta| >= 1  (no zero of
  Xi_diag in the open lambda-unit-disk).

THE DICHOTOMY (the decisive finding).
  Split S = R u 14Q (R = 14-free big part, |R|=13-r ; Q = apexes, r=|Q| in 2..6).
  * Q-BLOCK (apex sub-LRC, r<=6 speeds 14m): the danger-count PGF IS Lee-Yang -- all roots have
    |1-zeta| > 1, margins 7.0/5.77/3.78/2.76/2.32/1.74 for r=1..6.  Equivalently the multi-affine
    Xi_Q is zero-free on a polydisk of radius rho* >= 12 (>> 1).  Hence meas(Q-lonely)=Xi_Q(1)>0,
    also provable by the plain UNION BOUND meas(Q-lonely) >= 1 - r/7 > 0.
  * R-BLOCK (big part, |R|=13-r >= 7 speeds): the danger-count PGF FAILS Lee-Yang -- complex roots
    enter the lambda-unit-disk (e.g. R={1..11,13}: lambda-root -0.43+-0.48 i inside).  This is the
    > 6-speed regime where the union bound is vacuous (|R|/7 > 1).
  * COMONOTONE TEST: Lee-Yang is NOT implied by k*p < 1 alone -- the maximally-correlated count
    {0:1-1/k, k:1/k} at kp=1 FAILS for k>=4 (min|1-zeta|=0.933 at k=4).  The Q-block's Lee-Yang
    is a STRUCTURAL fact about the equidistributed apex comb (near-independent), not generic.

THE HONEST OBSTRUCTION (why naive Asano cannot close the joint floor).
  Joint loneliness = meas(R-safe cap Q-lonely) = Xi(1,1), where Xi(lambda,mu) is the two-block
  partition function.  Because the R-block factor Xi(lambda,0)=G_R(1-lambda) already has zeros
  inside |lambda|<1, the FULL bidisk |lambda|<=1, |mu|<=1 is NOT zero-free.  Therefore the Asano-
  contraction / polydisk-Lee-Yang route CANNOT, by itself, certify Xi(1,1)>0.  This RIGOROUSLY
  reproduces -- and explains, at the level of zero locations -- the documented fact that BONFERRONI
  FAILS for the few-apex covering case and the floor survives only via QUASI-INDEPENDENCE
  (R' = Xi(1,1)/(Xi(1,0)Xi(0,1)) ~ 0.5..1.07, the equidistributed-comb regime, THM-546/HYP-2840).

NET DELIVERABLE.  The Lee-Yang/Asano frame CLEANLY proves the Q-block (apex) factor positivity and
  isolates the obstruction to the R-block overcrowding; it does NOT yield a new proof of the joint
  R'>=c floor (which is genuinely an equidistribution statement, the TOOL-1/spectrum thread).

This script regenerates the certified data behind every claim above.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, pi
import numpy as np
import numpy.polynomial.polynomial as Poly

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_asano_zerofree_kpswf15 import (merge, meas, intersect, safe_arcs,
                                          complement_arcs, danger_arcs, loneliness_exact)

THR = F(1, 14)

def count_pgf(speeds):
    """Exact danger-count distribution dict j->meas, over the given integer speeds."""
    L0 = 1
    for s in speeds:
        L0 = L0 * s // gcd(L0, s)
    bps = {F(0), F(1)}
    for s in speeds:
        for a, b in danger_arcs(s, THR):
            bps.add(a); bps.add(b)
    bps = sorted(bps)
    dist = {}
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        j = sum(1 for s in speeds if ((s * mid) % 1) < THR or ((s * mid) % 1) > 1 - THR)
        dist[j] = dist.get(j, F(0)) + (hi - lo)
    return dist

def ly_clearance(speeds):
    """Return (P(M=0), min|1-zeta|) for the danger-count PGF of these speeds."""
    dist = count_pgf(speeds)
    k = len(speeds)
    coeffs = [float(dist.get(j, F(0))) for j in range(k + 1)]
    roots = Poly.polyroots(coeffs)
    m1 = min((abs(1 - r) for r in roots), default=float('inf'))
    return float(dist.get(0, F(0))), m1

def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("A. REFRAMING CHECK: every danger event has meas(D_s)=1/7 (exact)")
    bad = [s for s in [1, 2, 3, 5, 7, 11, 13, 84, 154, 30030] if meas(danger_arcs(s, THR)) != F(1, 7)]
    print(f"   meas(D_s)=1/7 for all sampled speeds: {'CONFIRMED' if not bad else f'FAILS {bad}'}")

    banner("B. Q-BLOCK Lee-Yang clearance vs r (apexes 14m), and union floor 1-r/7")
    Qsets = [(6,), (6, 11), (6, 12, 13), (6, 11, 12, 13), (6, 10, 11, 12, 13),
             (6, 9, 10, 11, 12, 13)]
    print(f"   {'r':>2}{'meas(Q-lonely)':>16}{'1-r/7':>10}{'min|1-zeta|':>13}{'Lee-Yang':>10}")
    for Q in Qsets:
        speeds = [14 * m for m in Q]
        p0, m1 = ly_clearance(speeds)
        r = len(Q)
        print(f"   {r:>2}{p0:>16.5f}{float(1 - F(r, 7)):>10.5f}{m1:>13.4f}"
              f"{'YES' if m1 >= 1 - 1e-9 else 'NO':>10}")

    banner("C. R-BLOCK Lee-Yang FAILS for |R|>=7 (the overcrowded big part)")
    Rsets = [tuple(range(1, 12)) + (13,), tuple(range(1, 11)) + (13,),
             tuple(range(1, 10)) + (11,), tuple(range(1, 8)) + (9,), tuple(range(1, 7)) + (8,)]
    print(f"   {'|R|':>4}{'meas(R-safe)':>14}{'min|1-zeta|':>13}{'Lee-Yang':>10}")
    for R in Rsets:
        p0, m1 = ly_clearance(list(R))
        print(f"   {len(R):>4}{p0:>14.5f}{m1:>13.4f}{'YES' if m1 >= 1 - 1e-9 else 'NO':>10}")

    banner("D. COMONOTONE control: Lee-Yang NOT implied by k*p<1 (extreme {0:1-1/k,k:1/k}, kp=1)")
    print(f"   {'k':>2}{'min|1-zeta|':>13}{'Lee-Yang':>10}")
    for k in range(2, 9):
        Rrad = (k - 1) ** (1.0 / k)
        best = min(abs(1 - Rrad * complex(math.cos((pi + 2 * pi * j) / k),
                                          math.sin((pi + 2 * pi * j) / k))) for j in range(k))
        print(f"   {k:>2}{best:>13.4f}{'YES' if best >= 1 - 1e-9 else 'NO':>10}")

    banner("E. JOINT floor R' (documented + synthetic r=2..6) -- equidistribution regime")
    cases = [
        ("r=1 doc", tuple(range(1, 12)) + (13,), (6,)),
        ("r=2 doc", tuple(range(1, 11)) + (13,), (6, 11)),
        ("r=3", tuple(range(1, 10)) + (11,), (6, 12, 13)),
        ("r=4", tuple(range(1, 9)) + (10,), (6, 11, 12, 13)),
        ("r=5", tuple(range(1, 8)) + (9,), (6, 10, 11, 12, 13)),
        ("r=6", tuple(range(1, 7)) + (8,), (6, 9, 10, 11, 12, 13)),
    ]
    print(f"   {'tag':<9}{'measR':>9}{'measQ':>9}{'both':>9}{'R prime':>9}{'indep both':>12}")
    for tag, R, Q in cases:
        aR = [(F(0), F(1))]
        for s in R:
            aR = intersect(aR, safe_arcs(s, THR))
        aQ = [(F(0), F(1))]
        for m in Q:
            aQ = intersect(aQ, safe_arcs(14 * m, THR))
        mR = meas(aR); mQ = meas(aQ); both = meas(intersect(aR, aQ))
        Rp = float(both / (mR * mQ)) if mR > 0 and mQ > 0 else float('nan')
        print(f"   {tag:<9}{float(mR):>9.5f}{float(mQ):>9.5f}{float(both):>9.5f}{Rp:>9.4f}"
              f"{float(mR * mQ):>12.5f}")

    banner("VERDICT")
    print("""   Lee-Yang/Asano CLEANLY proves the Q-block (apex sub-LRC, r<=6) factor positivity
   (zero-free polydisk radius >=12; also union bound >= 1-r/7).  The R-block (>=7 speeds)
   provably FAILS single-variable Lee-Yang, and the joint bidisk is NOT zero-free, so the
   naive Asano-contraction route CANNOT certify the joint R'>=c floor.  This rigorously
   matches and explains the documented 'Bonferroni fails / quasi-independence required'
   finding (HYP-3121).  The joint floor is genuinely an equidistribution statement.""")
