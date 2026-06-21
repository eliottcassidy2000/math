#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C (kps-2026-06-21-S?-wf9): the ADJACENT-DOUBLET = SINGLE-FAR + DRIFT decomposition.

THE ONE OPEN LEG of LRC(14) (OPEN-Q-108 leg C, genuine-wide): prove
    p0(genuine-wide) < cap_k  with margin >= 0.16.
The binding genuine-wide maximizer is the ADJACENT DOUBLET
    E_M = consec_{k-1} U {M, M+1}    (opus HYP-2794, exact k=8..13).
THM-563 closes SINGLE-far (consec_{k-1} U {M}) because M*Delta_M is EXACTLY periodic
in M (period P = 7*lcm(consec_{k-1})). The doublet is NOT exactly periodic
(2nd-diff ~0.01-0.03, decaying) -- it is ALMOST periodic.

CREATIVE REDUCTION (Thread C): write
    p0(E_M) = p0(S_M) + D(M),    S_M = consec_{k-1} U {M}  (single-far, THM-563-closed)
where D(M) = the DRIFT coupling = effect of adding the locked partner runner M+1.
Equivalently D(M) is itself a "one-far addition" but relative to the base S_M which
already contains M. The joint phase of the pair {M,M+1} is (frac(Mx), frac(Mx + x)):
a SINGLE fast runner at M plus a SLOW drift x. So D(M) measures the slow drift.

This engine computes EXACTLY (rationals):
  (1) p0(E_M), p0(S_M), and D(M) = p0(E_M) - p0(S_M)  over M in [15, Mmax].
  (2) the single-far term: confirm M*(p0(S_M) - Phi_1) is EXACTLY periodic (THM-563),
      report its period-max  PM1 = sup_M M*(p0(S_M)-Phi_1).
  (3) the DRIFT term D(M): is M*D(M) bounded? is D(M) -> 0 (and how fast: M*D, M^2*D)?
      report sup_M |M*D(M)| and sup_M |M^2*D(M)| over the window.
  (4) the SIGN of D(M): adding a runner can only COVER MORE sectors, so p0 (lonely
      measure / all-6-covered measure) ... determine the sign empirically.
  (5) ASSEMBLE: p0(E_M) = Phi_1 + (periodic_1(M) + M*D(M))/M. If M*D(M) is bounded by
      C_D, then p0(E_M) <= Phi_1 + (PM1 + C_D)/M <= Phi_1 + (PM1+C_D)/15.
      Compare to cap_k. (Phi_1 here is the SINGLE-far plateau, not the doublet plateau.)

Exact rationals. p0_fast matches repo p0 (cross-checked elsewhere).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from math import gcd, lcm
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL, MARGIN


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def consec(k):
    return tuple(range(k - 1))  # base speeds 0..k-2 (size k-1)


def single_far(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M])))


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


@lru_cache(maxsize=None)
def p0c(E):
    return p0_fast(E)


def single_far_plateau(k, w0, P):
    """Phi_1 from THM-563: S_w = w*p0(S_w) is w*Phi_1 + periodic(w). Phi_1 = slope over period.
    Returns (Phi_1, period_max_PM1, periodic_ok)."""
    S = {}
    for w in range(w0, w0 + 2 * P + 1):
        S[w] = w * p0c(single_far(k, w))
    slopes = set((S[w + P] - S[w]) for w in range(w0, w0 + P))
    ok = (len(slopes) == 1)
    Phi1 = next(iter(slopes)) / P if ok else None
    PM1 = None
    if Phi1 is not None:
        PM1 = max(S[w] - w * Phi1 for w in range(w0, w0 + P))
    return Phi1, PM1, ok


def analyze(k, Mmax, w0=15):
    base = list(range(1, k - 1))  # nonzero base speeds
    P = 7 * lcm_list(base)
    # single-far plateau + period-max (THM-563)
    Phi1, PM1, ok1 = single_far_plateau(k, w0, P)
    # scan
    rows = []
    sup_p0E = F(-1); argE = None
    supMD = F(-1); argMD = None
    supM2D = F(-1); argM2D = None
    Dsigns = set()
    for M in range(w0, Mmax + 1):
        pE = p0c(doublet(k, M))
        pS = p0c(single_far(k, M))
        D = pE - pS
        rows.append((M, pE, pS, D))
        if pE > sup_p0E:
            sup_p0E, argE = pE, M
        md = abs(M * D)
        if md > supMD:
            supMD, argMD = md, M
        m2d = abs(M * M * D)
        if m2d > supM2D:
            supM2D, argM2D = m2d, M
        if D > 0:
            Dsigns.add('+')
        elif D < 0:
            Dsigns.add('-')
        else:
            Dsigns.add('0')
    cap = CAP[k]
    return dict(k=k, P=P, Phi1=Phi1, PM1=PM1, ok1=ok1, cap=cap,
                sup_p0E=sup_p0E, argE=argE, supMD=supMD, argMD=argMD,
                supM2D=supM2D, argM2D=argM2D, Dsigns=Dsigns, rows=rows)


def main():
    print("=" * 88)
    print("THREAD C: ADJACENT-DOUBLET = SINGLE-FAR + DRIFT  (kps-wf9)")
    print("E_M = consec_{k-1} U {M,M+1} = S_M U {M+1},  S_M = consec_{k-1} U {M} (single-far)")
    print("D(M) = p0(E_M) - p0(S_M) = the DRIFT coupling from the locked partner M+1")
    print("=" * 88)
    Mmax = 700
    for k in range(8, 13):
        r = analyze(k, Mmax)
        cap = r['cap']; Phi1 = r['Phi1']; PM1 = r['PM1']
        print(f"\nk={k}  P=7*lcm={r['P']}  cap={cap}={float(cap):.6f}")
        print(f"  SINGLE-FAR (THM-563): periodic={r['ok1']}  Phi_1={Phi1}={float(Phi1):.6f}  "
              f"period-max PM1={PM1}={float(PM1):.5f}")
        print(f"  DRIFT sign(s) of D(M)={sorted(r['Dsigns'])}")
        print(f"  sup_M |M*D(M)|   = {float(r['supMD']):.5f} at M={r['argMD']}   (D(M) bounded*M?)")
        print(f"  sup_M |M^2*D(M)| = {float(r['supM2D']):.5f} at M={r['argM2D']}  (D(M) ~ 1/M^2?)")
        # ASSEMBLY: p0(E_M) <= Phi_1 + (PM1 + sup|M*D|)/M  for M>=15
        C = PM1 + r['supMD']  # but PM1 may be negative-ish; use signed worst:
        # the actual bound we want: p0(E_M)-Phi_1 = (S_M_periodic + M*D)/M <= (PM1 + sup M*D)/M
        bound15 = Phi1 + C / 15
        print(f"  ASSEMBLE: p0(E_M) <= Phi_1 + (PM1 + sup|M*D|)/M ; at M=15: "
              f"{float(bound15):.6f}  vs cap {float(cap):.6f}  margin {float(cap-bound15):+.6f}")
        print(f"  (empirical sup p0(E_M) over [15,{Mmax}] = {float(r['sup_p0E']):.6f} at M={r['argE']})")
    print("\n" + "=" * 88)
    print("KEY QUESTIONS: (a) is sup|M*D| bounded (drift is a bounded coupling)?  ")
    print("(b) does the assembled bound clear cap with margin >= 0.16?")
    print("If M*D bounded -> doublet error = single-far(periodic, THM-563) + bounded drift -> CLOSES.")


if __name__ == "__main__":
    main()
