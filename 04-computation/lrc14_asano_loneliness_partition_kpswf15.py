#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_loneliness_partition_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2: Asano/Lee-Yang)

GOAL: frame the LRC(14) loneliness measure as a PARTITION FUNCTION and probe whether it
has a Lee-Yang / Asano-contraction zero-free structure that PROVES the floor positivity.

SETUP (the r=2..6 multi-far covering core, kps HYP-3121):
  A covering set is S = R u 14Q.  After the substitution u=14t the multiples 14Q become a
  sub-LRC for Q; the 14-free "small part" R has 13-r >= 7 speeds.
  Loneliness of the lifted problem = meas(R-safe cap Q-lonely), where
      R-safe   = { t : ||r t|| >= 1/14 for all r in R }
      Q-lonely = { t : ||14 m t|| >= 1/14 } (the sub-LRC for Q)
  We focus on the R-safe FLOOR object: each safe-arc indicator
      phi_s(t) = 1[ ||s t|| >= 1/14 ]   = psi(s t),   psi(u)=1[1/14 <= {u} <= 13/14]
  is an arc of measure 6/7 on the circle.  The R-safe loneliness measure is
      L(R) = integral_0^1 prod_{s in R} phi_s(t) dt.
  Positivity L(R) > 0 (and the sharper floor) is a ZERO-FREENESS statement.

THE PARTITION-FUNCTION FRAME.
  Fix a common period.  Write z = e(t) = exp(2 pi i t).  Then phi_s(t) = psi(s t) where psi is
  a fixed circle function.  The integral L = [n=0 Fourier coefficient of prod_s phi_s] is the
  partition function of a 1-D "gas" of phase-constraints (one constraint per speed s in R).

  RIESZ-PRODUCT VIEW (HYP-2606 singular series): psi has Fourier series
      psi(u) = sum_m c_m e(m u),  c_0 = 6/7,  c_m = (e(-m/14) - e(-13 m/14))/(2 pi i m)  (m!=0)
            = sin(6 pi m /7)/(pi m)             [real, even in m]
  so phi_s(t) = sum_m c_m e(m s t) has frequencies only in s Z.  Then
      L = integral prod_s ( sum_{m_s} c_{m_s} e(m_s s t) ) dt
        = sum over (m_s)_{s in R} with sum_s m_s s = 0  of  prod_s c_{m_s}.
  This is EXACTLY a singular-series / lattice sum over the relation lattice {n : sum n_s s = 0}.

THE MULTI-AFFINE / LEE-YANG QUESTION.
  Asano contraction needs a MULTI-AFFINE polynomial (each variable degree <= 1), zero-free on a
  polydisk, and combines two variables preserving zero-freeness.  We test several encodings:

  (A) FINITE 14N-MODEL.  Replace the integral by the average over the 14*L_period-th roots of
      unity (exact rational arithmetic).  psi becomes a 0/1 vector on residues; L is a finite
      average.  We build the multi-affine generating polynomial Z(z_1,...,z_k) whose constant
      term (in the relation directions) is L, then locate its zeros.

  (B) PER-SPEED VARIABLES z_s.  Treat each speed's phase as an independent circle variable and
      build the multilinear "subset" extension; check Asano contractibility.

VALIDATION: every L is computed two independent ways:
  (i)  EXACT arc-sweep (Fraction):   L = meas( intersect_s {||s t||>=1/14} ).
  (ii) Riesz lattice sum / finite-root average  ->  must match (i).

This is route selection + a positivity engine, NOT yet a closed LRC proof.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from math import gcd, pi, sin
from functools import reduce

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

# ======================================================================== exact arcs
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out

def meas(arcs):
    return sum((b - a for a, b in arcs), F(0))

def intersect(A, B):
    out = []; i = j = 0
    A = merge(A); B = merge(B)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: out.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return out

def safe_arcs(s, thr=F(1, 14)):
    """{ t in [0,1) : ||s t|| >= thr }.  s positive integer."""
    # ||s t|| >= thr  <=>  frac(s t) in [thr, 1-thr].  Pull back: t in union over j of
    # [(j+thr)/s, (j+1-thr)/s], j=0..s-1.
    out = []
    for j in range(s):
        lo = F(j) / s + thr / s
        hi = F(j + 1) / s - thr / s
        if lo < hi:
            out.append((lo, hi))
    return merge(out)

def loneliness_exact(R, thr=F(1, 14)):
    """L(R) = meas( intersect_{s in R} {||s t|| >= thr} )  as an exact Fraction."""
    arcs = [(F(0), F(1))]
    for s in R:
        arcs = intersect(arcs, safe_arcs(s, thr))
    return meas(arcs), arcs

# ======================================================================== Riesz / psi Fourier
def psi_coeff(m):
    """Fourier coeff c_m of psi(u)=1[1/14<={u}<=13/14].  Real.  c_0=6/7."""
    if m == 0:
        return 6.0 / 7.0
    # c_m = integral_{1/14}^{13/14} e(-m u) du = (e(-m/14)-e(-13m/14))/(2 pi i m)
    # = sin(2 pi m * 6/14)/(pi m)  -- real, even.  6/14 = 3/7.
    return math.sin(2 * pi * m * F(6, 14)) / (pi * m)

def loneliness_riesz(R, M):
    """L via the relation-lattice sum truncated to |m_s|<=M.  Returns float approx."""
    total = 0.0
    ranges = [range(-M, M + 1) for _ in R]
    for combo in itertools.product(*ranges):
        if sum(mc * s for mc, s in zip(combo, R)) == 0:
            term = 1.0
            for mc in combo:
                term *= psi_coeff(mc)
            total += term
    return total

# ======================================================================== finite 14N-root model
def loneliness_finite(R, thr=F(1, 14), refine=1):
    """
    Exact finite-average model.  All breakpoints of prod phi_s lie in (1/(L*) ) Z where
    L* = lcm of the speeds times 14.  Sample at midpoints of the breakpoint grid so the 0/1
    value is exact, and average.  Returns exact Fraction == loneliness_exact (sanity model
    for the 'partition function = finite average' framing).
    """
    L0 = 1
    for s in R:
        L0 = L0 * s // gcd(L0, s)
    D = 14 * L0 * refine
    # breakpoints are j/D; sample midpoints (2j+1)/(2D)
    cnt = 0
    for j in range(D):
        t = F(2 * j + 1, 2 * D)
        ok = all(((s * t) % 1) >= thr and ((s * t) % 1) <= 1 - thr for s in R)
        if ok:
            cnt += 1
    return F(cnt, D)

# ======================================================================== Q-lonely (sub-LRC)
def qlonely_arcs(Q, thr=F(1, 14)):
    """
    Q-lonely = { t : ||14 m t|| >= 1/14 for all m in Q }.  Each m contributes safe_arcs(14*m).
    """
    arcs = [(F(0), F(1))]
    for m in Q:
        arcs = intersect(arcs, safe_arcs(14 * m, thr))
    return arcs

def covering_floor(R, Q, thr=F(1, 14)):
    """
    The quasi-independence ratio R' = meas(both)/(meas(R-safe)*meas(Q-lonely)).
    Returns (measR, measQ, measBoth, Rprime) exact where possible.
    """
    arcsR = [(F(0), F(1))]
    for s in R:
        arcsR = intersect(arcsR, safe_arcs(s, thr))
    arcsQ = qlonely_arcs(Q, thr)
    measR = meas(arcsR); measQ = meas(arcsQ)
    both = intersect(arcsR, arcsQ); measB = meas(both)
    if measR == 0 or measQ == 0:
        return measR, measQ, measB, None
    Rp = measB / (measR * measQ)
    return measR, measQ, measB, Rp


# ======================================================================== main probes
def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("LRC(14) loneliness as a partition function -- exact validation + Riesz sum")

    # ---- small multi-far R-safe sets (the 13-r small part R; here just probe small R) ----
    test_R = [
        (1, 2),
        (1, 2, 3),
        (2, 3),
        (1, 3),
        (1, 2, 3, 4),
        (3, 5, 7),
        (2, 5),
    ]
    print(f"\n{'R':<22}{'L_exact':>16}{'L_finite':>16}{'L_riesz(M=40)':>18}")
    for R in test_R:
        Lex, _ = loneliness_exact(R)
        Lfin = loneliness_finite(R)
        Lrz = loneliness_riesz(R, 40)
        match = "OK" if Lex == Lfin else "MISMATCH"
        print(f"{str(R):<22}{str(Lex):>16}{str(Lfin):>16}{Lrz:>18.6f}  [{match}]")

    banner("Covering floor R' for genuine few-apex covering sets S = R u 14Q")
    # r = |14Z cap S| ; e.g. R = {1..11,13}, Q = {6} -> S has the apex 84=14*6
    # Use the documented few-apex examples and small multi-far variants.
    cases = [
        (tuple(range(1, 12)) + (13,), (6,)),          # r=1, M=7/89 documented
        (tuple(range(1, 11)) + (13,), (6, 11)),       # r=2, M=14/155 documented
        ((1, 2, 3), (1,)),                            # tiny
        ((2, 3, 5), (1, 2)),                          # tiny r=2
    ]
    print(f"\n{'R':<28}{'Q':<12}{'measR':>12}{'measQ':>12}{'both':>14}{'Rprime':>12}")
    for R, Q in cases:
        mR, mQ, mB, Rp = covering_floor(R, Q)
        Rps = f"{float(Rp):.5f}" if Rp is not None else "n/a"
        print(f"{str(R):<28}{str(Q):<12}{float(mR):>12.6f}{float(mQ):>12.6f}{float(mB):>14.6f}{Rps:>12}")
