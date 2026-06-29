#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 covering FLOOR (CRUX-1): a clean closed-form sufficient condition via
the coefficient of variation of the 14-sheet count.

Instance: mac-mini-2026-06-29-S2.

The one remaining theorem (kps HYP-3415/3129): for covering S = R u 14Q,
  R' := meas(lonely(S)) / [ meas(lonely(R)) * meas(lonely(Q)) ] > 0   uniformly.
Equivalently |SPEC| < product, SPEC = sum_{n!=0} chat(n) conj(ghat(n)),
chat = Fourier of 1_{R-safe} (support gcd(R)Z), ghat = Fourier of 1_{Q-lonely}
pulled back by u=14t (support 14Z).  Certified R' >= 0.642 (HYP-3129, exact-low
+ Cauchy-Schwarz tail); open as a closed-form uniform bound.

NEW (this script): a single Cauchy-Schwarz on ALL of SPEC, using that ghat lives
ONLY on 14Z, gives the CLOSED-FORM bound
  |SPEC| <= sqrt( sum_{N!=0}|chat(14N)|^2 ) * sqrt( var(Q-lonely) ).
The projection of R-safe onto 14Z-frequencies is exactly (1/14)*N_R(t), the
14-SHEET COUNT N_R(t) = #{a in 0..13 : t+a/14 is R-safe} (HYP-3140 object).  So
  sum_{N!=0}|chat(14N)|^2 = Var(N_R)/196,   var(Q-lonely) = m_Q - m_Q^2,
  product = m_R * m_Q = (E[N_R]/14)*m_Q,
hence
  R' >= 1 - CV(N_R) * sqrt( (1 - m_Q)/m_Q ),   CV(N_R)=sqrt(Var N_R)/E[N_R].

This DERIVES the sheet-count / fiber-PGF framing and reduces the floor to a
coefficient-of-variation bound (an R-only quantity on the finite 14-grid) versus
m_Q (a Q-only quantity, >= C(14-r,2)/91 by THM-576).
"""
from __future__ import annotations
import functools, math
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd

W = F(1, 14)


def danger_intervals(p):
    ivs = []
    for k in range(p + 1):
        lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
        if hi > lo:
            ivs.append((lo, hi))
    return ivs


def union_intervals(lists):
    ivs = sorted(iv for lst in lists for iv in lst)
    if not ivs:
        return []
    out = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo > out[-1][1]:
            out.append([lo, hi])
        else:
            out[-1][1] = max(out[-1][1], hi)
    return [(a, b) for a, b in out]


def complement(intervals):
    """[0,1) minus the union (intervals assumed sorted/disjoint)."""
    out = []
    cur = F(0)
    for lo, hi in intervals:
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        out.append((cur, F(1)))
    return out


def lonely_set(P):
    """complement of union of danger combs = {x: ||p x||>=1/14 all p}."""
    u = union_intervals([danger_intervals(p) for p in P])
    return complement(u)


def measure(intervals):
    return sum((hi - lo for lo, hi in intervals), F(0))


def intersect(a, b):
    out = []
    i = j = 0
    for (lo1, hi1) in a:
        for (lo2, hi2) in b:
            lo = max(lo1, lo2); hi = min(hi1, hi2)
            if hi > lo:
                out.append((lo, hi))
    return out


def shift(intervals, s):
    """shift intervals by s (mod 1), returning sorted disjoint intervals in [0,1)."""
    res = []
    for lo, hi in intervals:
        length = hi - lo
        a = lo + s
        fa = a - F(math.floor(a))      # start in [0,1)
        end = fa + length
        if end <= 1:
            res.append((fa, end))
        else:
            res.append((fa, F(1)))
            res.append((F(0), end - 1))
    return union_intervals([res])


def autocorr(L, s):
    """meas( L  cap  (L shifted by s) )."""
    return measure(intersect(L, shift(L, s)))


def floor_row(R, Q):
    Rt = tuple(sorted(R))
    Qt = tuple(sorted(Q))
    S = tuple(sorted(set(Rt) | set(14 * m for m in Qt)))
    Lr = lonely_set(Rt)
    Lq = lonely_set(Qt)            # m_Q = meas(lonely(Q)) on standard circle
    Ls = lonely_set(S)
    m_R = measure(Lr); m_Q = measure(Lq); m_S = measure(Ls)
    Rprime = m_S / (m_R * m_Q) if m_R * m_Q != 0 else None
    # Var(N_R) via autocorrelations of lonely(R) at shifts d/14
    EN = 14 * m_R
    EN2 = F(0)
    for d in range(-13, 14):
        c = autocorr(Lr, F(d, 14))
        EN2 += (14 - abs(d)) * c       # multiplicity of shift d among 14x14 sheet pairs
    VarN = EN2 - EN * EN
    CV2 = VarN / (EN * EN)
    bound = 1 - math.sqrt(float(CV2)) * math.sqrt(float((1 - m_Q) / m_Q))
    return dict(R=Rt, Q=Qt, m_R=m_R, m_Q=m_Q, Rprime=Rprime,
                CV2=CV2, bound=bound, cond=float(CV2) < float(m_Q / (1 - m_Q)))


def main():
    print("=" * 80)
    print("LRC14 covering floor: R' >= 1 - CV(N_R)*sqrt((1-m_Q)/m_Q)   (mac-mini-S2)")
    print("=" * 80)
    rows = [
        (list(range(1, 13)), [1, 2]),     # r=2 binding row (agent: actual R'=0.7015)
        (list(range(1, 12)), [1, 2, 3]),  # r=3
        (list(range(1, 11)), [1, 2, 3, 4]),
        (list(range(1, 10)), [1, 2, 3, 4, 5]),
        (list(range(1, 9)),  [1, 2, 3, 4, 5, 6]),  # r=6
        # a few even-heavy R (the 2-adic binding worry, S259)
        ([2, 4, 6, 8, 10, 12, 1, 3], [1, 2]),
        ([1, 2, 4, 6, 8, 10, 12], [1, 3]),
    ]
    print(f"\n{'R':<26}{'Q':<12}{'m_R':>9}{'m_Q':>8}{'R_actual':>10}{'CV^2':>8}{'bound':>9}{'>0?':>5}")
    for R, Q in rows:
        d = floor_row(R, Q)
        Rp = float(d['Rprime'])
        print(f"{str(d['R']):<26}{str(d['Q']):<12}{float(d['m_R']):>9.4f}{float(d['m_Q']):>8.4f}"
              f"{Rp:>10.4f}{float(d['CV2']):>8.4f}{d['bound']:>9.4f}{'YES' if d['bound']>0 else 'no':>5}")
    print("\nInterpretation:")
    print(" - bound = the closed-form C-S lower bound on R'; actual R' is exact.")
    print(" - If bound>0 for a row, the floor is PROVED for it by the clean CV criterion,")
    print("   with NO exact-low computation -- purely CV(N_R) (R-only) and m_Q (Q-only).")
    print(" - CV(N_R)^2 < m_Q/(1-m_Q) is the sufficient condition; m_Q>=C(14-r,2)/91 (THM-576).")
    print("=" * 80)


if __name__ == "__main__":
    main()
