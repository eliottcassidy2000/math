#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 floor: is the 14-sheet count SUB-BINOMIAL?  (mac-mini-2026-06-29-S3)

THM-579: covering floor holds if CV(N_R)^2 < m_Q/(1-m_Q),
N_R(t) = #{a in 0..13 : t+a/14 R-safe}, a sum of 14 correlated Bernoulli(m_R).

CREATIVE ANGLE: if the sheets are NEGATIVELY ASSOCIATED on average
   Var(N_R) <= 14 m_R (1-m_R)            (sub-binomial / independent-or-better)
then CV(N_R)^2 <= (1-m_R)/(14 m_R), and the floor criterion
   CV^2 < m_Q/(1-m_Q)  is IMPLIED by  (1-m_R)/(14 m_R) < m_Q/(1-m_Q),
i.e. by the ELEMENTARY two-variable inequality
   (1-m_R)(1-m_Q) < 14 m_R m_Q   <=>   1 < m_R + m_Q + 13 m_R m_Q,
which holds whenever m_R, m_Q exceed their THM-576 caps.  That would reduce the
WHOLE covering floor to (a) a sub-binomial sheet lemma and (b) cap lower bounds.

This script measures Var(N_R) vs the binomial value 14 m_R (1-m_R) for many R
(consec / even-heavy / 7-heavy / random / large-speed), reporting the ratio
rho = Var_real / Var_binom.  rho<=1 everywhere => sub-binomial holds.  Where it
FAILS (rho>1, expected for even/7-heavy R that cluster the sheets), we still
check the actual criterion and the clean inequality.
"""
from __future__ import annotations
import functools, math, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations

# import shared machinery
src = open('04-computation/lrc14_floor_CV_uniformity_sweep_macmini_20260629.py').read()
G = {'__name__': 'x'}
exec(compile(src.replace('if __name__', 'if False and __name__'), 'm', 'exec'), G)
lonely_set = G['lonely_set']; measure = G['measure']; intersect = G['intersect']; shift = G['shift']


def sheet_stats(R):
    Lr = lonely_set(tuple(sorted(R)))
    m_R = measure(Lr)
    if m_R == 0:
        return None
    EN = 14 * m_R
    EN2 = F(0)
    for d in range(-13, 14):
        EN2 += (14 - abs(d)) * measure(intersect(Lr, shift(Lr, F(d, 14))))
    Var = EN2 - EN * EN
    Var_binom = 14 * m_R * (1 - m_R)
    return m_R, Var, Var_binom


def cap_for_size(j, universe=15):
    best = None
    for Q in combinations(range(1, universe + 1), j):
        m = measure(lonely_set(Q))
        if best is None or m < best:
            best = m
    return best


def main():
    print("=" * 90)
    print("LRC14 sheet-count sub-binomial test (mac-mini-S3)")
    print("=" * 90)
    rng = random.Random(2024)

    # caps for the small parts (THM-576 pattern); compute once
    print("\nCaps C-(14-j,2)/91 vs exact min meas(lonely) for |P|=j:")
    caps = {}
    for j in range(2, 12):
        cap = cap_for_size(j, universe=15)
        caps[j] = cap
        tri = F((14 - j) * (13 - j), 2) / 91 if j <= 13 else F(0)
        print(f"  j={j:2d}: cap={cap} ({float(cap):.4f})  C(14-j,2)/91={float(tri):.4f}")

    print("\nSub-binomial test: rho = Var(N_R)/[14 m_R(1-m_R)].  rho<=1 => sub-binomial.")
    print(f"{'R':<34}{'m_R':>8}{'Var':>9}{'Var_bin':>9}{'rho':>7}{'CV^2':>8}{'(1-mR)/14mR':>12}")
    families = {
        'consec_7': list(range(1, 8)),
        'consec_8': list(range(1, 9)),
        'consec_11': list(range(1, 12)),
        'even-heavy_7': [2, 4, 6, 8, 10, 12, 1],
        'even-heavy_8': [2, 4, 6, 8, 10, 12, 1, 3],
        '7-heavy': [7, 14 - 7, 1, 2, 3, 4, 5],   # contains 7
        'mult-of-7-ish': [7, 21, 1, 2, 3, 4, 5],
        'large-coprime_7': [1, 3, 5, 9, 11, 13, 15],
    }
    maxrho = 0.0
    for name, R in families.items():
        st = sheet_stats(R)
        if st is None:
            continue
        m_R, Var, Var_binom = st
        rho = float(Var / Var_binom) if Var_binom != 0 else float('nan')
        cv2 = float(Var / (14 * m_R) ** 2)
        cv2_bin = float((1 - m_R) / (14 * m_R))
        maxrho = max(maxrho, rho)
        flag = '' if rho <= 1.0001 else '  <-- super-binomial'
        print(f"{name:<14}{str(tuple(sorted(R))):<20}{float(m_R):>8.4f}{float(Var):>9.4f}"
              f"{float(Var_binom):>9.4f}{rho:>7.3f}{cv2:>8.4f}{cv2_bin:>12.4f}{flag}")

    # random sweep for max rho per size
    print("\nMax rho over random 14-free R (does sub-binomial ever fail badly?):")
    for size in (7, 8, 9, 10, 11):
        mr = 0.0; argR = None
        for _ in range(250):
            R = rng.sample([x for x in range(1, 60) if x % 14], size)
            st = sheet_stats(R)
            if st is None:
                continue
            m_R, Var, Var_binom = st
            if Var_binom == 0:
                continue
            rho = float(Var / Var_binom)
            if rho > mr:
                mr = rho; argR = tuple(sorted(R))
        print(f"  |R|={size}: max rho={mr:.3f} at {argR}")

    # the clean elementary inequality, using cap lower bounds for m_R (|R|=13-r), m_Q (|Q|=r)
    print("\nClean 2-var inequality  1 < m_R + m_Q + 13 m_R m_Q  at WORST measures (caps):")
    for r in range(2, 7):
        mR = caps[13 - r]      # worst m_R for |R|=13-r
        mQ = caps[r]           # worst m_Q for |Q|=r
        lhs = float(mR + mQ + 13 * mR * mQ)
        # what the clean ineq needs vs the binomial-CV criterion threshold
        print(f"  r={r}: m_R>={float(mR):.4f} (|R|={13-r}), m_Q>={float(mQ):.4f} (|Q|={r})  "
              f"=> m_R+m_Q+13 m_R m_Q >= {lhs:.4f}  {'> 1 OK' if lhs>1 else '<=1 FAIL'}")

    print("\n" + "=" * 90)
    print(f"max rho over named families = {maxrho:.3f}")
    print("If sub-binomial holds (rho<=1) -> floor reduces to the clean inequality (all r OK).")
    print("If rho>1 for even/7-heavy R -> need the actual Var, but the clean ineq margin may absorb it.")
    print("=" * 90)


if __name__ == "__main__":
    main()
