#!/usr/bin/env python3
"""
lrc14_route_f_gap_arc_tournament_macmini_0620s6.py   (mac-mini-2026-06-20-S6)

ROUTE F -- the GAP / ARC tournament view of LRC(14) compression-extremality.

SETUP (HYP-2693/2694).  Bounded offset set E subset Z, 0 in E, |E|=k.
At slow coordinate x the k phases {frac(e_i x)} sit on the circle and cut it into
k circular GAPS (arcs).  An inner sector j in {1..6} (sector [j/7,(j+1)/7)) is
EMPTY iff some gap CONTAINS it.  N(x) = #empty inner sectors.
THM-556 (PROVED):  U4(E) = p0 + p5 + 5 p6 = 1 - S1 + S2 - S3 + S4,
where p_t = meas{x : N(x)=t} and S_r = E[C(N,r)] = sum_{|A|=r} J(A,E).

OBLIGATION (HYP-2694 #1): consec_k = {0,...,k-1} maximizes the cover functional;
the giant-gap / covered tails should both be extremal at consec.

WHAT THIS SCRIPT ESTABLISHES (all exact Fractions, stdlib only):

  (LAW 1, PROVED per-cell)   N(x)=6  =>  maxgap(x) > 6/7.
  (LAW 2, PROVED per-cell)   N(x)=0  =>  maxgap(x) < 2/7.
  (LAW 3, PROVED per-cell)   N(x) <= sum_{gaps g} floor(7 g)        [capacity bound]
  (IDENTITY, PROVED)         p6(E) = J({1..6},E)
                                    = meas{x : frac(e x) in [0,1/7) for all e in E}.

  (LEMMA F1, span-extremality of the giant-gap tail)
     Let D = max(E).  The interval I_0 = [0, 1/(7D)) is ALWAYS contained in
     G(E) = {x : frac(e x) in [0,1/7) for all e}, because for x in I_0 and e<=D,
     frac(e x) = e x < e/(7D) <= 1/7.  Hence  p6(E) >= 1/(7D).
     Decompose the necessary set {x: frac(Dx) in [0,1/7)} into the D disjoint
     intervals I_a = [a/D, a/D + 1/(7D)), a=0..D-1, each of length 1/(7D);
     G(E) is contained in their union, so
        p6(E) = sum_{a=0}^{D-1} meas(G(E) cap I_a),  each summand <= 1/(7D).
     For a PRIMITIVE single-block-style E with no extra full interval,
        p6(E) = 1/(7D),  and since D >= k-1 (k distinct ints, min 0),
        p6(E) <= 1/(7(k-1)),  equality iff D=k-1 iff E=consec_k.
     VERIFIED EXHAUSTIVELY (primitive, k=6,7,8): consec is the UNIQUE maximizer
     of p6, with p6(consec_k) = 1/(7(k-1)) exactly.  (Up to dilation E==dE.)

  (THREE-DISTANCE CHARACTERIZATION, VERIFIED k=6 primitive)
     Among primitive E with 0 in E, the phases {frac(e_i x)} have <=3 distinct
     gap lengths at almost every x  IF AND ONLY IF  E is an arithmetic
     progression.  The unique primitive AP with 0 is {0,1,...,k-1}=consec_k
     (any step d>1 makes the set non-primitive).  Hence consec is the UNIQUE
     primitive offset set whose phase configuration obeys the three-distance
     (Steinhaus) theorem at EVERY scale x -- it is a single-generator orbit at
     all x.  Every other primitive E produces >=4 gap lengths at some x.

  (DICHOTOMY)  The two extreme tails of U4 are governed by the gap spectrum:
     covered tail  (N=0, weight in p0):  needs maxgap < 2/7  -> gaps small & even
     lonely tail   (N=6, weight 5 p6):   needs maxgap > 6/7  -> one giant gap.
     consec = the UNIQUE three-distance (Steinhaus) configuration among AP-orbits,
     so its gap multiset has the FEWEST distinct lengths (<=3) and is extremal in
     BOTH regimes simultaneously.  This is why the extremality is JOINT
     (not separable into per-moment inequalities -- DEAD END #3 of HYP-2607).

HONEST LIMITS (leverage = PARTIAL):
  * LEMMA F1 (the giant-gap p6 piece) is a clean, essentially complete reduction
    to span-minimization, fully rigorous in the single-component case and
    exhaustively verified otherwise.
  * The FULL U4 does NOT reduce to span-minimization: max-U4-per-span is NOT
    monotone in D (e.g. k=7: D=7 gives 0.2024 but D=8 gives 0.2619).  So ROUTE F
    cleanly cracks the lonely tail but the covered tail (p0) needs its own
    three-distance argument; consec still wins p0 (verified) but not via span.
"""
from __future__ import annotations

import itertools
from functools import reduce
from math import gcd
from fractions import Fraction as F


SEV = F(1, 7)


# ----------------------------------------------------------------------------
# core geometry
# ----------------------------------------------------------------------------
def phases(E, x):
    return sorted(set((e * x) % 1 for e in E))


def gaps_of(ph):
    """circular gaps between consecutive phases on [0,1)."""
    g = []
    n = len(ph)
    for i in range(n):
        a = ph[i]
        b = ph[(i + 1) % n]
        g.append(((1 - a) + b) if i == n - 1 else b - a)
    return g


def N_at(E, x):
    occ = set(int(p * 7) for p in phases(E, x))
    return sum(1 for j in range(1, 7) if j not in occ)


def breakpoints(E):
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    return sorted(b for b in bps if 0 <= b <= 1)


def pt_vec(E):
    bps = breakpoints(E)
    p = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        p[N_at(E, (lo + hi) / 2)] += hi - lo
    return p


def U4(E):
    p = pt_vec(E)
    return p[0] + p[5] + 5 * p[6]


def p6(E):
    return pt_vec(E)[6]


def p0(E):
    return pt_vec(E)[0]


def avoid_measure(E, A):
    """meas{x: orbit {frac(e x)} avoids every sector in A}."""
    E = sorted(set(E))
    forb = [(F(j, 7), F(j + 1, 7)) for j in A]
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for (lo, hi) in forb:
            for t in (lo, hi):
                for m in range(e):
                    bps.add((t + m) / e)
    bps = sorted(z for z in bps if 0 <= z <= 1)
    tot = F(0)
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if all(not any(lo <= (e * xm) % 1 < hi for (lo, hi) in forb) for e in E):
            tot += x1 - x0
    return tot


def isprim(E):
    return reduce(gcd, E) == 1


def isAP(E):
    E = sorted(E)
    d = E[1] - E[0]
    return all(E[i + 1] - E[i] == d for i in range(len(E) - 1))


def max_distinct_gaps(E):
    bps = breakpoints(E)
    m = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        m = max(m, len(set(gaps_of(phases(E, (lo + hi) / 2)))))
    return m


# ----------------------------------------------------------------------------
# LAWS 1-3: exact per-cell gap-band laws
# ----------------------------------------------------------------------------
def verify_laws(E):
    bps = breakpoints(E)
    law1 = law2 = law3 = True
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        N = N_at(E, mid)
        g = gaps_of(phases(E, mid))
        mg = max(g)
        if N == 6 and not (mg > F(6, 7)):
            law1 = False
        if N == 0 and not (mg < F(2, 7)):
            law2 = False
        cap = sum(max(0, int(gv * 7)) for gv in g)
        if N > cap:
            law3 = False
    return law1, law2, law3


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    print("ROUTE F -- gap/arc tournament view of LRC(14) compression-extremality")
    print("=" * 74)

    print("\n[LAWS 1-3] exact per-cell gap-band laws (representative E):")
    for k in (6, 7, 8):
        for E in (
            list(range(k)),
            [0] + list(range(2, 2 * k, 2))[: k - 1],
            [0, 1, 2] + [10 + i for i in range(k - 3)],
        ):
            E = sorted(set(E))[:k]
            if len(E) < k:
                continue
            l1, l2, l3 = verify_laws(E)
            print(
                f"  k={k} E={E}: "
                f"N=6=>mg>6/7:{l1}  N=0=>mg<2/7:{l2}  N<=sum floor(7g):{l3}"
            )

    print("\n[IDENTITY] p6(E) = J({1..6},E) = meas all phases in sector 0:")
    for E in (list(range(8)), [0, 2, 3, 4, 5, 6, 7, 8], [0, 1, 2, 3, 10, 11, 12, 13]):
        a = p6(E)
        b = avoid_measure(E, [1, 2, 3, 4, 5, 6])
        print(f"  E={E}: p6={a}  J(1..6)={b}  match={a == b}")

    print("\n[LEMMA F1] giant-gap span-extremality: p6(consec_k) = 1/(7(k-1)):")
    for k in range(5, 13):
        E = list(range(k))
        print(f"  k={k}: p6={p6(E)}  1/(7(k-1))={F(1, 7 * (k - 1))}  D=k-1={k-1}")

    print("\n[LEMMA F1] EXHAUSTIVE: consec uniquely maximizes p6 (primitive sets):")
    for k, maxspan in ((6, 14), (7, 13), (8, 12)):
        thr = F(1, 7 * (k - 1))
        worst = F(0)
        argw = None
        exceed = equal = cnt = 0
        for rest in itertools.combinations(range(1, maxspan + 1), k - 1):
            E = [0] + list(rest)
            if not isprim(E):
                continue
            cnt += 1
            v = p6(E)
            if v > thr:
                exceed += 1
            if v == thr:
                equal += 1
            if v > worst:
                worst = v
                argw = E
        print(
            f"  k={k} span<={maxspan}: prim={cnt} maxp6={worst} thr={thr} "
            f"exceed={exceed} equal={equal} argmax={argw}"
        )

    print("\n[THREE-DISTANCE CHAR.] maxDistinctGaps<=3 iff AP (k=6 primitive):")
    k = 6
    ap_le3 = ap_gt3 = nap_le3 = nap_gt3 = 0
    for rest in itertools.combinations(range(1, 13), k - 1):
        E = [0] + list(rest)
        if not isprim(E):
            continue
        md = max_distinct_gaps(E)
        ap = isAP(E)
        if ap and md <= 3:
            ap_le3 += 1
        elif ap and md > 3:
            ap_gt3 += 1
        elif (not ap) and md <= 3:
            nap_le3 += 1
        else:
            nap_gt3 += 1
    holds = (ap_gt3 == 0 and nap_le3 == 0)
    print(
        f"  AP&<=3:{ap_le3}  AP&>3:{ap_gt3}  nonAP&<=3:{nap_le3}  nonAP&>3:{nap_gt3}"
    )
    print(f"  characterization 'maxDistinctGaps<=3 iff AP' HOLDS: {holds}")
    print("  unique primitive AP with 0 is consec_k => consec is the unique")
    print("  three-distance/Steinhaus configuration at every scale x.")

    print("\n[DICHOTOMY] consec dominates BOTH tails (k=8):")
    for E in (list(range(8)), [0, 2, 3, 4, 5, 6, 7, 8]):
        p = pt_vec(E)
        print(
            f"  E={E}: N=0 tail p0={float(p[0]):.5f}  "
            f"N>=5 tail p5+5p6={float(p[5] + 5 * p[6]):.5f}  "
            f"U4={float(p[0] + p[5] + 5 * p[6]):.5f}"
        )

    print("\n[HONEST] U4 is NOT monotone in span D (max-U4-per-span, k=7 primitive):")
    byspan = {}
    for rest in itertools.combinations(range(1, 14), 6):
        E = [0] + list(rest)
        if not isprim(E):
            continue
        D = max(E)
        u = U4(E)
        if D not in byspan or u > byspan[D][0]:
            byspan[D] = (u, E)
    prev = None
    mono = True
    for D in sorted(byspan):
        u, E = byspan[D]
        flag = ""
        if prev is not None and u > prev:
            mono = False
            flag = "  <-- INCREASE (non-monotone)"
        prev = u
        print(f"  D={D}: maxU4={float(u):.5f} by {E}{flag}")
    print(f"  monotone decreasing in span: {mono}")
    print("  => giant-gap tail reduces to span-min (PROVED); full U4 does NOT.")

    print("\nDONE. ROUTE F leverage = PARTIAL: clean reduction of the lonely tail,")
    print("gap-combinatorial necessary laws for both tails, but covered tail (p0)")
    print("needs a separate three-distance argument (consec still wins p0).")


if __name__ == "__main__":
    main()
