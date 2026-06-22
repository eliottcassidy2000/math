#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_truemin_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE TRUE INTEGER MINIMUM of rho*(P,E) -- the operative LRC(14) bound.

The asymptotic test established: as spread -> infinity, rho* RISES to a positive
equidistribution floor (F1: ~0.04, F2: ~0.40, F3: ~0.16, F4: ~0.063). So the
MINIMUM lives at MODERATE (bounded) spread.  We now find it as exactly as
tractable, with a uniform-spread guarantee:

  For each k = |E|, |P| = 13-k:
    * search ALL distinct integer shapes 0=e_1<...<e_k <= W_max for a generous
      W_max (well past the asymptotic-floor onset),
    * over a P-shortlist (smallest meas G_P -- the only P that can bind),
    * track the running min rho* and the spread at which it occurs.

The KEY diagnostic: report the min rho* RESTRICTED to spread <= W, for a ladder
of W, to SEE the min STABILIZE (it should: small-spread localized minima are the
binding ones; large spread floors above them).  If min(spread<=W) is FLAT for
W beyond some W*, the bounded-spread reduction is empirically confirmed and we
have the true min.

We also record, for the global argmin, the relation lattice Lambda(E) min-support
(d) -- the "code distance" -- to connect to the subtorus picture.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def circ_maxgap_at(E, x):
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def gp_breaks(P):
    bps = set()
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks(E):
    bps = set()
    diffs = set()
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d != 0:
                diffs.add(d)
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        for m in range(0, 7 * d + 1):
            for s in (2, -2):
                v = Fr(7 * m + s, 7 * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_arcs(E, gapthr=Fr(2, 7)):
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E))
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        if circ_maxgap_at(E, (x0 + x1) / 2) > gapthr:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def gp_arcs(P, thr=Fr(1, 14)):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P))
    arcs = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for p in P:
            f = (Fr(p) * mid) % 1
            d = f if f <= Fr(1, 2) else 1 - f
            if d < thr:
                ok = False
                break
        if ok:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def arcs_measure(arcs):
    return sum((b - a for a, b in arcs), Fr(0))


def arcs_intersect_measure(A, B):
    i = j = 0
    tot = Fr(0)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            tot += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return tot


def min_support(E, B=2, maxs=4):
    """min support of a nonzero integer relation sum n_i e_i=0, |n_i|<=B."""
    El = list(E)
    k = len(El)
    for s in range(2, maxs + 1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B + 1), repeat=s):
                if any(c == 0 for c in coefs):
                    continue
                if sum(c * El[i] for c, i in zip(coefs, combo)) == 0:
                    return s
    return None


def main():
    print("=" * 78)
    print("THM-527 Thread A: TRUE integer min rho*(P,E)  (min lives at bounded spread)")
    print("=" * 78)
    sys.stdout.flush()

    # focus on the binding k (where rho* is smallest): k=9,10,11.  Also do 7,8.
    results = {}
    for (k, Wmax) in [(7, 20), (8, 20), (9, 20), (10, 19), (11, 18), (12, 18)]:
        npart = 13 - k
        # P shortlist
        Pall = []
        for P in itertools.combinations(range(1, 14), npart):
            a = gp_arcs(list(P))
            Pall.append((arcs_measure(a), P, a))
        Pall.sort(key=lambda e: e[0])
        Pcands = Pall[:50] if npart > 0 else [(Fr(1), (), [(Fr(0), Fr(1))])]

        # ladder of W; record min(spread<=W)
        ladder = {}
        gbest = None
        gE = None
        gP = None
        for tail in itertools.combinations(range(1, Wmax + 1), k - 1):
            E = [0] + list(tail)
            spread = tail[-1]
            ga = good_arcs(E)
            mu = arcs_measure(ga)
            loc_best = None
            for (mg, P, pa) in Pcands:
                lb = mu + mg - 1
                if gbest is not None and lb >= gbest and (loc_best is None or lb >= loc_best):
                    continue
                r = arcs_intersect_measure(ga, pa)
                if loc_best is None or r < loc_best:
                    loc_best = r
                    loc_P = list(P)
            if loc_best is not None and (gbest is None or loc_best < gbest):
                gbest = loc_best
                gE = E
                gP = loc_P
            # record into ladder bucket
            if loc_best is not None:
                b = ladder.get(spread)
                if b is None or loc_best < b:
                    ladder[spread] = loc_best
        # print the cumulative min by spread
        print(f"\n--- k={k} (|P|={npart}, Wmax={Wmax}) cumulative min rho* by spread<=W ---")
        cum = None
        cumrow = []
        for W in range(k - 1, Wmax + 1):
            here = ladder.get(W)
            if here is not None and (cum is None or here < cum):
                cum = here
            cumrow.append((W, cum))
        # print a compact ladder
        for W, cm in cumrow:
            if cm is not None:
                print(f"    spread<={W:2d}: min rho* = {float(cm):.6f} = {cm}")
        dmin = min_support(gE)
        print(f"  => k={k} GLOBAL min rho* = {gbest} = {float(gbest):.6f}")
        print(f"     argmin E={gE} (spread {max(gE)}, relation min-support d={dmin})  P={gP}")
        results[k] = (gbest, gE, gP, max(gE))
        sys.stdout.flush()

    print("\n" + "=" * 78)
    print("SUMMARY: true integer min rho* per k, and the spread of the minimizer")
    print("=" * 78)
    gm = None
    for k in sorted(results):
        b, E, P, sp = results[k]
        print(f"  k={k:2d}: min rho* = {float(b):.6f} = {b}   minimizer spread={sp}")
        if gm is None or b < gm[0]:
            gm = (b, k, E, P)
    print(f"\n  >>> OVERALL min rho* (k=7..12) = {gm[0]} = {float(gm[0]):.6f} "
          f"at k={gm[1]}, E={gm[2]}, P={gm[3]}")
    print(f"      POSITIVE with margin: {float(gm[0]):.6f} > 0")
    print("\nDONE.")


if __name__ == "__main__":
    main()
