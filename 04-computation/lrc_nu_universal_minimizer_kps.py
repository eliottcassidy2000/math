#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_nu_universal_minimizer_kps.py  (kind-pasteur 2026-06-22, THM-527 witness route)

DECISIVE TEST for the ELEMENTARY (compactness-free) witness-floor closure.

The Bonferroni route (lrc_rhoglob_closedform_kpswf10.py):
    rho*_glob(P,E) = meas(GOOD_{1/7}(E) cap G_P)
                   >= nu(E) + meas(G_P) - 1,
    nu(E) := meas{x in [0,1): circular maxgap of {frac(e*x): e in E} > 1/7}.
If nu(E) >= nu_min(k) UNIVERSALLY (all primitive k-shapes E, ANY spread) and
meas(G_P) >= cap_k, then rho*_glob >= nu_min(k) + cap_k - 1 > 0 with NO
compactness / NO bounded-spread reduction.

KEY STRUCTURAL FACT (proved here in comment): nu is SCALE-INVARIANT,
    nu(c*E) = nu(E)  for every integer c >= 1,
because x -> c*x mod 1 is measure-preserving and the integrand is 1-periodic.
So nu depends ONLY on the primitive shape E/gcd(E); "wide via scaling" does NOT
lower nu.  The ONLY danger is a wide PRIMITIVE shape (large spread, gcd 1).

THIS SCRIPT tests whether consecutive {0,1,...,k-1} minimizes nu over ALL
primitive k-shapes, INCLUDING wide ones (spread up to W).  If consec is the
global min, nu_min(k) = nu_consec(k) is an exact rational and the closure is
elementary.  If some wide primitive shape beats consec, we record it (then the
bounded-spread reduction is genuinely needed).
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


def good_breaks(E, thr_den=7):
    """Breakpoints where the maxgap>1/7 indicator can change: rationals t/d and
    (7m+-1)/(7d) for every pairwise difference d=|e_i-e_j|."""
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
        for m in range(0, thr_den * d + 1):
            for s in (1, -1):
                v = Fr(thr_den * m + s, thr_den * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def nu(E, gapthr=Fr(1, 7)):
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        if circ_maxgap_at(E, (x0 + x1) / 2) > gapthr:
            tot += (x1 - x0)
    return tot


def primitive(E):
    g = reduce(gcd, [abs(e) for e in E if e != 0], 0)
    return g == 1


def main():
    print("=" * 78)
    print("THM-527 witness route: is nu MINIMIZED at consecutive over ALL")
    print("primitive k-shapes (including WIDE ones)?  (=> elementary closure)")
    print("=" * 78)
    sys.stdout.flush()

    cap = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91),
           11: Fr(66, 91), 12: Fr(6, 7), 13: Fr(1)}

    # nu is scale-invariant: verify once.
    Etest = [0, 1, 2, 3, 4, 5, 6, 7]
    print("\n[scale-invariance check] nu(consec_8) vs nu(2*consec_8) vs nu(3*..):")
    for c in (1, 2, 3, 5):
        Ec = [c * e for e in Etest]
        print(f"    c={c}: nu = {nu(Ec)} = {float(nu(Ec)):.6f}")
    sys.stdout.flush()

    # Main: minimize nu over primitive k-shapes with spread up to W.
    # We scan ALL tails for moderate W exhaustively; for the wide stress test we
    # additionally scan shapes with one or two FAR elements appended to a small
    # base (the genuine-wide structure), spread up to WIDE.
    print("\n--- nu_min(k) over primitive shapes, exhaustive bounded + wide stress ---")
    results = {}
    for k in range(8, 14):
        consec_nu = nu(list(range(k)))
        best = consec_nu
        bestE = list(range(k))
        nscan = 0
        # (1) exhaustive moderate window W = k+6
        W1 = k + 6
        for tail in itertools.combinations(range(1, W1 + 1), k - 1):
            E = [0] + list(tail)
            if not primitive(E):
                continue
            v = nu(E)
            nscan += 1
            if v < best:
                best = v
                bestE = E
            if nscan > 300000:
                break
        # (2) WIDE stress: base = consec_{k-1} or perforated, plus a FAR element
        #     up to spread WIDE (the wide primitive shapes scaling can't reach).
        WIDE = 80
        for far in range(k, WIDE + 1):
            for base in (list(range(k - 1)),
                         [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12][:k - 1],
                         list(range(0, 2 * (k - 1), 2))):
                E = sorted(set(base + [far]))
                if len(E) != k or not primitive(E):
                    continue
                v = nu(E)
                nscan += 1
                if v < best:
                    best = v
                    bestE = E
        # (3) WIDE stress: TWO far elements (doublet structure)
        for f1 in range(k, 40):
            for f2 in range(f1 + 1, f1 + 40):
                base = list(range(k - 2))
                E = sorted(set(base + [f1, f2]))
                if len(E) != k or not primitive(E):
                    continue
                v = nu(E)
                nscan += 1
                if v < best:
                    best = v
                    bestE = E
        results[k] = (best, bestE, consec_nu)
        beats = "" if best == consec_nu else "  <<< NON-CONSEC MIN!"
        print(f"  k={k:2d} ({nscan} shapes): nu_min={best}={float(best):.5f}  "
              f"consec={float(consec_nu):.5f}{beats}")
        if best != consec_nu:
            print(f"      minimizer E = {bestE}")
        sys.stdout.flush()

    print("\n--- ELEMENTARY closed-form floor rho*_glob >= nu_min(k)+cap_k-1 ---")
    print("    k   nu_min(k)     cap_k       nu_min+cap-1    >0?")
    allpos = True
    worst = None
    for k in range(8, 14):
        nm = results[k][0]
        ub = nm + cap[k] - 1
        pos = ub > 0
        allpos = allpos and pos
        if worst is None or ub < worst:
            worst = ub
        print(f"    {k:2d}  {str(nm):>12s}  {float(cap[k]):.5f}   {float(ub):+.5f}     "
              f"{'YES' if pos else 'NO'}")
    print(f"\n  k<=7: nu_k = 1 (pigeonhole 1/k>=1/7), rho*_glob = meas(G_P) > 0.")
    print(f"  ALL k=8..13 positive?  {allpos}")
    print(f"  worst closed-form floor = {worst} = {float(worst):.5f}")
    print("=" * 78)
    if allpos and all(results[k][0] == results[k][2] for k in range(8, 14)):
        print("  *** ELEMENTARY UNIVERSAL CLOSURE CONFIRMED (computationally). ***")
        print("  Consecutive minimizes nu over ALL tested primitive shapes (incl. wide),")
        print("  so nu_min(k)=nu_consec(k) exact rational; nu_min+cap-1>0 => rho*_glob>0")
        print("  for ALL P,E with NO compactness/bounded-spread reduction. Remaining:")
        print("  PROVE consec globally minimizes nu (the one pure three-distance lemma).")
    else:
        print("  Either a wide shape beats consec (bounded-spread reduction needed),")
        print("  or some k fails positivity.  See minimizer above.")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
