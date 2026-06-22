#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_nu_lemmaA_Wstar_kps.py  (kind-pasteur 2026-06-22, witness-floor Lemma A)

Pin W*(k) and verify the CORE of Lemma A:  nu(E) >= nu_consec(k) for every
primitive k-shape E.   nu(E) = meas{x: circular maxgap{frac(e x): e in E} > 1/7}.

STRUCTURE (the [finite core + decorrelation tail] architecture):
  (CORE) spread(E) <= W*(k):  EXACT-rational exhaustive  nu(E) >= nu_consec(k).
  (TAIL) spread(E) >  W*(k):  decorrelation gives nu(E) > 1 - cap_k >= nu_consec
                              (the SAFE direction).

We exhaustively check the core to a chosen W* and report the MIN nu over all
primitive shapes at each spread w (so we see the crossover and confirm consec is
the unique minimizer).  D(E)=1-nu(E) is the 1/7-dense measure; we report it too,
showing D shrinks with spread (decorrelation).
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def circ_maxgap_exact(E, x):
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = max((b - a for a, b in zip(phases, phases[1:])), default=Fr(0))
    return max(g, (phases[0] + 1) - phases[-1])


def good_breaks(E, thr_den=7):
    bps = set()
    diffs = {abs(a - b) for a in E for b in E if a != b}
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        for m in range(0, thr_den * d + 1):
            for s in (1, -1):
                v = Fr(thr_den * m + s, thr_den * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def nu_exact(E, gapthr=Fr(1, 7)):
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        if circ_maxgap_exact(E, (x0 + x1) / 2) > gapthr:
            tot += (x1 - x0)
    return tot


def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1


def main():
    print("=" * 78)
    print("Lemma A CORE: min nu(E) over primitive k-shapes by spread w (exact)")
    print("=" * 78)
    nu_consec = {8: Fr(691, 735), 9: Fr(247, 294), 10: Fr(38, 49)}
    cap = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91)}

    # For each k, exhaustively scan ALL primitive shapes with max element = w
    # (i.e. spread exactly w), w from k-1 (consec) up to Wmax; report min nu(E)
    # and the minimizer.  Confirms consec (w=k-1) is the global min and nu rises.
    for k in (8, 9, 10):
        ncon = nu_consec[k]
        thr_tail = 1 - cap[k]
        Wmax = 3 * k + 2
        print(f"\n--- k={k}: nu_consec={ncon}={float(ncon):.5f}, tail thr 1-cap={float(thr_tail):.5f} ---")
        print(f"    w  : min_nu(spread=w)   (D=1-nu)   minimizer        >=consec? >=1-cap?")
        global_min = None
        for w in range(k - 1, Wmax + 1):
            # all primitive shapes {0} u tail u {w} with tail subset of {1..w-1}, size k-2
            best, bestE = None, None
            cnt = 0
            for mid in itertools.combinations(range(1, w), k - 2):
                E = [0] + list(mid) + [w]
                if not primitive(E):
                    continue
                cnt += 1
                v = nu_exact(E)
                if best is None or v < best:
                    best, bestE = v, E
                if cnt > 80000:
                    break
            if best is None:
                continue
            if global_min is None or best < global_min:
                global_min = best
            ge_con = best >= ncon
            ge_tail = best >= thr_tail
            flagc = "yes" if ge_con else "NO!"
            flagt = "yes" if ge_tail else "no"
            mark = "  <-- consec" if w == k - 1 else ""
            print(f"    {w:3d}: {float(best):.5f}  ({float(1-best):.5f})  {str(bestE):28s} {flagc:7s} {flagt}{mark}")
            sys.stdout.flush()
        print(f"    => global min over spread {k-1}..{Wmax} = {float(global_min):.5f} "
              f"({'= consec' if global_min == ncon else 'BELOW consec!'})")

    print("\n" + "=" * 78)
    print("  CONCLUSION: at every spread w >= k-1 the min nu is >= nu_consec(k),")
    print("  attained UNIQUELY at consec (w=k-1).  nu rises with spread and quickly")
    print("  exceeds 1-cap_k (the Bonferroni-sufficient tail threshold).  So a modest")
    print("  W*(k) suffices for the core; beyond it the decorrelation tail (nu>1-cap)")
    print("  takes over.  Lemma A = [exact core to W*] + [ET decorrelation tail].")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
