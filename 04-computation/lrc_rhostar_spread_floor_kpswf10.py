#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_spread_floor_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE DECISIVE QUESTION for the compactness proof:

   Does  inf_E rho*(P,E)  STAY POSITIVE as the spread -> infinity,
   or does it SHRINK to 0 (compactness FAILS)?

THM-527 part D ASSERTS the minimizer has bounded spread (<= ~30) because
"increasing spread RAISES mu". But the previous compact-min search showed the
min rho* DROPPING with the spread cap W (1/84 at consec-k9 -> 1/210 at k11,W16).
We must determine whether this is (a) a genuine floor that we have not yet
reached, or (b) rho* -> 0 along some spread-growing family => NO compactness.

We test, AT FIXED k, the min of rho*(P,E) as W = spread cap grows, AND we track
the WORST family explicitly to see WHY it shrinks (the thin-window mechanism vs
the subtorus-equidistribution floor).

Mechanism to watch (grounding point 2):
  rho* = meas( G_P  cap  {maxgap{e_i x} > 2/7} ).
  As spread grows with offsets in GENERAL POSITION (relation-free), the orbit
  {e_i x} -> iid on the full k-torus => mu -> F(k) > 0; intersect G_P => a
  positive product-like floor.
  BUT structured spread (short relations / arithmetic progressions among the
  e_i) keeps the orbit on a thin subtorus and can localize the GOOD set onto a
  measure that shrinks (the 0-anchored window ~ 1/maxE). The min thus probes
  whether STRUCTURED spread can drive rho* -> 0.

OUTPUT: per k, the min rho* as W grows; the argmin family; whether it plateaus
(floor) or keeps dropping (-> investigate as the genuine crux).
"""
import itertools
import sys
from fractions import Fraction as Fr


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


def main():
    print("=" * 78)
    print("THM-527 Thread A: DOES inf rho* STAY POSITIVE as spread grows?")
    print("=" * 78)
    sys.stdout.flush()

    # fix a representative k and its npart; sweep the spread cap W upward.
    # k=9 (npart=4) is the canonical worst region; also do k=7 (npart=6).
    for (k, npart) in [(7, 6), (9, 4), (11, 2)]:
        print(f"\n========== k={k}, |P|={npart} ==========")
        # P shortlist: smallest meas(G_P)
        Pall = []
        for P in itertools.combinations(range(1, 14), npart):
            a = gp_arcs(list(P))
            Pall.append((arcs_measure(a), P, a))
        Pall.sort(key=lambda e: e[0])
        Pcands = Pall[:40]
        print(f"  (P shortlist: {len(Pcands)} smallest-meas G_P; "
              f"min meas={float(Pcands[0][0]):.4f})")
        sys.stdout.flush()

        running = {}     # W -> (min rho*, argE, argP)
        for W in range(k - 1, k + 14):
            best = None
            bestE = None
            bestP = None
            ns = 0
            for tail in itertools.combinations(range(1, W + 1), k - 1):
                if tail[-1] != W:
                    continue    # only NEW shapes that USE the new max (spread==W)
                E = [0] + list(tail)
                ga = good_arcs(E)
                mu = arcs_measure(ga)
                ns += 1
                for (mg, P, pa) in Pcands:
                    lb = mu + mg - 1
                    if best is not None and lb >= best:
                        continue
                    r = arcs_intersect_measure(ga, pa)
                    if best is None or r < best:
                        best = r
                        bestE = E
                        bestP = list(P)
                if ns > 80000:
                    break
            # carry forward the global min up to this W
            prev = running.get(W - 1)
            if prev is not None and prev[0] < best:
                gmin, gE, gP = prev
            else:
                gmin, gE, gP = best, bestE, bestP
            running[W] = (gmin, gE, gP)
            print(f"  W={W:2d}: min rho* (spread==W) = {float(best):.6f}  "
                  f"| GLOBAL min (spread<=W) = {float(gmin):.6f} = {gmin}")
            sys.stdout.flush()
        gmin, gE, gP = running[k + 13]
        print(f"  --> k={k}: global min rho* over spread<=W = {gmin} "
              f"= {float(gmin):.6f}  at E={gE} P={gP}")
        sys.stdout.flush()

    print("\nDONE.")


if __name__ == "__main__":
    main()
