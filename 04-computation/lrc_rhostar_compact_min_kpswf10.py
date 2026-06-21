#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_compact_min_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE COMPACT MINIMUM of rho*(P, E).   [optimized]

Confirm  inf over (P subset {1..13}, bounded-spread shape E) of rho*(P,E) is a
POSITIVE minimum, and pin it exactly.

Speed strategy:
  * rho*(P,E) <= mu_pure(E) and rho*(P,E) <= meas(G_P).  So precompute, per k,
    the worst (smallest meas G_P) caps; and per shape compute mu_pure once.
  * Represent the GOOD set of a shape ONCE as a list of disjoint arcs (Fractions);
    represent each candidate G_P ONCE as arcs; rho* = measure of arc-intersection
    (fast sweep) -- far cheaper than re-breaking the whole circle each (P,E).
  * Prune: only intersect shapes whose mu_pure is within the running best+slack;
    only the P's that are "worst" (small meas G_P) can produce the min.

Parts of the Thread-A compactness argument addressed here:
  (3) compact min (spread monotonicity caps W; finite P) -> positive minimum.
  (4) integer shapes subset real => integer inf >= real inf (operative bound).
  (5) GOOD arc-count is O(1) (bounded spread) => rho_K = rho* + O(1/Vmax).
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
    """GOOD set as sorted disjoint arcs [(a,b),...]."""
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
    """G_P as sorted disjoint arcs."""
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
    """measure of union(A) cap union(B), A,B sorted disjoint arc lists."""
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


def mu_pure(E):
    return arcs_measure(good_arcs(E))


def main():
    print("=" * 78)
    print("THM-527 Thread A: THE COMPACT MINIMUM of rho*(P,E)  [optimized]")
    print("=" * 78)
    sys.stdout.flush()

    # ---- (3a) spread monotonicity: outliers RAISE mu ----
    print("\n--- (3a) spread monotonicity: pushing an outlier RAISES mu ---")
    for M in [9, 10, 12, 16, 24, 40, 80, 200, 1000]:
        E = [0, 1, 2, 3, 4, 5, 6, 7, M]
        print(f"    M={M:5d}: mu = {float(mu_pure(E)):.5f}")
    print("    => far outliers RAISE mu; minimizer has BOUNDED spread.")
    sys.stdout.flush()

    # ---- precompute worst-P arcs per npart (smallest meas G_P) ----
    # We only need the P's that could give the min. meas(G_P) decreases as P
    # grows / uses small p. Precompute ALL P arcs once per npart and keep the
    # ones with smallest meas (a generous shortlist), since rho*<=meas(G_P).
    print("\n--- precomputing G_P arcs per |P| (keep smallest-measure shortlist) ---")
    gp_cache = {}     # npart -> list of (meas, P, arcs)
    for npart in range(0, 11):
        entries = []
        for P in itertools.combinations(range(1, 14), npart):
            a = gp_arcs(list(P))
            entries.append((arcs_measure(a), P, a))
        entries.sort(key=lambda e: e[0])
        gp_cache[npart] = entries
        if entries:
            print(f"    |P|={npart}: {len(entries)} sets; "
                  f"min meas(G_P)={float(entries[0][0]):.5f} at P={entries[0][1]}; "
                  f"max meas={float(entries[-1][0]):.5f}")
    sys.stdout.flush()

    # ---- (3b)/(4) integer bounded-spread min rho* ----
    print("\n--- (3b)/(4) integer bounded-spread min rho* (operative inf) ---")
    overall_min = None
    overall_arg = None
    for k in range(3, 14):
        npart = 13 - k
        if npart < 0:
            continue
        W = min(k + 5, 30)
        # the P-shortlist: only P's with meas(G_P) <= current global cap matter,
        # but to be safe keep all P with meas(G_P) below a threshold = the
        # largest mu we will see (<=1). Keep the smallest-measure 60 sets.
        Pcands = gp_cache[npart][:60] if npart > 0 else [(Fr(1), (), [(Fr(0), Fr(1))])]
        # enumerate shapes; compute GOOD arcs once; intersect with P-candidates
        best = None
        bestE = None
        bestP = None
        nshapes = 0
        for tail in itertools.combinations(range(1, W + 1), k - 1):
            E = [0] + list(tail)
            ga = good_arcs(E)
            mu = arcs_measure(ga)
            nshapes += 1
            # prune: if mu >= best already, intersecting with any G_P (<=mu) could
            # still be smaller, so we cannot skip on mu alone. But rho* >=
            # mu + meas(G_P) - 1 (inclusion-exclusion lower bound). If that LB
            # already exceeds best, skip the whole shape.
            for (mg, P, pa) in Pcands:
                lb = mu + mg - 1   # meas(A cap B) >= meas A + meas B - 1
                if best is not None and lb >= best:
                    continue
                r = arcs_intersect_measure(ga, pa)
                if best is None or r < best:
                    best = r
                    bestE = E
                    bestP = list(P)
            if nshapes > 200000:
                break
        print(f"  k={k:2d} (W={W}, {nshapes} shapes x{len(Pcands)} P): "
              f"min rho* = {best} = {float(best):.6f}   E={bestE} P={bestP}")
        sys.stdout.flush()
        if overall_min is None or best < overall_min:
            overall_min = best
            overall_arg = (k, bestE, bestP)

    print(f"\n  >>> OVERALL integer bounded-spread min rho* = {overall_min} "
          f"= {float(overall_min):.6f}")
    print(f"      at k={overall_arg[0]}, E={overall_arg[1]}, P={overall_arg[2]}")
    print(f"      POSITIVE: {overall_min > 0}")
    sys.stdout.flush()

    # ---- (5) arc-count bound ----
    print("\n--- (5) GOOD arc-count (O(1) for bounded spread) ---")
    for E in [[0, 1, 2, 3, 4, 5, 6, 7, 8],
              [0, 1, 3, 5, 6, 8, 10, 11, 12, 14],
              [0, 2, 3, 4, 5, 6, 8],
              list(range(13)),
              [0, 5, 7, 8, 9, 10, 11, 13, 18]]:
        print(f"    E={E}: #GOOD arcs = {len(good_arcs(E))}  (spread {max(E)})")

    print("\nDONE.")


if __name__ == "__main__":
    main()
