#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_good_arc_count_kps.py  (kind-pasteur 2026-06-22, THM-527 Part A residual)

Count #arcs = connected components of GOOD(E) = {x in [0,1): maxgap{frac(e x)} > 1/7}.

This is the quantity in THM-527 Part A's finite-Vmax correction
rho_K = rho* + O(#arcs / Vmax).  The convergence-node question (kps-S30 -> mac-mini):
is #arcs BOUNDED for the BINDING clusters (consec, single-far), so #arcs/Vmax -> 0
and Part A closes effectively from the uniform floor delta_k?

THM-527 Part C predicts GOOD(consec) = neighborhoods of {a/b : b<=6}, a BOUNDED set
of arcs (independent of k beyond the b<=6 count).  We verify #arcs vs k for consec,
single-far, and a wide cluster, to see the scaling.
"""
import sys
from fractions import Fraction as Fr


def circ_maxgap(E, x):
    ph = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(ph) <= 1:
        return Fr(1)
    g = max((b - a for a, b in zip(ph, ph[1:])), default=Fr(0))
    return max(g, (ph[0] + 1) - ph[-1])


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


def good_arc_count_and_measure(E, gapthr=Fr(1, 7)):
    """Return (#arcs, measure) of GOOD(E)={maxgap>1/7}, exact."""
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E))
    meas = Fr(0)
    arcs = 0
    prev_good = False
    # we also need to handle wrap (0 and 1 identified): treat [0,1) linearly but
    # merge first and last interval if both good (circular component).
    first_good = None
    last_good = None
    seg_goods = []
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0:
            continue
        g = circ_maxgap(E, (x0 + x1) / 2) > gapthr
        seg_goods.append(g)
        if g:
            meas += (x1 - x0)
    # count maximal runs of True in seg_goods, with circular wrap merge
    n = len(seg_goods)
    if n == 0:
        return 0, Fr(0)
    if all(seg_goods):
        return 1, meas   # whole circle
    # count linear runs
    runs = 0
    for i in range(n):
        if seg_goods[i] and (i == 0 or not seg_goods[i - 1]):
            runs += 1
    # circular merge: if first and last segment both good, they are one arc
    if seg_goods[0] and seg_goods[-1]:
        runs -= 1
    return runs, meas


def main():
    print("=" * 72)
    print("THM-527 Part A: #arcs of GOOD(E)={maxgap>1/7} vs cluster type")
    print("=" * 72)
    print("\n--- CONSEC clusters E={0,1,...,k-1} (the binding family) ---")
    print("   k   #arcs   measure(GOOD)=nu   (Part C predicts bounded #arcs)")
    for k in range(8, 14):
        E = list(range(k))
        na, m = good_arc_count_and_measure(E)
        print(f"  {k:3d}   {na:4d}    {float(m):.5f} = {m}")
        sys.stdout.flush()

    print("\n--- SINGLE-FAR clusters E={0,1,...,k-2, F} (tight base + one far) ---")
    print("   k   F    #arcs   measure(GOOD)")
    for k in (9, 10, 11):
        for F in (k + 5, k + 20, k + 60):
            E = list(range(k - 1)) + [F]
            na, m = good_arc_count_and_measure(E)
            print(f"  {k:3d} {F:4d}  {na:4d}    {float(m):.5f}")
        sys.stdout.flush()

    print("\n--- WIDE clusters (spread ~ several*k) ---")
    print("   E                                  #arcs   measure(GOOD)")
    wides = [
        [0, 3, 7, 12, 18, 25, 33, 42],          # k=8 increasing gaps
        [0, 5, 11, 17, 23, 29, 35, 41, 47],     # k=9 AP-ish wide
        [0, 2, 9, 16, 23, 30, 37, 44, 51, 58],  # k=10 wide
    ]
    for E in wides:
        na, m = good_arc_count_and_measure(E)
        print(f"  {str(E):34s} {na:4d}    {float(m):.5f}")
        sys.stdout.flush()

    print("\n" + "=" * 72)
    print("  READING: if #arcs is BOUNDED (small, ~const) for consec & single-far")
    print("  (the binding clusters), then #arcs/Vmax -> 0 and THM-527 Part A's")
    print("  finite-Vmax correction is uniform from the floor delta_k => Part A")
    print("  closes effectively for the binding cases. Wide clusters have larger")
    print("  #arcs but also larger nu (safe), so the threshold #arcs/delta is the")
    print("  quantity to watch.")
    print("=" * 72)
    print("\nDONE.")


if __name__ == "__main__":
    main()
