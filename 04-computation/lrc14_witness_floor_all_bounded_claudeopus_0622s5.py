#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22-S5: WITNESS FLOOR for ALL BOUNDED E (span<=14), k=8..13.

HYP-2830: Extends the kpswf10 admissible probe (consecutive E only) to check
G2(P,E) = meas{x in G_P : maxgap{frac(e*x): e in E} > 1/7} >= m_P for ALL
BOUNDED PRIMITIVE k-speed sets E (span<=14) at k=8..13.

This is the key missing computation: the kpswf10 probe and S4 (gw doublets)
together suggest that consecutive E minimizes G2, but this has not been verified
over all bounded E.

PHASE 1: For each k=8..13, find the worst-witness P (min G2 over all P subsets
         of {1..13} of size 13-k) for CONSECUTIVE E={0..k-1}.

PHASE 2: For each k=8..13, fix the worst-witness P from Phase 1. Check G2>=m_P
         for ALL PRIMITIVE BOUNDED E (span<=14) with |E|=k.

m_P = 14249/252252 ~ 0.0565 (proved admissible witness floor, THM-530, kps-S29).

If PASS: witness floor G2>=m_P holds for all (worst-witness-P, bounded E) pairs,
strongly supporting that consecutive E is the binding witness case and closing
the computational gap for the bounded E regime.
"""
from __future__ import annotations
import itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import time, sys

P7 = 7
T1 = Fr(1, 7)    # witness threshold
HALF = Fr(1, 14)
m_P = Fr(14249, 252252)  # proved admissible floor (THM-530)


def primitive(E):
    if len(E) < 2:
        return len(E) == 1 and E[0] != 0
    diffs = [E[i] - E[0] for i in range(1, len(E))]
    return reduce(gcd, diffs) == 1


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(ph):
    if not ph:
        return Fr(1)
    if len(ph) == 1:
        return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, ph[0] + 1 - ph[-1])


def in_GP(Pset, x):
    for p in Pset:
        r = (int(p) * x) % 1
        if min(r, 1 - r) < HALF:
            return False
    return True


def grid(E, Pset):
    """All breakpoints for the G2 computation."""
    bp = {Fr(0), Fr(1)}
    for e in list(E) + list(Pset):
        e = int(e)
        if e == 0:
            continue
        # G_P constraint breakpoints: ||e*x|| = 1/14
        for ki in range(0, e + 1):
            for s in (HALF, -HALF):
                v = (Fr(ki) + s) / e
                if Fr(0) <= v <= Fr(1):
                    bp.add(v)
        # Phase collision at multiples of 1/7
        for t in range(0, P7 * e + 1):
            bp.add(Fr(t, P7 * e))
    # Cross-phase gap = 1/7 breakpoints
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d == 0:
                continue
            for m in range(0, d + 1):
                for s in (T1, -T1):
                    v = Fr(m, d) + s / d
                    if Fr(0) <= v <= Fr(1):
                        bp.add(v)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def measure_G2(E, Pset):
    """Exact rational meas(G_P) and G2 = meas{G_P ∩ {maxgap > 1/7}}."""
    pts = grid(E, Pset)
    mGP = Fr(0)
    G2 = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        if not in_GP(Pset, mid):
            continue
        mGP += w
        if maxgap(phases_at(E, mid)) > T1:
            G2 += w
    return mGP, G2


def all_bounded_primitive_E(k, maxspeed=14):
    """All primitive k-speed subsets of {0,...,maxspeed} containing 0."""
    for combo in itertools.combinations(range(1, maxspeed + 1), k - 1):
        E = (0,) + combo
        if primitive(E):
            yield E


def consec_E(k):
    return list(range(k))


# ============================================================
def phase1_find_worst_witness_P(k, verbose=False):
    """Find P ⊂ {1..13}, |P|=13-k, minimizing G2 for consecutive E."""
    psize = 13 - k
    E = consec_E(k)
    if psize == 0:
        mGP, G2 = measure_G2(E, [])
        return [], G2
    min_wit = (Fr(2), None)
    n = 0
    for Pset in itertools.combinations(range(1, 14), psize):
        _, G2 = measure_G2(E, list(Pset))
        if G2 < min_wit[0]:
            min_wit = (G2, list(Pset))
        n += 1
        if verbose and n % 100 == 0:
            print(f"    P1 k={k} {n} P-sets scanned, current worst G2={float(min_wit[0]):.5f}", flush=True)
    return min_wit[1], min_wit[0]


def phase2_check_all_bounded(k, worst_P, verbose=True):
    """Check G2(worst_P, E) >= m_P for all bounded primitive E with |E|=k."""
    n_configs = 0
    n_fail = 0
    min_G2 = (Fr(1), None)
    t0 = time.time()

    for E in all_bounded_primitive_E(k):
        n_configs += 1
        _, G2 = measure_G2(list(E), worst_P)
        if G2 < min_G2[0]:
            min_G2 = (G2, E)
        if G2 < m_P:
            n_fail += 1
            if verbose:
                print(f"  *** FAIL k={k} E={E} G2={float(G2):.6f} < m_P ***", flush=True)
        if verbose and n_configs % 500 == 0:
            dt = time.time() - t0
            print(f"  k={k} {n_configs} configs, min_G2={float(min_G2[0]):.5f}, {dt:.1f}s", flush=True)

    return n_configs, n_fail, min_G2


# ============================================================
def main():
    print("=" * 90)
    print("WITNESS FLOOR G2 >= m_P for ALL BOUNDED E (span<=14), k=8..13")
    print(f"  m_P = {m_P} ~ {float(m_P):.6f} (THM-530)")
    print(f"  G2 = meas{{x in G_P : maxgap > 1/7}}")
    print("=" * 90)

    # k=8..13: number of P-subsets
    # k=8: C(13,5)=1287, k=9: C(13,4)=715, k=10: C(13,3)=286,
    # k=11: C(13,2)=78, k=12: C(13,1)=13, k=13: 1
    p_counts = {8: 1287, 9: 715, 10: 286, 11: 78, 12: 13, 13: 1}

    global_min = (Fr(1), None)
    worst_Ps = {}

    print("\n--- PHASE 1: Find worst-witness P for consecutive E ---")
    for k in range(8, 14):
        t0 = time.time()
        print(f"k={k}: scanning {p_counts[k]} P-subsets ...", flush=True)
        wP, min_G2 = phase1_find_worst_witness_P(k)
        dt = time.time() - t0
        worst_Ps[k] = wP
        consec_G2 = min_G2
        print(f"  k={k} worst P={wP}  min witness(consec)={min_G2} ~ {float(min_G2):.6f}  [{dt:.1f}s]")
        if min_G2 < global_min[0]:
            global_min = (min_G2, (k, 'consec', wP))

    print(f"\n  Phase 1 global min G2 (consec E) = {global_min[0]} ~ {float(global_min[0]):.6f} @ {global_min[1]}")
    print(f"  G2 >= m_P for all consec E: {global_min[0] >= m_P}")

    print("\n--- PHASE 2: Check all bounded primitive E with worst-witness P ---")
    global_min2 = (Fr(1), None)

    for k in range(8, 14):
        wP = worst_Ps[k]
        t0 = time.time()
        print(f"\nk={k} worst P={wP}", flush=True)
        n_configs, n_fail, (min_G2, worst_E) = phase2_check_all_bounded(k, wP, verbose=True)
        dt = time.time() - t0
        status = "PASS" if n_fail == 0 else f"FAIL({n_fail})"
        print(f"  k={k}  {status}  configs={n_configs}  min G2={min_G2} ~ {float(min_G2):.6f}  [{dt:.1f}s]")
        print(f"  worst E={worst_E}  min G2 / m_P = {float(min_G2/m_P):.3f}x")
        if min_G2 < global_min2[0]:
            global_min2 = (min_G2, (k, worst_E))

    print("\n" + "=" * 90)
    print(f"PHASE 2 GLOBAL min G2 (all bounded E, worst-witness P) = {global_min2[0]} ~ {float(global_min2[0]):.6f}")
    print(f"  @ {global_min2[1]}")
    print(f"  min G2 / m_P = {float(global_min2[0]/m_P):.3f}x")
    print(f"  G2 >= m_P everywhere: {global_min2[0] >= m_P}")
    print()
    print("INTERPRETATION:")
    print("  If PASS: the witness floor G2>=m_P holds for ALL (worst-witness-P, bounded E)")
    print("  pairs at k=8..13.  Combined with k<=7 pigeonhole (G2=meas(G_P)>=m_P) and")
    print("  S4 gw-doublet check (witness~0.48>>m_P), this covers ALL LRC(14) config families.")
    print("  (Assumes worst-witness P does not change when E changes from consecutive.)")


if __name__ == "__main__":
    main()
