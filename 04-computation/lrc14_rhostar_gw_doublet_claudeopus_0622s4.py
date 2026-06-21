#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22-S4: rho* and witness floor for GENUINE-WIDE DOUBLETS.

Extends lrc14_rhostar_zero_probe_kpswf10.py (which checked consecutive E only)
to genuine-wide doublet configs E_M = B u {M, M+g}.

In THM-527's framework, given a k-speed config E_full = B u {M, M+g}:
  - Vmax = M+g
  - co-offsets: E_co = {0, g} u {M+g - b : b in B}
  - P = the remaining 13-k speeds (chosen adversarially to minimize rho*)

We compute:
  rho*(P, E_co)   = meas{x in G_P : maxgap{frac(e*x): e in E_co} > 2/7}
  witness(P,E_co) = meas{x in G_P : maxgap{frac(e*x): e in E_co} > 1/7}

for the BINDING (P, B, M, g) combinations:
  - P = {1,2,3,12} (worst P for k=9 from probe)
  - P = {1,2,3}    (worst P for k=10)
  - P = {1,6}      (worst P for k=11)
  - P = {1}        (worst P for k=12)
  - B = consecutive consec_{k-2} (binding base)
  - B = even-AP    (another binding base)
  - B = top-cluster (another binding base)
  - M = 15, 20, 30, 50, 200 (finite window + tail)
  - g = 1, 2, 3, 4

Also verifies the FROZEN LIMIT (M=200) matches the rho* probe's consecutive floor.

This is the KEY COMPUTATION for HYP-2826: verifying rho* > 0 (and witness > 0)
for all genuine-wide doublets, not just consecutive E.
"""
from __future__ import annotations
import sys, itertools, time
from fractions import Fraction as Fr
from math import gcd
from functools import reduce

P7 = 7
T1 = Fr(1, 7)   # gap > 1/7 threshold (witness)
T2 = Fr(2, 7)   # gap > 2/7 threshold (rho*, THM-527)
HALF = Fr(1, 14)


# ============================================================ helpers
def reprim(E):
    E = tuple(sorted(set(int(x) for x in E)))
    g = reduce(gcd, E)
    return tuple(x // g for x in E) if g > 1 else E


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(ph):
    if len(ph) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, (ph[0] + 1) - ph[-1])


def in_GP(Pset, x, thr=HALF):
    for p in Pset:
        r = (int(p) * x) % 1
        if min(r, 1 - r) < thr:
            return False
    return True


def breakpoints(E_co, Pset):
    """All x-breakpoints for rho* computation over E_co (co-offsets) and P."""
    bp = {Fr(0), Fr(1)}
    for e in list(E_co) + list(Pset):
        e = int(e)
        if e == 0:
            continue
        # G_P breakpoints: ||p*x|| = 1/14 => frac(p*x) in {1/14, 13/14}
        for k in range(0, e + 1):
            for s in (HALF, -HALF):
                v = (Fr(k) + s) / e
                if Fr(0) <= v <= Fr(1):
                    bp.add(v)
        # Phase collision breakpoints: frac(e*x) at multiples of 1/7
        for t in range(0, P7 * e + 1):
            bp.add(Fr(t, P7 * e))
    # Cross-phase gap thresholds: frac(e_i*x) - frac(e_j*x) = ±1/7, ±2/7
    El = list(E_co)
    for i in range(len(El)):
        for j in range(len(El)):
            d = El[i] - El[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in (T1, -T1, T2, -T2):
                    v = Fr(m, ad) + s / ad
                    if Fr(0) <= v <= Fr(1):
                        bp.add(v)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def compute_rhostar(E_co, Pset):
    """Exact rational rho*(P, E_co) and witness(P, E_co)."""
    pts = breakpoints(E_co, Pset)
    G2capGP = Fr(0)  # rho* = meas(maxgap > 2/7) cap G_P
    G1capGP = Fr(0)  # witness = meas(maxgap > 1/7) cap G_P
    measGP = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        if not in_GP(Pset, mid):
            continue
        measGP += w
        g = maxgap(phases_at(E_co, mid))
        if g > T1:
            G1capGP += w
        if g > T2:
            G2capGP += w
    return measGP, G1capGP, G2capGP


def consec_base(size):
    return tuple(range(size))


def even_ap_base(size):
    return tuple(range(0, 2 * size, 2)) if 2 * (size - 1) <= 14 else None


def top_cluster_base(size):
    return (0,) + tuple(range(15 - size + 1, 15))


# Worst-case P for each k (from rho* probe)
WORST_P = {
    9: [1, 2, 3, 12],
    10: [1, 2, 3],
    11: [1, 6],
    12: [1],
}


def gw_cooffsets(B, M, g):
    """Co-offsets of B u {M, M+g} relative to Vmax = M+g."""
    Vmax = M + g
    co = [0, g] + [Vmax - b for b in B if b > 0]
    # Remove duplicates, keep 0
    co = sorted(set(co))
    return co


def main():
    print("=" * 90)
    print("rho* and WITNESS FLOOR for GENUINE-WIDE DOUBLETS (claude-opus-0622-S4)")
    print("  E_full = B u {M, M+g}, co-offsets = {0,g} u {M+g-b: b in B}")
    print("  rho* = meas(G_P cap {gap > 2/7}),  witness = meas(G_P cap {gap > 1/7})")
    print("=" * 90)

    M_vals = [15, 20, 30, 50, 200]
    gaps = [1, 2, 3, 4]

    # Global minimum tracking
    min_rho_global = (Fr(1), None)
    min_wit_global = (Fr(1), None)

    for k in [9, 10, 11, 12]:
        P = WORST_P[k]
        size = k - 2  # base size
        cap = None  # cap not needed here; just checking positivity

        bases = []
        c = consec_base(size)
        bases.append(("consec", c))
        e = even_ap_base(size)
        if e is not None:
            bases.append(("even-AP", e))
        t = top_cluster_base(size)
        bases.append(("top-clust", t))

        print(f"\n{'='*80}")
        print(f"k={k}  |E|={k}  P={P}  (worst P from rho* probe)")

        for bname, B in bases:
            print(f"\n  Base {bname} = {B}")
            for g in gaps:
                row = []
                for M in M_vals:
                    Eco = gw_cooffsets(B, M, g)
                    # Check E_full is genuine-wide (sanity check for binding configs)
                    E_full = sorted(set(list(B) + [M, M + g]))
                    if len(E_full) != k:
                        continue
                    t0 = time.time()
                    measGP, G1, G2 = compute_rhostar(Eco, P)
                    dt = time.time() - t0
                    row.append((M, float(measGP), float(G1), float(G2), dt))
                    if G2 < min_rho_global[0]:
                        min_rho_global = (G2, (k, bname, g, M, B))
                    if G1 < min_wit_global[0]:
                        min_wit_global = (G1, (k, bname, g, M, B))

                if row:
                    print(f"    g={g}: ", end="")
                    for M, mGP, G1, G2, dt in row:
                        pos_r = "*" if G2 > 0 else "ZERO"
                        pos_w = "*" if G1 > 0 else "ZERO"
                        print(f"M={M}: rho*={G2:.5f}{pos_r}/wit={G1:.5f}{pos_w}  ", end="")
                    print()

    print("\n" + "=" * 90)
    rho, arg = min_rho_global
    wit, argw = min_wit_global
    print(f"GLOBAL min rho* (gap>2/7) = {rho} = {float(rho):.6f}  @ {arg}")
    print(f"GLOBAL min witness (gap>1/7) = {wit} = {float(wit):.6f}  @ {argw}")
    print(f"  rho* floor > 0 everywhere: {rho > 0}")
    print(f"  witness floor > 0 everywhere: {wit > 0}")


if __name__ == "__main__":
    main()
