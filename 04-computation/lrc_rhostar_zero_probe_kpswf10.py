#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_zero_probe_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

*** rho* = 0 FOUND -- the decisive probe. ***

The subtorus scan found  E=[0,3,4,5,7,8,9,10,11], P=[1,2,3,12]  with EXACT
rho* = 0.  This CONTRADICTS THM-527 part D's claim "rho* positive everywhere
tested / no rho*=0".  We must (a) CONFIRM it exactly, (b) understand WHY (the
cluster fills all gaps <=2/7 on ALL of G_P), (c) determine whether it is a real
obstruction to LRC(14) or excluded by the other reductions (k<=2, criterion C
is only SUFFICIENT, the offset shape may not be a genuine co-offset set, etc.).

This is HONEST DUE DILIGENCE: a rho*=0 shape means the Vmax-criterion-C route
gives NO witness for that S; LRC(14) for that S must then come from a DIFFERENT
v (not Vmax) or a different mechanism.  We report precisely what it means.
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


def in_GP(P, x, thr=Fr(1, 14)):
    for p in P:
        f = (Fr(p) * x) % 1
        dd = f if f <= Fr(1, 2) else 1 - f
        if dd < thr:
            return False
    return True


def rho_star_exact(P, E, verbose=False):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    total = Fr(0)
    goodcells = []
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid) and circ_maxgap_at(E, mid) > Fr(2, 7):
            total += (x1 - x0)
            goodcells.append((x0, x1))
    if verbose:
        return total, goodcells
    return total


def maxgap_on_GP_max(P, E):
    """On the cells of G_P, what is the MAX of maxgap? If <= 2/7 everywhere in
       G_P, then rho*=0. Report the max and where."""
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    best = Fr(0)
    bestx = None
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid):
            g = circ_maxgap_at(E, mid)
            if g > best:
                best = g
                bestx = mid
    return best, bestx


def main():
    print("=" * 78)
    print("THM-527 Thread A: rho* = 0 PROBE (the decisive due-diligence)")
    print("=" * 78)

    E = [0, 3, 4, 5, 7, 8, 9, 10, 11]
    P = [1, 2, 3, 12]
    print(f"\nE = {E}  (k={len(E)})")
    print(f"P = {P}  (|P|={len(P)})")

    # (a) confirm exactly
    r, cells = rho_star_exact(P, E, verbose=True)
    print(f"\n(a) EXACT rho* = {r}  ({'ZERO' if r == 0 else float(r)})")
    print(f"    # GOOD cells inside G_P: {len(cells)}")

    # (b) max maxgap over G_P
    g, gx = maxgap_on_GP_max(P, E)
    print(f"\n(b) max over x in G_P of maxgap{{e_i x}} = {g} = {float(g):.6f}")
    print(f"    attained at x = {gx} = {float(gx):.6f}")
    print(f"    2/7 = {Fr(2,7)} = {float(Fr(2,7)):.6f}")
    print(f"    => on ALL of G_P the cluster max-gap is <= 2/7: "
          f"{'YES (rho*=0 confirmed)' if g <= Fr(2,7) else 'NO'}")

    # show the phase pattern at the max-gap point
    ph = sorted(set((Fr(e) * gx) % 1 for e in E))
    print(f"    phases at that x: {[str(p) for p in ph]}")

    # (c) is E a genuine co-offset shape?  E = Vmax - u for the cluster u>13.
    # co-offsets are e_i = Vmax - u_i with u_i > 13 distinct, so e_i are
    # DISTINCT nonneg integers with 0 = e (the Vmax tooth) and all e_i < Vmax-13.
    # The shape [0,3,4,5,7,8,9,10,11] has spread 11 => Vmax - min(u) = 11, and the
    # cluster u_i = Vmax - e_i must ALL exceed 13.  That needs Vmax - 11 > 13 =>
    # Vmax > 24.  Possible. BUT the small part P and the cluster together must
    # be a PRIMITIVE COVERING 13-set in case S3.  Check the structural gates:
    print("\n(c) Is this shape an admissible S3 co-offset set?")
    print(f"    spread = {max(E)}; needs Vmax > 13 + spread = {13+max(E)} for u_i>13.")
    print(f"    |P|+k = {len(P)}+{len(E)} = {len(P)+len(E)} (must be 13 for a 13-set).")
    print(f"    => |S| = {len(P)+len(E)} {'== 13 OK' if len(P)+len(E)==13 else '!= 13'}")
    # gcd / primitivity of the would-be speed set is determined by actual values,
    # not offsets; the offset reformulation is blind to the absolute Vmax. The
    # criterion-C route is SUFFICIENT only: rho*=0 means NO Vmax-period witness.

    # (d) does a NEARBY shape have rho*>0? (robustness of the zero)
    print("\n(d) nearby shapes (perturb one offset) -- is rho*=0 isolated or a basin?")
    cnt0 = 0
    tested = 0
    for i in range(1, len(E)):
        for delta in (-2, -1, 1, 2):
            E2 = list(E)
            E2[i] += delta
            if len(set(E2)) != len(E2) or min(E2) < 0 or E2 != sorted(E2):
                continue
            r2 = rho_star_exact(P, E2)
            tested += 1
            if r2 == 0:
                cnt0 += 1
    print(f"    perturbations tested: {tested}; also rho*=0: {cnt0}")

    # (e) how many (P, E) with rho*=0 exist? scan moderate shapes for the worst P
    print("\n(e) census: how common is rho*=0 over k=9 shapes (spread<=14) x worst P?")
    worstP = [[1, 2, 3, 12], [1, 6, 7, 13], [1, 2, 3, 13]]
    z = 0
    tot = 0
    examples = []
    for tail in itertools.combinations(range(1, 15), 8):
        E2 = [0] + list(tail)
        for Pp in worstP:
            r2 = rho_star_exact(Pp, E2)
            tot += 1
            if r2 == 0:
                z += 1
                if len(examples) < 6:
                    examples.append((E2, Pp))
    print(f"    rho*=0 cases: {z} / {tot} (k=9, spread<=14, 3 worst P)")
    for E2, Pp in examples:
        # min-support
        print(f"      E={E2} P={Pp}")

    print("\n" + "=" * 78)
    print("VERDICT (honest):")
    print("  rho*=0 shapes EXIST.  THM-527 part D's 'no rho*=0' claim is FALSE.")
    print("  Meaning: for such (P,E) the Vmax-criterion-C route yields NO witness;")
    print("  the compactness-floor program (inf rho* > 0) DOES NOT hold as posed.")
    print("  LRC(14) for these S must come from a DIFFERENT certificate (another")
    print("  v != Vmax, or these shapes are excluded as primitive covering sets).")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
