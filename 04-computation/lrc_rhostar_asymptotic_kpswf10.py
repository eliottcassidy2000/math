#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_asymptotic_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE ASYMPTOTIC MECHANISM: along a spread -> infinity family, does rho* -> a
POSITIVE limit (a subtorus floor) or -> 0 (compactness FAILS)?

This is THE decisive test for the compactness proof.  Part D of THM-527 claims
bounded spread; the spread-floor search showed the min of rho* still creeping
down with the spread cap.  We now isolate the asymptotics directly.

THREE families at fixed k, with one or more offsets pushed to infinity:

  (F1) ONE outlier: E = {0,1,...,k-2, M}, M -> infinity.
       Relation-free in the limit (M huge, coprime to the rest): the phase
       e_{k}x = Mx decorrelates from the others -> the orbit fills the product
       (sub)torus.  rho* -> meas( GOOD_{k} on the product ) cap G_P -- a clean
       limit.

  (F2) TWO far blocks: E = {0,1,2} U {M, M+1, M+2}, M -> infinity.
       Two correlated triples, far apart: the relation lattice has the
       within-block relations only; the orbit fills a 4-dim subtorus
       (3 within + the block separation phase).  rho* -> subtorus integral.

  (F3) AP with growing step: E = {0, s, 2s, ..., (k-1)s}, s -> infinity.
       Pure scaling: meas{maxgap{j s x}>2/7} = meas{maxgap{j y}>2/7} = mu_k
       (substitute y = s x mod 1, measure-preserving) -- but G_P does NOT scale,
       so rho* = meas( consec-GOOD(y)/s-periodized cap G_P ).  As s->inf the
       GOOD set becomes s thin copies; intersect FIXED G_P -> mu_k * meas(G_P)
       (equidistribution of the fine comb against G_P).  POSITIVE floor.

For each family we compute rho* (EXACT) along increasing M (or s) and watch the
limit.  We also compute the PREDICTED product floor F(k)*cap and mu_k*cap to see
which floor the family approaches.

Conclusion target: identify whether ANY spread-growing family drives rho*->0
(=> the bounded-spread reduction is FALSE and compactness needs the subtorus
argument), or whether all families floor at a positive value (=> compactness
holds with an explicit floor).
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import comb


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


def rho_star(P, E, thr=Fr(1, 14), gapthr=Fr(2, 7)):
    bps = {Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E)
    pts = sorted(bps)
    total = Fr(0)
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for p in P:
            f = (Fr(p) * mid) % 1
            d = f if f <= Fr(1, 2) else 1 - f
            if d < thr:
                ok = False
                break
        if ok and circ_maxgap_at(E, mid) > gapthr:
            total += (x1 - x0)
    return total


def meas_GP(P, thr=Fr(1, 14)):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P))
    total = Fr(0)
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
            total += (x1 - x0)
    return total


def F_iid(k):
    """P(k uniform iid pts have maxgap > 2/7) = sum_j (-1)^{j+1} C(k,j)(1-2j/7)_+^{k-1}."""
    tot = Fr(0)
    for j in range(1, k + 1):
        base = 1 - Fr(2 * j, 7)
        if base <= 0:
            break
        tot += (-1) ** (j + 1) * comb(k, j) * base ** (k - 1)
    return tot


def main():
    print("=" * 78)
    print("THM-527 Thread A: ASYMPTOTIC -- does rho* -> positive limit or 0?")
    print("=" * 78)
    sys.stdout.flush()

    P = [1, 2, 3, 12]      # the canonical worst small-part (k=9 region)
    capP = meas_GP(P)
    print(f"\nfixed P={P}: meas(G_P)={capP}={float(capP):.5f}")
    print(f"F_iid(k) ceilings: " + ", ".join(f"k={k}:{float(F_iid(k)):.4f}" for k in range(7, 11)))

    # ---- F1: one outlier, M -> infinity (k=9: {0..7,M}) ----
    print("\n--- F1: E = {0,1,..,7, M}, M growing (k=9) ---")
    k = 9
    Fk = F_iid(k)
    pred_iid = Fk * capP
    print(f"    predicted product floor F({k})*cap = {float(Fk):.4f}*{float(capP):.4f} "
          f"= {float(pred_iid):.5f}")
    for M in [8, 9, 11, 13, 17, 23, 31, 47, 61, 101, 211, 401]:
        E = list(range(8)) + [M]
        r = rho_star(P, E)
        print(f"    M={M:4d}: rho* = {float(r):.6f}")
    sys.stdout.flush()

    # ---- F2: two far triples, M -> infinity (k=6: {0,1,2}U{M,M+1,M+2}) ----
    print("\n--- F2: E = {0,1,2} U {M,M+1,M+2}, M growing (k=6) ---")
    Pk = [1, 2, 3, 12]
    capk = meas_GP(Pk)
    for M in [3, 4, 6, 10, 16, 30, 61, 121, 241]:
        E = [0, 1, 2, M, M + 1, M + 2]
        r = rho_star(Pk, E)
        print(f"    M={M:4d}: rho* = {float(r):.6f}  (k=6)")
    print(f"    [F(6)*cap = {float(F_iid(6)*capk):.5f}]")
    sys.stdout.flush()

    # ---- F3: AP with growing step s (k=9: {0,s,..,8s}) ----
    print("\n--- F3: E = {0, s, 2s, ..., 8s}, s growing (k=9 AP) ---")
    print("    pure mu is CONSTANT (=mu_9) under scaling; rho* tests GOOD-comb vs G_P")
    mu9 = Fr(277, 980)
    pred_f3 = mu9 * capP
    print(f"    predicted floor mu_9*cap = {float(mu9):.4f}*{float(capP):.4f} = {float(pred_f3):.5f}")
    for s in [1, 2, 3, 5, 7, 11, 13, 17, 23, 31, 43]:
        E = [s * j for j in range(9)]
        r = rho_star(P, E)
        print(f"    s={s:3d}: rho* = {float(r):.6f}")
    sys.stdout.flush()

    # ---- F4: the WORST observed family extended: {0,1,3,5,6,7,9,13,14}-type ----
    # push the tail outliers of the k=9 minimizer further to see if rho* -> 0
    print("\n--- F4: extend the observed k=9 minimizer tail (does it -> 0?) ---")
    print("    base [0,1,3,5,6,7,9, a, a+1] with a growing (mimic the W-min family)")
    for a in [11, 13, 17, 23, 31, 47, 71, 101]:
        E = [0, 1, 3, 5, 6, 7, 9, a, a + 1]
        r = rho_star([1, 6, 7, 13], E)
        print(f"    a={a:4d}: rho* = {float(r):.6f}  P=[1,6,7,13]")

    print("\nNet read: if every family FLOORS at a positive value, the inf is")
    print("positive (the subtorus floors F(k)*cap, mu_k*cap, etc. are all > 0);")
    print("the only way compactness fails is a family with rho* -> 0.")
    print("\nDONE.")


if __name__ == "__main__":
    main()
