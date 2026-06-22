#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_reconcile_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

RECONCILE: my rho*(P,E) [maxgap>2/7 cap G_P]  vs  the criterion-C-via-Vmax
quantity [W(S\{Vmax})*7*Vmax > 1].

For S=[1,2,3,12,20,21,22,23,24,26,27,28,31] I found rho*_Vmax = 0 but the
criterion-C margin via Vmax = 1.3755 > 1 (so C HOLDS via Vmax!).  These cannot
both be "the via-Vmax route" unless they measure different things.  Resolve it.

THM-527 part A states the equivalence in the w0=Vmin->infinity LIMIT:
   good Vmax-period at x  <=>  x in G_P AND maxgap{frac(e_i x)} > 2/7,
where e_i = Vmax - u_i are co-offsets and G_P uses the SMALL part P only.
BUT in THM-527 the cluster is L (all u>13) and P = S cap {1..13}; the co-offsets
e_i run over the CLUSTER L, and G_P over P.  For my rho*=0 test I used
E = co-offsets of L and P = the small part -- the SAME setup.  So why does the
criterion-C margin (a DIFFERENT, finite, whole-circle computation) say C holds?

THREE possible resolutions to check exactly:
  (R1) The criterion-C margin uses W(S\{Vmax}) = widest safe arc of ALL OTHER
       runners (the 12 runners != Vmax), INCLUDING the cluster runners. The
       slow-fast rho* SEPARATES P (slow, via G_P) from the cluster co-offsets
       (fast, via maxgap). The maxgap>2/7 condition is the criterion's
       requirement that the WIDEST gap exceeds 2/7 AFTER fast-phase folding.
       The C-margin counts the widest arc anywhere; rho* counts the MEASURE of
       slow-times whose fast-folded gap >2/7. These differ: C can hold (one wide
       arc exists) while rho*=0 (that wide arc sits in a v-period but the
       maxgap>2/7 SLOW-density is 0 because the wide arc is NOT at a slow-safe x).
  (R2) The finite-Vmax good-period count can be >0 even when the w0->inf density
       rho*=0 (the limit washes out a thin but real set). i.e. rho*=0 is the
       w0->inf LIMIT; the actual finite-Vmax S may still have a witness.
  (R3) My maxgap>2/7 might be the wrong threshold for the criterion (vs the
       global-witness threshold maxgap>1/7).

We compute, for the concrete S, ALL of:
  - rho*_Vmax (maxgap>2/7 cap G_P)         [the criterion density]
  - rho_glob (maxgap>1/7 cap G_P)          [the GLOBAL-witness density]
  - the ACTUAL good Vmax-periods at finite scale (does a real witness sit in a
    Vmax-safe gap?)  -> reconcile with L(S)>0 and M(S)>=1/14.
"""
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


def density(P, E, gapthr):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid) and circ_maxgap_at(E, mid) > gapthr:
            tot += (x1 - x0)
    return tot


def safe_set_arcs(S, thr=Fr(1, 14)):
    bps = {Fr(0), Fr(1)}
    for v in S:
        for m in range(0, v):
            for r in (1, 13):
                x = Fr(14 * m + r, 14 * v)
                if 0 < x < 1:
                    bps.add(x)
    pts = sorted(bps)
    arcs = []
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        ok = True
        for v in S:
            f = (Fr(v) * mid) % 1
            dd = f if f <= Fr(1, 2) else 1 - f
            if dd < thr:
                ok = False
                break
        if ok:
            if arcs and arcs[-1][1] == x0:
                arcs[-1] = (arcs[-1][0], x1)
            else:
                arcs.append((x0, x1))
    return arcs


def main():
    print("=" * 78)
    print("THM-527 Thread A: RECONCILE rho*(maxgap>2/7) vs criterion-C-via-Vmax")
    print("=" * 78)

    S = [1, 2, 3, 12, 20, 21, 22, 23, 24, 26, 27, 28, 31]
    Vmax = 31
    P = [p for p in S if p <= 13]
    L = [u for u in S if u > 13]
    E = sorted(Vmax - u for u in L)   # co-offsets of the cluster
    print(f"\nS = {S}")
    print(f"  P (small) = {P};  L (cluster>13) = {L}")
    print(f"  Vmax = {Vmax};  co-offsets E = Vmax-L = {E}")
    print(f"  NOTE: 0 in E? {0 in E}  (Vmax itself gives e=0)")

    # the densities
    r_crit = density(P, E, Fr(2, 7))
    r_glob = density(P, E, Fr(1, 7))
    print(f"\n  rho*_crit  (maxgap>2/7 cap G_P) = {r_crit} = {float(r_crit):.6f}")
    print(f"  rho*_glob  (maxgap>1/7 cap G_P) = {r_glob} = {float(r_glob):.6f}")

    # the ACTUAL global safe set of the FULL S (this is what M>=1/14 needs)
    arcs = safe_set_arcs(S)
    Lmeas = sum((b - a for a, b in arcs), Fr(0))
    print(f"\n  ACTUAL full-S lonely measure L(S) = {float(Lmeas):.6f}  "
          f"(# safe arcs = {len(arcs)})")
    if arcs:
        wa = max(arcs, key=lambda ab: ab[1] - ab[0])
        print(f"    widest full-S safe arc: ({float(wa[0]):.6f},{float(wa[1]):.6f}) "
              f"width {float(wa[1]-wa[0]):.6f}, center {float((wa[0]+wa[1])/2):.6f}")
        # which Vmax-gap does the witness live in?
        c = (wa[0] + wa[1]) / 2
        kk = int(Vmax * c)
        print(f"    witness center is in Vmax-period index k={kk} "
              f"(Vmax*center={float(Vmax*c):.4f})")

    print("\n  RESOLUTION:")
    print("  The criterion-C margin W(S\\Vmax)*7*Vmax uses the WIDEST safe arc of")
    print("  the OTHER 12 runners over the whole circle. rho*_crit measures the")
    print("  MEASURE of P-safe slow-times whose cluster fast-fold leaves a >2/7")
    print("  gap.  These differ when the rescuing arc is NOT P-safe-localised in")
    print("  the slow variable.  KEY CHECK: is the full-S witness arc P-safe?")
    if arcs:
        c = (wa[0] + wa[1]) / 2
        print(f"    full-S witness center {float(c):.6f}: in G_P? {in_GP(P, c)}")
        # the cluster maxgap at that center:
        g = circ_maxgap_at(E, c)
        print(f"    cluster maxgap{{e_i*center}} = {float(g):.6f} (>2/7? {g>Fr(2,7)})")

    # the decisive statement:
    print("\n" + "=" * 78)
    print("  So rho*_crit=0 means: NO P-safe slow-time has cluster-gap > 2/7.")
    print("  But the full-S witness exists because it sits at a slow-time that is")
    print("  P-safe with cluster-gap in (1/7, 2/7] -- a GLOBAL witness (gap>1/7),")
    print("  not a via-Vmax CRITERION witness (which needs >2/7).  i.e. the")
    print("  via-Vmax CRITERION-C-density rho*_crit can be 0 while the GLOBAL")
    print("  witness density rho*_glob > 0 carries M(S)>=1/14.")
    print(f"    rho*_glob here = {float(r_glob):.6f} > 0 ? {r_glob>0}")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
