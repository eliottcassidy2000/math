#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhoglob_largek_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

DOES rho*_glob (gap>1/7, global-witness density) STAY POSITIVE for LARGE k?

At k=9 we found rho*_glob >= 419/1092 ~ 0.384, no zeros.  But the hard regime is
k >= 7 (where maxgap >= 1/k is NOT automatically > 1/7).  As k -> 13, 1/k -> 1/13
< 1/7, so gap>1/7 becomes genuinely binding.  We must confirm rho*_glob does not
collapse to 0 at k = 10,11,12,13.

We scan, per k, integer shapes (bounded spread + spread-growing families) x the
worst-P shortlist, looking for ANY rho*_glob = 0, and report the min rho*_glob.

We also separately verify the CONSECUTIVE-cluster pure value
   nu_k := meas{ maxgap{j x : j<k} > 1/7 }
which is the pure (no-G_P) global-witness measure -- the analogue of mu_k but at
the 1/7 threshold.  nu_k > 0 for all k<=13 would be the pure floor; then
intersect with G_P (which only removes a bounded amount).
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd, comb
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


def good_breaks(E, thr_den=7):
    """gap=+-1/7 crossings at x=(thr_den*m +- 1)/(thr_den*d); plus collisions."""
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
        # gap = +-1/7 crossings (for the 1/7 threshold)
        for m in range(0, thr_den * d + 1):
            for s in (1, -1):
                v = Fr(thr_den * m + s, thr_den * d)
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


def rho_glob(P, E, gapthr=Fr(1, 7)):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid) and circ_maxgap_at(E, mid) > gapthr:
            tot += (x1 - x0)
    return tot


def nu_pure(E, gapthr=Fr(1, 7)):
    bps = sorted({Fr(0), Fr(1)} | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if circ_maxgap_at(E, mid) > gapthr:
            tot += (x1 - x0)
    return tot


def gp_meas(P):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid):
            tot += (x1 - x0)
    return tot


def main():
    print("=" * 78)
    print("THM-527 Thread A: rho*_glob (gap>1/7) at LARGE k -- stays positive?")
    print("=" * 78)
    sys.stdout.flush()

    # pure consec nu_k at the 1/7 threshold
    print("\n--- pure consecutive nu_k = meas{maxgap{j x:j<k} > 1/7} ---")
    print("    (1/k>1/7 automatic for k<=7; binding for k>=8)")
    for k in range(3, 14):
        E = list(range(k))
        nu = nu_pure(E)
        print(f"    k={k:2d}: nu_k = {nu} = {float(nu):.6f}")
    sys.stdout.flush()

    # per-k min rho*_glob over bounded-spread shapes x worst-P
    print("\n--- per-k min rho*_glob (bounded spread + worst-P), looking for 0 ---")
    for k in range(7, 14):
        npart = 13 - k
        W = k + 4
        # worst-P shortlist
        Pall = []
        for P in (itertools.combinations(range(1, 14), npart) if npart > 0 else [()]):
            Pall.append((gp_meas(list(P)) if npart > 0 else Fr(1), list(P)))
        Pall.sort(key=lambda e: e[0])
        Pcands = [p for _, p in Pall[:30]] if npart > 0 else [[]]
        best = None
        bestarg = None
        zeros = 0
        ns = 0
        for tail in itertools.combinations(range(1, W + 1), k - 1):
            E = [0] + list(tail)
            ns += 1
            for P in Pcands:
                rg = rho_glob(P, E)
                if rg == 0:
                    zeros += 1
                if best is None or rg < best:
                    best = rg
                    bestarg = (E, P)
            if ns > 60000:
                break
        # admissibility of the argmin
        E, P = bestarg
        sp = max(E)
        cov = 0
        for Vmax in range(14 + sp, 14 + sp + 120):
            L = [Vmax - e for e in E]
            if len(set(L)) != len(L) or min(L) <= 13:
                continue
            S = sorted(set(P) | set(L))
            if len(S) != 13 or reduce(gcd, S) != 1:
                continue
            if all(any(v % q == 0 for v in S) for q in range(2, 15)):
                cov += 1
        print(f"  k={k:2d} (|P|={npart}, W={W}, {ns} shapes x{len(Pcands)} P): "
              f"min rho*_glob = {best} = {float(best):.6f}  zeros={zeros}")
        print(f"        argmin E={E} P={P}  (admissible covering S: {cov})")
        sys.stdout.flush()

    print("\n" + "=" * 78)
    print("If min rho*_glob > 0 with zeros=0 for ALL k=7..13, the global-witness")
    print("density floor is POSITIVE across the whole hard regime -> the corrected")
    print("compactness program (floor rho*_glob, not rho*_crit) is on solid ground.")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
