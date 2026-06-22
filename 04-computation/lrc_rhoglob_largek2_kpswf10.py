#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhoglob_largek2_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

Complete the rho*_glob floor at k=10,11,12,13 (the slow tail), with arc-based
intersection for speed. Confirm min rho*_glob > 0 and zeros=0.

Reuses: pure nu_13 = 477/1078 ~ 0.442 > 0 (from largek run). Here we intersect
with the worst-P to get the actual floor at the largest k.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
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
        for m in range(0, thr_den * d + 1):
            for s in (1, -1):
                v = Fr(thr_den * m + s, thr_den * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_arcs(E, gapthr=Fr(1, 7)):
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


def am(arcs):
    return sum((b - a for a, b in arcs), Fr(0))


def aim(A, B):
    i = j = 0
    t = Fr(0)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            t += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return t


def main():
    print("=" * 78)
    print("THM-527 Thread A: rho*_glob floor at k=10..13 (completing the hard tail)")
    print("=" * 78)
    sys.stdout.flush()
    for k in [10, 11, 12, 13]:
        npart = 13 - k
        W = k + 3
        Pall = []
        for P in (itertools.combinations(range(1, 14), npart) if npart > 0 else [()]):
            pa = gp_arcs(list(P))
            Pall.append((am(pa), list(P), pa))
        Pall.sort(key=lambda e: e[0])
        Pcands = Pall[:25] if npart > 0 else [(Fr(1), [], [(Fr(0), Fr(1))])]
        best = None
        bestarg = None
        zeros = 0
        ns = 0
        for tail in itertools.combinations(range(1, W + 1), k - 1):
            E = [0] + list(tail)
            ga = good_arcs(E)
            mu = am(ga)
            ns += 1
            for (mg, P, pa) in Pcands:
                lb = mu + mg - 1
                if best is not None and lb >= best:
                    continue
                rg = aim(ga, pa)
                if rg == 0:
                    zeros += 1
                if best is None or rg < best:
                    best = rg
                    bestarg = (E, P)
            if ns > 80000:
                break
        E, P = bestarg
        sp = max(E)
        cov = 0
        for Vmax in range(14 + sp, 14 + sp + 80):
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
        print(f"        argmin E={E} P={P}  admissible covering S: {cov}")
        sys.stdout.flush()
    print("\nDONE.")


if __name__ == "__main__":
    main()
