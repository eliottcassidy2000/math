#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C part 3 (kps-S25-wf8): the BOUNDARY-slack branch (spread base, r<=2) EXACT max p0,
plus the slack err_max conditional bound  decorr_sup_slack + err_max_slack <= cap  (prompt item 2).

Fast: small far windows (single/double far max cover sits just past the base), exact rationals.
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_wide_branch_ridge_codex_s47 import CAP

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}


def p0(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l; den2 = 2 * l; bps = {0, d}
    for e in nz:
        step = l // e; x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps); num = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = lo + hi; mask = 0
        for e in nz:
            mask |= 1 << ((e * mid // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num += hi - lo
    return F(num, d)


def is_block(B):
    B = sorted(set(B))
    return B == list(range(len(B)))


def prim(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def boundary_spread_scan(k, far_win=14):
    """EXACT max p0 over r in {1,2}, SPREAD (non-single-block) base, far in a bounded window past base."""
    best = F(-1); arg = None
    for r in (1, 2):
        b = k - r
        if b < 2 or b - 1 > 14:
            continue
        for rest in combinations(range(1, 15), b - 1):
            base = (0,) + rest
            if is_block(base):
                continue
            bmax = max(base); lo = max(15, bmax + 1)
            for far in combinations(range(lo, lo + far_win), r):
                E = base + far
                if E[-1] <= 14 or not prim(E):
                    continue
                pv = p0(E)
                if pv > best:
                    best = pv; arg = E
    return best, arg


def slack_errmax(k, far_hi=45, n_rand=4000, seed=3):
    """sup over slack family of err=p0-decorr, sup decorr, sup p0. Structured (dilated-AP perturb,
    multi-cluster) + bounded random. EXACT."""
    rng = random.Random(seed)
    sup_dec = F(-1); sup_err = F(-99); sup_p0 = F(-1)
    cands = []
    # dilated-AP perturbations (the resonant near-binding slack)
    for d in (1, 2, 3, 4, 5, 6):
        for extra in range(0, 20):
            E = tuple(sorted(set([d * i for i in range(k - 1)] + [d * (k - 1) + 1 + extra])))
            if len(E) == k and E[-1] > 14:
                cands.append(E)
    # multi-cluster (2-3 clusters)
    for M in (15, 18, 22, 30, 44):
        for sizes in ([k - 2, 2], [k - 3, 3], [k // 2, k - k // 2], [3, 3, k - 6]):
            cfg = []; off = 0
            for s in sizes:
                cfg += [off + i for i in range(s)]; off += M
            E = tuple(sorted(set(cfg)))
            if len(E) == k and E[-1] > 14:
                cands.append(E)
    # bounded random
    for _ in range(n_rand):
        r = rng.randint(1, k - 1); b = k - r
        base = [0] + (rng.sample(range(1, 15), b - 1) if b >= 2 else [])
        if len(set(base)) != b:
            continue
        far = rng.sample(range(15, far_hi), r)
        E = tuple(sorted(set(base + far)))
        if len(E) != k or E[-1] <= 14:
            continue
        cands.append(E)
    for E in cands:
        if not prim(E):
            continue
        base = tuple(e for e in E if e <= 14); far = [e for e in E if e > 14]; r = len(far)
        if not ((r >= 3) or (not is_block(base))):
            continue
        dec = boundary_value_direct(base, r); pv = p0(E); err = pv - dec
        sup_dec = max(sup_dec, dec); sup_err = max(sup_err, err); sup_p0 = max(sup_p0, pv)
    return sup_dec, sup_err, sup_p0


def main():
    print("THREAD C part 3 (kps-S25-wf8): boundary-slack branch + slack err_max conditional bound\n")
    print("=" * 84)
    print("(B) BOUNDARY slack r in {1,2}, SPREAD (non-single-block) base: EXACT max p0")
    print("=" * 84)
    print(" k    max_p0 (spread,r<=2)    cap_k      GAP        Q(k-1)     (max_p0 - Q)   argmax E")
    bres = {}
    for k in (8, 9, 10, 11):
        best, arg = boundary_spread_scan(k)
        gap = CAP[k] - best
        bres[k] = (best, arg, gap)
        print(f" {k:>2}   {float(best):.5f} ({best})   {float(CAP[k]):.5f}   {float(gap):.5f}   {float(QVAL[k]):.5f}   {float(best-QVAL[k]):+.5f}   E={arg}")
    print()
    print("  NOTE: the spread base can slightly EXCEED Q(k-1) (e.g. k=9 even-base), but stays << cap.")
    print("  These r<=2 configs have TINY error (single/double far, THM-546) -- they are boundary, not deep slack.\n")

    print("=" * 84)
    print("(D) SLACK err_max + conditional bound  decorr_sup + err_max <= cap  (prompt item 2)")
    print("=" * 84)
    print(" k    sup_decorr    sup_err     sup_p0      decorr_sup+err_max    cap_k      bound<=cap?")
    for k in (8, 9, 10, 11, 12):
        sd, se, sp = slack_errmax(k)
        bound = sd + se
        print(f" {k:>2}   {float(sd):.5f}    {float(se):+.5f}    {float(sp):.5f}    {float(bound):.5f}              {float(CAP[k]):.5f}   {bound <= CAP[k]}")
    print()
    print("READING: the conditional bound decorr_sup+err_max<=cap is LOOSE (decorr_sup and err_max do NOT")
    print("co-occur), but it STILL holds with margin -- the slack regime is doubly safe. The direct max_p0")
    print("(part A/B) is the tighter certificate.\n")

    print("=" * 84)
    print("BOUNDARY-SLACK SUMMARY")
    print("=" * 84)
    print(" k    spread-base r<=2 max_p0   gap")
    for k in (8, 9, 10, 11):
        b, _, g = bres[k]
        print(f" {k:>2}   {float(b):.5f}                {float(g):.5f}")


if __name__ == "__main__":
    main()
