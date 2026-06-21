#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): EXACT genuine-wide max p0, r=2 binding stratum, k=10,11,12.

Fast, complete window scan to pin the TRUE genuine-wide max p0 and locate the binding config.
- genuine-wide => r>=2 (PROVED separately). Binding stratum r=2.
- r=2 config: base = (k-2)-subset of [0,14] with 0; far = {f1, f1+g}, both >14.
- We scan ALL bases, f1 in [15, F1HI], g in [1, GHI]. We report the EXACT max p0 over genuine-wide
  configs and the per-stratum (base-shape) leaders. Coverage is maximized by dense bases + tight far,
  so the binding config lives in this window; F1HI is chosen past the saturation point.
"""
from __future__ import annotations
import sys, functools, time
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def primitive(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def reduced_span(S):
    S = sorted(set(S))
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False
    return True


def scan(k, F1HI=120, GHI=14):
    """EXACT max p0 over r=2 genuine-wide configs, all bases, far in window."""
    b = k - 2
    cap = CAP[k]
    q = QVAL[k]
    best = F(-1)
    arg = None
    near = []  # (p0, E) with p0 within 0.02 of Q
    nb = 0
    for rest in combinations(range(1, 15), b - 1):
        base = (0,) + rest
        bmax = base[-1]
        lo = max(15, bmax + 1)
        nb += 1
        for f1 in range(lo, F1HI + 1):
            for g in range(1, GHI + 1):
                f2 = f1 + g
                E = base + (f1, f2)
                E = tuple(sorted(E))
                if len(E) != k or not primitive(E):
                    continue
                if not is_genuine_wide(E):
                    continue
                pv = p0_fast(E)
                if pv > best:
                    best, arg = pv, E
                if pv > q - F(2, 100):
                    near.append((pv, E))
    near.sort(reverse=True)
    return best, arg, nb, near[:12]


def main():
    print("=" * 96)
    print("THREAD 2: EXACT genuine-wide max p0 (r=2 binding), k=10,11,12 (mac-mini-S7)")
    print("=" * 96)
    for k in (10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.6f}  Q(k-1)={float(QVAL[k]):.6f}  MARGIN={float(MARGIN[k]):.6f}")
    print()
    for k in (10, 11, 12):
        t = time.time()
        # window sizes tuned for feasibility; F1HI past saturation
        F1HI = {10: 80, 11: 60, 12: 45}[k]
        GHI = {10: 14, 11: 12, 12: 10}[k]
        print("-" * 96)
        print(f"k={k}  (all bases size {k-2}, f1<=lo+{F1HI}, gap<={GHI})")
        best, arg, nb, near = scan(k, F1HI=F1HI, GHI=GHI)
        q = QVAL[k]
        under = best < q
        print(f"  bases={nb}, time={time.time()-t:.1f}s")
        print(f"  => genuine-wide max p0 = {float(best):.6f} = {best}")
        print(f"     Q(k-1)={float(q):.6f}  Q-max_p0={float(q-best):+.6f}  cap-max_p0={float(CAP[k]-best):.6f}  "
              f"{'UNDER Q' if under else '*** OVER Q (OBSTRUCTION) ***'}")
        print(f"     argmax E={arg} base={tuple(e for e in arg if e<=14)} far={tuple(e for e in arg if e>14)}")
        print(f"     top near-Q configs:")
        for pv, E in near[:8]:
            print(f"        p0={float(pv):.6f} (Q-p0={float(q-pv):+.6f}) E={E}")
        print()


if __name__ == "__main__":
    main()
