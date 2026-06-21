#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): EXACT genuine-wide r=2 max p0 over DENSE bases, k=10,11,12.

Why dense bases: p0(B u F) <= sum_{t<=r} prof_t(B) (deterministic), and the decorrelated cover
boundary_value_direct(B,r) is monotone in the base coverage. The base maximizing coverage (hence
the slack max p0) is the densest; SPREAD bases put profile mass at HIGH missing-count t, killed by
c_t(t,r)=0 for t>r. So the binding genuine-wide config has a DENSE base. We scan ALL bases with
max consecutive-diff <= D and verify the leader stabilizes as D grows (D=1 consec, D=2, D=3).

For each dense base we sweep far pair (f1, f1+g), f1 in [15, F1HI], g in [1, GHI], EXACT p0,
restricted to genuine-wide. Report the EXACT max and the binding config.

Also reports, per base, the decorrelated 2-far floor boundary_value_direct(B,2) and the deterministic
P_{<=2}(B), for the tail/structural argument.
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

from lrc14_wide_branch_ridge_codex_s47 import CAP, missed_distribution
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


def dense_bases(b, D):
    """All bases (0 + (b-1) from [1,14]) with max consecutive diff <= D."""
    out = []
    for rest in combinations(range(1, 15), b - 1):
        base = (0,) + rest
        diffs = [y - x for x, y in zip(base, base[1:])]
        if not diffs or max(diffs) <= D:
            out.append(base)
    return out


def scan_dense(k, D, F1HI, GHI):
    b = k - 2
    q = QVAL[k]
    best = F(-1)
    arg = None
    near = []
    bases = dense_bases(b, D)
    for base in bases:
        bmax = base[-1]
        lo = max(15, bmax + 1)
        for f1 in range(lo, F1HI + 1):
            for g in range(1, GHI + 1):
                E = tuple(sorted(base + (f1, f1 + g)))
                if len(E) != k or not primitive(E):
                    continue
                if not is_genuine_wide(E):
                    continue
                pv = p0_fast(E)
                if pv > best:
                    best, arg = pv, E
                if pv > q - F(15, 1000):
                    near.append((pv, E))
    near.sort(reverse=True)
    return best, arg, len(bases), near[:10]


def main():
    print("=" * 96)
    print("THREAD 2: EXACT genuine-wide r=2 max p0 over DENSE bases, k=10,11,12 (mac-mini-S7)")
    print("=" * 96)
    for k in (10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.6f}  Q(k-1)={float(QVAL[k]):.6f}  MARGIN={float(MARGIN[k]):.6f}")
    print()
    for k in (10, 11, 12):
        print("-" * 96)
        print(f"k={k}: dense-base leader vs density D (stabilization check)")
        q = QVAL[k]
        F1HI = {10: 120, 11: 90, 12: 70}[k]
        GHI = {10: 24, 11: 20, 12: 16}[k]
        prev = None
        for D in (1, 2, 3):
            t = time.time()
            best, arg, nb, near = scan_dense(k, D, F1HI, GHI)
            stable = (prev is not None and best == prev)
            print(f"  D={D}: bases={nb} max_p0={float(best):.6f}={best} arg={arg} "
                  f"(Q-p0={float(q-best):+.6f}) {'[stable]' if stable else ''} t={time.time()-t:.0f}s")
            prev = best
        print(f"  Q(k-1)={float(q):.6f}. top near-Q (D=3):")
        for pv, E in near[:6]:
            print(f"     p0={float(pv):.6f} (Q-p0={float(q-pv):+.6f}) E={E} base={tuple(e for e in E if e<=14)} far={tuple(e for e in E if e>14)}")
        print()


if __name__ == "__main__":
    main()
