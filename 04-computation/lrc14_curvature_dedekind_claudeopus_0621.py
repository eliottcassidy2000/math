#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: EXACT characterization of the doublet joint curvature C(M)
as a SIGNED double-Dedekind sum — pins down the single object codex must bound (HYP-2796).

From the pointwise second difference of the cover indicator (derivation in HYP-2797):
  C(M) = p0(B+{M,M+1}) - p0(B+{M}) - p0(B+{M+1}) + p0(B)
       = meas{phi: B misses EXACTLY {j,j'},  {sector_M(phi),sector_{M+1}(phi)} = {j,j'}}   [+ part]
       - meas{phi: B misses EXACTLY {j},     sector_M(phi)=sector_{M+1}(phi)=j}             [- part]
The (+) part lives on the 2-MISS arcs of B (doublet usefully covers both missing sectors);
the (-) part lives on the 1-MISS arcs (doublet redundantly hits the same missing sector).
Both are double-sawtooth (Asano multiple Dedekind) sums in (M, M+1); C(M)->C_sat by
Dedekind/equidistribution at rate O(1/M).

This script VERIFIES the characterization EXACTLY (rationals): it computes C(M) directly
and via the (+)/(-) decomposition on the base miss-arcs, and confirms they match. It also
reports the 1-miss and 2-miss arc total measures (the supports), and the (+)/(-) split of C.
Base B = consec {0,...,k-3}.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast

ALL = frozenset(range(1, 7))


def lcm(a, b):
    return a // gcd(a, b) * b


def sectors_breakpoints(speeds):
    """Exact breakpoints of phi in [0,1) for the given integer speeds (sector = floor(7*frac(e*phi)))."""
    bps = {F(0), F(1)}
    for e in speeds:
        if e == 0:
            continue
        for j in range(7 * abs(e) + 1):
            bps.add(F(j, 7 * abs(e)))
    return sorted(b for b in bps if 0 <= b <= 1)


def sector_at(e, phi):
    # sector = floor(7*frac(e*phi)) in 0..6
    return int((e * phi) % 1 * 7) % 7 if e else 0


def miss_arcs(B, kmax=2):
    """Return dict: for r=1,2 the total measure where B misses EXACTLY r inner sectors,
    plus the per-arc (missing-set, lo, hi) list. Exact."""
    nz = [e for e in B if e]
    bps = sectors_breakpoints(nz)
    arcs = {1: [], 2: []}
    meas = {1: F(0), 2: F(0)}
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cov = set()
        for e in nz:
            s = sector_at(e, mid)
            if 1 <= s <= 6:
                cov.add(s)
        missing = ALL - cov
        r = len(missing)
        if r in (1, 2):
            arcs[r].append((frozenset(missing), lo, hi))
            meas[r] += hi - lo
    return arcs, meas


def curvature_via_arcs(B, M):
    """C(M) reconstructed from the (+)/(-) decomposition on base miss-arcs (exact)."""
    arcs, _ = miss_arcs(B)
    plus = F(0)   # 2-miss arcs, doublet covers both
    minus = F(0)  # 1-miss arcs, both doublet sectors = the single missing one
    # need to refine each base arc by the M, M+1 sector breakpoints
    for r, lst in arcs.items():
        for missing, lo, hi in lst:
            sub = sorted(set([lo, hi]) |
                         {b for b in sectors_breakpoints([M, M + 1]) if lo < b < hi})
            for a, c in zip(sub, sub[1:]):
                if c <= a:
                    continue
                mid = (a + c) / 2
                sM = sector_at(M, mid)
                sM1 = sector_at(M + 1, mid)
                if r == 2:
                    if {sM, sM1} == set(missing):
                        plus += c - a
                else:  # r == 1
                    j = next(iter(missing))
                    if sM == j and sM1 == j:
                        minus += c - a
    return plus - minus, plus, minus


def main():
    print("=" * 78)
    print("EXACT curvature C(M) = signed double-Dedekind sum  (claude-opus 2026-06-21)")
    print("C(M) = [2-miss arcs: doublet covers both] - [1-miss arcs: doublet redundant]")
    print("=" * 78)
    for k in (8, 9, 10):
        B = tuple(range(k - 2))
        arcs, meas = miss_arcs(B)
        p0B = p0_fast(B)
        print(f"\nk={k}  B=consec{B}")
        print(f"  base 1-miss arc measure = {float(meas[1]):.5f}  2-miss arc measure = {float(meas[2]):.5f}")
        ok = True
        for M in (15, 21, 40, 91):
            Cdirect = (p0_fast(tuple(sorted(B + (M, M + 1)))) - p0_fast(tuple(sorted(B + (M,))))
                       - p0_fast(tuple(sorted(B + (M + 1,)))) + p0B)
            Cviaarcs, plus, minus = curvature_via_arcs(B, M)
            match = (Cdirect == Cviaarcs)
            ok = ok and match
            print(f"   M={M:3d}: C_direct={float(Cdirect):+.6f}  C_via_arcs={float(Cviaarcs):+.6f} "
                  f"(+{float(plus):.5f} / -{float(minus):.5f})  match={match}")
        print(f"  => characterization {'VERIFIED exact' if ok else 'MISMATCH (check)'}")
    print("\n" + "=" * 78)
    print("If VERIFIED: codex's curvature-approach bound = bound a SIGNED double Dedekind sum")
    print("on the base 2-miss/1-miss arcs; Dedekind-Rademacher reciprocity gives the O(1/M) rate.")


if __name__ == "__main__":
    main()
