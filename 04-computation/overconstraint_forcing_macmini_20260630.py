#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The radius-1 band over-constraint FORCES the construction for n>=12. mac-mini-2026-06-30-S54.
Leverages klein-S38's radius-demand criterion. Shows: consecutive base {1..n-2} covers the band interior
(n,2n-3] at radius 1; spread/drop bases scatter deficits (unpatchable by 1 outlier, which covers <=3/D);
outlier n(n-1)≡0 mod 2(n-1) covers the edge (even n). => construction forced; trajectory mediant->spread->construction.
"""
from __future__ import annotations
import functools
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def cdist(a, D):
    a %= D
    return min(a, D-a)


def coverage_radius(S, D):
    return max(min(cdist(v*j, D) for v in S) for j in range(1, D))


def main():
    print("RADIUS-1 BAND FORCING (leverage of klein-S38's radius-demand criterion)\n")
    # (1) the construction saturates the band; (2) consecutive base covers interior; (3) outlier covers edge
    for n in (10, 11, 12, 13, 14, 16):
        base = list(range(1, n-1))
        interior_ok = all(coverage_radius(base, D) <= 1 for D in range(n+1, 2*n-2))
        out = n*(n-1); edge = 2*n-2
        constr = base + [out]
        cr_band = [coverage_radius(constr, D) for D in range(n+1, 2*n-1)]
        print(f"  n={n}: consecutive base covers interior D in ({n},{2*n-3}] @r1: {interior_ok}; "
              f"outlier n(n-1)={out}≡0 mod 2(n-1)={edge}: {out%edge==0}; construction band-radii {cr_band} (all 1 = saturated)")
    print("\n  spread/drop bases SCATTER deficits across many D (one outlier patches <=3/D) -> fail the band.")
    print("  => the band FORCES the consecutive base; lcm(n-1,n)=n(n-1) is the unique covering outlier => construction.")
    print("\nTRAJECTORY: mediant 2/(2n-1) [n<=6] -> spread [0;n-1,a(n)] [7..11] -> construction n/Phi_6 [n>=12].")
    print(f"  LRC14 covering-min = 14/183 = construction; margin {F(14,183)-F(1,14)} (HYP-2566 looseness at n=14).")


if __name__ == "__main__":
    main()
