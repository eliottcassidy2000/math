#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) wide-cover as an APEX-PRIME PARTITION FUNCTION (kps-2026-06-20-S21).
==========================================================================
The sole LRC sector residual is HYP-2675: span(E)>14 => p0(E)=meas(S7(E)) <= cap_k.
For wide E the runners DECORRELATE (Weyl), and the cover factors into a DECORRELATED
part (runners as independent particles on Z/7 -- the cut-space inclusion-exclusion)
plus the decorrelation error (the cycle-space relation-lattice correction, HYP-2606).

The decorrelated cover is a PARTITION FUNCTION over the CLUSTER SHAPE of E.  We show
(HYP-2694) its sup is the SINGLE COHERENT BLOCK -- the apex-prime twin of THM-027/555's
"max c3 = regular score".  A cluster {M..M+m-1} at slow time x has residues
frac((M+d)x) = phi + frac(d x), phi=frac(Mx) the anchor (decorrelated for separated
scales).  The decorrelated cover of a partition into blocks B_1..B_r is

    p0_decorr = mean_x [ mean_{phi_1..phi_r indep} 1( union of {phi_i + frac(d x)} hits all 6 inner) ].

Single block: the phi-integral is exact (piecewise-constant in phi, breakpoints at
(s/7 - frac(d x)) mod 1).  Multi-block: phi-grid.  EXACT rationals for the single block.
"""
import sys, itertools
from fractions import Fraction as F
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))


def single_block_decorr(m, Nx=1260):
    """Exact-in-phi decorrelated cover of one coherent block of m sweeping points."""
    tot = F(0)
    for ix in range(Nx):
        x = F(2 * ix + 1, 2 * Nx)
        r = [(j * x) % 1 for j in range(m)]
        bps = sorted({(F(s, 7) - rj) % 1 for rj in r for s in range(7)})
        bps.append(bps[0] + 1)
        good = F(0)
        for a, b in zip(bps, bps[1:]):
            mid = (a + b) / 2
            hit = {int(((mid + rj) % 1) * 7) for rj in r}
            if len(hit & INNER) == 6:
                good += b - a
        tot += good
    return tot / Nx


def multiblock_decorr(blocks, Nx=240, Nphi=28):
    """Decorrelated cover of a partition into blocks (phi grid)."""
    tot = 0.0
    grid = [(p + 0.5) / Nphi for p in range(Nphi)]
    for ix in range(Nx):
        x = (ix + 0.5) / Nx
        bres = [[(d * x) % 1 for d in B] for B in blocks]
        cnt = 0
        N = Nphi ** len(blocks)
        for phis in itertools.product(grid, repeat=len(blocks)):
            hit = set()
            for i, B in enumerate(bres):
                ph = phis[i]
                for rj in B:
                    hit.add(int(((ph + rj) % 1) * 7))
            if len(hit & INNER) == 6:
                cnt += 1
        tot += cnt / N
    return tot / Nx


def main():
    print(__doc__)
    print("=== SINGLE-BLOCK decorrelated cover (the wide extremal shape) vs cap_k (exact) ===")
    for k in range(8, 13):
        m = k - 1
        v = single_block_decorr(m, 1260)
        cap = CAPS[k]
        print(f"  k={k:2d} (m={m:2d}): p0_decorr={float(v):.4f}  cap={float(cap):.4f}  "
              f"margin={float(cap - v):.4f}  {'OK' if v < cap else 'EXCEEDS'}")
    print("\n=== Splitting strictly LOWERS coverage (single block is the sup) ===")
    splits = {9: [[(0, 1, 2, 3), (0, 1, 2, 3)], [(0, 1, 2, 3, 4, 5), (0, 1)],
                  [(0, 1, 2), (0, 1, 2), (0, 1)]],
              10: [[(0, 1, 2, 3, 4), (0, 1, 2, 3)], [(0, 1, 2), (0, 1, 2), (0, 1, 2)]]}
    for k in (9, 10):
        sb = multiblock_decorr([tuple(range(k - 1))], 200, 24)
        print(f"  k={k}: single-block ~ {sb:.4f}")
        for P in splits[k]:
            v = multiblock_decorr(P, 160, 20)
            print(f"      split {[len(b) for b in P]}: {v:.4f}  {'< single' if v < sb else '>= SINGLE'}")


if __name__ == "__main__":
    main()
