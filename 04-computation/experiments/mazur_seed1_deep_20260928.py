#!/usr/bin/env python3
"""The seed-1 harmonic mass beyond level 18 by the depth decomposition (Theorem A applied twice):

    H_(m+d)(1) = sum_(y in T_d(1)) 3^d 2^(-A_d(y)) H_m(y),   H_m(y) = 3^m mu_m(y mod 3^m),

with mu_m for m <= 18 read off the level-18 law (consistency: reduction mod 3^m gives mu_m) and the depth-d layer
T_d(1) of the tree of 1 enumerated by a weight-pruned walk (children (2^a y - 1)/3 of the right parity; a child is
expanded iff its weight 3^(d') 2^(-A) is at least THETA).  Validation: for m + d = 18 the identity must reproduce the
exact H_18(1) = 0.41046 up to the pruned weight; the deficit measures the pruning loss at each depth d.

Reports H_n(1) for n = 19 .. 18 + DMAX with the pruned-weight sums, n^(1/6) H_n, the ratios, and the class split
(mod 3) of the truncated layers T_d(1) (the recursion H_(d+1) = H^(1) + 2 H^(2)).
Session: opus, collatz-poset-dag-20260927 (S19), 2026-09-28.
Run: python 04-computation/experiments/mazur_seed1_deep_20260928.py   (about 20 GB, a few minutes)
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from mazur_harmonic_mass_deep_20260928 import level_up  # noqa: E402

TOP = 18
DMAX = 20
THETA = 1e-8

if __name__ == "__main__":
    t0 = time.time()
    mu = np.array([0.0, 1 / 3, 2 / 3])
    for n in range(2, TOP + 1):
        mu = level_up(mu, n)
    print(f"level-{TOP} law computed [{time.time() - t0:.0f}s]; H_{TOP}(1) = {3 ** TOP * mu[1]:.6f}")
    # reductions: mu_m(z) = sum over lifts; stored as arrays H_m(z) = 3^m mu_m(z) for m = 0..TOP
    Hm = {}
    red = mu.copy()
    for m in range(TOP, -1, -1):
        Hm[m] = (3 ** m) * red
        if m > 0:
            red = red.reshape(3, 3 ** (m - 1)).sum(axis=0)
    del red
    assert abs(Hm[0][0] - 1) < 1e-9
    # pruned walk over the tree of 1: layers[d] = list of (y, weight); pruned[d] = weight lost at depth d (children not expanded)
    layers = {0: [(1, 1.0)]}
    pruned = {d: 0.0 for d in range(DMAX + 1)}
    for d in range(DMAX):
        nxt = []
        lost = 0.0
        for y, w in layers[d]:
            r = y % 3
            if r == 0:
                continue
            eps = 1 if r == 2 else 0  # a odd iff y = 2 mod 3
            a = 2 - eps
            while True:
                cw = w * 3 * 2.0 ** (-a)
                if cw < THETA:
                    # geometric tail of the pruned children: cw (1 + 1/4 + ...) = cw 4/3
                    lost += cw * 4 / 3
                    break
                z = ((1 << a) * y - 1) // 3
                nxt.append((z, cw))
                a += 2
        layers[d + 1] = nxt
        pruned[d + 1] = pruned[d] + lost  # cumulative weight removed from the walk by depth d + 1
        print(f"   depth {d + 1}: {len(nxt)} nodes kept, cumulative pruned weight {pruned[d + 1]:.3e} [{time.time() - t0:.0f}s]")
    # validation at n = TOP: H_TOP(1) = sum_(y in T_d) w H_(TOP-d)(y) for every d
    print(f"== validation: reproduce H_{TOP}(1) = {Hm[TOP][1]:.6f} from the depth-d layers ==")
    for d in (1, 2, 5, 10, 15, DMAX - 2, DMAX):
        if d > TOP:
            continue
        m = TOP - d
        mod = 3 ** m
        val = sum(w * Hm[m][y % mod] for y, w in layers[d])
        print(f"   d={d:2d}, m={m:2d}: {val:.6f}  (deficit {Hm[TOP][1] - val:.2e}; pruned weight {pruned[d]:.2e})")
    # the extension: H_n(1) for n = TOP + 1 .. TOP + DMAX, from (d, m = TOP) and cross-checked with (d + 1, m = TOP - 1)
    print("== H_n(1) beyond the computed levels ==")
    modT = 3 ** TOP
    modT1 = 3 ** (TOP - 1)
    prev = Hm[TOP][1]
    for d in range(1, DMAX + 1):
        n = TOP + d
        val = sum(w * Hm[TOP][y % modT] for y, w in layers[d])
        alt = sum(w * Hm[TOP - 1][y % modT1] for y, w in layers[d + 1]) if d + 1 <= DMAX else float("nan")
        print(f"   n={n}: H_n(1) = {val:.5f} (via d={d}; via d={d + 1}: {alt:.5f}); n^(1/6) H_n = {n ** (1 / 6) * val:.4f}; ratio {val / prev:.4f}; pruned weight <= {pruned[d]:.1e}")
        prev = val
    # class split of the truncated layers (recursion check)
    print("== class split (mod 3) of the truncated layer weights of the tree of 1 ==")
    for d in range(1, DMAX + 1):
        tot = [0.0, 0.0, 0.0]
        for y, w in layers[d]:
            tot[y % 3] += w
        H = sum(tot)
        nxt = tot[1] + 2 * tot[2]
        print(f"   d={d:2d}: H_d = {H:.5f} (exact {Hm[d][1 % 3 ** d] if d <= TOP else float('nan'):.5f}); shares (0,1,2) = ({tot[0] / H:.3f}, {tot[1] / H:.3f}, {tot[2] / H:.3f}); H^(1) + 2H^(2) = {nxt:.5f} = next layer's mass")
    print("DONE")
