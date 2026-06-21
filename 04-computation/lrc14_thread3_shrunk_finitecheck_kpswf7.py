#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- ANGLE 3: is the n=14 finite check FEASIBLE after the repo's reductions?

Rosenfeld 2025 (arXiv:2604.23906) proves LRC(n) for n=8,9,10 by a FINITE check: the gap
function is determined by speeds bounded by ~n^{2n}, then a (huge but finite) enumeration. For
n=14 the raw bound n^{2n}=14^28 ~ 1.7e32 is hopeless.

The repo reductions that SHRINK the relevant space:
  (R1) WIDE handled analytically: span(E) > 14 => p0 <= Q(k-1) + tiny err < cap. So the finite
       check only needs BOUNDED-span configs (span <= 14, i.e. all speeds within a window of 14).
  (R2) FULL-RESIDUE stratum: only configs hitting all residues mod 7 matter (else measS7=0 trivially
       safe). Cuts the count.
  (R3) consec-max LAYERS 1-2 proven: the optimum has a specific shape (minimal leg profile +
       doubled identity residue). Only need to check configs "near" consec.
  (R4) commensurable far-pairs reduce to a FINITE resonant atlas (1-D torus curves).

GOAL: estimate the reduced enumeration count for n=14 and compare to feasible (~1e9-1e12).
The LRC reformulation: speeds are DIFFERENCES; the relevant object is the multiset of speeds
mod 14 within a bounded window. With WIDE removed, speeds live in a window [0, W]. The
quantity that matters for measS7 is the residue mod 7 AND the magnitude (for conductance).

This script COUNTS:
  - raw Rosenfeld bound 14^28
  - bounded-window configs: choose k speeds from {1..W} (W = window after WIDE reduction)
  - full-residue + dedup-by-dilation reductions
for k = 13 (the lonely runner with n=14 => 13 OTHER runners => speeds), W ranging.
"""
import sys
from math import comb, log10, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def main():
    print("=== ANGLE 3: feasibility of the n=14 finite check after repo reductions ===\n")
    n = 14; k = n - 1  # 13 nonzero speed-differences
    raw = 2 * n * log10(n)
    print(f"Rosenfeld raw bound n^(2n) = 14^28 ~ 1e{raw:.0f}  (hopeless)\n")

    print("AFTER WIDE reduction: only span(E) <= 14 matters. Speeds in a window of width W<=14.")
    print("Count = #{ size-k subsets of {1..W} } (the bounded-span stratum):")
    for W in (14, 21, 28):
        for kk in (12, 13):
            if W >= kk:
                c = comb(W, kk)
                print(f"  W={W}, k={kk}: C({W},{kk}) = {c:.3e}  (log10={log10(c):.1f})")
    print()

    print("AFTER full-residue + dilation-quotient + consec-locality:")
    print("  - dilation quotient: /6 (Z/7* orbits)")
    print("  - full-residue: only ~ (fraction hitting all 6 residues) survive")
    print("  - consec-locality (LAYERS 1-2): optimum within Hamming distance d of consec")
    # estimate: configs within distance d of consec (swap d speeds out of window)
    print("  configs within edit-distance d of consec_13 (swap d of the 13 speeds to window {1..14}):")
    for d in (1, 2, 3, 4):
        # choose d positions to change, replace with other window values
        cnt = comb(13, d) * comb(14 - 13 + d, d)  # rough: d new values from remaining window slots
        # better: number of ways to delete d of 13 consec speeds and add d new from {1..14}\consec
        avail = 14 - 13  # window {1..14} has 14 slots, consec uses 13
        cnt = comb(13, d) * comb(avail + 0, 0) if avail < d else comb(13, d) * comb(avail, d)
        # general: window must be enlarged; use W=20 window
        W = 20
        cnt = comb(13, d) * comb(W - 13, d)
        print(f"    d={d} (window W=20): ~{cnt:.2e} configs")
    print()
    print("VERDICT REASONING:")
    print("  Even the FULL bounded-span stratum C(20,13)~7.8e4 to C(28,13)~3.7e7 is TINY and feasible.")
    print("  The real cost is the PER-CONFIG measS7 (exact rational) + the WIDE-glue + commensurable")
    print("  far-pair atlas. The bounded-span finite check at n=14 is ENUMERABLE (<1e8 configs).")
    print("  BUT: it only closes the FINITE/bounded piece -- still needs the WIDE analytic bound +")
    print("  the consec-max PROOF (not just check) to glue. So this is an ENGINEERING-feasible")
    print("  verification of the bounded piece, NOT a standalone proof. Rating depends on whether")
    print("  the bounded enumeration + interval-arithmetic over the resonant atlas is rigorous.")

if __name__ == "__main__":
    main()
