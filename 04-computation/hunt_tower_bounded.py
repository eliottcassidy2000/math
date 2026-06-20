#!/usr/bin/env python3
"""
Tractable adversarial hunt: tower bit deleted, exhaustive over SMALL tails
(where danger arcs are widest and measure is smallest). Plus the floor argument
verification: show base measure M_B (k holes, before tails) times (6/7)^(#tails)
stays above thr2, so far tails cannot drag below threshold.
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

below = []
total = 0
TMAX = 40
for bit in [1,2,4,8]:
    cnt = 0
    for nh in range(1, 7):
        pool = [d for d in range(1,14) if d != bit]
        for other in itertools.combinations(pool, nh-1):
            holes = set(other) | {bit}
            kept = [d for d in range(1,14) if d not in holes]
            nt = 12 - len(kept)
            if nt < 0 or nt > 6: continue
            for tails in itertools.combinations(range(14, TMAX+1), nt):
                C = tuple(sorted(kept + list(tails)))
                if len(C) != 12: continue
                if reduce(gcd, C) != 1: continue
                total += 1; cnt += 1
                L = lonely_measure(C)
                if L < THR2:
                    below.append((L, bit, tuple(sorted(holes)), tails))
    print(f"  bit {bit}: scanned {cnt}, below={len([b for b in below if b[1]==bit])}", flush=True)

print(f"\nTOTAL scanned (tails<= {TMAX}): {total}; counterexamples below thr2: {len(below)}")
for L,bit,holes,tails in sorted(below)[:20]:
    print(f"  COUNTEREX meas={float(L):.6f}={L} bit={bit} holes={holes} tails={tails}")
if not below:
    print("  NONE below thr2.")

# ---- Floor argument: tails only DECREASE measure by factor >= 6/7 each (comb).
# Smallest base (bit deleted + max other holes, no tails) measure:
print("\n=== Floor: min base measure with bit deleted, varying #extra holes ===")
for bit in [1,2,4,8]:
    print(f"  delete {bit}:")
    for nh_extra in range(0, 5):
        best = None
        pool = [d for d in range(1,14) if d != bit]
        for other in itertools.combinations(pool, nh_extra):
            holes = set(other) | {bit}
            kept = tuple(d for d in range(1,14) if d not in holes)
            if len(kept) < 12 - 0:  # this is just the base set, may be <12; measure still meaningful
                pass
            M = lonely_measure(kept)
            if best is None or M < best[0]:
                best = (M, holes)
        # base set has 12 - nh_extra elements; adding nh_extra tails each cuts by >=6/7
        floor = best[0] * Fraction(6,7)**(nh_extra+0)  # tails count = nh_extra+1 (bit) ... approx
        print(f"    extra_holes={nh_extra}: min base M={float(best[0]):.5f} (holes {sorted(best[1])}), "
              f"(6/7)^? floor stays>thr2: {best[0]*Fraction(6,7)**(nh_extra+1) > THR2}")
