#!/usr/bin/env python3
"""
ADVERSARIAL counterexample hunt for THM-545 tower-deletion lemma.
Look for a PRIMITIVE AP-tail 12-core C with a tower bit (1,2,4, or 8) DELETED
and meas(G_C) < thr2 = 426/35035. If found, the lemma is REFUTED.

Higher layers k>=2 (the unproven part): multiple holes, multiple tails.
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, ap_tail_core, THR2

def core_from(holes, tails):
    C = tuple(sorted([d for d in range(1,14) if d not in holes] + list(tails)))
    return C

below_global = []
total = 0
# For each tower bit deleted, scan configs with that bit absent.
# Config: holes (incl. bit) of size nh in {1..13}, tails of size nt in [14, TMAX], len==12, primitive.
TMAX = 60
for bit in [1,2,4,8]:
    for nh in range(1, 6):  # 1..5 holes total (k = nh-1 extra beyond the bit... but bit IS a hole)
        # bit must be among the holes
        other_holes_pool = [d for d in range(1,14) if d != bit]
        for other in itertools.combinations(other_holes_pool, nh-1):
            holes = set(other) | {bit}
            kept = [d for d in range(1,14) if d not in holes]
            nt = 12 - len(kept)   # tails needed to reach size 12
            if nt < 0: continue
            for tails in itertools.combinations(range(14, TMAX+1), nt):
                C = tuple(sorted(kept + list(tails)))
                if len(C) != 12: continue
                if reduce(gcd, C) != 1: continue
                total += 1
                L = lonely_measure(C)
                if L < THR2:
                    below_global.append((L, bit, tuple(sorted(holes)), tails))
    print(f"  bit {bit}: scanned (cumulative total={total}), below-so-far={len(below_global)}", flush=True)

print(f"\nTOTAL primitive tower-deleted cores scanned (tails<= {TMAX}): {total}")
print(f"Cores with meas < thr2 (= COUNTEREXAMPLES): {len(below_global)}")
for L,bit,holes,tails in sorted(below_global)[:20]:
    print(f"  COUNTEREX meas={L}={float(L):.6f} bit={bit} holes={holes} tails={tails}")
if not below_global:
    print("  NONE. No tower-deleted core falls below thr2 in this exhaustive window.")
