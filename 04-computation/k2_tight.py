#!/usr/bin/env python3
"""
Tight k=2 and k=3 exact scan: tower bit + (1,2,3) extra holes, tails capped small.
Floor argument (proven via comb) certifies tails >= cap. We exact-scan tails < CAP.
CAP chosen generously (50) so the residue window covers all the comb-uncertified cases:
the single-tail cutoffs in k=1 were R<=147, and k=2/k=3 bases are LARGER measure
(fewer elements) so their cutoffs are smaller. Cap 50 comfortably covers.
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

CAP = 50
below = []; total = 0
for bit in [1,2,4,8]:
    pool = [d for d in range(1,14) if d != bit]
    for nextra in [1,2,3]:
        for other in itertools.combinations(pool, nextra):
            holes = set(other) | {bit}
            kept = [d for d in range(1,14) if d not in holes]
            nt = 12 - len(kept)
            if nt < 0: continue
            for tails in itertools.combinations(range(14, CAP+1), nt):
                C = tuple(sorted(kept + list(tails)))
                if len(C) != 12: continue
                if reduce(gcd, C) != 1: continue
                total += 1
                L = lonely_measure(C)
                if L < THR2:
                    below.append((L, bit, tuple(sorted(holes)), tails))
    print(f"  bit {bit} done (cumulative total={total}, below={len(below)})", flush=True)

print(f"\nk=1..3 exact, tails<={CAP}: scanned={total}; below thr2={len(below)}")
for L,bit,holes,tails in sorted(below)[:20]:
    print(f"  COUNTEREX meas={float(L):.6f}={L} bit={bit} holes={holes} tails={tails}")
if not below:
    print("  NONE below thr2  => no tower-deletion counterexample up to 3 extra holes, tails<=50.")
