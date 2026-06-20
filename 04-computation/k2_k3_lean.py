#!/usr/bin/env python3
"""
Lean exact scan of the unproven layers:
 - k=2 (tower bit + 2 holes, 2 tails) FULL up to tail 40
 - k=3 (tower bit + 3 holes, 3 tails) with tails up to 28 (small window where measure is lowest)
The comb floor argument certifies larger tails. This targets the danger zone (small tails).
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

below = []; total = 0
for bit in [1,2,4,8]:
    pool = [d for d in range(1,14) if d != bit]
    # k=2
    for other in itertools.combinations(pool, 2):
        holes = set(other) | {bit}
        kept = [d for d in range(1,14) if d not in holes]
        nt = 12 - len(kept)
        for tails in itertools.combinations(range(14, 41), nt):
            C = tuple(sorted(kept + list(tails)))
            if len(C)!=12 or reduce(gcd,C)!=1: continue
            total+=1; L=lonely_measure(C)
            if L<THR2: below.append((L,bit,tuple(sorted(holes)),tails,'k2'))
    # k=3 small window
    for other in itertools.combinations(pool, 3):
        holes = set(other) | {bit}
        kept = [d for d in range(1,14) if d not in holes]
        nt = 12 - len(kept)
        for tails in itertools.combinations(range(14, 29), nt):
            C = tuple(sorted(kept + list(tails)))
            if len(C)!=12 or reduce(gcd,C)!=1: continue
            total+=1; L=lonely_measure(C)
            if L<THR2: below.append((L,bit,tuple(sorted(holes)),tails,'k3'))
    print(f"  bit {bit} done (total={total}, below={len(below)})", flush=True)

print(f"\nk=2(tails<=40) + k=3(tails<=28) exact: scanned={total}; below thr2={len(below)}")
for L,bit,holes,tails,lay in sorted(below)[:20]:
    print(f"  COUNTEREX[{lay}] meas={float(L):.6f}={L} bit={bit} holes={holes} tails={tails}")
if not below:
    print("  NONE below thr2 => k=2/k=3 danger zone clean, no tower-deletion counterexample.")
