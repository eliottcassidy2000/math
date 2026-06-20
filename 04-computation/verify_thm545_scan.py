#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXHAUSTIVE adversarial scan for THM-545.
(1) Reproduce: exactly TWO cores strictly < thr2 among AP-tail 12-cores (|holes|<=3, tails<=30).
(2) Hunt a counterexample: any sub-thr2 core that BREAKS the dyadic tower {1,2,4,8}
    (i.e. carry[1] != 15). The Level-1 obstruction claims NONE exist.
"""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd
from verify_thm545_indep import lonely_measure, carry_profile, ap_tail_core, THR2

def tower_intact(C):
    """carry[1]==15 iff {1,2,4,8} all in C."""
    return carry_profile(C).get(1, 0) == 15

def scan(max_holes=3, tail_max=30):
    base = {1,2,3,4,5,6,7,8,9,10,11,12,13}
    below = []          # all cores strictly < thr2
    below_broken = []   # sub-thr2 cores with tower BROKEN (would refute Level-1)
    seen = set()
    count = 0
    tail_pool = list(range(14, tail_max + 1))
    for nh in range(0, max_holes + 1):
        for holes in itertools.combinations(range(1, 14), nh):
            kept = [x for x in range(1, 14) if x not in holes]
            n_tail = 12 - len(kept)
            if n_tail < 0:
                continue
            for tails in itertools.combinations(tail_pool, n_tail):
                C = tuple(sorted(kept + list(tails)))
                if len(C) != 12:
                    continue
                if reduce(gcd, C) != 1:
                    continue
                if C in seen:
                    continue
                seen.add(C)
                count += 1
                L = lonely_measure(C)
                if L < THR2:
                    below.append((L, C))
                    if not tower_intact(C):
                        below_broken.append((L, C))
    return count, below, below_broken

def main():
    print("=== EXHAUSTIVE SCAN (|holes|<=3, tails 14..30) ===")
    count, below, below_broken = scan(3, 30)
    below.sort()
    print(f"  cores scanned (distinct, gcd=1, len=12): {count}")
    print(f"  cores strictly < thr2 (426/35035): {len(below)}")
    for L, C in below:
        print(f"    meas={L}={float(L):.8f}  C={C}  carry={carry_profile(C)}  tower_intact={tower_intact(C)}")
    print(f"\n  *** COUNTEREXAMPLES (sub-thr2 with BROKEN tower): {len(below_broken)} ***")
    for L, C in below_broken:
        print(f"    meas={L}  C={C}  carry={carry_profile(C)}")
    if not below_broken:
        print("    NONE. Level-1 obstruction holds on this domain: meas<thr2 => tower intact.")

    # Wider tail range as extra adversarial pressure
    print("\n=== WIDER SCAN (|holes|<=3, tails 14..40) ===")
    count2, below2, below_broken2 = scan(3, 40)
    below2.sort()
    print(f"  cores scanned: {count2};  strictly < thr2: {len(below2)}")
    for L, C in below2:
        print(f"    meas={float(L):.8f}  C={C}  tower_intact={tower_intact(C)}")
    print(f"  COUNTEREXAMPLES (broken tower, sub-thr2): {len(below_broken2)}")
    for L, C in below_broken2:
        print(f"    meas={L} C={C} carry={carry_profile(C)}")

if __name__ == "__main__":
    main()
