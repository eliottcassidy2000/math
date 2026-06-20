#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Adversarial counterexample hunt: AP-tail cores with LARGER tails (beyond cap20).
Look for any sub-THR2 core that is TOWER-BROKEN (missing one of {1,2,4,8})."""
import itertools
from fractions import Fraction
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])

def ap_tail_core(holes, tails):
    return tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))

TMAX = 80
found_sub = []
counterexamples = []
scanned = 0
# k holes -> k tails, tails in [14,TMAX]
for k in range(1,5):
    for holes in itertools.combinations(range(1,14), k):
        for tails in itertools.combinations(range(14,TMAX+1), k):
            C = ap_tail_core(set(holes), tails)
            if len(C)!=12: continue
            if not primitive(C): continue
            scanned += 1
            L = lonely_measure(C)
            if L < THR2:
                tw = has_tower(C)
                found_sub.append((L, holes, tails, tw))
                if not tw:
                    counterexamples.append((L,holes,tails,C))

found_sub.sort()
print(f"scanned {scanned} AP-tail cores (k holes/tails, tails<= {TMAX})")
print(f"sub-THR2 AP-tail cores: {len(found_sub)}")
for L,h,t,tw in found_sub:
    print(f"  meas={float(L):.6f}={L}  holes={h} tails={t}  tower={tw}")
print(f"\nTOWER-BROKEN sub-THR2 counterexamples: {len(counterexamples)}")
for L,h,t,C in counterexamples:
    print(f"  *** meas={L} holes={h} tails={t} C={C} ***")
if not counterexamples:
    print("  NONE. HYP-2661 holds over AP-tail family with tails up to", TMAX)
