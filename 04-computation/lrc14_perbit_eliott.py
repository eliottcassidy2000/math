#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from fractions import Fraction
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])

drop6 = set(ap_tail_core({6},[]))  # (1,2,3,4,5,7,8,9,10,11,12,13)
print("drop6 =", tuple(sorted(drop6)), "has tower:", has_tower(drop6))
print("THR2 =", THR2, float(THR2))
print()
# Break exactly one tower bit: remove t from {1,2,4,8}, add a replacement to keep |C|=12 primitive.
# The claim's "breaking only X" cores. Let's just delete each tower bit and add 14 (smallest free).
for t in [1,2,4,8]:
    # replacement must keep set size 12 and primitive; pick smallest not in set
    repl = None
    for cand in range(14, 60):
        if cand not in drop6 and cand != t:
            repl = cand; break
    C = tuple(sorted((drop6 - {t}) | {repl}))
    L = lonely_measure(C)
    print(f"break {t} (->{repl}): C={C}")
    print(f"   meas={L}={float(L):.6f}  tower? {has_tower(C)}  >=THR2? {L>=THR2}")
