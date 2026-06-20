#!/usr/bin/env python3
"""Quick k=2 exact: tower bit + 2 holes + 2 tails, tails in [14,30]. Danger zone only."""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

below=[]; total=0
for bit in [1,2,4,8]:
    pool=[d for d in range(1,14) if d!=bit]
    for other in itertools.combinations(pool,2):
        holes=set(other)|{bit}
        kept=[d for d in range(1,14) if d not in holes]
        nt=12-len(kept)
        for tails in itertools.combinations(range(14,31),nt):
            C=tuple(sorted(kept+list(tails)))
            if len(C)!=12 or reduce(gcd,C)!=1: continue
            total+=1; L=lonely_measure(C)
            if L<THR2: below.append((L,bit,tuple(sorted(holes)),tails))
print(f"k=2 exact, tails<=30: scanned={total}; below thr2={len(below)}")
for L,bit,holes,tails in sorted(below)[:20]:
    print(f"  COUNTEREX meas={float(L):.6f}={L} bit={bit} holes={holes} tails={tails}")
if not below: print("  NONE below thr2 => k=2 danger zone clean.")
