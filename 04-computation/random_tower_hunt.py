#!/usr/bin/env python3
"""Fast randomized hunt: random primitive 12-cores MISSING a tower bit; check meas<thr2."""
from fractions import Fraction
import random
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

rng = random.Random(2026)
below = []; n = 0
TRIALS = 8000
for _ in range(TRIALS):
    bit = rng.choice([1,2,4,8])
    # build a 12-core: keep some of {1..13}\{bit}, fill rest with tails in [14, T]
    smalls = [d for d in range(1,14) if d != bit]
    nkeep = rng.randint(6, 12)
    nkeep = min(nkeep, 12)
    keep = sorted(rng.sample(smalls, min(nkeep, len(smalls))))
    ntail = 12 - len(keep)
    if ntail < 0:
        keep = keep[:12]; ntail = 0
    T = rng.choice([30,45,60,80])
    tails = sorted(rng.sample(range(14, T+1), ntail)) if ntail>0 else []
    C = tuple(sorted(keep + tails))
    if len(C) != 12: continue
    if bit in C: continue
    if reduce(gcd, C) != 1: continue
    n += 1
    L = lonely_measure(C)
    if L < THR2:
        below.append((L, bit, C))

print(f"random tower-missing 12-cores checked: {n}; below thr2: {len(below)}")
for L,bit,C in sorted(below)[:10]:
    print(f"  COUNTEREX meas={float(L):.6f}={L} missing_bit={bit} C={C}")
if not below:
    print("  NONE below thr2 (consistent with lemma).")
