#!/usr/bin/env python3
"""Verify base case k=0 measures and the floor-argument numbers from THM-545."""
from fractions import Fraction
import itertools
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

print("=== Base case k=0: single-hole bases {1..13}\{2^a}, no tail ===")
expect = {0:Fraction(1,14),1:Fraction(11,364),2:Fraction(97,4004),3:Fraction(950,21021)}
for a in [0,1,2,3]:
    bit=2**a
    B=tuple(d for d in range(1,14) if d!=bit)
    M=lonely_measure(B)
    print(f"  a={a} del {bit}: meas={M} (claim {expect[a]}) match={M==expect[a]}  >=thr2: {M>=THR2}")

print("\n=== Floor numbers for a=3 (delete 8): (6/7)^k * min k-hole base M ===")
# claim: 0.0452, 0.0457, 0.0569, 0.0678, 0.0775 for k=0..4
bit=8
for k in range(0,5):
    pool=[d for d in range(1,14) if d!=bit]
    best=None
    for other in itertools.combinations(pool,k):
        holes=set(other)|{bit}
        kept=tuple(d for d in range(1,14) if d not in holes)
        M=lonely_measure(kept)
        if best is None or M<best: best=M
    print(f"  k={k}: min base M={float(best):.4f}  >>thr2={float(THR2):.4f}")
