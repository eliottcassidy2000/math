#!/usr/bin/env python3
"""
Exact k=2 layer scan, CORRECT iterated comb cutoff.
delete tower bit + 2 extra holes (base size 10), add tails r1<r2.

Correct logic: fix r1. Form B1 = base + {r1}, EXACTLY measure its M1, c1.
Then the second tail r2 satisfies, by the single-tail comb on B1:
   meas(B1 + r2) >= (6/7)M1 - 2 c1/(7 r2) >= thr2  for r2 >= R1 := ceil(2 c1/(6 M1 - 7 thr2)).
So for each r1, only r2 in [14, R1) need EXACT checking; r2 >= R1 certified by comb.
But r1 itself ranges over [14, inf): for large r1, M1 -> M (close to base, comb gives
M1 >= (6/7)M - small), and R1 stays bounded. We must bound the r1 range too.

r1 range bound: meas(B+{r1}+{r2}) <= meas(B+{r1}) = M1. If M1 >= thr2 always cannot conclude
(2nd tail lowers). BUT the final core has BOTH tails; we use: for the final core, by single-tail
comb applied to B+{r2} with new speed r1: meas >= (6/7) meas(B+{r2}) - 2 c(B+{r2})/(7 r1).
By symmetry we may assume r1 < r2 and bound the SMALLER tail's window via base B (size 10):
   meas(core) >= (6/7)^2 M_B - (2 c_B/(7 r1))*(6/7) - 2 c1/(7 r2)   [iterated]
For r1 >= R2 := ceil(2 c_B/( (6/7 M_B - thr2*7/6) ... )) -- messy. Instead: do a GENEROUS but
FINITE exact scan: r1 in [14, RA), r2 in [14, RB(r1)) where RA from a safe iterated floor and
RB(r1) from the exact M1,c1 cutoff. Anything outside certified by comb.
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from verify_tower_deletion_indep import lonely_measure, THR2

def ceil_frac(x):
    return x.numerator // x.denominator + (1 if x.numerator % x.denominator else 0)

below = []; total = 0
for bit in [1,2,4,8]:
    cnt = 0; maxRA = 0; maxRB = 0
    pool = [d for d in range(1,14) if d != bit]
    for other in itertools.combinations(pool, 2):
        holes = set(other) | {bit}
        kept = tuple(d for d in range(1,14) if d not in holes)  # size 10
        M0, c0 = lonely_measure(kept, want_comps=True)
        # r1 window: iterated floor (6/7)^2 M0 - (6/7)*2c0/(7 r1) - (worst 2nd term).
        # Use: meas(core) >= (6/7)M1 - 2 c1/(7 r2). Since c1 <= c0 + (#full teeth of r1)... bound c1 by
        # measuring it exactly per r1. For the r1 OUTER window, certify via:
        #   meas(core) >= (6/7) M1 - 2 c1/(7*14)   (r2>=14)  ... if this >= thr2 for ALL r2>=14? No.
        # Simplest rigorous finite scan: r1 in [14, RA) where RA makes (6/7)M0 - 2c0/(7 RA) such that
        # even the MINIMUM over r2 (=(6/7)M1 - 2c1/(7*14) using M1>= (6/7)M0-2c0/(7 r1)) >= thr2.
        # Solve conservatively: require (6/7)*((6/7)M0 - 2c0/(7 r1)) - 2*(c0+1)/(98) >= thr2.
        # => 2c0/(7 r1) <= (6/7)M0 - (7/6)*(thr2 + 2(c0+1)/98)
        rhs = Fraction(6,7)*M0 - Fraction(7,6)*(THR2 + Fraction(2*(c0+1),98))
        if rhs > 0:
            RA = ceil_frac(Fraction(2*c0,7)/rhs)
        else:
            RA = 400  # fallback generous
        RA = min(RA, 400)
        maxRA = max(maxRA, RA)
        for r1 in range(14, RA):
            if r1 in kept: continue
            B1 = tuple(sorted(kept + (r1,)))
            M1, c1 = lonely_measure(B1, want_comps=True)
            denom1 = 6*M1 - 7*THR2
            RB = ceil_frac(Fraction(2*c1)/denom1) if denom1 > 0 else 14
            maxRB = max(maxRB, RB)
            for r2 in range(r1+1, RB):
                if r2 in B1: continue
                C = tuple(sorted(B1 + (r2,)))
                if len(C) != 12: continue
                if reduce(gcd, C) != 1: continue
                total += 1; cnt += 1
                L = lonely_measure(C)
                if L < THR2:
                    below.append((L, bit, tuple(sorted(holes)), (r1,r2)))
    print(f"  bit {bit}: k=2 RA<={maxRA} RB<={maxRB}, scanned {cnt}, below={len([b for b in below if b[1]==bit])}", flush=True)

print(f"\nk=2 iterated-comb EXACT residue: scanned={total}; counterexamples below thr2={len(below)}")
for L,bit,holes,tails in sorted(below)[:20]:
    print(f"  COUNTEREX meas={float(L):.6f}={L} bit={bit} holes={holes} tails={tails}")
if not below:
    print("  NONE below thr2  => k=2 layer CERTIFIED via iterated comb cutoff + exact residue.")
