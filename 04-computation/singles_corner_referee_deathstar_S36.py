#!/usr/bin/env python3
"""death-star-2026-07-17-S36 (HYP-7173): referee for LRCDeviationSingles.lean and the
QuadDenseCoreClosed corner closure.

SINGLES: N_i = q-1-bandSize under gcd(v_i,q)=1; D_i = -13/7 EXACT at 14|q;
-13/7 <= D_i <= 0 for all q >= 14 (gcd case).
CORNER: at j=10/11 the eps=0 top-block fee closes exactly when
  sum_u [2d/7 + 3/(7u)] < 2d, d = (13-j)/(14(j+1) w[j-1])."""
import random
from fractions import Fraction as Fr
from math import gcd

def band_size(q):
    return sum(1 for r in range(q) if q <= 14 * r <= 13 * q)

def N_single(vi, q):
    return sum(1 for p in range(1, q)
               if not (q <= 14 * ((vi * p) % q) <= 13 * q))

def referee_singles(trials=3000, seed=36):
    rnd = random.Random(seed)
    ok_exact = ok_const = ok_bounds = True
    for _ in range(trials):
        q = rnd.randint(14, 700)
        vi = rnd.randint(1, 10**6)
        if gcd(vi, q) != 1:
            continue
        Ni = N_single(vi, q)
        if Ni != (q - 1) - band_size(q):
            ok_exact = False
        D = Fr(Ni) - Fr(q - 1, 7)
        if not (Fr(-13, 7) <= D <= 0):
            ok_bounds = False
        if q % 14 == 0 and D != Fr(-13, 7):
            ok_const = False
    print(f"singles referee: exact-count {'PASS' if ok_exact else 'FAIL'}; "
          f"bounds [-13/7, 0] {'PASS' if ok_bounds else 'FAIL'}; "
          f"14|q constant {'PASS' if ok_const else 'FAIL'}")

def referee_corner(trials=40000, seed=136):
    rnd = random.Random(seed)
    n_close10 = n_core10 = n_close11 = n_core11 = 0
    for _ in range(trials):
        j = rnd.choice([10, 11])
        w = sorted(rnd.sample(range(1, 4000), 13))
        # force a dense pair at j and ladder/fee shape below is irrelevant for the FEE test
        d = Fr(13 - j, 14 * (j + 1) * w[j - 1])
        us = [w[10], w[11], w[12]] if j == 10 else [w[11], w[12]]
        fee = sum(2 * d / 7 + Fr(3, 7 * u) for u in us)
        if fee < 2 * d:
            if j == 10:
                n_close10 += 1
            else:
                n_close11 += 1
        else:
            if j == 10:
                n_core10 += 1
            else:
                n_core11 += 1
    print(f"corner fee referee: j=10 close {n_close10} / core {n_core10}; "
          f"j=11 close {n_close11} / core {n_core11}")
    print("  (the fee needs the top cluster to dwarf the cited base by the fat-mass")
    print("   margin; closures are the structurally-heavy top families, as designed)")

if __name__ == "__main__":
    referee_singles()
    referee_corner()
