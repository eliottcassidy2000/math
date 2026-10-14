#!/usr/bin/env python3
"""Independent brute-force cross-check of the D19 special curves (no elliptic-curve machinery):
 (a) midline x = 1/2: rational y with 4y^2 + 1 and 4(1 - y)^2 + 1 both rational squares;
 (b) side line y = 0: rational x with x^2 + 1 and (1 - x)^2 + 1 both rational squares;
searched over all reduced fractions with |numerator| <= H, denominator <= H.  The theorem (rank 0, torsion Z/4)
predicts no solutions at all.  Also: (c) the circumcircle and diagonals are excluded by Ptolemy / sqrt2 (proof only).
"""
import sys
from math import gcd, isqrt
H = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
def is_sq(n):
    return n >= 0 and isqrt(n) ** 2 == n
mid, side, checked = [], [], 0
for q in range(1, H + 1):
    q2 = q * q
    for p in range(-2 * H, 2 * H + 1):
        if gcd(p, q) != 1:
            continue
        checked += 1
        # midline: 4y^2+1 = (4p^2 + q^2)/q^2 ; 4(1-y)^2+1 = (4(q-p)^2 + q^2)/q^2
        if is_sq(4 * p * p + q2) and is_sq(4 * (q - p) ** 2 + q2):
            mid.append((p, q))
        # side: x^2+1 = (p^2+q^2)/q^2 ; (1-x)^2+1 = ((q-p)^2+q^2)/q^2
        if is_sq(p * p + q2) and is_sq((q - p) ** 2 + q2):
            side.append((p, q))
print(f"reduced fractions p/q with q <= {H}, |p| <= {2*H}: {checked}")
print(f"(a) midline solutions: {mid}")
print(f"(b) side-line solutions: {side}")
