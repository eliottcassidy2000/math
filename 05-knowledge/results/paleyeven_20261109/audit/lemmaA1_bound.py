#!/usr/bin/env python3
"""Audit P1: the Cauchy-Schwarz bound of Lemma A1, f(s) = (q - s)(s - sqrt(s))/2 >= (q - 1)/2 for 3 <= s <= q - 2,
checked exactly for every q = 3 mod 4 up to 20000 (using integer arithmetic: f(s) >= (q-1)/2 iff
(q - s)(s - sqrt s) >= q - 1 iff (q - s) s - (q - 1) >= (q - s) sqrt s, both sides compared after squaring when the
left side is nonnegative). Also reports the threshold of the s = 3 case and which q fail (q = 7 is the only one)."""
from math import isqrt, sqrt
def ok(q, s):
    L = (q - s) * s - (q - 1)
    if L < 0:
        return False
    return L * L >= (q - s) * (q - s) * s
fails = [(q, s) for q in range(7, 20001, 4) for s in range(3, q - 1) if not ok(q, s)]
print("pairs (q, s) with f(s) < (q-1)/2, q = 3 mod 4 in [7, 20000], 3 <= s <= q-2:", fails)
print("s = 3 threshold: q >= 1.401925/0.133975 =", round((3 * (3 - sqrt(3)) / 2 - 0.5) / ((3 - sqrt(3)) / 2 - 0.5), 4))
