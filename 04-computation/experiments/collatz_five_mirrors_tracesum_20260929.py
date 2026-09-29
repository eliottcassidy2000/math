#!/usr/bin/env python3
"""The trace of a rational cycle is reversal-invariant (S20 addendum, Proposition 6).

For a word w = (a_1..a_d) let rot_j(w) be its cyclic rotations and C_w(u, v) = sum_i u^(d-i) v^(a_1+...+a_(i-1)).
The cycle-sum polynomial S_w(u, v) = sum_j C_(rot_j w)(u, v) equals sum_(i=1)^d u^(d-i) B_(i-1)(w; v), where
B_m(w; v) = sum_j v^(sum of the m cyclically consecutive entries starting at j) depends only on the multiset of
cyclic block sums, which reversal preserves.  Hence S_(rev w) = S_w, and the elements of the cycle of w for the map
x -> (3x + k)/2^v, y_j = k C_(rot_j w)(3,2)/(2^A - 3^d), have a reversal-invariant sum.
Checks: the identity on random words and rational (u, v); the verified reversal pairs of 3x+13 and 3x+37; the
seven-cycle of -17 and its rational partner.
Run: python 04-computation/experiments/collatz_five_mirrors_tracesum_20260929.py
"""
from __future__ import annotations

import random
from fractions import Fraction


def carry(word, u, v):
    d = len(word)
    s, pref = 0, 0
    for i in range(1, d + 1):
        s += u ** (d - i) * v ** pref
        pref += word[i - 1]
    return s


def cycle_sum_poly(word, u, v):
    d = len(word)
    return sum(carry(word[j:] + word[:j], u, v) for j in range(d))


def cycle_points(word, k):
    d, A = len(word), sum(word)
    den = 2 ** A - 3 ** d
    return [Fraction(k * carry(word[j:] + word[:j], 3, 2), den) for j in range(len(word))]


if __name__ == "__main__":
    rng = random.Random(7)
    ok = True
    for _ in range(300):
        d = rng.randint(1, 9)
        w = [rng.randint(1, 6) for _ in range(d)]
        u = Fraction(rng.randint(1, 9), rng.randint(1, 9))
        v = Fraction(rng.randint(1, 9), rng.randint(1, 9))
        ok &= cycle_sum_poly(w, u, v) == cycle_sum_poly(w[::-1], u, v)
        # the sum of squares is not invariant in general
    print(f"   S_w(u,v) = S_rev(w)(u,v) on 300 random words and rational (u,v): {ok}")
    sq_inv = all(sum(c * c for c in [carry(w[j:] + w[:j], 3, 2) for j in range(len(w))]) == sum(c * c for c in [carry(w[::-1][j:] + w[::-1][:j], 3, 2) for j in range(len(w))]) for w in ([1, 2, 3], [1, 1, 2, 3], [2, 1, 1, 3, 1]))
    print(f"   sum of squares of the carries invariant on (1,2,3), (1,1,2,3), (2,1,1,3,1)? {sq_inv}")
    for k, w, name in ((13, [1, 1, 1, 2, 3], "3x+13, cycle of 227"), (37, [1, 2, 3], "3x+37, cycle of 23"), (1, [1, 1, 1, 2, 1, 1, 4], "3x+1, seven-cycle of -17")):
        pts = cycle_points(w, k)
        ptsr = cycle_points(w[::-1], k)
        print(f"   {name}: word {w}: points {[int(p) if p.denominator == 1 else str(p) for p in pts]} sum {sum(pts)}; reversed word {w[::-1]}: points {[int(p) if p.denominator == 1 else str(p) for p in ptsr]} sum {sum(ptsr)}")
    print("DONE")
