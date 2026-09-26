#!/usr/bin/env python3
"""Orchestrator audit of lane `rank` (forced charges on backward trees), written
from the note's statements; the lane's scripts were not read.

  1. The backward tree of -1 in Z_(2) has 2^(j-1) points at depth j (j = 1..12).
  2. Shadow violation: for the bank {-1: 0.9 kappa} (a = 1, h = 0), a good integer
     y at exact 2-adic depth D from -1 has R(Ty) - R(y) -> kappa - 0.9 kappa.
  3. Entry (seam) violation: for the bank {-1: 1.2 kappa} alone, the good integer
     at depth D from the uncharged preimage z = -4 (T^2 z = -1) has
     R(T^2 y) - R(y) ~ 1.2 kappa (D - 2), linear in D (Theorem A(ii)).
  4. Theorem A(i) on the cycle -5 -> -7 -> -10 (chi = log(9/8)/3): the bank
     {cycle points: 0.95 a chi} fails inside the shadow, while 1.05 a chi passes
     the one-period test (Proposition S direction).
"""
import math
from fractions import Fraction

KAPPA = math.log(1.5)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def T(x):
    """Collatz shortcut on Z_(2) (Fractions with odd denominator)."""
    x = Fraction(x)
    return (3 * x + 1) / 2 if x.numerator % 2 else x / 2


def preimages(y):
    y = Fraction(y)
    return [2 * y, (2 * y - 1) / 3]


# 1. backward tree of -1
level = {Fraction(-1)}
seen = {Fraction(-1)}
for j in range(1, 13):
    nxt = set()
    for y in level:
        for z in preimages(y):
            assert T(z) == y
            if z not in seen:
                nxt.add(z)
    level = nxt
    seen |= nxt
    assert len(level) == 2 ** (j - 1), (j, len(level))
check(True, "backward tree of -1 has 2^(j-1) new points at depth j for j = 1..12")


def v2_frac(x):
    x = Fraction(x)
    if x == 0:
        return 10 ** 9
    n, d = x.numerator, x.denominator
    v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def good_integer(z, D):
    """integer y in [2^(D+1), 2^(D+2)) with v2(y - z) = D exactly (z in Z_(2))."""
    z = Fraction(z)
    M = 2 ** (D + 1)
    r = (z.numerator * pow(z.denominator, -1, M)) % M      # z mod 2^(D+1)
    y = r ^ (1 << D)                                        # differ exactly at bit D
    y = y + M * ((2 ** (D + 1) - y + M - 1) // M)           # lift into [2^(D+1), 2^(D+2))
    assert v2_frac(Fraction(y) - z) == D
    return y


def Rank(n, bank, a=1.0):
    return a * math.log(n) + sum(c * v2_frac(Fraction(n) - b) for b, c in bank)


# 2. shadow violation at -1
bank = [(Fraction(-1), 0.9 * KAPPA)]
for D in (30, 60, 120):
    y = good_integer(-1, D)
    inc = Rank(T(y), bank) - Rank(y, bank)
    assert abs(inc - 0.1 * KAPPA) < 1e-6, (D, inc)
check(True, "bank {-1: 0.9 kappa}: every shadow step increases R by 0.1 kappa = %.5f (D = 30, 60, 120)" % (0.1 * KAPPA))

# 3. entry violation from the uncharged preimage -4 of -2 of -1
bank = [(Fraction(-1), 1.2 * KAPPA)]
incs = []
for D in (40, 80, 160):
    y = good_integer(-4, D)
    y2 = T(T(y))
    inc = Rank(y2, bank) - Rank(y, bank)
    incs.append(inc)
    assert abs(inc - (1.2 * KAPPA * (D - 2) + math.log(y2 / y))) < 1e-6
assert incs[1] - incs[0] > 1.2 * KAPPA * 39 and incs[2] - incs[1] > 1.2 * KAPPA * 79
check(True, "bank {-1: 1.2 kappa}: entering the shadow from the uncharged preimage -4 raises R by 1.2 kappa (D - 2) + O(1), linear in D")

# 4. forced charge a*chi on the cycle -5 -> -7 -> -10
cyc = [Fraction(-5), Fraction(-7), Fraction(-10)]
chi = math.log(9 / 8) / 3
for factor, should_fail in ((0.95, True), (1.05, False)):
    bank = [(x, factor * chi) for x in cyc]
    D = 200
    y = good_integer(-5, D)
    y3 = T(T(T(y)))
    inc = Rank(y3, bank) - Rank(y, bank)
    # one period: log(T^3 y / y) -> log(9/8); counters drop by 3 at charge factor*chi
    pred = math.log(9 / 8) - 3 * factor * chi
    assert abs(inc - pred) < 1e-6, (factor, inc, pred)
    assert (inc > 0) == should_fail
check(True, "cycle -5,-7,-10: charge 0.95 chi fails over a period (increment > 0), 1.05 chi passes (increment < 0)")
