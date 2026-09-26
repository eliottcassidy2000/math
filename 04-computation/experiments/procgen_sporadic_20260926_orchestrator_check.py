#!/usr/bin/env python3
"""Orchestrator audit of lane `sporadic`, written from the lane's report; the
lane's scripts were not read.
  1. The two long primitive cycles that resolve the Belaga-Mignotte off-by-one:
     d = 14303 with least element 101 (clock (2155, 1092)) and d = 17021 with
     least element 5 (clock (2140, 1088)), for T_d(y) = y/2, (3y+d)/2.
  2. The divisibility facts 2^27 - 3^17 = 355*14303 and
     2^65 - 3^41 = 17021*19*29*44835377399.
  3. Shift criterion: a mixed shape (p,a) of T_{q,d} is free (every word integral)
     iff (2^p - q^a) | d (exhaustive for small q, d, p).
  4. |2^p - q^a| = 1 with a >= 1 only for a = 1 (q = 2^p +- 1) or (q,p,a) = (3,3,2)
     (q <= 201 odd, p <= 300, a <= 60).
"""
import math
from itertools import combinations


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def cycle_of(y0, d, cap=10 ** 6):
    y, seen = y0, []
    for _ in range(cap):
        seen.append(y)
        y = (3 * y + d) // 2 if y % 2 else y // 2
        if y == y0:
            return seen
    return None


for d, least, p, a in ((14303, 101, 2155, 1092), (17021, 5, 2140, 1088)):
    cyc = cycle_of(least, d)
    assert cyc is not None, d
    odd = sum(1 for y in cyc if y % 2)
    g = 0
    for y in cyc:
        g = math.gcd(g, y)
    assert len(cyc) == p and odd == a and min(cyc) == least and math.gcd(g, d) == 1, (d, len(cyc), odd, min(cyc), g)
    print(f"   d={d}: cycle through {least}: period {len(cyc)}, odd steps {odd}, min {min(cyc)}, max {max(cyc)}, gcd(elements, d) = {math.gcd(g, d)}, a/p = {odd/len(cyc):.4f}")
check(True, "the two long primitive cycles exist with the stated least elements, clocks (2155,1092) and (2140,1088), and gcd 1 with d")

check(2 ** 27 - 3 ** 17 == 355 * 14303 and 2 ** 65 - 3 ** 41 == 17021 * 19 * 29 * 44835377399,
      "2^27 - 3^17 = 355*14303 and 2^65 - 3^41 = 17021*19*29*44835377399")


def c_word(bits, q):
    c = 0
    for j, b in enumerate(bits):
        if b:
            c = q * c + 2 ** j
    return c


for q in (3, 5, 7):
    for d in (1, 5, 7, 11, 13, 23, 139):
        if math.gcd(d, 2 * q) != 1:
            continue
        for p in range(2, 13):
            for a in range(1, p):
                D = 2 ** p - q ** a
                allint = all((d * c_word([1 if j in ones else 0 for j in range(p)], q)) % D == 0
                             for ones in combinations(range(p), a))
                assert allint == (d % D == 0), (q, d, p, a)
check(True, "shift criterion: mixed shape (p,a) of (qx+d)/2 is free iff (2^p - q^a) | d (q = 3,5,7; d <= 139; p <= 12)")

sols = []
for q in range(3, 202, 2):
    for p in range(1, 301):
        for a in range(1, 61):
            if abs(2 ** p - q ** a) == 1:
                sols.append((q, p, a))
assert all(a == 1 for q, p, a in sols if (q, p, a) != (3, 3, 2)) and (3, 3, 2) in sols
check(True, f"|2^p - q^a| = 1 (a >= 1) only for a = 1 or (q,p,a) = (3,3,2): {len(sols)} solutions, q <= 201")
