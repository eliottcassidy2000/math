#!/usr/bin/env python3
"""THM-4591 statement 3: E(L,k) = Lyn(L,k)/|2^L - 3^k| ranking over 12 <= L <= 1500 (ALL k), and the
rational-cycle denominators realised at the Ellison shapes (own enumeration)."""
from math import lgamma, log, gcd
from fractions import Fraction
from collections import Counter
from b3_lyndon import lyn_count, lyndon_words, cycle_value

# --- ranking (log-scale screen over all k, exact for the leaders)
cand = []
for L in range(12, 1501):
    lc = [lgamma(L + 1) - lgamma(k + 1) - lgamma(L - k + 1) for k in range(L + 1)]
    for k in range(1, L + 1):
        D = abs(2 ** L - 3 ** k)
        le = lc[k] - log(L) - log(D)
        cand.append((le, L, k))
cand.sort(reverse=True)
top = []
for le, L, k in cand[:40]:
    E = Fraction(lyn_count(L, k), abs(2 ** L - 3 ** k))
    top.append((float(E), L, k))
top.sort(reverse=True)
print("top shapes by exact E (L>=12, all k):")
for E, L, k in top[:10]:
    print(f"  ({L},{k})  E={E:.5f}  Lyn={lyn_count(L,k)}  gap 2^L-3^k={2**L-3**k}")
print("screen: largest log-E outside the 40 exact candidates:", cand[40][1:], "approx E", 2.718281828 ** cand[40][0])
print("sum of top five E:", sum(t[0] for t in top[:5]))
# Ellison exceptions check: |2^x - 3^y| < 2^x e^{-x/10}, x >= 12 (y nearest)
import math
exc = []
for x in range(12, 2000):
    for y in (int(x / math.log2(3)), int(x / math.log2(3)) + 1):
        d = abs(2 ** x - 3 ** y)
        if d and math.log(d) < x * math.log(2) - x / 10:
            exc.append((x, y))
print("Ellison inequality fails (x>=12, x<2000) at:", exc)

import sys
if len(sys.argv) < 2: sys.exit(0)
# --- rational cycles at the shapes
for (L, k) in [(11, 7), (12, 7), (13, 8), (14, 9), (16, 10), (19, 12), (27, 17)]:
    W = lyndon_words(L, k)
    assert len(W) == lyn_count(L, k)
    D = 2 ** L - 3 ** k
    cnt = Counter()
    ex = {}
    for w in W:
        d, DD = cycle_value(w)
        g = gcd(d, abs(DD))
        red = abs(DD) // g
        cnt[red] += 1
        ex.setdefault(red, []).append(Fraction(d, DD))
    small = {r: c for r, c in sorted(cnt.items()) if r < abs(D)}
    print(f"shape ({L},{k}) gap {D}: Lyn={len(W)}  reduced denominators < |gap| (count): {small}   full-gap words: {cnt[abs(D)]}")
    if (L, k) == (27, 17):
        for x in ex[5]:
            n = x * 5
            assert n.denominator == 1
            n = int(n)
            # 3x+5 cycle containing n; minimum element
            orb = [n]; y = n
            while True:
                y = y // 2 if y % 2 == 0 else (3 * y + 5) // 2
                if y == n:
                    break
                orb.append(y)
            print(f"    3x+5 cycle from rational {x}: min {min(orb)}, period {len(orb)}, odd {sum(1 for t in orb if t % 2)}")
    if (L, k) == (19, 12):
        mins = []
        for x in ex[23]:
            n = int(x * 23)
            orb = [n]; y = n
            while True:
                y = y // 2 if y % 2 == 0 else (3 * y + 23) // 2
                if y == n:
                    break
                orb.append(y)
            assert len(orb) == 19
            mins.append(max(orb))  # closest to zero (negative cycle)
        print(f"    3x+23 cycles (on NEGATIVE integers) at (19,12): elements closest to 0 = {sorted(mins)}")
