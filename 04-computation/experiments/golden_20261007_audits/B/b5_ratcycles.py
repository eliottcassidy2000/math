#!/usr/bin/env python3
"""THM-4591 statement 3: E(L,k) = Lyn(L,k)/|2^L - 3^k| ranking over 12 <= L <= 1500 (ALL k), and the
rational-cycle denominators realised at the Ellison shapes (own enumeration)."""
from math import lgamma, log, gcd
from fractions import Fraction
from collections import Counter
from b3_lyndon import lyn_count, lyndon_words, cycle_value

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

