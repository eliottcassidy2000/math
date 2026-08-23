#!/usr/bin/env python3
"""Independent bounded corroboration of the THM-3743 -> THM-3793 lane.

The source universe is the exact support-two Graver ratio atlas
  1 <= a < b, gcd(a,b)=1, a+b <= 356.
The candidate address keeps only sums whose prime factors are 2 mod 3 and
occur to exponent at most two, then maps (a,b) to a^3+b^3.

This script verifies injectivity and full positive two-cube singleton fibres
on the bounded atlas, plus a larger nonprimitive control through sum 1000.
Its test puts the exponent cap on the whole pair sum, so it audits a strict
subcase of proved THM-3793, whose cap is only on the primitive quotient.
Nothing here preserves loneliness or supplies a new LRC consequence.
"""

from collections import defaultdict
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        out.append((n, 1))
    return out


def admissible_sum(n):
    fs = factor(n)
    return bool(fs) and all(p % 3 == 2 and e <= 2 for p, e in fs)


def icbrt(n):
    lo, hi = 0, 1
    while hi**3 <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def positive_pair_fibres(max_value):
    max_y = icbrt(max_value)
    fibres = defaultdict(list)
    for y in range(2, max_y + 1):
        y3 = y**3
        for x in range(1, y):
            value = x**3 + y3
            if value <= max_value:
                fibres[value].append((x, y))
    return fibres


def bounded_nonprimitive_control(sum_cap):
    candidates = []
    max_value = 0
    for d in range(3, sum_cap + 1):
        if not admissible_sum(d):
            continue
        for x in range(1, (d + 1) // 2):
            y = d - x
            if x < y:
                value = x**3 + y**3
                candidates.append((value, x, y, d))
                max_value = max(max_value, value)
    fibres = positive_pair_fibres(max_value)
    failures = [(v, x, y, d, fibres[v]) for v, x, y, d in candidates if len(fibres[v]) != 1]
    return candidates, failures


def main():
    ratios = []
    admissible = []
    admissible_sums = set()
    for a in range(1, 356):
        for b in range(a + 1, 357 - a):
            if gcd(a, b) != 1:
                continue
            ratios.append((a, b))
            if admissible_sum(a + b):
                admissible.append((a, b))
                admissible_sums.add(a + b)

    require(len(ratios) == 19314, f"unexpected Graver ratio count {len(ratios)}")
    require(len(admissible_sums) == 94, f"unexpected admissible sum count {len(admissible_sums)}")
    require(len(admissible) == 5855, f"unexpected admissible address count {len(admissible)}")

    values = [a**3 + b**3 for a, b in admissible]
    require(len(values) == len(set(values)), "admissible address map collided internally")
    max_value = max(values)
    fibres = positive_pair_fibres(max_value)
    bad = [(a, b, fibres[a**3 + b**3]) for a, b in admissible if len(fibres[a**3 + b**3]) != 1]
    require(not bad, f"bounded LRC address has nonsingleton fibres: {bad[:3]}")

    all_ratio_fibres = defaultdict(list)
    for a, b in ratios:
        all_ratio_fibres[a**3 + b**3].append((a, b))
    collision_values = {m: ps for m, ps in all_ratio_fibres.items() if len(ps) > 1}
    require(1729 in collision_values, "taxicab hostile 1729 missing")
    require(set(collision_values[1729]) == {(1, 12), (9, 10)}, "1729 fibre changed")

    large_candidates, large_failures = bounded_nonprimitive_control(1000)
    require(not large_failures, f"nonprimitive control failed: {large_failures[:3]}")

    exp3 = 54**3 + 71**3
    require(exp3 == 15**3 + 80**3 == 515375, "exponent-three hostile arithmetic failed")
    require(factor(54 + 71) == [(5, 3)], "exponent-three hostile sum factor changed")
    split = 1**3 + 12**3
    require(split == 9**3 + 10**3 == 1729, "split-prime hostile arithmetic failed")
    require(factor(1 + 12) == [(13, 1)] and factor(9 + 10) == [(19, 1)], "split hostile factors changed")

    print("LRC SUPPORT-TWO / TWO-CUBE ADDRESS AUDIT")
    print(f"THM-3743 ratio universe: {len(ratios)}")
    print(f"admissible inert exponent<=2 sums: {len(admissible_sums)}")
    print(f"admissible ratio addresses: {len(admissible)}")
    print(f"labelled coordinate-pair addresses: {len(admissible) * 78}")
    print(f"largest audited address value: {max_value}")
    print(f"all admissible address fibres singleton in complete positive-pair universe: yes")
    print(f"all-ratio collision values inside cap: {len(collision_values)}")
    print(f"first canonical collision: 1729 -> {sorted(collision_values[1729])}")
    print(f"nonprimitive admissible pairs through sum 1000: {len(large_candidates)}")
    print("nonprimitive singleton failures through sum 1000: 0")
    print("sharp boundary: 515375=(54,71)=(15,80), with 54+71=5^3")
    print("split-prime hostile: 1729=(1,12)=(9,10), sums 13 and 19")
    print("VERDICT: exact bounded reversible address; loneliness, owner, phase, and the other eleven speeds are lost.")
    print("STATUS FIREWALL: bounded corroboration of a strict subcase of proved THM-3793; no LRC consequence.")


if __name__ == "__main__":
    main()
