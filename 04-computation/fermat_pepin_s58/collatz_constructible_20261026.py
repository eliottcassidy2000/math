#!/usr/bin/env python3
"""Collatz meets the constructible polygons (opus-2026-10-26-S58).

n >= 3 is constructible (Gauss-Wantzel) iff n = 2^a * (product of distinct Fermat primes) (OEIS A003401: 3, 4, 5, 6, 8,
10, 12, 15, 16, 17, 20, 24, 30, 32, 34, ...).  Odd constructible numbers below 2^(2^33) are exactly the 32 divisors of
2^32 - 1 = 3*5*17*257*65537 (F_5..F_32 are composite; any further Fermat prime exceeds 2^(2^33)).
Checks:
 (T) the Collatz trunk (4^k - 1)/3 (the odd numbers whose 3n+1 is a power of 2) meets the constructible numbers
     exactly at k in {1, 2, 4, 8, 16}: (4^k-1)/3 = F_1 ... F_(j-1) for k = 2^(j-1), unconditionally (F_5 is composite);
 (I) initial-segment principle: 2^k - 1 is constructible iff k in {1,2,4,8,16,32}; 2^k + 1 iff it is a Fermat prime;
 (C) the odd constructible n (dividing 2^32 - 1) with 3n+1 constructible, and the Collatz-closed constructible set:
     n whose whole Collatz orbit stays constructible.
Usage: python collatz_constructible_20261026.py
"""
import json, os, math
from sympy import factorint

HERE = os.path.dirname(os.path.abspath(__file__))
FP = [3, 5, 17, 257, 65537]
M32 = 2 ** 32 - 1


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def constructible_below_bound(n):
    """n < 2^(2^33): constructible iff odd part divides 2^32 - 1 (and n >= 3 or n in {1, 2})."""
    return M32 % odd_part(n) == 0


if __name__ == '__main__':
    out = {}
    # (T)
    trunk_hits = [k for k in range(1, 4097) if constructible_below_bound((4 ** k - 1) // 3)]
    assert trunk_hits == [1, 2, 4, 8, 16], trunk_hits
    vals = [(4 ** k - 1) // 3 for k in trunk_hits]
    prods = [math.prod(FP[1:j]) for j in range(1, 6)]
    assert vals == prods, (vals, prods)
    out['trunk_constructible_k'] = trunk_hits
    out['trunk_constructible_values'] = vals
    # (I)
    out['2^k-1_constructible_k<=4096'] = [k for k in range(1, 4097) if constructible_below_bound(2 ** k - 1)]
    out['2^k+1_constructible_k<=4096'] = [k for k in range(1, 4097) if constructible_below_bound(2 ** k + 1)]
    # (C)
    odd_c = sorted(math.prod(FP[i] for i in range(5) if s >> i & 1) for s in range(32))
    step_ok = [n for n in odd_c if constructible_below_bound(3 * n + 1)]
    out['odd_constructible_with_3n+1_constructible'] = step_ok
    # Collatz-closed: iterate the Collatz map from n; closed iff every term constructible (terms stay < 2^(2^33))
    closed = []
    for n in odd_c:
        x, ok, path = n, True, [n]
        while x != 1:
            x = x // 2 if x % 2 == 0 else 3 * x + 1
            path.append(x)
            if not constructible_below_bound(x):
                ok = False
                break
        if ok:
            closed.append((n, path))
    out['collatz_closed_odd_roots'] = [n for n, _ in closed]
    out['collatz_closed_paths'] = {n: p for n, p in closed}
    for n, p in closed:
        print(n, '->', p if len(p) < 12 else p[:6] + ['...'] + p[-3:])
    print('trunk constructible k:', trunk_hits, vals)
    print('2^k-1 constructible:', out['2^k-1_constructible_k<=4096'], ' 2^k+1 constructible:', out['2^k+1_constructible_k<=4096'])
    print('odd constructible n with 3n+1 constructible:', step_ok)
    json.dump(out, open(os.path.join(HERE, 'collatz_constructible.json'), 'w'), indent=1)
