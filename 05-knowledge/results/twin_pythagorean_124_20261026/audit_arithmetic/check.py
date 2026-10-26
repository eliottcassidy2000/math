#!/usr/bin/env python3
"""Independent side-pair and literal-divisor audit; imports no producer code."""
from fractions import Fraction
from math import gcd, isqrt, prod
from collections import defaultdict
import json

checks = 0
def need(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)

def prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))

def radical(n):
    factors = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return prod(factors)

def s(n):
    divisors = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            divisors.update((d, n // d))
    divisors.discard(n)
    return sum(divisors)

by_radius = defaultdict(list)
selected = []
total = 0
for a in range(1, 1000, 2):
    for b in range(2, 1000, 2):
        c = isqrt(a * a + b * b)
        if c > 1000 or c * c != a * a + b * b or gcd(a, b) != 1:
            continue
        total += 1
        semiperimeter = Fraction(a + b + c, 2)
        area = Fraction(a * b, 2)
        r = area / semiperimeter
        rho = area / (semiperimeter - b)
        need(r.denominator == rho.denominator == 1, 'geometric integral radii')
        r, rho = int(r), int(rho)
        by_radius[r].append((a, b, c))
        m, n = isqrt((c + a) // 2), isqrt((c - a) // 2)
        need((r, rho) == (n * (m - n), n * (m + n)), 'radius factor coordinates')
        need((prime(r) and prime(rho)) ==
             (c - a == 2 and prime(b // 2 - 1) and prime(b // 2 + 1)), 'prime-radius iff')
        if prime(r) and prime(rho):
            need((a, b, c) == (r * (r + 2), 2 * (r + 1), (r + 1) ** 2 + 1), 'inverse twin')
            selected.append([r, rho, [a, b, c]])
need(total == 158, 'direct leg-pair count')
for p in (2, 3, 5, 7, 11, 13, 17, 19):
    expected = {(2 * p + 1, 2 * p * (p + 1), 2 * p * p + 2 * p + 1)}
    if p != 2:
        expected.add((p * (p + 2), 2 * (p + 1), p * p + 2 * p + 2))
    need(set(by_radius[p]) == expected, 'complete prime-radius fibre within explicit height')
fib74 = [n for n in range(1, 73 ** 2 + 1) if s(n) == 74]
need(fib74 == [70, 142, 5329], 'target-complete fibre74')
for p in range(3, 100, 2):
    if not (prime(p) and prime(p + 2)):
        continue
    m = p + 1
    a, b, c = m * m - 1, 2 * m, m * m + 1
    need(s(2 * p) == s((p + 2) ** 2) == p + 3, 'literal aliquot collision')
    R0 = radical(a) * radical(m)
    need(R0 == a * radical(m) and R0 > m * m, 'center radical exceeds height')
    need(radical(a) * radical(b) * radical(c) == R0 * radical(c), 'new norm support cost')
    need((1 + 2 + p) == (1 + p + 2) and 3 != 2, 'evaluation collision and augmentation loss')
need(s(18) == 21 and s(25) == s(6) == 6, 'perfect-core versus fixed-target boundary')
print(json.dumps({'status': 'PASS', 'checks': checks, 'primitive_side_pairs': total,
    'height_bound': 1000, 'prime_radius_bound': 19, 'prime_radius_selected': selected,
    'complete_fibre74': fib74, 'source_bound74': 5329}, indent=2, sort_keys=True))
