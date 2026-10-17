#!/usr/bin/env python3
"""Why 463 = 1 + 462 is prime: Euler's idoneal criterion, checked on all 65 idoneal numbers (S40).

Euler: if n is idoneal and an odd m > 1 prime to n has exactly one proper representation m = x^2 + n y^2 (x, y >= 0,
gcd(x, ny) = 1), then m is a prime or a prime power (twice one of these if m is even). For m = n + 1 the
representation 1^2 + n 1^2 is the only proper one (y = 1 forces x = 1; y = 0 is improper unless m = 1). Hence n + 1 is
a prime power or twice a prime power for every idoneal n; for n = 462 this is the prime 463, the norm of the marked
ideal (1 + sqrt(-462)) that generates the 1848 conductor clock (prime463_conductor_clock_20261017e.md). The script
checks this for all 65 idoneal numbers, checks that the proper representation of n + 1 is unique, contrasts with the
non-idoneal n <= 2000, and replays Euler's primality certificate 18518809 = 197^2 + 1848 * 100^2.
Usage: python idoneal_euler_criterion_463_20261017.py
Session opus-2026-10-17-S40.
"""
import math
from sympy import factorint, isprime

IDONEAL = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57, 58,
           60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120, 130, 133, 165, 168, 177, 190, 210, 232, 240, 253, 273, 280,
           312, 330, 345, 357, 385, 408, 462, 520, 760, 840, 1320, 1365, 1848]


def proper_reps(m, n):
    out = []
    y = 0
    while n * y * y <= m:
        r = m - n * y * y
        x = math.isqrt(r)
        if x * x == r and math.gcd(x, n * y) == 1:
            out.append((x, y))
        y += 1
    return out


def kind(m):
    f = factorint(m)
    if len(f) == 1:
        return 'prime' if list(f.values())[0] == 1 else 'prime power'
    if m % 2 == 0 and len(f) == 2 and f[2] == 1:
        return 'twice a prime' if list(factorint(m // 2).values())[0] == 1 else 'twice a prime power'
    return 'composite'


def main():
    bad = []
    rows = []
    for n in IDONEAL:
        m = n + 1
        reps = proper_reps(m, n)
        k = kind(m) if m > 1 else 'unit'
        rows.append((n, m, k, reps))
        if k == 'composite' or len(reps) != 1:
            bad.append((n, m, k, reps))
    print('idoneal n with n + 1 not a prime power / twice one, or with a non-unique proper representation:', bad)
    print('n + 1 for the idoneal numbers that are 2 x odd squarefree (clock_duality P1 family):')
    print('   ', [(n, n + 1, kind(n + 1)) for n in IDONEAL if n % 4 == 2 and all(e == 1 for e in factorint(n // 2).values())])
    print('n + 1 square:', [(n, n + 1, math.isqrt(n + 1)) for n in IDONEAL if math.isqrt(n + 1) ** 2 == n + 1])
    nonid = [n for n in range(1, 2001) if n not in IDONEAL]
    fail = [n for n in nonid if kind(n + 1) == 'composite']
    print('non-idoneal n <= 2000: %d; n + 1 composite (not a prime power or twice one): %d (the criterion detects these at m = n + 1)'
          % (len(nonid), len(fail)))
    m = 197 ** 2 + 1848 * 100 ** 2
    print("Euler's certificate: 197^2 + 1848*100^2 = %d; proper representations by x^2 + 1848 y^2: %s; prime: %s"
          % (m, proper_reps(m, 1848), isprime(m)))


if __name__ == '__main__':
    main()
