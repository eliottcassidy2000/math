#!/usr/bin/env python3
"""Exact finite controls for the all-poles/core resummation."""
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial, gcd, lcm, isqrt
from pathlib import Path
import json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def partitions(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n + 1):
        for tail in partitions(n - d, d):
            yield ((d,) + tail)

def oddparts(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n + 1, 2):
        for tail in oddparts(n - d, d):
            yield ((d,) + tail)

def centralizer(mu):
    z = 1
    for d, c in Counter(mu).items():
        z *= d ** c * factorial(c)
    return z

def phi(n):
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            result -= result // p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result -= result // n
    return result

def divisors(n):
    ans = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ans.append(d)
            if d * d != n:
                ans.append(n // d)
    return ans

def edge_orbits_cyclic_average(parts):
    order = lcm(*parts) if parts else 1
    total = sum((phi(order // d) * comb(sum((x for x in parts if d % x == 0)), 2) for d in divisors(order)))
    require(total % order == 0, 'total % order == 0')
    return total // order

def normalized_full(n):
    return sum((F(factorial(n) * 2 ** edge_orbits_cyclic_average(mu), centralizer(mu) * 2 ** comb(n, 2)) for mu in oddparts(n)), F(0))

def core(mu):
    m = sum(mu)
    c = len(mu)
    a = m - c
    e = sum(((x - 1) // 2 for x in mu)) + sum((gcd(x, y) for i, x in enumerate(mu) for y in mu[i + 1:]))
    loss = comb(m, 2) - e
    require(a % 2 == 0 and 3 * a >= 2 * m and (3 * loss >= m * (m - 1)), 'a % 2 == 0 and 3 * a >= 2 * m and (3 * loss >= m * (m - 1))')
    return (m, a, F(factorial(m), centralizer(mu) * 2 ** loss))

def core_coefficient(mu, n):
    m, a, C = core(mu)
    return C * comb(n, m) * F(1, 2 ** (a * (n - m))) if n >= m else F(0)
N = 35
cores = [mu for m in range(N + 1) for mu in oddparts(m, 3)]
for n in range(N + 1):
    full = normalized_full(n)
    extracted = sum((core_coefficient(mu, n) for mu in cores if sum(mu) <= n), F(0))
    require(full == extracted, 'full == extracted')
    raw = full * F(2 ** comb(n, 2), factorial(n))
    require(raw.denominator == 1, 'raw.denominator == 1')
banks = []
for k in range(11):
    ds = list(partitions(k))
    bank = [tuple((2 * d + 1 for d in p)) for p in ds]
    exhaustive = [mu for m in range(3 * k + 1) for mu in oddparts(m, 3) if core(mu)[1] == 2 * k]
    require(set(bank) == set(exhaustive), 'set(bank) == set(exhaustive)')
    require([mu for mu in bank if sum(mu) == 3 * k] == [(3,) * k], '[mu for mu in bank if sum(mu) == 3 * k] == [(3,) * k]')
    top = (3,) * k
    m, a, C = core(top)
    amplitude = C * (4 ** k) ** m
    expected = F(factorial(3 * k), 3 ** k * factorial(k)) * 2 ** (3 * k * k + k)
    require(amplitude == expected, 'amplitude == expected')
    banks.append({'k': k, 'pole': 4 ** k, 'order': 3 * k + 1, 'core_types': len(bank), 'leading_amplitude': str(amplitude)})
require(core((3,)) == (3, 2, F(1, 2)), 'core((3,)) == (3, 2, F(1, 2))')
require(core((5,)) == (5, 4, F(3, 32)), 'core((5,)) == (5, 4, F(3, 32))')
require(core((3, 3)) == (6, 4, F(5, 128)), 'core((3, 3)) == (6, 4, F(5, 128))')
for n in range(N + 1):
    require(core_coefficient((3,), n) == F(32 * comb(n, 3), 4 ** n) if n >= 3 else core_coefficient((3,), n) == 0, 'core_coefficient((3,), n) == F(32 * comb(n, 3), 4 ** n) if n >= 3 else core_coefficient((3,), n) == 0')
    actual = core_coefficient((5,), n) + core_coefficient((3, 3), n)
    expected = F(98304 * (comb(n, 5) if n >= 5 else 0) + 655360 * (comb(n, 6) if n >= 6 else 0), 16 ** n)
    require(actual == expected, 'actual == expected')
report = {'status': 'PASS', 'coefficient_n_max': N, 'core_count': len(cores), 'independent_orbit_count': 'cyclic subgroup Burnside vs pairwise gcd', 'banks': banks, 'claims_not_proved_by_finite_checks': ['normal convergence', 'infinite pole set', 'non-P-recursiveness']}
Path(__file__).with_name('all_poles_verification.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
