#!/usr/bin/env python3
"""Finite exact controls for power-depth-rough-pairs.md; no onset claims."""
from fractions import Fraction
from math import gcd, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def primes(n):
    a = bytearray(b'\1') * (n + 1)
    a[:2] = b'\0\0'
    for p in range(2, int(n**0.5) + 1):
        if a[p]:
            a[p*p:n+1:p] = b'\0' * ((n-p*p)//p + 1)
    return [p for p in range(2, n+1) if a[p]]


def rho(p, m):
    return 2 if m % p == 0 else 3 if (9*m*m-1) % p == 0 else 4


def forms(a, m):
    return (6*a-1, 6*a+1, 6*(m-a)-1, 6*(m-a)+1)


def root(a, m, d):
    return prod(forms(a, m)) % d == 0


def small_clear(a, m, q):
    return gcd(prod(forms(a, m)), q) == 1


def factor_squarefree(d):
    out = []
    for p in primes(d):
        if d % p == 0:
            out.append(p)
            assert d % (p*p) != 0
    assert prod(out) == d
    return out


def main():
    local_cases = 0
    for p in primes(97):
        if p < 5:
            continue
        for m in range(p):
            actual = sum(root(a, m, p) for a in range(p))
            assert actual == rho(p, m)
            local_cases += 1

    q = 35  # illustrative fixed pre-sieve, not an asserted asymptotic cutoff
    ds = [1, 11, 13, 17, 11*13, 11*17, 13*17, 11*13*17]
    factors = {d: factor_squarefree(d) for d in ds}
    crt_cases = 0
    for m in range(1, 25):
        r0 = (5-rho(5, m)) * (7-rho(7, m))
        assert sum(small_clear(a, m, q) for a in range(q)) == r0
        for d in ds:
            rd = prod(rho(p, m) for p in factors[d])
            actual = sum(small_clear(a, m, q) and root(a, m, d)
                         for a in range(q*d))
            assert actual == r0 * rd
            crt_cases += 1

    discrepancy_cases = 0
    max_relative_error = Fraction(0)
    for m in range(31, 250):
        r0 = (5-rho(5, m)) * (7-rho(7, m))
        aa = range(m//3 + 1, (2*m)//3 + 1)
        for a in aa:
            assert all(1 < n < 6*m for n in forms(a, m))
        for d in ds:
            rd = prod(rho(p, m) for p in factors[d])
            actual = sum(small_clear(a, m, q) and root(a, m, d) for a in aa)
            expected = Fraction(m*r0*rd, 3*q*d)
            error = abs(actual-expected)
            assert error <= r0 * rd
            max_relative_error = max(max_relative_error, error/(r0*rd))
            discrepancy_cases += 1

    # The stronger printed dimension factor needed for the primary-source
    # specialization: (1-g)/(1-2g) <= (p-4)/(p-8).
    factor_cases = 0
    for p in primes(997):
        if p < 11:
            continue
        for r in (2, 3, 4):
            g = Fraction(r, p)
            assert 0 < g < Fraction(1, 2)
            assert (1-g)/(1-2*g) <= Fraction(p-4, p-8)
            factor_cases += 1

    # Check the finite Euler identity used to translate the published lemma.
    pp = (11, 13, 17, 19)
    divisors = [prod(pp[i] for i in range(len(pp)) if mask >> i & 1)
                for mask in range(1 << len(pp))]
    lam = {d: Fraction((7*i % 11)-5, 5) for i, d in enumerate(divisors)}
    gs = {p: Fraction(4, p) for p in pp}
    g = {d: prod(gs[p] for p in pp if d % p == 0) for d in divisors}
    h = {d: prod(gs[p]/(1-gs[p]) for p in pp if d % p == 0)
         for d in divisors}
    theta = {b: sum(lam[d] for d in divisors if b % d == 0) for b in divisors}
    v = prod(1-gs[p] for p in pp)
    assert sum(lam[d]*g[d] for d in divisors) == v*sum(theta[b]*h[b] for b in divisors)

    out = {
        'status': 'all finite exact controls passed',
        'local_root_cases': local_cases,
        'mixed_crt_cases': crt_cases,
        'half_open_interval_discrepancy_cases': discrepancy_cases,
        'maximum_observed_error_divided_by_r0_rho': str(max_relative_error),
        'stronger_dimension_factor_cases': factor_cases,
        'finite_euler_identity': 'exact rational equality passed',
        'not_verified_by_this_script': [
            'the cited analytic fundamental lemma',
            'an explicit asymptotic cutoff w0 or m0',
            'prime representations or a four-prime conjecture',
        ],
    }
    text = '\n'.join(f'{k}: {v}' for k, v in out.items()) + '\n'
    (HERE/'power_depth_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    (HERE/'power_depth_checks.txt').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
