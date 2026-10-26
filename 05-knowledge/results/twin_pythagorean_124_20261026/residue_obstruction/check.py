#!/usr/bin/env python3
"""Explicit CRT hostile family; no imports from the producer lanes."""
import json
from math import gcd, isqrt

checks = 0


def require(test, message):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def triple(p):
    return (p * (p + 2), 2 * (p + 1), (p + 1) ** 2 + 1)


def hostile(p, modulus):
    candidates = []
    q = 3
    while len(candidates) < 2:
        if prime(q) and (2 * modulus) % q:
            candidates.append(q)
        q += 2
    q, r = candidates
    # x=p+2M*t, first solve t modulo q and r.
    tq = -p * pow(2 * modulus, -1, q) % q
    tr = -(p + 2) * pow(2 * modulus, -1, r) % r
    t = tq + q * ((tr - tq) * pow(q, -1, r) % r)
    x = p + 2 * modulus * (t + q * r)
    return x, q, r


twins = [p for p in range(3, 501, 2) if prime(p) and prime(p + 2)]
examples = []
for p in twins:
    for modulus in range(1, 121):
        x, q, r = hostile(p, modulus)
        a, b, c = triple(x)
        require(x > p and x > q and x + 2 > r, 'proper nontrivial factors')
        require(x % q == 0 and (x + 2) % r == 0, 'two composite certificates')
        require(x % 2 == 1 and (x - p) % (2 * modulus) == 0, 'odd CRT lift')
        require(((x - 1) // 2 - (p - 1) // 2) % modulus == 0, 'depth retained')
        require(a * a + b * b == c * c and gcd(a, b) == 1, 'primitive Euclid lift')
        require(tuple(z % modulus for z in triple(p)) ==
                tuple(z % modulus for z in (a, b, c)), 'coordinate reduction retained')
        require((a + b - c) // 2 == x and (a - b + c) // 2 == x * (x + 1),
                'inradius and other exradius are distinguished')
        require((-a + b + c) // 2 == x + 2, 'even-leg exradius')
        if p in (3, 5, 11) and modulus in (1, 7, 14, 30):
            examples.append({'p': p, 'M': modulus, 'x': x, 'q': q, 'r': r})

x, q, r = hostile(3, 7)
require((x, q, r) == (213, 3, 5), 'named CRT witness')
require(triple(x) == (45795, 428, 45797), 'named triple')
require(tuple(z % 7 for z in triple(x)) == (1, 1, 3), 'named conic reduction')
require(pow(x + 1, -1, 7) == pow(4, -1, 7) == 2, 'same marked slope')
print(json.dumps({'status': 'PASS', 'checks': checks,
    'universe': {'twin_lower_prime_max': 500, 'twin_pairs': len(twins),
                 'moduli': [1, 120], 'CRT_instances': len(twins) * 120},
    'examples': examples, 'F7_hostile': {'p': 3, 'x': x, 'factors': [q, r],
        'twin_triple': triple(3), 'composite_triple': triple(x),
        'common_reduction': [1, 1, 3], 'common_slope': 2},
    'scope': 'Fixed-modulus loss only; all-order theorem uses CRT proof.'},
    indent=2, sort_keys=True))
