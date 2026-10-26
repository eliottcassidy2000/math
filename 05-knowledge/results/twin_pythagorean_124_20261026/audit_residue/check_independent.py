#!/usr/bin/env python3
"""Independent CRT audit by literal search, rather than modular inverses."""
from math import gcd, isqrt
import json

checks = 0


def require(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)


def prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def T(x):
    m = x + 1
    return m * m - 1, 2 * m, m * m + 1


twins = [p for p in range(3, 102, 2) if prime(p) and prime(p + 2)]
instances = 0
for M in range(1, 41):
    qs = [q for q in range(3, 200, 2) if prime(q) and gcd(q, 2 * M) == 1]
    q, r = qs[:2]
    for p in twins:
        candidates = [p + 2 * M * k for k in range(q * r)
                      if (p + 2 * M * k) % q == 0
                      and (p + 2 * M * k + 2) % r == 0]
        require(len(candidates) == 1, 'unique literal CRT residue')
        base = candidates[0]
        for lift in (1, 2, 9):
            x = base + 2 * M * q * r * lift
            instances += 1
            require(x > q and x + 2 > r and x % q == (x + 2) % r == 0,
                    'proper composite factors persist under lifts')
            a, b, c = T(x)
            require((x - p) % (2 * M) == 0 and x % 2 == 1,
                    'odd residue and depth modulus')
            require(a * a + b * b == c * c and gcd(a, b) == 1,
                    'independent primitive polynomial lift')
            require(tuple(z % M for z in T(p)) == tuple(z % M for z in T(x)),
                    'triple reductions agree')
            area, semiperimeter = a * b // 2, (a + b + c) // 2
            require(area == x * semiperimeter,
                    'inradius via area/semiperimeter')
            require(area == (x + 2) * (semiperimeter - b),
                    'even-leg exradius via area formula')
            require(((x - 1) // 2 - (p - 1) // 2) % M == 0,
                    'D depth observer agrees')

require(T(213) == (45795, 428, 45797), 'literal named hostile triple')
require((213 - 3) % 14 == 0 and 213 == 3 * 71 and 215 == 5 * 43,
        'literal named factors')
require(tuple(v % 7 for v in T(213)) == tuple(v % 7 for v in T(3)) == (1, 1, 3),
        'literal F7 equality')
print(json.dumps({'status': 'PASS', 'checks': checks,
                  'universe': {'moduli': [1, 40], 'twin_lower_prime_max': 101,
                               'twins': len(twins), 'three_lift_instances': instances},
                  'method': 'Literal search of each CRT class; independent area-derived radii.'},
                 indent=2, sort_keys=True))
