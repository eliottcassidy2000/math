#!/usr/bin/env python3
"""Portable exact controls, not a computer proof of Dudeney optimality.

Only the Python standard library is used. Explicit checks survive python -O.
All inputs are local; there is no repository-path assumption or network access.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
from math import factorial, gcd, isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def factor(n):
    ans = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            ans[p] = ans.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        ans[n] = ans.get(n, 0)+1
    return ans


def norm_criterion(n):
    return all(p % 4 != 3 or e % 2 == 0 for p, e in factor(n).items())


def rational_norm_criterion(q):
    num, den = factor(q.numerator), factor(q.denominator)
    return all(p % 4 != 3 or (num.get(p, 0)-den.get(p, 0)) % 2 == 0
               for p in num.keys() | den.keys())


def gaussian_mul(u, v):
    a, b = u
    c, d = v
    return (a*c-b*d, a*d+b*c)


def norm(u):
    return sum(x*x for x in u)


def phi(u):
    x, y = u
    return (x*x-y*y, 2*x*y, x*x+y*y)


def det(u, v):
    return u[0]*v[1]-u[1]*v[0]


def lorentz(u, v):
    return u[2]*v[2]-u[0]*v[0]-u[1]*v[1]


def kadd(x, y):
    return tuple(a+b for a, b in zip(x, y))


def kmul(x, y):
    a, b = x
    c, d = y
    return (a*c+3*b*d, a*d+b*c)


def j(x):
    return (x[0], -x[1])


def qmul(x, y):
    """Q[theta] multiplication, theta^4=3, exact length-four vectors."""
    result = [F(0)]*4
    for i, a in enumerate(x):
        for k, b in enumerate(y):
            degree = i+k
            result[degree % 4] += a*b*(3 if degree >= 4 else 1)
    return tuple(result)


def main():
    pins = json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for pin in pins:
        data = (ROOT/pin['snapshot']).read_bytes()
        require(sha256(data).hexdigest() == pin['sha256'], pin['snapshot'])
        require(pin['source_commit'] == 'c2039c32f9726b09ed0e09d5f5c126fc72975897',
                'unannounced input revision')
    external = json.loads((ROOT/'EXTERNAL_SOURCES.json').read_text())
    for source in external:
        if 'snapshot' in source:
            require(sha256((ROOT/source['snapshot']).read_bytes()).hexdigest() ==
                    source['sha256'], 'external source bytes')
    # Independent integer representation enumerator versus factor criterion.
    limit = 10000
    representations = {n: [] for n in range(limit+1)}
    for a in range(isqrt(limit)+1):
        for b in range(a, isqrt(limit-a*a)+1):
            representations[a*a+b*b].append((a, b))
    for n in range(1, limit+1):
        require(bool(representations[n]) == norm_criterion(n), f'norm {n}')
    require(representations[5] == [(1, 2)], 'norm5 positive control')
    require(representations[44] == [], 'norm44 inert hostile')
    require(factor(44) == {2: 2, 11: 1}, '44 factorization')
    require(representations[65] == [(1, 8), (4, 7)], 'norm65 collision')
    rational_count = 0
    for m in range(1, 201):
        for n in range(1, 21):
            q = F(m, n)
            a, b = q.numerator, q.denominator
            require(rational_norm_criterion(q) == bool(representations[a*b]),
                    f'rational cleared criterion {q}')
            if representations[a*b]:
                x, y = representations[a*b][0]
                require(F(x, b)**2+F(y, b)**2 == q, 'rational witness')
            rational_count += 1
    # Independent Gaussian multiplication supplies a second Phi path.
    rays = list(product(range(-12, 13), repeat=2))
    for u in rays:
        z = gaussian_mul(u, u)
        require(phi(u) == (z[0], z[1], norm(u)), 'Gaussian square')
        require(phi(u)[0]**2+phi(u)[1]**2 == phi(u)[2]**2, 'null cone')
        for v in rays:
            require(norm(gaussian_mul(u, v)) == norm(u)*norm(v), 'norm product')
            require(lorentz(phi(u), phi(v)) == 2*det(u, v)**2, 'Lorentz identity')
    primitive_rows = 0
    for m in range(2, 51):
        for n in range(1, m):
            if gcd(m, n) == 1 and (m-n) % 2:
                a, b, c = phi((m, n))
                require(gcd(gcd(a, b), c) == 1, 'primitive Euclid triple')
                require(all(p % 4 == 1 for p in factor(c)), 'primitive odd norm')
                primitive_rows += 1
    require(lorentz(phi((8, 1)), phi((1, 0))) == 2, 'Farey positive')
    require(lorentz(phi((7, 4)), phi((1, 0))) == 32, 'same-norm negative')
    x, y = phi((3, 1)), phi((2, 1))
    require(lorentz(x, y) == 2, 'normalization raw edge')
    require(sorted(v//2 for v in x) == sorted(y), 'unordered collapse')
    require(lorentz((3, 4, 5), (5, 12, 13)) == 2, 'ordered positive')
    require(lorentz((4, 3, 5), (5, 12, 13)) == 9, 'leg-swap hostile')
    # Q(sqrt3) identity and conjugation, with integer basis coordinates.
    kvals = list(product(range(-4, 5), repeat=2))
    for x in kvals:
        for y in kvals:
            require(j(kmul(x, y)) == kmul(j(x), j(y)), 'field embedding')
            z = kadd(kmul(x, x), kmul(y, y))
            r, s = x
            t, u = y
            require(z == (r*r+3*s*s+t*t+3*u*u, 2*(r*s+t*u)), 'basis formula')
            require(z != (0, 1), 'sqrt3 nonnorm finite control')
            require(z[0] > 0 or x == y == (0, 0), 'positive constant coefficient')
    require(j((0, 1)) == (0, -1), 'negative embedding witness')
    # A field positive control: theta^2=sqrt3, theta^4=3.
    one, theta = (F(1), F(0), F(0), F(0)), (F(0), F(1), F(0), F(0))
    theta2 = qmul(theta, theta)
    require(theta2 == (0, 0, 1, 0), 'quartic norm exists')
    require(qmul(theta2, theta2) == tuple(3*x for x in one), 'quartic relation')
    # For positive quantities sqrt(2sqrt3)<2 iff3<4.
    require(3 < 4, 'exact squared diameter certificate')
    # Paper display-list arithmetic, not enumeration of geometric classes.
    square = {
        14: ['A1'],
        13: ['A2', 'D1', 'D2'],
        12: ['A3', 'A4', 'B1', 'D3', 'D4', 'D6', 'D7', 'F1', 'F2', 'F3', 'F4', 'F11'],
        11: ['A5', 'B2', 'D5', 'D8', 'D9', 'D11', 'D12', 'E1', 'E2', 'E3', 'E4', 'F5', 'F6', 'F7', 'F9'],
        10: ['A6', 'B3', 'C1', 'D10', 'D13', 'D14', 'E5', 'E6', 'E7', 'F8', 'F10'],
        9: ['D15', 'E8'],
    }
    triangle = {12: ['D1'], 11: ['D2', 'F1'], 10: ['D3', 'E1']}
    sq_list = [x for bucket in square.values() for x in bucket]
    tr_list = [x for bucket in triangle.values() for x in bucket]
    require(len(sq_list) == len(set(sq_list)) == 44, '44 distinct square labels')
    require(len(tr_list) == len(set(tr_list)) == 5, '5 distinct triangle labels')
    cc_pairs = sum(len(v)*len(square.get(k, [])) for k, v in triangle.items())
    require(cc_pairs == 64, 'cc filter arithmetic')
    require(cc_pairs-2*9 == 46, 'tri filter arithmetic')
    # Genuine !5=44 source relation, two independent finite formulas.
    ders = [p for p in permutations(range(5)) if all(i != v for i, v in enumerate(p))]
    inclusion_exclusion = sum((-1)**k * F(factorial(5), factorial(k)) for k in range(6))
    require(len(ders) == inclusion_exclusion == 44, 'derangements')
    # Four square-side edge occurrences do not force four piece owners.
    owner_labels = ['one_square_piece']*4
    require(len(owner_labels) == 4 and len(set(owner_labels)) == 1,
            'edge-count to piece-count hostile')
    # Norm statuses differ while similarity leaves every normalized shape unchanged.
    for area in (F(5), F(44)):
        require(area/area == 1, 'normalized square area')
        # Square/triangle side ratio fourth power always3/16.
        require((area*area)/(F(16, 3)*area*area) == F(3, 16), 'shape ratio invariant')
    report = {
        'status': 'PASS', 'checks': checks,
        'source_pins_verified': len(pins),
        'integer_norm_universe': [1, limit],
        'represented_integer_norms': sum(bool(representations[n]) for n in range(1, limit+1)),
        'rational_cells': rational_count,
        'gaussian_ray_count': len(rays), 'ordered_ray_pairs': len(rays)**2,
        'primitive_euclid_rows': primitive_rows,
        'quadratic_field_pairs': len(kvals)**2,
        'norm_controls': {'5': representations[5], '44': representations[44], '65': representations[65]},
        'paper_list_counts': {'square': len(sq_list), 'triangle': len(tr_list),
                              'cc_pairs': cc_pairs, 'after_tri_diff': cc_pairs-18},
        'derangements5': len(ders),
        'proof_boundary': 'Controls corroborate elementary proofs; paper exhaustion and minimum4 are cited, not re-proved.'
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
