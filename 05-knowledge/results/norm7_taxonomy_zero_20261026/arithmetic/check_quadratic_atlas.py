#!/usr/bin/env python3
"""Exact controls for the quadratic residue-field exceptional-conic atlas.

Standard library only; no producer imports. The global splitting theorem
is cited in the proof, not inferred from these finite computations.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt
import json

GATES = 0


def check(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)


def factor(n):
    n = abs(n)
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def squareclass(value):
    f = Fraction(value)
    if not f:
        raise ValueError('twists must be nonzero')
    n = f.numerator*f.denominator
    out = 1 if n > 0 else -1
    for p, e in factor(n).items():
        if e % 2:
            out *= p
    return out


def rational_square(x):
    x = Fraction(x)
    return (x >= 0 and isqrt(x.numerator)**2 == x.numerator
            and isqrt(x.denominator)**2 == x.denominator)


def ramified(delta):
    q = squareclass(delta)
    u = q//2 if q % 2 == 0 else q
    answer = {p for p in factor(q) if p % 4 == 3}
    if u % 4 == 3:
        answer.add(2)
    if q < 0:
        answer.add('infinity')
    return answer


def split_place(d, p):
    d = squareclass(d)
    if p == 'infinity':
        return d > 0
    if p == 2:
        return d % 8 == 1
    return d % p != 0 and pow(d % p, (p-1)//2, p) == 1


def atlas(D, delta):
    d = squareclass(D)
    R = ramified(delta)
    obstruction = R if d == 1 else {p for p in R if split_place(d, p)}
    degree = (1 if d == 1 else 2)*(2 if obstruction else 1)
    return obstruction, degree


def primitive_conic_mod(delta, p):
    """Independent exact local test on squarefree integral delta.

    Modulus8 at2, p² otherwise. Normalize a unit coordinate to1;
    A/B symmetry leaves only two charts. Every returned solution is
    primitive; the text explains why these moduli suffice in this scope.
    """
    modulus = 8 if p == 2 else p*p
    sq = [a*a % modulus for a in range(modulus)]
    # C is a unit, normalized to1.
    for a, b in product(range(modulus), repeat=2):
        if (sq[a]+sq[b]-delta) % modulus == 0:
            return (a, b, 1)
    # A or B is a unit. Symmetry lets us normalize A=1.
    for b, c in product(range(modulus), repeat=2):
        if (1+sq[b]-delta*sq[c]) % modulus == 0:
            return (1, b, c)
    return None


def quadratic_square_mod(d, p):
    modulus = 8 if p == 2 else p*p
    return any((x*x-d) % modulus == 0 for x in range(modulus))


def field_add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def field_mul(a, b, d):
    return (a[0]*b[0]+d*a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def conic_witness(d, delta, A, B, C):
    left = field_add(field_mul(A, A, d), field_mul(B, B, d))
    c2 = field_mul(C, C, d)
    right = (delta*c2[0], delta*c2[1])
    return left == right and any(x for V in (A, B, C) for x in V)


def sorted_places(values):
    return sorted(values, key=lambda p: (p == 'infinity', str(p) if p == 'infinity' else p))


def main():
    # Signed rational squareclass normalization and scaling invariance.
    rational_cases = 0
    for n in range(-30, 31):
        if not n:
            continue
        for denominator in range(1, 16):
            x = Fraction(n, denominator)
            q = squareclass(x)
            check(rational_square(x/q), 'squareclass retains exact rational square')
            check(squareclass(x*Fraction(25, 49)) == q, 'rational-square invariance')
            check(len(ramified(x)) % 2 == 0, 'Hilbert reciprocity parity')
            rational_cases += 1

    deltas = [x for x in range(-40, 41) if x and squareclass(x) == x]
    local_cases = 0
    local_rejections = Counter()
    for delta in deltas:
        for p in (2, 3, 5, 7, 11):
            witness = primitive_conic_mod(delta, p)
            check((witness is None) == (p in ramified(delta)), 'primitive local conic versus Hilbert symbol')
            local_cases += 1
            if witness is None:
                local_rejections[p] += 1

    # Independent polynomial-root test for the quadratic splitting places.
    ds = [d for d in range(-100, 101) if d and squareclass(d) == d]
    primes = (2, 3, 5, 7, 11, 19, 23, 31)
    for d in ds:
        for p in primes:
            check(split_place(d, p) == quadratic_square_mod(d, p), 'local squareclass root test')

    # Atlas consequence object, with an independent route through local roots.
    census = Counter()
    delta_three_pass = []
    for d in ds:
        for delta in deltas:
            bad, degree = atlas(d, delta)
            independent_bad = {p for p in ramified(delta)
                               if (d > 0 if p == 'infinity' else quadratic_square_mod(d, p))}
            check(bad == independent_bad, 'atlas obstruction places by roots')
            check(degree in ((1, 2) if d == 1 else (2, 4)), 'possible minimum degrees')
            census[(d == 1, degree)] += 1
        if d != 1:
            bad, degree = atlas(d, 3)
            congruence = d % 8 != 1 and d % 3 != 1
            check((degree == 2) == congruence, 'delta-three closed congruence rule')
            if congruence:
                delta_three_pass.append(d)

    # Exact number-field point witnesses, including a projective infinity point.
    Q = Fraction
    witnesses = [
        (2, 3, (1, 0), (0, 1), (1, 0)),
        (3, 3, (0, 1), (0, 0), (1, 0)),
        (5, 3, (Q(1, 2), Q(1, 2)), (Q(1, 2), Q(-1, 2)), (1, 0)),
        (-1, 3, (1, 0), (0, 1), (0, 0)),
        (13, 2, (1, 0), (1, 0), (1, 0)),
        (-1, -1, (1, 0), (0, 1), (0, 0)),
    ]
    for d, delta, A, B, C in witnesses:
        check(conic_witness(d, delta, A, B, C), 'exact quadratic-field witness')
        check(atlas(d, delta)[1] == 2, 'witness matches complete-fibre minimum')

    expected = {(2, 3): (set(), 2), (3, 3): (set(), 2), (5, 3): (set(), 2),
                (13, 3): ({3}, 4), (17, 3): ({2}, 4), (73, 3): ({2, 3}, 4),
                (-1, 3): (set(), 2), (-7, 3): ({2}, 4),
                (2, -1): ({'infinity'}, 4), (1, 3): ({2, 3}, 2),
                (1, 2): (set(), 1)}
    for params, value in expected.items():
        check(atlas(*params) == value, 'positive or hostile complete-fibre control')
    check(ramified(-1) == {2, 'infinity'}, 'MISTAKE-506 Hamilton control')
    check(primitive_conic_mod(3, 2) is None, 'dyadic hostile has no primitive modular point')
    check(primitive_conic_mod(3, 3) is None, 'ternary hostile has no primitive modular point')
    # Deliberately dropping2 would misclassify D17; dropping infinity D2,-1.
    check(not any(split_place(17, p) for p in ramified(3)-{2}), 'odd-prime-only mutant misses hostile')
    check(not any(split_place(2, p) for p in ramified(-1)-{'infinity'}), 'finite-only mutant misses real hostile')

    rows = [{'D': d, 'delta': delta, 'ramified_places': sorted_places(ramified(delta)),
             'remaining_obstruction_places': sorted_places(atlas(d, delta)[0]),
             'minimum_point_field_degree': atlas(d, delta)[1]}
            for d, delta in expected]
    print(json.dumps({'status': 'FINITE-EXACT controls for proved quadratic-field atlas',
        'gates': GATES, 'passed': True,
        'universe': {'rational_delta': 'numerator -30..30 nonzero; denominator1..15',
                     'squarefree_delta': '-40..40 nonzero', 'squarefree_D': '-100..100 nonzero',
                     'local_conic_primes': [2, 3, 5, 7, 11], 'quadratic_root_primes': primes},
        'rational_squareclass_cases': rational_cases, 'primitive_local_cases': local_cases,
        'primitive_local_rejections': dict(sorted(local_rejections.items())),
        'atlas_degree_census': [{'D_square': flag, 'minimum_degree': degree, 'cases': count}
                               for (flag, degree), count in sorted(census.items())],
        'delta_three_pass_D': delta_three_pass, 'consequence_controls': rows,
        'scope': 'Specified complete reduced resolved fibre only. Classical global splitting theorem is cited, not proved computationally. No whole-surface rational-point or F1 claim.'},
        indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
