#!/usr/bin/env python3
"""Exact finite controls for guarded observational equivalence (stdlib only).
Universal assertions are proved in GUARDED_OBSERVERS.md, not inferred here.
"""
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod
import json


def require(condition, message):
    """Checks remain active under python -O as well as normal execution."""
    if not condition:
        raise AssertionError(message)


def factors(n):
    require(n >= 1, 'verification failed: n >= 1')
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = 1
    return out


def data(n):
    f = factors(n)
    return {
        'tau': prod(e + 1 for e in f.values()),
        'omega': len(f),
        'sigma': prod(sum(p**j for j in range(e + 1)) for p, e in f.items()),
        'rad': prod(f),
        'squarefree': all(e == 1 for e in f.values()),
        'prime': len(f) == 1 and next(iter(f.values())) == 1,
    }


def fus(n):
    require(n >= 2, 'verification failed: n >= 2')
    d = data(n)
    return (d['tau'] - 2,
            d['omega'] - int(d['prime']),
            2**d['omega'] - 1 - int(d['squarefree']))


def literal_fus(n):
    ds = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ds.update((d, n // d))
    ds.discard(1)
    ds.discard(n)
    return len(ds), sum(data(d)['prime'] for d in ds), sum(data(d)['squarefree'] for d in ds)


def delta(n):
    f, u, s = fus(n)
    return f - u - s


def aliquot(n):
    return data(n)['sigma'] - n


def v(n):
    d = data(n)
    return (d['tau'], 2**d['omega'], d['omega'], 1)


def det(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    out = Fraction(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return 0
        if j != i:
            a[i], a[j] = a[j], a[i]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i + 1, len(a)):
            ratio = a[j][i] / pivot
            for k in range(i + 1, len(a)):
                a[j][k] -= ratio * a[i][k]
    require(out.denominator == 1, 'verification failed: out.denominator == 1')
    return out.numerator


def main():
    limit = 5000
    balanced = []
    for n in range(2, limit + 1):
        require(fus(n) == literal_fus(n), n)
        exps = sorted(factors(n).values())
        expected = exps in ([1], [3], [1, 1, 2])
        require((delta(n) == 0) == expected, n)
        if expected:
            balanced.append(n)
        if not data(n)['squarefree']:
            t, q, w, _ = v(n)
            require(delta(n) == t - q - w - 1, 'verification failed: delta(n) == t - q - w - 1')

    coprime_pairs = nonsquarefree_pairs = 0
    for a, b in product(range(1, 121), repeat=2):
        if gcd(a, b) != 1:
            continue
        coprime_pairs += 1
        require(aliquot(a*b) == data(a)['sigma']*data(b)['sigma'] - a*b, "verification failed: aliquot(a*b) == data(a)['sigma']*data(b)['sigma'] - a*b")
        if not data(a)['squarefree']:
            nonsquarefree_pairs += 1
            ta, qa, wa, _ = v(a)
            tb, qb, wb, _ = v(b)
            require(delta(a*b) == ta*tb - qa*qb - wa - wb - 1, 'verification failed: delta(a*b) == ta*tb - qa*qb - wa - wb - 1')
            require(v(a*b) == (ta*tb, qa*qb, wa+wb, 1), 'verification failed: v(a*b) == (ta*tb, qa*qb, wa+wb, 1)')

    # A genuine coprime rectangle, including the unit convention s(1)=0.
    a2, b2 = [1, 2], [1, 3]
    k2 = [[aliquot(a*b) for b in b2] for a in a2]
    require(k2 == [[0, 1], [1, 6]] and det(k2) == -1, 'verification failed: k2 == [[0, 1], [1, 6]] and det(k2) == -1')
    # Prime blocks {2,3,5} and {7,11,13}; all rows and products nonsquarefree.
    a4 = [4, 8, 12, 60]
    b4 = [49, 343, 539, 7007]
    require(all(gcd(a, b) == 1 for a in a4 for b in b4), 'verification failed: all(gcd(a, b) == 1 for a in a4 for b in b4)')
    k4 = [[delta(a*b) for b in b4] for a in a4]
    expected = [[2,5,6,15], [5,9,12,27], [6,12,15,34], [15,27,34,73]]
    require(k4 == expected and det(k4) == 4, 'verification failed: k4 == expected and det(k4) == 4')
    require(det([v(a) for a in a4]) == det([v(b) for b in b4]) == 2, 'verification failed: det([v(a) for a in a4]) == det([v(b) for b in b4]) == 2')
    c4 = [[1,0,0,0], [0,-1,0,0], [0,0,0,-1], [0,0,-1,-1]]
    require(det(c4) == 1, 'verification failed: det(c4) == 1')

    # Same current observation, different common legal continuation.
    require((delta(4), delta(12), delta(20), delta(60)) == (-1,-1,-1,0), 'verification failed: (delta(4), delta(12), delta(20), delta(60)) == (-1,-1,-1,0)')
    require((delta(8), delta(60), delta(56), delta(420)) == (0,0,1,3), 'verification failed: (delta(8), delta(60), delta(56), delta(420)) == (0,0,1,3)')
    require((aliquot(2), aliquot(3), aliquot(10), aliquot(15)) == (1,1,8,9), 'verification failed: (aliquot(2), aliquot(3), aliquot(10), aliquot(15)) == (1,1,8,9)')
    # Equal full arithmetic features still cannot decide named coprimality guards.
    require(v(4) == v(9) and (gcd(4,2), gcd(9,2)) == (2,1), 'verification failed: v(4) == v(9) and (gcd(4,2), gcd(9,2)) == (2,1)')
    # (rad,tau) deliberately identifies 12 and 18 under coprime multiplication.
    require((data(12)['rad'], data(12)['tau']) == (data(18)['rad'], data(18)['tau']), "verification failed: (data(12)['rad'], data(12)['tau']) == (data(18)['rad'], data(18)['tau'])")
    for b in range(1, 121):
        require((gcd(12,b) == 1) == (gcd(18,b) == 1), 'verification failed: (gcd(12,b) == 1) == (gcd(18,b) == 1)')
        if gcd(12,b) == 1:
            require(delta(12*b) == delta(18*b), 'verification failed: delta(12*b) == delta(18*b)')
    # Enlarging the operation family to overlapping multiplication breaks it.
    require((delta(24), delta(36)) == (1,2), 'verification failed: (delta(24), delta(36)) == (1,2)')
    # A fresh prime plus the current aliquot output reconstructs the source.
    fresh_checks = 0
    for n in range(1, 501):
        for r in (2,3,5,7,11):
            if n % r:
                require(aliquot(n*r) - (r+1)*aliquot(n) == n, 'verification failed: aliquot(n*r) - (r+1)*aliquot(n) == n')
                fresh_checks += 1

    # Full-domain scalar defect: one fresh-prime continuation detects boundaries.
    boundary_checks = 0
    for n in range(2, 501):
        d = data(n)
        for r in (2,3,5,7,11):
            if n % r:
                c = 2*delta(n) + d['omega'] - delta(n*r)
                require(c == int(d['squarefree']) + 2*int(d['prime']), n)
                require(c in (0,1,3), n)
                boundary_checks += 1

    # Enlarging the context language with named exponent promotions is injective.
    promotion_checks = 0
    for n in range(2, 501):
        d = data(n)
        if d['squarefree']:
            continue
        for p, e in factors(n).items():
            increment = delta(p*n) - delta(n)
            require(increment == d['tau'] // (e+1), "verification failed: increment == d['tau'] // (e+1)")
            require(d['tau'] // increment - 1 == e, "verification failed: d['tau'] // increment - 1 == e")
            promotion_checks += 1
    require((data(972)['rad'], data(972)['tau'], delta(972)) == (6,18,11), "verification failed: (data(972)['rad'], data(972)['tau'], delta(972)) == (6,18,11)")
    require((data(288)['rad'], data(288)['tau'], delta(288)) == (6,18,11), "verification failed: (data(288)['rad'], data(288)['tau'], delta(288)) == (6,18,11)")
    require((delta(1944), delta(576)) == (17,14), 'verification failed: (delta(1944), delta(576)) == (17,14)')

    result = {
        'status': 'PASS',
        'literal_FUS_and_balanced_classification_n': [2, limit],
        'balanced_count_in_window': len(balanced),
        'coprime_pairs_a_b_le_120': coprime_pairs,
        'nonsquarefree_left_pairs_a_b_le_120': nonsquarefree_pairs,
        'fresh_prime_reconstruction_checks': fresh_checks,
        'full_domain_FUS_boundary_classifier_checks': boundary_checks,
        'prime_promotion_exponent_reconstruction_checks': promotion_checks,
        'aliquot_rank_2_minor': {'rows':a2, 'columns':b2, 'matrix':k2, 'determinant':det(k2)},
        'nonsquarefree_FUS_rank_4_minor': {'rows':a4, 'columns':b4, 'matrix':k4, 'determinant':det(k4)},
        'scope': 'Finite exact controls; universal proofs are in GUARDED_OBSERVERS.md.'
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
