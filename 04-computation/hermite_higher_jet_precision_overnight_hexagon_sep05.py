#!/usr/bin/env python3
"""Exact arbitrary-multiplicity inverse and dyadic three-jet boundary gates.

Universe: complete zero-translated three-node head of diameter<=80;
all four unit classes at scales0..16; deterministic signed/unit-affine
controls; full inverse/Smith controls on prescribed heterogeneous and
uniform observers. Uses exact rational arithmetic throughout. Tests fail
explicitly, including under python -O. This is a proof companion, not a
finite extrapolation to arbitrary node sets or jet orders.
"""

from fractions import Fraction as F
import hashlib
import itertools
import math
import random

from sympy import Matrix, ZZ, binomial
from sympy.matrices.normalforms import smith_normal_form


GATES = 0
DIGEST = hashlib.sha256()


def check(condition, data):
    global GATES
    if not condition:
        raise RuntimeError(repr(data))
    GATES += 1
    DIGEST.update(repr(data).encode())
    DIGEST.update(b'\n')


def vp(a, p):
    if not a:
        return math.inf
    if isinstance(a, F):
        return vp(a.numerator, p) - vp(a.denominator, p)
    out = 0
    while a % p == 0:
        a //= p
        out += 1
    return out


def reciprocal_jet(nodes, multiplicities, i):
    """Product of exact negative-binomial series; no matrix inversion."""
    bound = multiplicities[i]
    coeff = [F(1)] + [F(0)] * (bound - 1)
    for j, y in enumerate(nodes):
        if j == i:
            continue
        d = nodes[i] - y
        m = multiplicities[j]
        factor = [F((-1)**r * math.comb(m+r-1, r), d**(m+r))
                  for r in range(bound)]
        coeff = [sum(coeff[r]*factor[l-r] for r in range(l+1))
                 for l in range(bound)]
    return coeff


def precision(nodes, multiplicities, p):
    return max(-vp(a, p) for i in range(len(nodes))
               for a in reciprocal_jet(nodes, multiplicities, i))


def dyadic_three_prediction(nodes):
    pairs = [(vp(nodes[j]-nodes[i], 2), i, j)
             for i, j in itertools.combinations(range(3), 2)]
    e = min(row[0] for row in pairs)
    depth, i, j = max(pairs)
    d = depth - e
    if d >= 2:
        return 8*e+5*d-1
    z = nodes[3-i-j]
    ratio = F(2*(z-nodes[i]), nodes[j]-nodes[i])
    t = ratio.numerator * pow(ratio.denominator, -1, 16) % 16
    gamma = 3 if t % 4 == 3 else 2 if t % 8 == 5 else 1 if t == 1 else 0
    return max(7*e+4, 8*e+gamma)


def literal_observer(nodes, multiplicities):
    size = sum(multiplicities)
    return Matrix([[binomial(j, r)*x**(j-r) if j >= r else 0
                    for j in range(size)]
                   for x, m in zip(nodes, multiplicities) for r in range(m)])


def audit_matrix(nodes, multiplicities):
    matrix = literal_observer(nodes, multiplicities)
    inverse = matrix.inv()
    smith = smith_normal_form(matrix, domain=ZZ)
    denominator = 1
    column = 0
    for i, m in enumerate(multiplicities):
        jet = reciprocal_jet(nodes, multiplicities, i)
        for r in range(m):
            coefficient = F(inverse[-1, column])
            check(coefficient == jet[m-r-1],
                  ('attainment', nodes, multiplicities, i, r, coefficient))
            column += 1
        for a in jet:
            denominator = math.lcm(denominator, a.denominator)
    inverse_denominator = math.lcm(*(int(a.q) for a in inverse))
    largest = abs(int(smith[-1, -1]))
    check(denominator == inverse_denominator == largest,
          ('exact-global-denominator', nodes, multiplicities, largest))
    for p in [2, 3, 5, 7]:
        check(precision(nodes, multiplicities, p) == vp(largest, p),
              ('prime-denominator', nodes, multiplicities, p, vp(largest, p)))
    return tuple(vp(abs(int(smith[i, i])), 2) for i in range(matrix.rows))


def main():
    head = 0
    for a, b in itertools.combinations(range(1, 81), 2):
        nodes = (0, a, b)
        actual = precision(nodes, (3, 3, 3), 2)
        predicted = dyadic_three_prediction(nodes)
        check(actual == predicted, ('complete-head', nodes, actual, predicted))
        head += 1
    print('Complete three-node diameter<=80 head:', head)
    for e, t in itertools.product(range(17), [1, 3, 5, 9]):
        nodes = tuple(2**e*x for x in (0, 2, t))
        gamma = {1: 1, 3: 3, 5: 2, 9: 0}[t]
        expected = max(7*e+4, 8*e+gamma)
        check(precision(nodes, (3, 3, 3), 2) == expected,
              ('four-branches', e, t, expected))
    rng = random.Random(444039)
    for trial in range(1000):
        nodes = tuple(rng.sample(range(-2000, 2001), 3))
        shift = rng.randrange(-1000, 1001)
        unit = 2*rng.randrange(-20, 21)+1
        lifted = tuple(unit*x+shift for x in nodes)
        prediction = dyadic_three_prediction(nodes)
        check(precision(nodes, (3, 3, 3), 2) == prediction,
              ('signed-probe', trial, nodes, prediction))
        check(dyadic_three_prediction(lifted) == prediction
              == precision(lifted, (3, 3, 3), 2),
              ('unit-affine', trial, lifted, prediction))
    print('All four unit branches at17 scales;1000 signed/unit-affine pairs PASS')

    matrices = 0
    for e in range(6):
        results = []
        for base in [(0, 1, 2), (0, 1, 3)]:
            nodes = tuple(2**e*x for x in base)
            results.append(audit_matrix(nodes, (3, 3, 3)))
            matrices += 1
        check(results[0][-1] == max(7*e+4, 8*e+1)
              and results[1][-1] == max(7*e+4, 8*e+3),
              ('uniform-hostile', e, results))
        print('Uniform hostile scale', e, 'Smith', results)
        mixed = []
        for base in [(0, 2, 1), (1, 3, 0)]:
            nodes = tuple(2**e*x for x in base)
            mixed.append(audit_matrix(nodes, (2, 2, 1)))
            matrices += 1
        check(mixed[0][-1] == max(3*e+2, 4*e+1)
              and mixed[1][-1] == max(3*e+2, 4*e),
              ('mixed-hostile', e, mixed))
        print('Mixed hostile scale', e, 'Smith', mixed)

    for n in range(1, 6):
        for trial in range(5):
            nodes = tuple(rng.sample(range(-12, 13), n))
            multiplicities = tuple(rng.randrange(1, 5) for _ in nodes)
            audit_matrix(nodes, multiplicities)
            matrices += 1
    for k, e, p in itertools.product(range(1, 8), range(4), [2, 3, 5, 7]):
        gap = p**e
        predicted = max((k+l)*e-vp(math.comb(k+l-1, l), p)
                        for l in range(k))
        actual = precision((0, gap), (k, k), p)
        check(actual == predicted, ('two-node', k, e, p, actual))
    print('Independent literal inverse plus Smith matrices:', matrices)
    print('Two-node positive controls:', 7*4*4)
    print('GATES', GATES, 'SEMANTIC', DIGEST.hexdigest())
    print('PASS exact inverse and dyadic boundary controls; general theorem is analytic.')


if __name__ == '__main__':
    main()
