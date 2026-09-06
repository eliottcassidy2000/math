#!/usr/bin/env python3
"""Independent literal Smith audit of higher/mixed Hasse-jet boundaries.

No primary compiler imports. Every named case compares the entire largest
integer Smith factor against exact truncated reciprocal-product denominators.
The dyadic cases additionally test the closed valuation formulas.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import comb, lcm
from random import Random

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0
CASES = 0
DIGEST = sha256()


def need(value, payload):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(payload)


def v2(value):
    value = abs(int(value))
    need(value != 0, "nonsingular Smith entry")
    result = 0
    while value % 2 == 0:
        value //= 2
        result += 1
    return result


def reciprocal_denominator(nodes, multiplicities):
    result = 1
    for i, x in enumerate(nodes):
        length = multiplicities[i]
        series = [Q(1)] + [Q(0)] * (length - 1)
        for j, y in enumerate(nodes):
            if i == j:
                continue
            power = multiplicities[j]
            factor = [Q((-1)**ell * comb(power + ell - 1, ell),
                        (x - y)**(power + ell)) for ell in range(length)]
            series = [sum((series[r] * factor[ell - r] for r in range(ell + 1)), Q())
                      for ell in range(length)]
        for coefficient in series:
            result = lcm(result, coefficient.denominator)
    return result


def audit(nodes, multiplicities):
    global CASES
    degree = sum(multiplicities)
    need(len(set(nodes)) == len(nodes), "distinct nodes")
    matrix = Matrix([
        [comb(j, r) * x**(j - r) if j >= r else 0 for j in range(degree)]
        for x, multiplicity in zip(nodes, multiplicities)
        for r in range(multiplicity)
    ])
    diagonal = smith_normal_form(matrix, domain=ZZ)
    factors = tuple(abs(int(diagonal[i, i])) for i in range(degree))
    expected = reciprocal_denominator(nodes, multiplicities)
    need(factors[-1] == expected, ("full integer factor", nodes, multiplicities, factors, expected))
    CASES += 1
    DIGEST.update((repr((nodes, multiplicities, factors)) + "\n").encode())
    return tuple(v2(factor) for factor in factors)


def main():
    for exponent in range(6):
        result = []
        for base, constants in (((0, 1, 2), (3, 4, 1)), ((0, 1, 3), (3, 4, 3))):
            nodes = tuple(2**exponent * x for x in base)
            spectrum = audit(nodes, (3, 3, 3))
            expected = max((6 + ell) * exponent + constants[ell] for ell in range(3))
            need(spectrum[-1] == expected, ("uniform twins", exponent, base, spectrum))
            result.append(spectrum)
        print("UNIFORM3 e", exponent, "A", result[0], "B", result[1])
    for exponent in range(6):
        first = audit(tuple(2**exponent * x for x in (0, 2, 1)), (2, 2, 1))
        second = audit(tuple(2**exponent * x for x in (1, 3, 0)), (2, 2, 1))
        need(first[-1] == max(3 * exponent + 2, 4 * exponent + 1), ("mixed A", exponent, first))
        need(second[-1] == max(3 * exponent + 2, 4 * exponent), ("mixed B", exponent, second))
        print("MIXED221 e", exponent, "A", first, "B", second)
    unit_cases = 0
    for exponent in range(5):
        for unit, gamma in ((1, 1), (3, 3), (5, 2), (9, 0)):
            spectrum = audit(tuple(2**exponent * x for x in (0, 2, unit)), (3, 3, 3))
            need(spectrum[-1] == max(7 * exponent + 4, 8 * exponent + gamma),
                 ("unit residue table", exponent, unit, spectrum))
            unit_cases += 1
        for gap in range(2, 5):
            for unit in (1, 3, 5):
                spectrum = audit(tuple(2**exponent * x for x in (0, 2**gap, unit)), (3, 3, 3))
                need(spectrum[-1] == 8 * exponent + 5 * gap - 1,
                     ("deep gap", exponent, gap, unit, spectrum))
                unit_cases += 1
    need(unit_cases == 65, "complete named dyadic controls")
    for multiplicity in range(1, 6):
        need(audit((-3,), (multiplicity,))[-1] == 0, "singleton jet has no loss")
    rng = Random(4440)
    for _ in range(32):
        n = rng.randrange(1, 5)
        nodes = tuple(sorted(rng.sample(range(-10, 16), n)))
        multiplicities = tuple(rng.randrange(1, 4) for _ in nodes)
        audit(nodes, multiplicities)
    need(CASES == 126, ("literal matrix count", CASES))
    print("NAMED DYADIC UNIT/GAP CONTROLS", unit_cases)
    print("LITERAL INTEGER SMITH MATRICES", CASES, "EXPLICIT GATES", CHECKS)
    print("FULL INTEGER FACTOR/RECIPROCAL SEMANTIC SHA256", DIGEST.hexdigest())
    print("PASS uniform3 metric hostile; mixed221 earlier type boundary; arbitrary-multiplicity inverse denominators")


if __name__ == "__main__":
    main()
