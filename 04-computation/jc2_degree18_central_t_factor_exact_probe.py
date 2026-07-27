#!/usr/bin/env python3
"""Exact rational reconstruction of the degree-18 central weight-30 factor.

On D=25*B^2/126 the branch discriminant of the THM-2297 spectral cubic
has a universal y^2 factor.  Its degree-ten residual has a repeated-root
resultant of weight 90.  This script reconstructs that resultant exactly on
the B=1 weighted-projective chart by a 31 x 19 rational tensor grid and
checks the full coefficient identity

    Res_y(delta/y^2, d_y(delta/y^2)) = constant * R^6 * S^3 * T.

Weighted homogeneity then lifts the identity to (B,C,W).  No CAS package is
used.
"""

from fractions import Fraction
from math import factorial, gcd


T_MONOMIALS = [
    (0, 0, 6),
    (0, 5, 3),
    (0, 10, 0),
    (1, 1, 5),
    (1, 6, 2),
    (2, 2, 4),
    (2, 7, 1),
    (3, 3, 3),
    (3, 8, 0),
    (4, 4, 2),
    (5, 0, 4),
    (5, 5, 1),
    (6, 1, 3),
    (6, 6, 0),
    (7, 2, 2),
    (8, 3, 1),
    (9, 4, 0),
    (10, 0, 2),
    (11, 1, 1),
    (12, 2, 0),
    (15, 0, 0),
]

T_COEFFICIENTS = [
    1857584148891020160556640625,
    -720893165611681749052800000,
    187775146320922282243129344,
    6184255303140290720859375000,
    -1511892036737051357606400000,
    8180469698416791719343750000,
    -1038920491424993478174720000,
    5592912523479745046100000000,
    -224453394609952422813696000,
    2217787332403775132700000000,
    12754325524369432500000000,
    550162720634152586284800000,
    31485168419060160000000000,
    73571138739313427520000000,
    25958526361650288000000000,
    7389030553435238400000000,
    206537825241216000000000,
    23073416987040000000000,
    40577760518400000000000,
    17795598912000000000000,
    256000000000000000,
]


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add(first, second):
    out = [Fraction(0)] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] += value
    for index, value in enumerate(second):
        out[index] += value
    return trim(out)


def scale(poly, scalar):
    return trim([scalar * value for value in poly])


def mul(first, second):
    out = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            out[i + j] += left * right
    return trim(out)


def power(poly, exponent):
    out = [Fraction(1)]
    base = poly[:]
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        base = mul(base, base)
        exponent //= 2
    return out


def derivative(poly):
    if len(poly) == 1:
        return [Fraction(0)]
    return trim([index * poly[index] for index in range(1, len(poly))])


def divmod_poly(numerator, denominator):
    numerator = trim(numerator[:])
    denominator = trim(denominator[:])
    if denominator == [0]:
        raise ZeroDivisionError
    if len(numerator) < len(denominator):
        return [Fraction(0)], numerator
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] -= coefficient * value
        trim(numerator)
    return trim(quotient), trim(numerator)


def resultant(first, second):
    first = trim(first[:])
    second = trim(second[:])
    m = len(first) - 1
    n = len(second) - 1
    if second == [0]:
        return Fraction(0)
    if n == 0:
        return second[0] ** m
    _, remainder = divmod_poly(first, second)
    if remainder == [0]:
        return Fraction(0)
    r = len(remainder) - 1
    sign = -1 if (m * n) & 1 else 1
    return (
        sign
        * second[-1] ** (m - r)
        * resultant(second, remainder)
    )


def interpolate_consecutive(values):
    """Exact monomial coefficients from values at 0,1,...,n."""
    differences = [Fraction(value) for value in values]
    out = [Fraction(0)]
    falling = [Fraction(1)]
    for order in range(len(values)):
        coefficient = differences[0] / factorial(order)
        out = add(out, scale(falling, coefficient))
        differences = [
            differences[index + 1] - differences[index]
            for index in range(len(differences) - 1)
        ]
        falling = mul(falling, [Fraction(-order), Fraction(1)])
    return trim(out)


def central_residual_branch(c_value, w_value):
    b_value = 1
    d_value = Fraction(25, 126)
    a = Fraction(-26040609)
    bpoly = [Fraction(49601160), 0, Fraction(1607445)]
    cpoly = [
        Fraction(-20995200) - Fraction(52907904) * d_value,
        0,
        Fraction(-2857680),
        0,
        Fraction(-138915),
    ]
    dpoly = [
        Fraction(33592320) * d_value,
        Fraction(-5878656 * w_value - 5598720 * c_value),
        Fraction(777600) + Fraction(1959552) * d_value,
        Fraction(-435456 * c_value),
        Fraction(78120),
        0,
        Fraction(1127),
    ]
    terms = [
        mul(power(bpoly, 2), power(cpoly, 2)),
        scale(power(cpoly, 3), -4 * a),
        scale(mul(power(bpoly, 3), dpoly), -4),
        scale(power(dpoly, 2), -27 * a * a),
        scale(mul(mul(bpoly, cpoly), dpoly), 18 * a),
    ]
    branch = [Fraction(0)]
    for term in terms:
        branch = add(branch, term)
    if len(branch) != 13 or branch[0] != 0 or branch[1] != 0:
        raise RuntimeError("central branch lost the universal y^2 factor")
    residual = trim(branch[2:])
    if len(residual) != 11:
        raise RuntimeError("central residual lost degree ten")
    return residual


def residual_resultant(c_value, w_value):
    residual = central_residual_branch(c_value, w_value)
    return resultant(residual, derivative(residual))


def bivariate_mul(first, second):
    out = {}
    for (j1, k1), left in first.items():
        for (j2, k2), right in second.items():
            key = (j1 + j2, k1 + k2)
            out[key] = out.get(key, 0) + left * right
    return {key: value for key, value in out.items() if value}


def bivariate_power(poly, exponent):
    out = {(0, 0): 1}
    base = dict(poly)
    while exponent:
        if exponent & 1:
            out = bivariate_mul(out, base)
        base = bivariate_mul(base, base)
        exponent //= 2
    return out


def expected_factor():
    r_factor = {(1, 0): 20, (0, 1): 21}
    s_factor = {
        (0, 0): 2888,
        (2, 0): 108864,
        (1, 1): 571536,
        (0, 2): 750141,
    }
    t_factor = {
        (j, k): coefficient
        for (_i, j, k), coefficient in zip(
            T_MONOMIALS, T_COEFFICIENTS, strict=True
        )
    }
    return bivariate_mul(
        bivariate_mul(
            bivariate_power(r_factor, 6),
            bivariate_power(s_factor, 3),
        ),
        t_factor,
    )


def reconstructed_resultant():
    # First interpolate in W at each C, then in C coefficientwise.
    w_coefficients_by_c = []
    for c_value in range(31):
        values = [
            residual_resultant(c_value, w_value)
            for w_value in range(19)
        ]
        coefficients = interpolate_consecutive(values)
        coefficients.extend(
            [Fraction(0)] * (19 - len(coefficients))
        )
        w_coefficients_by_c.append(coefficients)

    out = {}
    for w_power in range(19):
        c_coefficients = interpolate_consecutive(
            [
                w_coefficients_by_c[c_value][w_power]
                for c_value in range(31)
            ]
        )
        for c_power, coefficient in enumerate(c_coefficients):
            if coefficient:
                out[(c_power, w_power)] = coefficient
    return out


def main():
    if any(
        2 * i + 3 * j + 5 * k != 30
        for i, j, k in T_MONOMIALS
    ):
        raise RuntimeError("T support lost weight thirty")
    coefficient_gcd = 0
    for coefficient in T_COEFFICIENTS:
        coefficient_gcd = gcd(coefficient_gcd, abs(coefficient))
    if coefficient_gcd != 1:
        raise RuntimeError("T coefficients are not primitive")

    found = reconstructed_resultant()
    expected = expected_factor()
    common_key = (0, 18)
    if common_key not in found or common_key not in expected:
        raise RuntimeError("missing leading W coefficient")
    scalar = found[common_key] / expected[common_key]
    scaled_expected = {
        key: Fraction(value) * scalar for key, value in expected.items()
    }
    if found != scaled_expected:
        missing = set(found).symmetric_difference(scaled_expected)
        raise RuntimeError(f"central factor identity failed: {sorted(missing)[:5]}")

    print("JC2 DEGREE-18 CENTRAL WEIGHT-30 FACTOR EXACT AUDIT")
    print("universal_branch_factor=y^2")
    print("residual_branch_degree=10")
    print("residual_resultant_weight=90")
    print("factorization=constant*R^6*S^3*T")
    print("R=20*B*C+21*W")
    print("S=2888*B^5+108864*B^2*C^2+571536*B*C*W+750141*W^2")
    print("T_monomial_count=", len(T_MONOMIALS), sep="")
    print("T_support_and_coefficients=")
    for monomial, coefficient in zip(
        T_MONOMIALS, T_COEFFICIENTS, strict=True
    ):
        print(monomial, coefficient)
    print("resultant_scalar=", scalar, sep="")
    print("tensor_grid=31x19 exact rational points on B=1")
    print("primitive_coefficient_gcd=1")
    print("PASS")


if __name__ == "__main__":
    main()
