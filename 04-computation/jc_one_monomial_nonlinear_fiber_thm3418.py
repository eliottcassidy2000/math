#!/usr/bin/env python3
"""Exact standard-library audit for THM-3418.

For

    P=f(x)+g(x) z^d,       d>=2,

the script checks the gradient gate, derives the coefficient recurrence in
``Jac(P,Q)=kappa``, and verifies the nonterminating residue-one tower by an
independent bivariate-polynomial Jacobian calculation.  It also checks the
complete constant-``g`` normal form and its displayed inverse.  Arithmetic is
over ``Fraction``; no floating point and no Python ``assert`` statement is
used.
"""

from fractions import Fraction
from hashlib import sha256
import json


EXPECTED_SEMANTIC_SHA256 = "90d08019ad87638fecc375d76258e282bc6a4de17f53e01ca87802c8f05c5c92"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def trim(poly):
    values = [Fraction(value) for value in poly]
    while values and not values[-1]:
        values.pop()
    return tuple(values)


ZERO = ()
ONE = (Fraction(1),)
X = (Fraction(0), Fraction(1))


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ])


def poly_scale(poly, scalar):
    scalar = Fraction(scalar)
    return trim([scalar * value for value in poly])


def poly_sub(left, right):
    return poly_add(left, poly_scale(right, -1))


def poly_mul(left, right):
    if not left or not right:
        return ZERO
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return trim(answer)


def poly_pow(poly, exponent):
    require(exponent >= 0, ("negative polynomial exponent", exponent))
    answer = ONE
    base = trim(poly)
    current = exponent
    while current:
        if current & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        current //= 2
    return answer


def poly_deriv(poly):
    return trim([index * poly[index] for index in range(1, len(poly))])


def poly_divmod(dividend, divisor):
    divisor = trim(divisor)
    require(divisor, "division by zero polynomial")
    remainder = list(trim(dividend))
    quotient = [Fraction(0)] * max(1, len(remainder) - len(divisor) + 1)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = list(trim(remainder))
    return trim(quotient), trim(remainder)


def poly_mod(poly, modulus):
    return poly_divmod(poly, modulus)[1]


def poly_degree(poly):
    return len(trim(poly)) - 1


def poly_lc(poly):
    poly = trim(poly)
    require(poly, "zero polynomial has no leading coefficient")
    return poly[-1]


def constant(value):
    return trim([Fraction(value)])


def z_clean(poly):
    return {int(exponent): trim(coefficient)
            for exponent, coefficient in poly.items() if trim(coefficient)}


def z_add(left, right):
    answer = dict(z_clean(left))
    for exponent, coefficient in z_clean(right).items():
        answer[exponent] = poly_add(answer.get(exponent, ZERO), coefficient)
    return z_clean(answer)


def z_scale(poly, scalar):
    return z_clean({exponent: poly_scale(coefficient, scalar)
                    for exponent, coefficient in poly.items()})


def z_sub(left, right):
    return z_add(left, z_scale(right, -1))


def z_mul(left, right):
    answer = {}
    for left_exponent, left_coefficient in z_clean(left).items():
        for right_exponent, right_coefficient in z_clean(right).items():
            exponent = left_exponent + right_exponent
            product = poly_mul(left_coefficient, right_coefficient)
            answer[exponent] = poly_add(answer.get(exponent, ZERO), product)
    return z_clean(answer)


def z_pow(poly, exponent):
    require(exponent >= 0, ("negative bivariate exponent", exponent))
    answer = {0: ONE}
    base = z_clean(poly)
    current = exponent
    while current:
        if current & 1:
            answer = z_mul(answer, base)
        base = z_mul(base, base)
        current //= 2
    return answer


def z_dx(poly):
    return z_clean({exponent: poly_deriv(coefficient)
                    for exponent, coefficient in poly.items()})


def z_dz(poly):
    return z_clean({exponent - 1: poly_scale(coefficient, exponent)
                    for exponent, coefficient in poly.items() if exponent})


def jacobian(first, second):
    return z_sub(z_mul(z_dx(first), z_dz(second)),
                 z_mul(z_dz(first), z_dx(second)))


def compose_univariate(poly, argument):
    answer = {}
    for coefficient in reversed(trim(poly)):
        answer = z_add(z_mul(answer, argument), {0: constant(coefficient)})
    return z_clean(answer)


def gradient_gate(f, g):
    """Algebraic form of the geometric gradient-unimodularity boundary."""
    derivative_f = poly_deriv(f)
    if poly_degree(derivative_f) != 0:
        return False
    if not g:
        return True
    derivative_g = poly_deriv(g)
    return not poly_mod(poly_mul(derivative_g, derivative_g), g)


def residue_one_tower(slope, d, g, kappa, steps):
    slope = Fraction(slope)
    kappa = Fraction(kappa)
    require(slope and d >= 2 and poly_degree(g) >= 1,
            ("bad tower inputs", slope, d, g))
    derivative_g = poly_deriv(g)
    tower = [constant(kappa / slope)]
    for index in range(steps):
        n = 1 + index * d
        current = tower[-1]
        numerator = poly_sub(poly_scale(poly_mul(g, poly_deriv(current)), d),
                              poly_scale(poly_mul(derivative_g, current), n))
        tower.append(poly_scale(numerator, Fraction(1, slope * (n + d))))
    return tower


def response_operator(slope, d, g, n, q):
    """The sector transition L_n from the cokernel presentation."""
    slope = Fraction(slope)
    require(slope and d >= 2 and n >= 1, ("bad response inputs", slope, d, n))
    numerator = poly_sub(poly_scale(poly_mul(g, poly_deriv(q)), d),
                          poly_scale(poly_mul(poly_deriv(g), q), n))
    return poly_scale(numerator, Fraction(1, slope * n))


def expected_tower_lc(slope, d, g, kappa, index):
    slope = Fraction(slope)
    value = Fraction(kappa, 1) / slope
    leading_g = poly_lc(g)
    degree_g = poly_degree(g)
    for step in range(index):
        value *= (-leading_g * Fraction(d * step + degree_g,
                                        slope * (1 + (step + 1) * d)))
    return value


def check_tower_case(slope, d, g, kappa, depth, label):
    tower = residue_one_tower(slope, d, g, kappa, depth + 1)
    degree_g = poly_degree(g)
    for index, coefficient in enumerate(tower):
        expected_degree = index * (degree_g - 1)
        require(poly_degree(coefficient) == expected_degree,
                (label, "degree", index, poly_degree(coefficient), expected_degree))
        expected_leading = expected_tower_lc(slope, d, g, kappa, index)
        require(poly_lc(coefficient) == expected_leading,
                (label, "leading coefficient", index,
                 poly_lc(coefficient), expected_leading))

    # Independently form P and the truncated forced residue-one part of Q.
    p = {0: trim([2, slope]), d: g}
    q = {1 + index * d: tower[index] for index in range(depth + 1)}
    actual = jacobian(p, q)
    next_n = 1 + (depth + 1) * d
    expected = {
        0: constant(kappa),
        (depth + 1) * d: poly_scale(tower[depth + 1], -Fraction(slope) * next_n),
    }
    require(actual == z_clean(expected),
            (label, "telescoping defect", actual, z_clean(expected)))
    return len(tower)


def check_constant_normal_form(slope, intercept, g_constant, d, kappa, h, label):
    slope = Fraction(slope)
    intercept = Fraction(intercept)
    g_constant = Fraction(g_constant)
    kappa = Fraction(kappa)
    require(slope and kappa and d >= 2, (label, "bad normal-form inputs"))

    x_poly = {0: X}
    z_poly = {1: ONE}
    p = {
        0: trim([intercept, slope]),
        d: constant(g_constant),
    }
    h_of_p = compose_univariate(h, p)
    q = z_add(z_scale(z_poly, kappa / slope), h_of_p)
    require(jacobian(p, q) == {0: constant(kappa)},
            (label, "Jacobian", jacobian(p, q)))

    z_back = z_scale(z_sub(q, h_of_p), slope / kappa)
    require(z_back == z_poly, (label, "z inverse", z_back))
    x_back = z_scale(
        z_sub(z_sub(p, {0: constant(intercept)}),
              z_scale(z_pow(z_back, d), g_constant)),
        Fraction(1, slope),
    )
    require(x_back == x_poly, (label, "x inverse", x_back))


def check_response_sectors(polynomial_bank):
    checks = 0
    test_q = trim([2, -3, 1, 2])
    test_h = trim([-1, 2, 0, 1])
    for d in range(2, 9):
        for bank_index, g in enumerate(polynomial_bank):
            slope = (Fraction(1), Fraction(-2), Fraction(3, 2))[bank_index % 3]
            p = {0: trim([2, slope]), d: g}
            derivative_g = poly_deriv(g)
            for s in range(1, d):
                for level in range(3):
                    n = s + level * d
                    actual = jacobian(p, {n: test_q})
                    tail = poly_sub(poly_scale(poly_mul(derivative_g, test_q), n),
                                    poly_scale(poly_mul(g, poly_deriv(test_q)), d))
                    expected = {
                        n - 1: poly_scale(test_q, slope * n),
                        n + d - 1: tail,
                    }
                    require(actual == z_clean(expected),
                            ("response sector differential", d, bank_index,
                             s, level, actual, z_clean(expected)))
                    transition = response_operator(slope, d, g, n, test_q)
                    require(tail == poly_scale(transition, -slope * n),
                            ("response sector relation", d, bank_index,
                             s, level, tail, transition))
                    checks += 1

            for level in range(3):
                exponent = level + 1
                n = exponent * d
                source = poly_mul(poly_pow(g, exponent), test_h)
                actual = response_operator(slope, d, g, n, source)
                expected = poly_scale(
                    poly_mul(poly_pow(g, exponent + 1), poly_deriv(test_h)),
                    Fraction(1, slope * exponent),
                )
                require(actual == expected,
                        ("wrap quotient", d, bank_index, level, actual, expected))
                checks += 1
    return checks


def main():
    x2 = poly_pow(X, 2)
    nonsplit_repeated = poly_pow(poly_add(x2, ONE), 2)
    gradient_controls = [
        ("f=x^2, g=x^2", x2, x2, False),
        ("f=x, g=x", X, X, False),
        ("f=x, g=x^2", X, x2, True),
        ("f=x, g=(x^2+1)^2", X, nonsplit_repeated, True),
        ("f=x, g=3", X, constant(3), True),
        ("f=x, g=0", X, ZERO, True),
    ]
    gradient_rows = []
    for label, f, g, expected in gradient_controls:
        actual = gradient_gate(f, g)
        require(actual == expected, (label, "gradient gate", actual, expected))
        gradient_rows.append((label, actual))

    polynomial_bank = [
        X,
        x2,
        trim([2, -1, 0, 3]),
        trim([-3, 0, 2, 0, -2]),
        nonsplit_repeated,
    ]
    recurrence_cases = 0
    coefficient_checks = 0
    for d in range(2, 9):
        for bank_index, g in enumerate(polynomial_bank):
            slope = (Fraction(1), Fraction(-2), Fraction(3, 2))[bank_index % 3]
            kappa = (Fraction(1), Fraction(5, 3))[bank_index % 2]
            depth = 5
            coefficient_checks += check_tower_case(
                slope, d, g, kappa, depth,
                ("tower", d, bank_index, slope, kappa),
            )
            recurrence_cases += 1

    response_sector_checks = check_response_sectors(polynomial_bank)

    normal_form_cases = 0
    h_bank = [trim([2, -3, 1]), trim([-1, 0, 2, 1])]
    for d in range(2, 9):
        for index, g_constant in enumerate((0, 3, Fraction(-5, 2))):
            check_constant_normal_form(
                slope=(1, -2, Fraction(3, 2))[index],
                intercept=(0, 4, -3)[index],
                g_constant=g_constant,
                d=d,
                kappa=(2, Fraction(-7, 3), 5)[index],
                h=h_bank[index % 2],
                label=("normal form", d, index),
            )
            normal_form_cases += 1

    semantic_packet = {
        "theorem": "THM-3418",
        "characteristic": 0,
        "fiber_degree_minimum": 2,
        "gradient_gate": "f' is a nonzero constant and g divides (g')^2",
        "gradient_controls": gradient_rows,
        "recurrence": "a(n+d)q[n+d]=d*g*q[n]'-n*g'*q[n]",
        "residue_one_initial": "q[1]=kappa/a",
        "recurrence_cases": recurrence_cases,
        "coefficient_checks": coefficient_checks,
        "response_sector_checks": response_sector_checks,
        "constant_normal_form_cases": normal_form_cases,
        "classification": "g constant and Q=(kappa/a)z+H(P)",
        "scope": "one top fiber monomial only; JC(2) remains open",
    }
    semantic = sha256(json.dumps(
        semantic_packet, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("THM-3418 ONE-MONOMIAL NONLINEAR-FIBER KELLER AUDIT")
    print("arithmetic: exact Fraction polynomial algebra; assertion-independent")
    print("gradient controls:")
    for label, value in gradient_rows:
        print(f"  {label}: {'PASS' if value else 'FAIL'}")
    print(f"nonconstant recurrence cases: {recurrence_cases}")
    print(f"nonzero degree/leading/telescoping coefficients checked: {coefficient_checks}")
    print("sharp hostile: P=x+x^2*z^2 has unimodular gradient and an infinite residue-one tail")
    print("nonsplit hostile: g=(x^2+1)^2 has the same nonterminating tail over Q")
    print(f"Hamiltonian sector transition/wrap checks: {response_sector_checks}")
    print(f"constant-g tame normal forms and inverses checked: {normal_form_cases}")
    print("classification: f=a*x+b, g=c, Q=(kappa/a)*z+H(P)")
    print("scope: d>=2 one-monomial fiber stratum only; JC(2) remains open")
    print(f"semantic_sha256={semantic}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
