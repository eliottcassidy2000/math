#!/usr/bin/env python3
"""Exact audit for THM-3362 (ordered-simplex odd-profile-pair FC(3) rays)."""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import sympy as sp


x, y, z, t = sp.symbols("x y z t")
T = x + y + z
U = x - y - z
V = x + y - z


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factorial_moment(poly: sp.Expr) -> sp.Expr:
    total = sp.Integer(0)
    for (i, j, k), coefficient in sp.Poly(sp.expand(poly), x, y, z).terms():
        total += coefficient * factorial(i) * factorial(j) * factorial(k)
    return sp.factor(total)


def homogenize(coefficients: dict[int, int], variable: sp.Expr) -> sp.Expr:
    degree = max(coefficients)
    require(degree % 2 == 1, "odd profile must have odd degree")
    require(all(power % 2 == 1 for power in coefficients), "profile has an even term")
    return sp.expand(
        sum(sp.Integer(coefficient) * variable**power * T ** (degree - power)
            for power, coefficient in coefficients.items())
    )


def profile_polynomial(coefficients: dict[int, int]) -> sp.Expr:
    return sum(sp.Integer(coefficient) * t**power for power, coefficient in coefficients.items())


def reduce_mod(poly: sp.Expr, modulus: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    return sp.factor(
        sp.rem(sp.Poly(sp.expand(poly), variable, domain=sp.QQ),
               sp.Poly(modulus, variable, domain=sp.QQ)).as_expr()
    )


Pair = tuple[Fraction, Fraction]
PairPolynomial = dict[tuple[int, int, int], Pair]


def pair_add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def pair_multiply(left: Pair, right: Pair) -> Pair:
    # a^2=-(42a+35)/15=-(14/5)a-7/3.
    cross = left[1] * right[1]
    return (
        left[0] * right[0] - Fraction(7, 3) * cross,
        left[0] * right[1] + left[1] * right[0] - Fraction(14, 5) * cross,
    )


def pair_polynomial_multiply(left: PairPolynomial, right: PairPolynomial) -> PairPolynomial:
    out: PairPolynomial = {}
    for (i, j, k), coefficient_left in left.items():
        for (r, s, u), coefficient_right in right.items():
            exponent = i + r, j + s, k + u
            out[exponent] = pair_add(out.get(exponent, (Fraction(0), Fraction(0))),
                                     pair_multiply(coefficient_left, coefficient_right))
    return {exponent: coefficient for exponent, coefficient in out.items() if coefficient != (0, 0)}


def pair_factorial_moment(poly: PairPolynomial) -> Pair:
    total = Fraction(0), Fraction(0)
    for (i, j, k), coefficient in poly.items():
        scale = factorial(i) * factorial(j) * factorial(k)
        total = pair_add(total, (scale * coefficient[0], scale * coefficient[1]))
    return total


def audit_tensor_formula() -> None:
    profile_pairs = [
        ("t|t", {1: 1}, {1: 1}),
        ("t+t^3|t+t^3", {1: 1, 3: 1}, {1: 1, 3: 1}),
        ("t-2t^3|t-2t^3", {1: 1, 3: -2}, {1: 1, 3: -2}),
        ("t|t^3", {1: 1}, {3: 1}),
    ]
    cells = 0
    heads: list[tuple[str, int, int, sp.Expr]] = []
    for name, h_coefficients, k_coefficients in profile_pairs:
        h = profile_polynomial(h_coefficients)
        k = profile_polynomial(k_coefficients)
        h_degree = max(h_coefficients)
        k_degree = max(k_coefficients)
        hu = homogenize(h_coefficients, U)
        kv = homogenize(k_coefficients, V)
        for radial_exponent in range(3):
            generator = sp.expand(T**radial_exponent * hu * kv)
            D = h_degree + k_degree + radial_exponent
            require(sp.Poly(generator, x, y, z).total_degree() == D, "degree mismatch")
            require(generator != 0, "generator collapsed")
            for r in range(5):
                direct = factorial_moment(generator**r)
                h_interval = sp.integrate(h**r, (t, -1, 1))
                k_interval = sp.integrate(k**r, (t, -1, 1))
                closed = sp.factor(
                    sp.Integer(factorial(D * r + 2)) * h_interval * k_interval / 8
                )
                require(sp.simplify(direct - closed) == 0, f"moment mismatch {name},e={radial_exponent},r={r}")
                if r in (2, 4) and radial_exponent == 0:
                    heads.append((name, D, r, direct))
                cells += 1
    unequal_prefix = []
    unequal_generator = sp.expand(U * V**3)
    for r in range(5):
        unequal_prefix.append(factorial_moment(unequal_generator**r))
    require(
        unequal_prefix == [1, 0, 86400, 0, 49249028505600],
        "unequal-profile prefix changed",
    )
    print("TENSOR_PAIR_FORMULA", cells, "direct cells PASS")
    print("TENSOR_PAIR_HEAD", heads)
    print("UNEQUAL_ODD_PREFIX h=t,k=t^3", unequal_prefix)


def audit_mixed_parity_hostile() -> None:
    u, v = sp.symbols("u v")
    ordered_angular = sp.integrate(sp.integrate(v, (v, u, 1)), (u, -1, 1))
    simplex_angular = ordered_angular / 4
    direct = factorial_moment(V)
    radial_reconstruction = factorial(3) * simplex_angular
    false_half_square = (
        sp.Integer(factorial(3))
        * sp.integrate(sp.Integer(1), (t, -1, 1))
        * sp.integrate(t, (t, -1, 1))
        / 8
    )
    require(ordered_angular == sp.Rational(2, 3), "mixed-parity triangle changed")
    require(simplex_angular == sp.Rational(1, 6), "mixed-parity Jacobian changed")
    require(direct == radial_reconstruction == 1, "mixed-parity direct moment changed")
    require(false_half_square == 0, "mixed-parity false product is nonzero")
    print(
        "MIXED_PARITY_HOSTILE h=1,k=t L(V)=1 false_half_square=0",
        "ordered_triangle=2/3 PASS",
    )


def audit_separate_pencil_hostile() -> None:
    left = sp.Matrix([[0, 1], [-1, 2]])
    middle = sp.Matrix([[0, 1], [1, 2]])
    right = sp.Matrix([[1, 0], [2, 1]])
    transfer_left = sp.Matrix([[1, -2, 2], [2, -1, 2], [2, -2, 3]])
    transfer_middle = sp.Matrix([[1, 2, 2], [2, 1, 2], [2, 2, 3]])
    transfer_right = sp.Matrix([[-1, 2, 2], [-2, 1, 2], [-2, 2, 3]])
    q_at_one = sp.det(left + middle + right)
    symmetric_square_det = q_at_one**3
    p_at_one = sp.det(transfer_left + transfer_middle + transfer_right)
    require(q_at_one == 1, "spinor pencil value changed")
    require(symmetric_square_det == 1, "symmetric-square determinant changed")
    require(p_at_one == -3, "triple-transfer determinant changed")
    require(symmetric_square_det != p_at_one, "separate pencils unexpectedly agree")
    print("SEPARATE_PENCIL_HOSTILE q=1 sym2_det=1 transfer_det=-3 PASS")


def audit_elimination() -> None:
    A, B, C, p, q, w = sp.symbols("A B C p q w")
    m1 = A + p * C
    m2 = A**2 + 2 * p * A * C + p * B**2 + q * C**2
    m3 = (
        A**3 + 3 * p * A**2 * C + 3 * p * A * B**2
        + 3 * q * A * C**2 + 3 * q * B**2 * C + w * C**3
    )
    m2_reduced = sp.factor(m2.subs(A, -p * C))
    m3_reduced = sp.factor(m3.subs(A, -p * C))
    require(m1.subs(A, -p * C) == 0, "first-moment elimination failed")
    require(sp.simplify(m2_reduced - (p * B**2 + (q - p**2) * C**2)) == 0, "second-moment formula failed")
    require(
        sp.simplify(m3_reduced - (3 * (q - p**2) * B**2 * C + (w - 3 * p * q + 2 * p**3) * C**3)) == 0,
        "third-moment formula failed",
    )
    b_squared = -(q - p**2) * C**2 / p
    eliminated = sp.factor(m3_reduced.subs(B**2, b_squared))
    J = p * w + 3 * p**2 * q - p**4 - 3 * q**2
    require(sp.simplify(eliminated - C**3 * J / p) == 0, "final eliminant failed")
    print("SYMBOLIC_ELIMINATION J=", J, "PASS")


def audit_factorial_gap() -> None:
    for D in range(2, 101):
        ratio = Fraction(factorial(2 * D + 2) * factorial(6 * D + 2), factorial(4 * D + 2) ** 2)
        product = Fraction(1)
        for j in range(1, 2 * D + 1):
            product *= Fraction(4 * D + 2 + j, 2 * D + 2 + j)
        require(ratio == product, f"factorial product mismatch D={D}")
        require(ratio >= Fraction(7, 5) ** (2 * D), f"factor bound failed D={D}")
        require(ratio > 3, f"factor gap failed D={D}")
    print("FACTORIAL_GAP D=2..100 PASS")
    print("UNIVERSAL_LOWER_BOUND", Fraction(7, 5) ** 4)
    print("FACTORIAL_RATIO_HEAD", [(D, Fraction(factorial(2 * D + 2) * factorial(6 * D + 2), factorial(4 * D + 2) ** 2)) for D in (2, 3, 4)])


def closed_p_moment(r: int) -> int:
    if r % 2:
        return 0
    return factorial(3 * r + 2) // (2 * (r + 1) ** 2)


def audit_determinant_specialization() -> None:
    for d in range(1, 100, 2):
        D = 3 * d
        for r in range(7):
            interval = sp.Integer(0) if r % 2 else sp.Rational(2, d * r + 1)
            tensor = sp.Integer(factorial(D * r + 2)) * interval**2 / 8
            require(tensor == closed_p_moment(d * r), f"P^d specialization failed d={d},r={r}")
    print("DETERMINANT_RAYS odd d=1..99, r=0..6 PASS")


def audit_complex_hostile() -> None:
    a = sp.symbols("a")
    modulus = 15 * a**2 + 42 * a + 35
    h = t + a * t**3
    generator_expr = sp.expand((U * T**2 + a * U**3) * (V * T**2 + a * V**3))
    generator: PairPolynomial = {}
    for exponent, coefficient in sp.Poly(generator_expr, x, y, z).terms():
        reduced = reduce_mod(coefficient, modulus, a)
        coefficient_poly = sp.Poly(reduced, a, domain=sp.QQ)
        generator[exponent] = (
            Fraction(int(coefficient_poly.nth(0).p), int(coefficient_poly.nth(0).q)),
            Fraction(int(coefficient_poly.nth(1).p), int(coefficient_poly.nth(1).q)),
        )
    power: PairPolynomial = {(0, 0, 0): (Fraction(1), Fraction(0))}
    remainders: list[Pair] = []
    for r in range(1, 5):
        power = pair_polynomial_multiply(power, generator)
        direct = pair_factorial_moment(power)
        interval = reduce_mod(sp.integrate(h**r, (t, -1, 1)), modulus, a)
        interval_poly = sp.Poly(interval, a, domain=sp.QQ)
        interval_pair = (
            Fraction(int(interval_poly.nth(0).p), int(interval_poly.nth(0).q)),
            Fraction(int(interval_poly.nth(1).p), int(interval_poly.nth(1).q)),
        )
        closed = pair_multiply(interval_pair, interval_pair)
        closed = (Fraction(factorial(6 * r + 2), 8) * closed[0],
                  Fraction(factorial(6 * r + 2), 8) * closed[1])
        require(direct == closed, f"complex formula failed r={r}")
        remainders.append(direct)
    require(remainders[:3] == [(0, 0), (0, 0), (0, 0)], "complex hostile does not kill moments 1--3")
    require(remainders[3] != (0, 0), "complex hostile also killed moment four")
    i2 = sp.integrate(h**2, (t, -1, 1))
    i4 = sp.integrate(h**4, (t, -1, 1))
    require(sp.factor(i2 - 2 * modulus / 105) == 0, "I2 formula failed")
    i4_remainder = reduce_mod(i4, modulus, a)
    expected_i4 = sp.Rational(128, 1126125) * (417 * a + 560)
    require(sp.simplify(i4_remainder - expected_i4) == 0, "I4 remainder failed")
    require(sp.gcd(sp.Poly(modulus, a), sp.Poly(417 * a + 560, a)).degree() == 0, "hostile gcd failed")
    print("COMPLEX_PROFILE moments_1_2_3_zero moment_4_nonzero PASS")
    print("COMPLEX_I4_REMAINDER", i4_remainder)
    print("COMPLEX_L4_REMAINDER_PAIR", remainders[3])


def audit_first_two_hostile_example() -> None:
    # h=k=t,e=0 gives R=UV, D=2 and especially small exact scalars.
    p = factorial(6) // 18
    q = factorial(10) // 50
    w = factorial(14) // 98
    require((p, q, w) == (40, 72576, 889574400), "small profile moments changed")
    b_squared = Fraction(-(q - p * p), p)
    J = p * w + 3 * p * p * q - p**4 - 3 * q * q
    require(b_squared == Fraction(-8872, 5), "small hostile changed")
    require(J > 0, "small third-moment obstruction failed")
    print("FIRST_TWO_HOSTILE h=k=t,e=0 B^2", b_squared, "J", J)


def main() -> None:
    audit_tensor_formula()
    audit_mixed_parity_hostile()
    audit_separate_pencil_hostile()
    audit_elimination()
    audit_factorial_gap()
    audit_determinant_specialization()
    audit_complex_hostile()
    audit_first_two_hostile_example()
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
