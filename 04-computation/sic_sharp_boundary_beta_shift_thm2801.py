#!/usr/bin/env python3
"""Exact controls for THM-2801.

The script checks the two displayed Special Image Conjecture witnesses by
two different exact routes:

1. literal sparse-polynomial application of E_n; and
2. the coefficient/Beta identities used in the all-m proofs.

All truth-bearing checks use ``require`` so ``python -O`` runs the same
audit.
"""

import ast
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def odd_double_factorial(index):
    require(index >= 1 and index % 2 == 1, "illegal odd double factorial")
    value = 1
    for term in range(1, index + 1, 2):
        value *= term
    return value


def image_operator(poly, coefficient_variables, position_variables):
    """Apply E_n(P(w)Q(z))=P(partial_z)Q(z) exactly."""
    variables = tuple(coefficient_variables) + tuple(position_variables)
    count = len(coefficient_variables)
    result = sp.Integer(0)
    for powers, coefficient in sp.Poly(sp.expand(poly), *variables).terms():
        alpha = powers[:count]
        beta = powers[count:]
        if any(a > b for a, b in zip(alpha, beta)):
            continue
        term = coefficient
        for a, b, variable in zip(alpha, beta, position_variables):
            term *= factorial(b) // factorial(b - a)
            term *= variable ** (b - a)
        result += term
    return sp.expand(result)


def two_pair_literal_audit():
    xi1, xi2, z1, z2 = sp.symbols("xi1 xi2 z1 z2")
    R = xi1 * z1 + xi2 * z2
    Z = xi1 * z2
    W = 2 * xi2 * z1
    T = xi1 * z1 - xi2 * z2
    F = sp.expand(
        (R + Z) * (R**2 * W - sp.Rational(1, 2) * (2 * R + Z) * T**2)
    )

    require(sp.expand(R**2 - T**2 - 2 * Z * W) == 0,
            "rank-one quadratic relation")
    require(len(F.as_ordered_terms()) == 16, "unexpected two-pair support")

    power = sp.Integer(1)
    checked = []
    for m in range(1, 9):
        power = sp.expand(power * F)
        pure = image_operator(power, (xi1, xi2), (z1, z2))
        marked = image_operator(Z * power, (xi1, xi2), (z1, z2))
        expected = (
            factorial(4 * m + 2)
            * factorial(m)
            // odd_double_factorial(2 * m + 1)
        )
        coordinate_marked = image_operator(
            z2 * power, (xi1, xi2), (z1, z2)
        )
        require(pure == 0, f"two-pair pure moment m={m}")
        require(marked == expected, f"two-pair marked moment m={m}")
        require(coordinate_marked == expected * z1,
                f"coordinate-only marked output m={m}")
        for s in range(2, min(4, m + 1)):
            coordinate_power = image_operator(
                z2**s * power, (xi1, xi2), (z1, z2)
            )
            scalar_power = image_operator(
                Z**s * power, (xi1, xi2), (z1, z2)
            )
            expected_scalar = (
                factorial(4 * m + s + 1)
                * factorial(m)
                * comb(m - 1, s - 1)
                // odd_double_factorial(2 * m + 1)
            )
            require(scalar_power == expected_scalar,
                    f"bilinear observer ladder m={m},s={s}")
            require(
                sp.Poly(coordinate_power, z1, z2).coeff_monomial(z1**s)
                == sp.Rational(expected_scalar, factorial(s)),
                f"coordinate observer ladder m={m},s={s}",
            )
        checked.append(expected)

    # A hostile control for the literal E_2 implementation.
    radial_control = image_operator(R**4, (xi1, xi2), (z1, z2))
    require(radial_control == factorial(5), "radial positive control")
    return F, checked


def two_pair_beta_audit():
    u, x, y = sp.symbols("u x y")
    checked = 0
    ladder_checked = 0
    for m in range(1, 31):
        antiderivative = sp.expand(
            sum(
                (-1) ** j * comb(m, j) * x ** (2 * j + 1) / sp.Integer(2 * j + 1)
                for j in range(m + 1)
            )
        )
        require(sp.expand(sp.diff(antiderivative, x) - (1 - x**2) ** m) == 0,
                f"Beta primitive m={m}")

        shifted = sp.Poly(
            sp.expand((1 + u) ** (m - 1) * antiderivative.subs(x, 1 + u)),
            u,
        )
        pure_coefficient = shifted.coeff_monomial(u**m)
        marked_coefficient = shifted.coeff_monomial(u ** (m - 1))
        beta_value = sp.simplify(antiderivative.subs(x, 1))
        expected_beta = sp.Rational(
            2**m * factorial(m), odd_double_factorial(2 * m + 1)
        )

        require(pure_coefficient == 0, f"missing Taylor coefficient m={m}")
        require(marked_coefficient == beta_value,
                f"one-step Taylor shift m={m}")
        require(beta_value == expected_beta, f"Beta evaluation m={m}")
        require(
            sp.Rational(1, 2**m) * marked_coefficient
            == sp.Rational(factorial(m), odd_double_factorial(2 * m + 1)),
            f"angular marked moment m={m}",
        )
        for s in range(1, min(7, m + 3)):
            coefficient_index = m - s
            shifted_coefficient = (
                shifted.coeff_monomial(u**coefficient_index)
                if coefficient_index >= 0
                else 0
            )
            expected_angular = (
                sp.Rational(factorial(m), odd_double_factorial(2 * m + 1))
                * comb(m - 1, s - 1)
                if m >= s
                else 0
            )
            require(
                sp.Rational(1, 2**m) * shifted_coefficient == expected_angular,
                f"Z-power observer ladder m={m},s={s}",
            )
            ladder_checked += 1
        checked += 1
    return checked, ladder_checked


def three_pair_polynomial(r, tau, w, v, t, z, y):
    return sp.expand(
        tau**r * (t + z) * (w * t**r - v * y * (t + y) ** (r - 1))
    )


def expected_three_pair_marked(r, m, t, z, y):
    result = factorial(r * m + 1) * factorial(m) * t
    if m == 1:
        result += factorial(r) * z
    if m == 2:
        result -= 4 * (r - 1) * factorial(2 * r) * y
    return result


def three_pair_literal_audit():
    tau, w, v, t, z, y = sp.symbols("tau w v t z y")
    f = three_pair_polynomial(3, tau, w, v, t, z, y)
    require(len(f.as_ordered_terms()) == 8, "unexpected three-pair support")

    power = sp.Integer(1)
    checked = []
    for m in range(1, 9):
        power = sp.expand(power * f)
        pure = image_operator(power, (tau, w, v), (t, z, y))
        marked = image_operator(z * power, (tau, w, v), (t, z, y))
        expected_t = factorial(3 * m + 1) * factorial(m)
        require(pure == 0, f"three-pair pure moment m={m}")
        require(marked == expected_three_pair_marked(3, m, t, z, y),
                f"three-pair full marked output m={m}")
        checked.append(expected_t)

    family_checks = 0
    for r in range(1, 5):
        family = three_pair_polynomial(r, tau, w, v, t, z, y)
        family_power = sp.Integer(1)
        for m in range(1, 5):
            family_power = sp.expand(family_power * family)
            pure = image_operator(family_power, (tau, w, v), (t, z, y))
            marked = image_operator(z * family_power, (tau, w, v), (t, z, y))
            require(pure == 0, f"general three-pair pure r={r},m={m}")
            require(marked == expected_three_pair_marked(r, m, t, z, y),
                    f"general three-pair marked r={r},m={m}")
            family_checks += 1

    f_one = three_pair_polynomial(1, tau, w, v, t, z, y)
    require(len(f_one.as_ordered_terms()) == 4, "r=1 is not four-term")
    require(
        image_operator(z * f_one, (tau, w, v), (t, z, y)) == 2 * t + z,
        "four-term witness first marked output",
    )
    return f, checked, family_checks


def three_pair_binomial_audit():
    checked = 0
    for m in range(1, 101):
        pure_sum = sum((-1) ** (m - k) * comb(m, k) for k in range(m + 1))
        marked_t_sum = sum(
            (-1) ** (m - k) * comb(m, k - 1)
            for k in range(1, m + 1)
        )
        marked_z_sum = sum(
            (-1) ** (m - k) * (k + 1) * comb(m, k)
            for k in range(m + 1)
        )
        marked_y_base = sum(
            (-1) ** (m - k) * k * (m - k) * comb(m, k)
            for k in range(1, m)
        )
        require(pure_sum == 0, f"binomial cancellation m={m}")
        require(marked_t_sum == 1, f"shifted binomial sum m={m}")
        require(marked_z_sum == (1 if m == 1 else 0),
                f"linear finite difference m={m}")
        require(marked_y_base == (-2 if m == 2 else 0),
                f"quadratic finite difference m={m}")
        for r in range(1, 9):
            t_coefficient = factorial(r * m + 1) * factorial(m) * marked_t_sum
            z_coefficient = factorial(r * m) * factorial(m) * marked_z_sum
            y_coefficient = (
                factorial(r * m)
                * factorial(m)
                * (r - 1)
                * marked_y_base
            )
            require(
                t_coefficient == factorial(r * m + 1) * factorial(m),
                f"general t coefficient r={r},m={m}",
            )
            require(z_coefficient == (factorial(r) if m == 1 else 0),
                    f"general z coefficient r={r},m={m}")
            require(
                y_coefficient
                == (-4 * (r - 1) * factorial(2 * r) if m == 2 else 0),
                f"general y coefficient r={r},m={m}",
            )
        checked += 1
    return checked


def main():
    script = Path(__file__)
    require(not has_asserts(script), "truth-bearing Python assert found")

    _, two_values = two_pair_literal_audit()
    beta_count, ladder_count = two_pair_beta_audit()
    _, three_values, family_count = three_pair_literal_audit()
    binomial_count = three_pair_binomial_audit()

    print("THM-2801 SHARP SPECIAL IMAGE BOUNDARY EXACT AUDIT")
    print("two-pair literal powers checked: 8")
    print(f"two-pair first marked value: {two_values[0]}")
    print(f"two-pair eighth marked value: {two_values[-1]}")
    print(f"angular Taylor/Beta identities checked: {beta_count}")
    print(f"Z-power observer-ladder identities checked: {ladder_count}")
    print("three-pair user-witness literal powers checked: 8")
    print(f"three-pair first t coefficient: {three_values[0]}")
    print(f"three-pair eighth t coefficient: {three_values[-1]}")
    print(f"general f_r literal pairs checked: {family_count}")
    print("r=1 four-term bidegree-(2,2) witness: PASS")
    print(f"three-pair binomial identities checked: {binomial_count}")
    print("radial Gamma(2), low-order side terms, and hostile control: PASS")


if __name__ == "__main__":
    main()
