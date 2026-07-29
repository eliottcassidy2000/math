#!/usr/bin/env python3
"""Exact referee for the Gamma-holotopy of the THM-2846 transverse hostile."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ALPHA, X, Y = sp.symbols("alpha X Y")


def rising(order: int) -> sp.Expr:
    return sp.prod(ALPHA + offset for offset in range(order))


def adjacent_tensor(indices: tuple[int, ...]) -> sp.Expr:
    """Return L_alpha(prod d_n) exactly."""
    denominator = sp.prod(rising(index + 1) for index in indices)
    numerator = sp.Integer(0)
    for bits in product((0, 1), repeat=len(indices)):
        term = (-1) ** (len(indices) - sum(bits))
        term *= rising(sum(index + bit for index, bit in zip(indices, bits)))
        for index, bit in zip(indices, bits):
            if bit == 0:
                term *= ALPHA + index
        numerator += term
    return sp.cancel(numerator / denominator)


TENSOR_CACHE: dict[tuple[int, ...], sp.Expr] = {}


def tensor(indices: tuple[int, ...]) -> sp.Expr:
    key = tuple(sorted(indices))
    if key not in TENSOR_CACHE:
        TENSOR_CACHE[key] = adjacent_tensor(key)
    return TENSOR_CACHE[key]


def multilinear(vectors: tuple[dict[int, sp.Expr], ...]) -> sp.Expr:
    return sp.factor(
        sum(
            sp.prod(coefficient for _, coefficient in choices)
            * tensor(tuple(index for index, _ in choices))
            for choices in product(*(tuple(vector.items()) for vector in vectors))
        )
    )


def poly_dictionary(expression: sp.Expr) -> dict[tuple[int, int, int], Fraction]:
    polynomial = sp.Poly(expression, ALPHA, X, Y)
    require(polynomial.domain == sp.ZZ, "certificate polynomial is not integral")
    return {
        powers: Fraction(int(coefficient), 1)
        for powers, coefficient in polynomial.terms()
    }


def multiply_sparse(
    left: dict[tuple[int, ...], Fraction],
    right: dict[tuple[int, ...], Fraction],
) -> dict[tuple[int, ...], Fraction]:
    result: dict[tuple[int, ...], Fraction] = {}
    for left_power, left_coefficient in left.items():
        for right_power, right_coefficient in right.items():
            power = tuple(a + b for a, b in zip(left_power, right_power))
            result[power] = (
                result.get(power, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        power: coefficient
        for power, coefficient in result.items()
        if coefficient
    }


def linear_powers(
    degree: int,
    constant: Fraction,
    first: Fraction,
    second: Fraction,
) -> tuple[dict[tuple[int, int], Fraction], ...]:
    """Powers of constant + first*A + second*U, indexed by (A,U)."""
    powers: list[dict[tuple[int, int], Fraction]] = [{(0, 0): Fraction(1)}]
    linear = {
        (0, 0): constant,
        (1, 0): first,
        (0, 1): second,
    }
    for _ in range(degree):
        powers.append(multiply_sparse(powers[-1], linear))
    return tuple(powers)


def pull_back_to_unit_cube(
    expression: sp.Expr,
) -> tuple[list[list[list[Fraction]]], tuple[int, int, int]]:
    """Pull back along the rational moving tube, as exact power coefficients.

    alpha = 9/10 + t/5,
    X = 26362/100000 + 73/1000*(alpha-1) + (2u-1)/10000,
    Y = 23420/1000000 + 7/500*(alpha-1) + (2v-1)/20000.
    """
    original = poly_dictionary(expression)
    degree_x = max(power[1] for power in original)
    degree_y = max(power[2] for power in original)

    x_powers = linear_powers(
        degree_x,
        Fraction(26362, 100000) - Fraction(73, 1000) - Fraction(1, 10000),
        Fraction(73, 1000),
        Fraction(2, 10000),
    )
    y_powers = linear_powers(
        degree_y,
        Fraction(23420, 1000000) - Fraction(7, 500) - Fraction(1, 20000),
        Fraction(7, 500),
        Fraction(2, 20000),
    )

    alpha_u_v: dict[tuple[int, int, int], Fraction] = {}
    for (power_alpha, power_x, power_y), coefficient in original.items():
        for (extra_alpha_x, power_u), coefficient_x in x_powers[power_x].items():
            for (extra_alpha_y, power_v), coefficient_y in y_powers[power_y].items():
                power = (
                    power_alpha + extra_alpha_x + extra_alpha_y,
                    power_u,
                    power_v,
                )
                alpha_u_v[power] = (
                    alpha_u_v.get(power, Fraction(0))
                    + coefficient * coefficient_x * coefficient_y
                )

    max_alpha = max(power[0] for power in alpha_u_v)
    alpha_powers: list[dict[int, Fraction]] = [{0: Fraction(1)}]
    for _ in range(max_alpha):
        previous = alpha_powers[-1]
        current: dict[int, Fraction] = {}
        for degree_t, coefficient in previous.items():
            current[degree_t] = (
                current.get(degree_t, Fraction(0))
                + coefficient * Fraction(9, 10)
            )
            current[degree_t + 1] = (
                current.get(degree_t + 1, Fraction(0))
                + coefficient * Fraction(1, 5)
            )
        alpha_powers.append(current)

    pulled: dict[tuple[int, int, int], Fraction] = {}
    for (power_alpha, power_u, power_v), coefficient in alpha_u_v.items():
        for power_t, alpha_coefficient in alpha_powers[power_alpha].items():
            power = (power_t, power_u, power_v)
            pulled[power] = (
                pulled.get(power, Fraction(0))
                + coefficient * alpha_coefficient
            )

    degrees = tuple(
        max(power[axis] for power in pulled) for axis in range(3)
    )
    degree_t, degree_u, degree_v = degrees
    array = [
        [
            [Fraction(0) for _ in range(degree_v + 1)]
            for _ in range(degree_u + 1)
        ]
        for _ in range(degree_t + 1)
    ]
    for (power_t, power_u, power_v), coefficient in pulled.items():
        array[power_t][power_u][power_v] = coefficient
    return array, degrees


def power_to_bernstein(
    power: list[list[list[Fraction]]],
    degrees: tuple[int, int, int],
) -> list[list[list[Fraction]]]:
    """Exact separable power-to-Bernstein conversion on [0,1]^3."""
    degree_t, degree_u, degree_v = degrees
    after_t = [
        [
            [Fraction(0) for _ in range(degree_v + 1)]
            for _ in range(degree_u + 1)
        ]
        for _ in range(degree_t + 1)
    ]
    for index_t in range(degree_t + 1):
        for power_u in range(degree_u + 1):
            for power_v in range(degree_v + 1):
                after_t[index_t][power_u][power_v] = sum(
                    power[power_t][power_u][power_v]
                    * Fraction(comb(index_t, power_t), comb(degree_t, power_t))
                    for power_t in range(index_t + 1)
                )

    after_u = [
        [
            [Fraction(0) for _ in range(degree_v + 1)]
            for _ in range(degree_u + 1)
        ]
        for _ in range(degree_t + 1)
    ]
    for index_t in range(degree_t + 1):
        for index_u in range(degree_u + 1):
            for power_v in range(degree_v + 1):
                after_u[index_t][index_u][power_v] = sum(
                    after_t[index_t][power_u][power_v]
                    * Fraction(comb(index_u, power_u), comb(degree_u, power_u))
                    for power_u in range(index_u + 1)
                )

    result = [
        [
            [Fraction(0) for _ in range(degree_v + 1)]
            for _ in range(degree_u + 1)
        ]
        for _ in range(degree_t + 1)
    ]
    for index_t in range(degree_t + 1):
        for index_u in range(degree_u + 1):
            for index_v in range(degree_v + 1):
                result[index_t][index_u][index_v] = sum(
                    after_u[index_t][index_u][power_v]
                    * Fraction(comb(index_v, power_v), comb(degree_v, power_v))
                    for power_v in range(index_v + 1)
                )
    return result


def all_coefficients(
    array: list[list[list[Fraction]]],
) -> tuple[Fraction, ...]:
    return tuple(
        coefficient
        for plane in array
        for row in plane
        for coefficient in row
    )


def face_coefficients(
    array: list[list[list[Fraction]]],
    degrees: tuple[int, int, int],
    axis: int,
    side: int,
) -> tuple[Fraction, ...]:
    degree_t, degree_u, degree_v = degrees
    index = (degree_t, degree_u, degree_v)[axis] if side else 0
    if axis == 0:
        return tuple(
            array[index][u][v]
            for u in range(degree_u + 1)
            for v in range(degree_v + 1)
        )
    if axis == 1:
        return tuple(
            array[t][index][v]
            for t in range(degree_t + 1)
            for v in range(degree_v + 1)
        )
    return tuple(
        array[t][u][index]
        for t in range(degree_t + 1)
        for u in range(degree_u + 1)
    )


def evaluate_power(
    array: list[list[list[Fraction]]],
    point: tuple[Fraction, Fraction, Fraction],
) -> Fraction:
    return sum(
        coefficient
        * point[0] ** power_t
        * point[1] ** power_u
        * point[2] ** power_v
        for power_t, plane in enumerate(array)
        for power_u, row in enumerate(plane)
        for power_v, coefficient in enumerate(row)
    )


def evaluate_bernstein(
    array: list[list[list[Fraction]]],
    degrees: tuple[int, int, int],
    point: tuple[Fraction, Fraction, Fraction],
) -> Fraction:
    degree_t, degree_u, degree_v = degrees
    parameter_t, parameter_u, parameter_v = point
    basis_t = tuple(
        comb(degree_t, index)
        * parameter_t**index
        * (1 - parameter_t) ** (degree_t - index)
        for index in range(degree_t + 1)
    )
    basis_u = tuple(
        comb(degree_u, index)
        * parameter_u**index
        * (1 - parameter_u) ** (degree_u - index)
        for index in range(degree_u + 1)
    )
    basis_v = tuple(
        comb(degree_v, index)
        * parameter_v**index
        * (1 - parameter_v) ** (degree_v - index)
        for index in range(degree_v + 1)
    )
    return sum(
        array[index_t][index_u][index_v]
        * basis_t[index_t]
        * basis_u[index_u]
        * basis_v[index_v]
        for index_t in range(degree_t + 1)
        for index_u in range(degree_u + 1)
        for index_v in range(degree_v + 1)
    )


def moving_tube_point(
    point: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    parameter_t, parameter_u, parameter_v = point
    alpha = Fraction(9, 10) + Fraction(1, 5) * parameter_t
    x = (
        Fraction(26362, 100000)
        + Fraction(73, 1000) * (alpha - 1)
        + (2 * parameter_u - 1) / 10000
    )
    y = (
        Fraction(23420, 1000000)
        + Fraction(7, 500) * (alpha - 1)
        + (2 * parameter_v - 1) / 20000
    )
    return alpha, x, y


def check_pullback(
    expression: sp.Expr,
    power: list[list[list[Fraction]]],
    point: tuple[Fraction, Fraction, Fraction],
) -> None:
    alpha, x, y = moving_tube_point(point)
    expected = evaluate_power(power, point)
    actual = sum(
        coefficient * alpha**power_alpha * x**power_x * y**power_y
        for (power_alpha, power_x, power_y), coefficient in poly_dictionary(
            expression
        ).items()
    )
    require(actual == expected, "sparse affine pullback mismatch")


def main() -> None:
    lower = {1: sp.Integer(1), 3: X}
    upper = {2: sp.Integer(1), 3: Y}

    # Independent normalization controls at alpha=1.
    for order in range(2, 5):
        for indices in product(range(4), repeat=order):
            expected = sp.Integer(0)
            for bits in product((0, 1), repeat=order):
                degrees = tuple(
                    index + bit for index, bit in zip(indices, bits)
                )
                term = (
                    (-1) ** (order - sum(bits))
                    * sp.factorial(sum(degrees))
                )
                for degree in degrees:
                    term /= sp.factorial(degree)
                expected += term
            require(
                sp.cancel(tensor(indices).subs(ALPHA, 1) - expected) == 0,
                "alpha=1 Pascal tensor mismatch",
            )

    g11 = multilinear((lower, lower))
    g12 = multilinear((lower, upper))
    g22 = multilinear((upper, upper))
    t111 = multilinear((lower, lower, lower))
    t112 = multilinear((lower, lower, upper))
    t122 = multilinear((lower, upper, upper))
    t222 = multilinear((upper, upper, upper))

    invariant_one = sp.factor(
        3 * t112 * g11 * g22
        - t222 * g11**2
        - 2 * t111 * g12 * g22
    )
    invariant_two = sp.factor(
        3 * t122 * g11 * g22
        - 2 * t222 * g12 * g11
        - t111 * g22**2
    )
    numerator_one, denominator_one = sp.together(
        invariant_one
    ).as_numer_denom()
    numerator_two, denominator_two = sp.together(
        invariant_two
    ).as_numer_denom()
    common_denominator = (
        ALPHA**4
        * (ALPHA + 1) ** 4
        * (ALPHA + 2) ** 4
        * (ALPHA + 3) ** 3
    )
    require(
        sp.expand(denominator_one - common_denominator) == 0
        and sp.expand(denominator_two - common_denominator) == 0,
        "common positive denominator changed",
    )

    transformed = sp.expand(
        30 * numerator_two - (41 - 6 * ALPHA) * numerator_one
    )
    jacobian = sp.expand(
        sp.diff(numerator_one, X) * sp.diff(transformed, Y)
        - sp.diff(numerator_one, Y) * sp.diff(transformed, X)
    )

    power_f, degrees_f = pull_back_to_unit_cube(numerator_one)
    power_g, degrees_g = pull_back_to_unit_cube(transformed)
    power_minus_fx, degrees_fx = pull_back_to_unit_cube(
        -sp.diff(numerator_one, X)
    )
    power_j, degrees_j = pull_back_to_unit_cube(jacobian)
    audit_point = (Fraction(2, 5), Fraction(1, 3), Fraction(3, 7))
    for expression, power in (
        (numerator_one, power_f),
        (transformed, power_g),
        (-sp.diff(numerator_one, X), power_minus_fx),
        (jacobian, power_j),
    ):
        check_pullback(expression, power, audit_point)

    bernstein_f = power_to_bernstein(power_f, degrees_f)
    bernstein_g = power_to_bernstein(power_g, degrees_g)
    bernstein_minus_fx = power_to_bernstein(power_minus_fx, degrees_fx)
    bernstein_j = power_to_bernstein(power_j, degrees_j)
    for power, bernstein_coefficients, degrees in (
        (power_f, bernstein_f, degrees_f),
        (power_g, bernstein_g, degrees_g),
        (power_minus_fx, bernstein_minus_fx, degrees_fx),
        (power_j, bernstein_j, degrees_j),
    ):
        require(
            evaluate_bernstein(bernstein_coefficients, degrees, audit_point)
            == evaluate_power(power, audit_point),
            "power-to-Bernstein conversion mismatch",
        )

    require(
        all(
            coefficient > 0
            for coefficient in face_coefficients(
                bernstein_f, degrees_f, 1, 0
            )
        ),
        "F is not positive on the left moving face",
    )
    require(
        all(
            coefficient < 0
            for coefficient in face_coefficients(
                bernstein_f, degrees_f, 1, 1
            )
        ),
        "F is not negative on the right moving face",
    )
    require(
        all(
            coefficient > 0
            for coefficient in face_coefficients(
                bernstein_g, degrees_g, 2, 0
            )
        ),
        "G is not positive on the bottom moving face",
    )
    require(
        all(
            coefficient < 0
            for coefficient in face_coefficients(
                bernstein_g, degrees_g, 2, 1
            )
        ),
        "G is not negative on the top moving face",
    )
    require(
        all(
            coefficient > 0
            for coefficient in all_coefficients(bernstein_minus_fx)
        ),
        "F_x lost its negative sign in the moving tube",
    )
    require(
        all(
            coefficient > 0
            for coefficient in all_coefficients(bernstein_j)
        ),
        "Jacobian lost its positive sign in the moving tube",
    )

    # Reduce the fourth moment modulo the Gram quadratic.  The displayed
    # expression is g22^3 times the coefficient of z in the remainder.
    q4 = tuple(
        multilinear(
            tuple([lower] * (4 - upper_count) + [upper] * upper_count)
        )
        for upper_count in range(5)
    )
    quartic_coefficients = tuple(
        comb(4, upper_count) * coefficient
        for upper_count, coefficient in enumerate(q4)
    )
    q0, q1, q2 = g11, 2 * g12, g22
    _, k1, k2, k3, k4 = quartic_coefficients
    quartic_linear = sp.cancel(
        k1 * q2**3
        - k2 * q1 * q2**2
        + k3 * (q1**2 - q0 * q2) * q2
        + k4 * (2 * q0 * q1 * q2 - q1**3)
    )
    quartic_numerator, quartic_denominator = sp.together(
        quartic_linear
    ).as_numer_denom()
    expected_quartic_denominator = (
        ALPHA**6
        * (ALPHA + 1) ** 6
        * (ALPHA + 2) ** 6
        * (ALPHA + 3) ** 4
    )
    require(
        sp.expand(quartic_denominator - expected_quartic_denominator) == 0,
        "quartic positive denominator changed",
    )
    power_quartic, degrees_quartic = pull_back_to_unit_cube(
        quartic_numerator
    )
    check_pullback(quartic_numerator, power_quartic, audit_point)
    bernstein_quartic = power_to_bernstein(
        power_quartic, degrees_quartic
    )
    require(
        evaluate_bernstein(
            bernstein_quartic, degrees_quartic, audit_point
        )
        == evaluate_power(power_quartic, audit_point),
        "quartic power-to-Bernstein conversion mismatch",
    )
    require(
        all(
            coefficient > 0
            for coefficient in all_coefficients(bernstein_quartic)
        ),
        "quartic imaginary exit lost positivity",
    )

    alpha_left, alpha_right = Fraction(9, 10), Fraction(11, 10)
    x_min = (
        Fraction(26362, 100000)
        + Fraction(73, 1000) * (alpha_left - 1)
        - Fraction(1, 10000)
    )
    y_min = (
        Fraction(23420, 1000000)
        + Fraction(7, 500) * (alpha_left - 1)
        - Fraction(1, 20000)
    )
    require(
        x_min > 0 and y_min > 0 and alpha_right > alpha_left > 0,
        "positive tube failed",
    )

    print("GAMMA TRANSVERSE HOSTILE HOLOTOPY -- exact Bernstein referee")
    print("alpha_interval=[9/10,11/10]")
    print(
        "moving_x_center=26362/100000+(73/1000)(alpha-1);"
        " halfwidth=1/10000"
    )
    print(
        "moving_y_center=23420/1000000+(7/500)(alpha-1);"
        " halfwidth=1/20000"
    )
    print("row_operation=G=30*I2-(41-6*alpha)*I1; determinant=30")
    print(f"F_degrees={degrees_f} G_degrees={degrees_g}")
    print(
        f"minus_Fx_degrees={degrees_fx}"
        f" Jacobian_degrees={degrees_j}"
    )
    print(f"quartic_exit_degrees={degrees_quartic}")
    print(
        "face_orientation=PASS"
        " unique_transverse_zero_for_every_alpha=PASS"
    )
    print("positive_cones=PASS quartic_imaginary_exit=PASS")
    print(
        "consequence=Gamma moments 1..3 vanish along a unique"
        " real-analytic positive branch"
    )
    print(
        "scope=uniform deformation hostile;"
        " only alpha=1 is the complex-Gaussian radial law"
    )


if __name__ == "__main__":
    main()
