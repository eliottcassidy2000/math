#!/usr/bin/env python3
"""Exact referee for the shifted positive-cone transverse family."""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


n, x, y, t = sp.symbols("n x y t")
a, b, A_coordinate, D_coordinate = sp.symbols(
    "a b A_coordinate D_coordinate"
)
sqrt_two = sp.sqrt(2)


def factorial_ratio(base: sp.Expr, shift: int) -> sp.Expr:
    """Return (base+shift)!/base! as a rational product."""
    if shift >= 0:
        return sp.prod(base + index for index in range(1, shift + 1))
    return 1 / sp.prod(
        base - index + 1 for index in range(1, -shift + 1)
    )


def normalized_bilinear_tensor(left: int, right: int) -> sp.Expr:
    """L(d_(n+left)d_(n+right))/L(d_(n+2)^2)."""
    return sp.cancel(
        factorial_ratio(2 * n + 4, left + right - 4)
        * factorial_ratio(n + left, 2 - left)
        * factorial_ratio(n + right, 2 - right)
    )


def cubic_shape(left: int, middle: int, right: int) -> sp.Expr:
    indices = (n + left, n + middle, n + right)
    total = sum(indices)
    second = sum(
        indices[first] * indices[second_index]
        for first in range(3)
        for second_index in range(first + 1, 3)
    )
    third = sp.prod(indices)
    return sp.cancel(
        (
            2 * (total + 1) ** 2
            + total * second
            - third
        )
        / sp.prod(index + 1 for index in indices)
    )


cubic_shape_222 = cubic_shape(2, 2, 2)


def normalized_trilinear_tensor(
    left: int, middle: int, right: int
) -> sp.Expr:
    """L(d_(n+left)d_(n+middle)d_(n+right))/L(d_(n+2)^3)."""
    return sp.cancel(
        factorial_ratio(3 * n + 6, left + middle + right - 6)
        * factorial_ratio(n + left, 2 - left)
        * factorial_ratio(n + middle, 2 - middle)
        * factorial_ratio(n + right, 2 - right)
        * cubic_shape(left, middle, right)
        / cubic_shape_222
    )


def direct_adjacent_tensor(indices: tuple[int, ...]) -> sp.Integer:
    """Direct factorial expansion, used only as an independent control."""
    value = sp.Integer(0)
    order = len(indices)
    for mask in range(1 << order):
        degrees = []
        sign = 1
        for position, index in enumerate(indices):
            if mask & (1 << position):
                degrees.append(index + 1)
            else:
                degrees.append(index)
                sign = -sign
        term = sp.factorial(sum(degrees))
        for degree in degrees:
            term /= sp.factorial(degree)
        value += sign * term
    require(value.q == 1, "direct adjacent tensor lost integrality")
    return sp.Integer(value)


def closed_cubic_tensor(left: int, middle: int, right: int) -> sp.Integer:
    total = left + middle + right
    second = left * middle + left * right + middle * right
    third = left * middle * right
    multinomial = factorial(total) // (
        factorial(left) * factorial(middle) * factorial(right)
    )
    value = sp.Rational(
        multinomial
        * (2 * (total + 1) ** 2 + total * second - third),
        (left + 1) * (middle + 1) * (right + 1),
    )
    require(value.q == 1, "closed cubic tensor lost integrality")
    return sp.Integer(value)


lower = {0: sp.Integer(1), 2: x}
upper = {1: sp.Integer(1), 2: y}


def bilinear(
    first: dict[int, sp.Expr], second: dict[int, sp.Expr]
) -> sp.Expr:
    return sp.factor(
        sum(
            first_value
            * second_value
            * normalized_bilinear_tensor(first_index, second_index)
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
        )
    )


def trilinear(
    first: dict[int, sp.Expr],
    second: dict[int, sp.Expr],
    third: dict[int, sp.Expr],
) -> sp.Expr:
    return sp.factor(
        sum(
            first_value
            * second_value
            * third_value
            * normalized_trilinear_tensor(
                first_index, second_index, third_index
            )
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
            for third_index, third_value in third.items()
        )
    )


g11 = bilinear(lower, lower)
g12 = bilinear(lower, upper)
g22 = bilinear(upper, upper)
t111 = trilinear(lower, lower, lower)
t112 = trilinear(lower, lower, upper)
t122 = trilinear(lower, upper, upper)
t222 = trilinear(upper, upper, upper)

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


def sign_variations(values: list[sp.Expr]) -> int:
    nonzero = [value for value in values if value != 0]
    return sum(
        bool(left > 0) != bool(right > 0)
        for left, right in zip(nonzero, nonzero[1:])
    )


def all_coefficients_have_sign(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...], sign: int
) -> bool:
    coefficients = sp.Poly(expression, *variables).coeffs()
    if sign > 0:
        return all(bool(coefficient > 0) for coefficient in coefficients)
    return all(bool(coefficient < 0) for coefficient in coefficients)


def normalized_four_tensor(indices: tuple[int, int, int, int]) -> sp.Expr:
    """Fourth adjacent tensor divided by the positive base multinomial."""
    value = sp.Integer(0)
    for mask in range(16):
        increments = tuple((mask >> position) & 1 for position in range(4))
        increment_sum = sum(increments)
        term = (-1) ** (4 - increment_sum)
        term *= factorial_ratio(
            4 * n + 8, sum(indices) + increment_sum - 8
        )
        for index, increment in zip(indices, increments):
            term *= factorial_ratio(
                n + index + increment, 2 - index - increment
            )
        value += term
    return sp.cancel(value)


def four_linear(
    first: dict[int, sp.Expr],
    second: dict[int, sp.Expr],
    third: dict[int, sp.Expr],
    fourth: dict[int, sp.Expr],
) -> sp.Expr:
    return sp.factor(
        sum(
            first_value
            * second_value
            * third_value
            * fourth_value
            * normalized_four_tensor(
                (first_index, second_index, third_index, fourth_index)
            )
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
            for third_index, third_value in third.items()
            for fourth_index, fourth_value in fourth.items()
        )
    )


def main() -> None:
    # Two independent tensor controls at several literal depths.
    for depth in range(1, 4):
        for left in range(3):
            for right in range(3):
                direct_two = direct_adjacent_tensor(
                    (depth + left, depth + right)
                )
                closed_two = sp.binomial(
                    2 * depth + left + right, depth + left
                )
                require(direct_two == closed_two, "bilinear control failed")
                for third in range(3):
                    direct_three = direct_adjacent_tensor(
                        (depth + left, depth + right, depth + third)
                    )
                    closed_three = closed_cubic_tensor(
                        depth + left, depth + right, depth + third
                    )
                    require(
                        direct_three == closed_three,
                        "trilinear control failed",
                    )

    # Every denominator used in the cleared invariant system is positive
    # on n>=1.  Its leading x coefficient is also never lost for y>0.
    numerator_one, denominator_one = sp.fraction(
        sp.together(invariant_one)
    )
    numerator_two, denominator_two = sp.fraction(
        sp.together(invariant_two)
    )
    for denominator in (denominator_one, denominator_two):
        require(
            all_coefficients_have_sign(denominator, (n,), 1),
            "a normalized-invariant denominator lost positivity",
        )
    for numerator in (numerator_one, numerator_two):
        leading = sp.Poly(numerator, x).LC()
        require(
            all_coefficients_have_sign(leading, (n, y), -1),
            "an x-leading coefficient can vanish in the positive quadrant",
        )

    # Exact limiting tensors and the singular limiting line.
    limiting_one = sp.factor(
        sp.limit(invariant_one.subs(n, 1 / t), t, 0)
    )
    limiting_two = sp.factor(
        sp.limit(invariant_two.subs(n, 1 / t), t, 0)
    )
    expected_one = -(
        (4 * x + 1)
        * (6 * x - 5 * y - 1) ** 2
        * (108 * x * y + 48 * x + 17 * y + 7)
        / sp.Integer(186624)
    )
    expected_two = -(
        (2 * y + 1)
        * (6 * x - 5 * y - 1) ** 2
        * (54 * x * y + 21 * x + 11 * y + 4)
        / sp.Integer(46656)
    )
    require(
        sp.factor(limiting_one - expected_one) == 0,
        "first limiting invariant changed",
    )
    require(
        sp.factor(limiting_two - expected_two) == 0,
        "second limiting invariant changed",
    )

    limiting_line = (6 * x - 1) / 5
    line_one = sp.factor(
        sp.limit(
            invariant_one.subs({n: 1 / t, y: limiting_line}) / t,
            t,
            0,
        )
    )
    line_two = sp.factor(
        sp.limit(
            invariant_two.subs({n: 1 / t, y: limiting_line}) / t,
            t,
            0,
        )
    )
    common_line_factor = (
        (4 * x + 1) ** 2
        * (9 * x + 1)
        * (18 * x**2 - 36 * x - 7)
    )
    require(
        sp.factor(
            line_one
            - common_line_factor
            / sp.Integer(3888000)
        )
        == 0,
        "first line coefficient changed",
    )
    require(
        sp.factor(
            line_two
            - common_line_factor
            / sp.Integer(3240000)
        )
        == 0,
        "second line coefficient changed",
    )

    x_zero = 1 + 5 * sqrt_two / 6
    y_zero = 1 + sqrt_two
    a_zero = -(74184 + 52463 * sqrt_two) / sp.Integer(1296)
    b_zero = -(14808 + 10495 * sqrt_two) / sp.Integer(216)
    d_zero = -sp.Rational(2, 3) + sqrt_two / 18

    prefactor_one = 25 * (-99 * sqrt_two - 140) / sp.Integer(5184)
    prefactor_two = sp.Rational(6, 5) * prefactor_one

    first_blowup = sp.cancel(
        invariant_one.subs(
            {
                n: 1 / t,
                x: x_zero + a * t,
                y: y_zero + b * t,
            }
        )
    )
    second_blowup = sp.cancel(
        invariant_two.subs(
            {
                n: 1 / t,
                x: x_zero + a * t,
                y: y_zero + b * t,
            }
        )
    )
    F1 = sp.factor(
        sp.limit(first_blowup / t**2, t, 0) / prefactor_one,
        extension=sqrt_two,
    )
    F2 = sp.factor(
        sp.limit(second_blowup / t**2, t, 0) / prefactor_two,
        extension=sqrt_two,
    )

    expected_F1 = (
        D_coordinate**2
        + (12 + 38 * sqrt_two / 5) * D_coordinate
        - 2 * sqrt_two * A_coordinate / 5
        - 829 * sqrt_two / 45
        - sp.Rational(13867, 540)
    ) / 36
    expected_F2 = (
        D_coordinate**2
        + (12 + 23 * sqrt_two / 3) * D_coordinate
        - 2 * sqrt_two * A_coordinate / 5
        - 827 * sqrt_two / 45
        - sp.Rational(13871, 540)
    ) / 36
    coordinate_substitution = {
        a: A_coordinate,
        b: (6 * A_coordinate - D_coordinate) / 5,
    }
    require(
        sp.factor(
            F1.subs(coordinate_substitution) - expected_F1,
            extension=sqrt_two,
        )
        == 0,
        "first blown-up equation changed",
    )
    require(
        sp.factor(
            F2.subs(coordinate_substitution) - expected_F2,
            extension=sqrt_two,
        )
        == 0,
        "second blown-up equation changed",
    )
    require(
        sp.factor(
            expected_F1
            - expected_F2
            + sqrt_two * (D_coordinate - d_zero) / 540,
            extension=sqrt_two,
        )
        == 0,
        "triangular difference changed",
    )
    require(
        sp.factor(expected_F1.subs(
            {A_coordinate: a_zero, D_coordinate: d_zero}
        ), extension=sqrt_two)
        == 0
        and sp.factor(expected_F2.subs(
            {A_coordinate: a_zero, D_coordinate: d_zero}
        ), extension=sqrt_two)
        == 0
        and sp.factor(6 * a_zero - 5 * b_zero - d_zero) == 0,
        "blown-up base point changed",
    )
    blowup_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(F1, a), sp.diff(F1, b)],
                    [sp.diff(F2, a), sp.diff(F2, b)],
                ]
            )
        ).subs({a: a_zero, b: b_zero}),
        extension=sqrt_two,
    )
    require(
        blowup_jacobian == sp.Rational(1, 4860),
        "blown-up Jacobian changed",
    )

    # Exact all-integer elimination.  The last subresultant is the
    # resultant; the preceding linear subresultant supplies x.
    subresultants = sp.subresultants(
        numerator_one, numerator_two, x
    )
    require(
        [sp.degree(item, x) for item in subresultants]
        == [4, 3, 2, 1, 0],
        "subresultant degree profile changed",
    )
    resultant = subresultants[-1]
    resultant_factors = sp.factor_list(resultant)[1]
    y_factors = [
        (factor, exponent)
        for factor, exponent in resultant_factors
        if sp.degree(factor, y) > 0
    ]
    require(
        sorted((sp.degree(factor, y), exponent) for factor, exponent in y_factors)
        == [(2, 2), (15, 1)],
        "resultant y-factor profile changed",
    )
    for factor, _ in resultant_factors:
        if sp.degree(factor, y) == 0:
            oriented = factor
            if sp.Poly(oriented, n).LC() < 0:
                oriented = -oriented
            require(
                all_coefficients_have_sign(oriented, (n,), 1),
                "an n-only resultant factor can vanish for n>=1",
            )

    quadratic_factor = next(
        factor for factor, _ in y_factors if sp.degree(factor, y) == 2
    )
    positive_gram_factor = (
        n + 2
        + 2 * (2 * n + 3) * y
        + 2 * (2 * n + 3) * y**2
    )
    gram_ratio = sp.cancel(
        sp.Poly(quadratic_factor, y).LC()
        / sp.Poly(positive_gram_factor, y).LC()
    )
    require(
        sp.expand(quadratic_factor - gram_ratio * positive_gram_factor)
        == 0,
        "quadratic resultant factor is not the positive Gram factor",
    )

    elimination_polynomial = next(
        factor for factor, _ in y_factors if sp.degree(factor, y) == 15
    )
    if sp.Poly(elimination_polynomial, y).LC().subs(n, 1) < 0:
        elimination_polynomial = -elimination_polynomial
    elimination_coefficients = sp.Poly(
        elimination_polynomial, y
    ).all_coeffs()
    require(
        len(elimination_coefficients) == 16,
        "degree-fifteen elimination factor changed",
    )
    for coefficient in elimination_coefficients[:3]:
        require(
            all_coefficients_have_sign(coefficient, (n,), 1),
            "a permanent-positive elimination coefficient changed",
        )

    first_negative_depths = (
        102,
        43,
        25,
        16,
        11,
        8,
        6,
        5,
        4,
        3,
        2,
        2,
        1,
    )
    require(
        all(
            left >= right
            for left, right in zip(
                first_negative_depths, first_negative_depths[1:]
            )
        ),
        "coefficient threshold staircase changed",
    )
    for coefficient, threshold in zip(
        elimination_coefficients[3:], first_negative_depths
    ):
        coefficient_polynomial = sp.Poly(coefficient, n)
        require(
            sign_variations(coefficient_polynomial.all_coeffs()) == 1,
            "an elimination coefficient lost its unique positive root",
        )
        require(
            coefficient.subs(n, threshold - 1) > 0
            and coefficient.subs(n, threshold) < 0,
            "an elimination coefficient threshold changed",
        )

    linear_subresultant = subresultants[-2]
    linear_coefficient = sp.Poly(linear_subresultant, x).coeff_monomial(x)
    constant_coefficient = sp.Poly(linear_subresultant, x).coeff_monomial(1)
    common_content = sp.gcd(linear_coefficient, constant_coefficient)
    for content_factor, _ in sp.factor_list(common_content)[1]:
        oriented_content_factor = content_factor
        if sp.Poly(oriented_content_factor, n, y).LC() < 0:
            oriented_content_factor = -oriented_content_factor
        require(
            all_coefficients_have_sign(
                oriented_content_factor, (n, y), 1
            ),
            "linear-subresultant content can vanish in the positive quadrant",
        )
    linear_coefficient = sp.cancel(linear_coefficient / common_content)
    constant_coefficient = sp.cancel(constant_coefficient / common_content)
    require(
        sp.degree(linear_coefficient, y) == 10
        and sp.degree(constant_coefficient, y) == 10
        and all_coefficients_have_sign(
            linear_coefficient, (n, y), 1
        )
        and all_coefficients_have_sign(
            constant_coefficient, (n, y), -1
        ),
        "positive linear subresultant selector changed",
    )

    # The fourth-moment remainder has a positive limiting imaginary
    # coefficient on the IFT branch.  This proves eventual, not all-depth,
    # moment-four escape.
    fourth_moments = tuple(
        four_linear(
            *([lower] * (4 - upper_count) + [upper] * upper_count)
        )
        for upper_count in range(5)
    )
    quartic_coefficients = tuple(
        sp.binomial(4, upper_count) * fourth_moments[upper_count]
        for upper_count in range(5)
    )
    q0, q1, q2 = g11, 2 * g12, g22
    _, k1, k2, k3, k4 = quartic_coefficients
    quartic_linear_numerator = sp.factor(
        k1 * q2**3
        - k2 * q1 * q2**2
        + k3 * (q1**2 - q0 * q2) * q2
        + k4 * (2 * q0 * q1 * q2 - q1**3)
    )
    branch_substitution = {
        n: 1 / t,
        x: x_zero + a_zero * t,
        y: y_zero + b_zero * t,
    }
    quartic_limit = sp.factor(
        sp.limit(
            sp.cancel(quartic_linear_numerator.subs(branch_substitution)),
            t,
            0,
        ),
        extension=sqrt_two,
    )
    expected_quartic_limit = (
        3 * (208885 + 147704 * sqrt_two) / sp.Integer(262144)
    )
    require(
        sp.radsimp(quartic_limit - expected_quartic_limit) == 0
        and expected_quartic_limit > 0,
        "eventual quartic escape changed",
    )
    gram_determinant = sp.factor(g11 * g22 - g12**2)
    gram_limit = sp.factor(
        sp.limit(
            sp.cancel(gram_determinant.subs(branch_substitution) / t),
            t,
            0,
        ),
        extension=sqrt_two,
    )
    require(
        sp.radsimp(gram_limit - (17 + 12 * sqrt_two) / 1152) == 0
        and (17 + 12 * sqrt_two) / 1152 > 0,
        "branch Gram-separation asymptotic changed",
    )

    # Literal n=1 hostile control: the unique positive y root lies in the
    # already-audited THM-2846 bracket, and its selected x is positive.
    depth_one_eliminant = sp.Poly(
        elimination_polynomial.subs(n, 1), y
    )
    bracket_left = sp.Rational(2341, 100000)
    bracket_right = sp.Rational(2342, 100000)
    require(
        depth_one_eliminant.eval(bracket_left) < 0
        < depth_one_eliminant.eval(bracket_right)
        and depth_one_eliminant.count_roots(
            bracket_left, bracket_right
        )
        == 1,
        "depth-one hostile control changed",
    )

    print("SHIFTED POSITIVE-CONE TRANSVERSE FAMILY -- exact referee")
    print("status=PROOF-COMPLETE-CANDIDATE+VERIFIED-EXACT")
    print("tensor_controls=bilinear+cubic direct expansion depths1..3")
    print("normalized_limit_line=6x-5y-1")
    print("x_limit=1+5sqrt2/6 y_limit=1+sqrt2")
    print("blowup_jacobian=1/4860")
    print("x_first_correction=-(74184+52463sqrt2)/1296")
    print("y_first_correction=-(14808+10495sqrt2)/216")
    print("resultant_profile=positive_factors*Gram(y)^2*P15(y)")
    print("P15_positive_root_count=1 for every integer n>=1")
    print("linear_subresultant=A*x+C with A>0 C<0")
    print("positive_transverse_solution_count=1 for every integer n>=1")
    print(
        "quartic_remainder_limit="
        "3(208885+147704sqrt2)/262144>0"
    )
    print("moment4_escape=proved for all sufficiently large n only")
    print("scope=moment-three hostile family; not a GMC2 counterexample")


if __name__ == "__main__":
    main()
