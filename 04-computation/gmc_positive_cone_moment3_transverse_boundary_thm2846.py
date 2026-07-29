#!/usr/bin/env python3
"""Exact referee for the arbitrary-positive-cone moment-three boundary."""

from __future__ import annotations

from math import comb, factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


Y, X = sp.symbols("Y X")

P = (
    933319220872000 * Y**15
    + 5884528273360000 * Y**14
    + 17194209549700000 * Y**13
    + 30835666750660000 * Y**12
    + 37890363245864000 * Y**11
    + 33723003604779000 * Y**10
    + 22402533775665500 * Y**9
    + 11275117818448500 * Y**8
    + 4315433784681000 * Y**7
    + 1247681141032500 * Y**6
    + 267329168334540 * Y**5
    + 40861458466200 * Y**4
    + 4132565626725 * Y**3
    + 231008559525 * Y**2
    + 2764978335 * Y
    - 258946389
)

A = 2 * (
    46665961043600 * Y**10
    + 517246931163000 * Y**9
    + 1734100192751000 * Y**8
    + 2973460269848000 * Y**7
    + 3076524208583550 * Y**6
    + 2056143278933070 * Y**5
    + 909053232722550 * Y**4
    + 264045781471725 * Y**3
    + 48371683260150 * Y**2
    + 5052446879130 * Y
    + 228165400071
)

N = (
    485253235285200 * Y**10
    + 2203659617816000 * Y**9
    + 4482965815747000 * Y**8
    + 5358371582968500 * Y**7
    + 4154724741540300 * Y**6
    + 2179627177491390 * Y**5
    + 783069490441350 * Y**4
    + 190343609071200 * Y**3
    + 29999951982225 * Y**2
    + 2773709682480 * Y
    + 114503538087
)


def quadratic(left: int, right: int) -> int:
    return comb(left + right, left)


def cubic(left: int, middle: int, right: int) -> int:
    total = left + middle + right
    e2 = left * middle + left * right + middle * right
    e3 = left * middle * right
    multinomial = factorial(total) // (
        factorial(left) * factorial(middle) * factorial(right)
    )
    numerator = 2 * (total + 1) ** 2 + total * e2 - e3
    denominator = (left + 1) * (middle + 1) * (right + 1)
    value = sp.Rational(multinomial * numerator, denominator)
    require(value.q == 1, "cubic tensor lost integrality")
    return int(value)


def adjacent_tensor(indices: tuple[int, ...]) -> int:
    """Direct 2^k expansion of L(prod d_i), independent of closed formulas."""
    order = len(indices)
    value = sp.Integer(0)
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
    require(value.q == 1, "adjacent tensor lost integrality")
    return int(value)


def bilinear(first: dict[int, sp.Expr], second: dict[int, sp.Expr]) -> sp.Expr:
    return sp.expand(
        sum(
            first_value * second_value * quadratic(first_index, second_index)
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
        )
    )


def trilinear(
    first: dict[int, sp.Expr],
    second: dict[int, sp.Expr],
    third: dict[int, sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            first_value
            * second_value
            * third_value
            * cubic(first_index, second_index, third_index)
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
            for third_index, third_value in third.items()
        )
    )


def four_linear(
    first: dict[int, sp.Expr],
    second: dict[int, sp.Expr],
    third: dict[int, sp.Expr],
    fourth: dict[int, sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            first_value
            * second_value
            * third_value
            * fourth_value
            * adjacent_tensor(
                (first_index, second_index, third_index, fourth_index)
            )
            for first_index, first_value in first.items()
            for second_index, second_value in second.items()
            for third_index, third_value in third.items()
            for fourth_index, fourth_value in fourth.items()
        )
    )


def main() -> None:
    lower = {1: sp.Integer(1), 3: X}
    upper = {2: sp.Integer(1), 3: Y}

    for left in range(5):
        for middle in range(5):
            require(
                quadratic(left, middle) == adjacent_tensor((left, middle)),
                "quadratic/direct tensor mismatch",
            )
            for right in range(5):
                require(
                    cubic(left, middle, right)
                    == adjacent_tensor((left, middle, right)),
                    "cubic/direct tensor mismatch",
                )

    g11 = bilinear(lower, lower)
    g12 = bilinear(lower, upper)
    g22 = bilinear(upper, upper)
    t111 = trilinear(lower, lower, lower)
    t112 = trilinear(lower, lower, upper)
    t122 = trilinear(lower, upper, upper)
    t222 = trilinear(upper, upper, upper)

    invariant_one = sp.expand(
        3 * t112 * g11 * g22
        - t222 * g11**2
        - 2 * t111 * g12 * g22
    )
    invariant_two = sp.expand(
        3 * t122 * g11 * g22
        - 2 * t222 * g12 * g11
        - t111 * g22**2
    )

    # The claimed algebraic substitution annihilates both invariants.
    numerator_one = sp.together(invariant_one.subs(X, N / A)).as_numer_denom()[0]
    numerator_two = sp.together(invariant_two.subs(X, N / A)).as_numer_denom()[0]
    require(sp.rem(numerator_one, P, domain=sp.QQ) == 0, "I1 remainder nonzero")
    require(sp.rem(numerator_two, P, domain=sp.QQ) == 0, "I2 remainder nonzero")

    # Recover the displayed positive-root factor from exact elimination.
    resultant = sp.resultant(invariant_one, invariant_two, X)
    quotient, remainder = sp.div(resultant, P, domain=sp.QQ)
    require(remainder == 0 and quotient != 0, "resultant factor mismatch")

    polynomial = sp.Poly(P, Y)
    coefficients = polynomial.all_coeffs()
    require(
        all(coefficient > 0 for coefficient in coefficients[:-1])
        and coefficients[-1] < 0,
        "Descartes sign pattern changed",
    )
    left = sp.Rational(2341, 100000)
    right = sp.Rational(2342, 100000)
    require(P.subs(Y, left) < 0 < P.subs(Y, right), "root bracket changed")
    require(polynomial.count_roots(left, right) == 1, "root isolation changed")
    require(
        all(coefficient > 0 for coefficient in sp.Poly(A, Y).all_coeffs())
        and all(coefficient > 0 for coefficient in sp.Poly(N, Y).all_coeffs()),
        "positive x certificate changed",
    )

    # Independent rational-box proof.  Positivity of all Bernstein
    # coefficients bounds a polynomial on the whole interval/rectangle.
    transformed = sp.expand(6 * invariant_two - 7 * invariant_one)
    x_left = sp.Rational(2636, 10000)
    x_right = sp.Rational(2637, 10000)
    y_bottom = sp.Rational(23418, 1000000)
    y_top = sp.Rational(23421, 1000000)

    def univariate_bernstein_coefficients(
        expression: sp.Expr,
        variable: sp.Symbol,
        interval_left: sp.Rational,
        interval_right: sp.Rational,
    ) -> tuple[sp.Expr, ...]:
        parameter = sp.symbols("parameter")
        pulled_back = sp.Poly(
            sp.expand(
                expression.subs(
                    variable,
                    interval_left + (interval_right - interval_left) * parameter,
                )
            ),
            parameter,
        )
        degree = pulled_back.degree()
        power_coefficients = tuple(
            pulled_back.nth(index) for index in range(degree + 1)
        )
        return tuple(
            sp.expand(
                sum(
                    power_coefficients[index]
                    * sp.Rational(comb(bernstein_index, index), comb(degree, index))
                    for index in range(bernstein_index + 1)
                )
            )
            for bernstein_index in range(degree + 1)
        )

    face_coefficients = (
        univariate_bernstein_coefficients(
            invariant_one.subs(X, x_left), Y, y_bottom, y_top
        ),
        univariate_bernstein_coefficients(
            -invariant_one.subs(X, x_right), Y, y_bottom, y_top
        ),
        univariate_bernstein_coefficients(
            transformed.subs(Y, y_bottom), X, x_left, x_right
        ),
        univariate_bernstein_coefficients(
            -transformed.subs(Y, y_top), X, x_left, x_right
        ),
    )
    require(
        all(coefficient > 0 for face in face_coefficients for coefficient in face),
        "rectangle face orientation changed",
    )

    def bivariate_bernstein_coefficients(
        expression: sp.Expr,
    ) -> tuple[sp.Expr, ...]:
        parameter_x, parameter_y = sp.symbols("parameter_x parameter_y")
        pulled_back = sp.Poly(
            sp.expand(
                expression.subs(
                    {
                        X: x_left + (x_right - x_left) * parameter_x,
                        Y: y_bottom + (y_top - y_bottom) * parameter_y,
                    }
                )
            ),
            parameter_x,
            parameter_y,
        )
        degree_x = pulled_back.degree(parameter_x)
        degree_y = pulled_back.degree(parameter_y)
        coefficients = []
        for bernstein_x in range(degree_x + 1):
            for bernstein_y in range(degree_y + 1):
                coefficients.append(
                    sp.expand(
                        sum(
                            pulled_back.coeff_monomial(
                                parameter_x**power_x * parameter_y**power_y
                            )
                            * sp.Rational(
                                comb(bernstein_x, power_x),
                                comb(degree_x, power_x),
                            )
                            * sp.Rational(
                                comb(bernstein_y, power_y),
                                comb(degree_y, power_y),
                            )
                            for power_x in range(bernstein_x + 1)
                            for power_y in range(bernstein_y + 1)
                        )
                    )
                )
        return tuple(coefficients)

    jacobian = sp.expand(
        sp.diff(invariant_one, X) * sp.diff(transformed, Y)
        - sp.diff(invariant_one, Y) * sp.diff(transformed, X)
    )
    signed_cell_polynomials = (
        -sp.diff(invariant_one, X),
        sp.diff(invariant_one, Y),
        -sp.diff(transformed, X),
        -sp.diff(transformed, Y),
        jacobian,
    )
    require(
        all(
            coefficient > 0
            for expression in signed_cell_polynomials
            for coefficient in bivariate_bernstein_coefficients(expression)
        ),
        "rectangle transversality changed",
    )

    # The two exact remainder equations force cubic divisibility.
    z = sp.symbols("z")
    q_form = g11 + 2 * g12 * z + g22 * z**2
    c_form = t111 + 3 * t112 * z + 3 * t122 * z**2 + t222 * z**3
    division_remainder = sp.Poly(c_form, z).rem(sp.Poly(q_form, z)).as_expr()
    restored_remainder = sp.together(division_remainder.subs(X, N / A))
    restored_numerator = restored_remainder.as_numer_denom()[0]
    for coefficient in sp.Poly(restored_numerator, z).all_coeffs():
        require(sp.rem(coefficient, P, domain=sp.QQ) == 0, "division did not vanish")

    # Exact positive control: neither quadratic root also kills moment four.
    # Here q4[j]=L(U^(4-j)V^j).
    displayed_q4 = (
        20790000 * X**4 + 4132800 * X**3 + 407400 * X**2 + 24840 * X + 864,
        20790000 * X**3 * Y
        + 5197500 * X**3
        + 3099600 * X**2 * Y
        + 887670 * X**2
        + 203700 * X * Y
        + 68460 * X
        + 6210 * Y
        + 2550,
        20790000 * X**2 * Y**2
        + 10395000 * X**2 * Y
        + 1388100 * X**2
        + 2066400 * X * Y**2
        + 1183560 * X * Y
        + 183120 * X
        + 67900 * Y**2
        + 45640 * Y
        + 8440,
        20790000 * X * Y**3
        + 15592500 * X * Y**2
        + 4164300 * X * Y
        + 398160 * X
        + 1033200 * Y**3
        + 887670 * Y**2
        + 274680 * Y
        + 30870,
        20790000 * Y**4
        + 20790000 * Y**3
        + 8328600 * Y**2
        + 1592640 * Y
        + 123480,
    )
    derived_q4 = tuple(
        four_linear(
            *([lower] * (4 - upper_count) + [upper] * upper_count)
        )
        for upper_count in range(5)
    )
    require(
        all(
            sp.expand(derived - displayed) == 0
            for derived, displayed in zip(derived_q4, displayed_q4)
        ),
        "quartic tensor reconstruction changed",
    )
    quartic_coefficients = [
        comb(4, upper_count) * expression
        for upper_count, expression in enumerate(derived_q4)
    ]
    q0 = g11
    q1 = 2 * g12
    q2 = g22

    # Reduce the quartic modulo q2*z^2+q1*z+q0.  Multiplication by q2^3
    # clears every division.  Only the z coefficient is needed: at the
    # upper-half-plane quadratic root it is the sign of Im L(H^4).
    k0, k1, k2, k3, k4 = quartic_coefficients
    del k0
    quartic_linear_numerator = sp.expand(
        k1 * q2**3
        - k2 * q1 * q2**2
        + k3 * (q1**2 - q0 * q2) * q2
        + k4 * (2 * q0 * q1 * q2 - q1**3)
    )
    require(
        all(
            coefficient > 0
            for coefficient in bivariate_bernstein_coefficients(
                quartic_linear_numerator
            )
        ),
        "quartic imaginary transversality changed",
    )

    positive_root = next(
        root
        for root in sp.nroots(P, n=40, maxsteps=200)
        if abs(sp.im(root)) < sp.Rational(1, 10) ** 30 and sp.re(root) > 0
    )
    x_value = sp.N((N / A).subs(Y, positive_root), 30)

    print("ARBITRARY POSITIVE PASCAL CONE MOMENT-THREE BOUNDARY -- exact referee")
    print("status=SCRATCH THEOREM; NO CANON OR ID RESERVATION")
    print("supports=U:d1+d3; V:d2+d3; both divisible by s")
    print("positive_root_count=1")
    print(f"y_bracket=({left},{right})")
    print(f"y_decimal={sp.N(positive_root, 24)}")
    print(f"x_decimal={x_value}")
    print("I1_mod_p=0 I2_mod_p=0 cubic_remainder_mod_p=0")
    print("rational_box_unique_transverse_zero=1")
    print("quartic_remainder_imaginary_positive=1 Gaussian_moment8_nonzero=1")
    print("consequence=factorial moments 1..3 and Gaussian moments 1..6 can vanish")
    print("scope=finite-moment hostile only; not a GMC2 counterexample")


if __name__ == "__main__":
    main()
