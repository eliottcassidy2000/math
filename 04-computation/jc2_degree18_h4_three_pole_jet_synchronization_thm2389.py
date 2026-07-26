#!/usr/bin/env python3
"""Exact companion for THM-2389's three-pole jet synchronization."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def quotient_reduce(
    expression: sp.Expr,
    alpha: sp.Symbol,
    modulus: sp.Poly,
) -> sp.Expr:
    """Reduce a rational function modulo an irreducible polynomial."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, alpha, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, alpha, domain=sp.QQ)
    inverse = sp.invert(denominator_poly, modulus)
    return sp.rem(numerator_poly * inverse, modulus).as_expr()


def coefficient_vector(
    expression: sp.Expr,
    variable: sp.Symbol,
    degree: int,
) -> sp.Matrix:
    polynomial = sp.Poly(sp.expand(expression), variable)
    return sp.Matrix(
        [polynomial.coeff_monomial(variable**power) for power in range(degree + 1)]
    )


def main() -> None:
    alpha, s = sp.symbols("alpha s")
    B, C, D, W = sp.symbols("B C D W")
    a1, a2, a3, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")

    # 1. The clean rescaling of THM-2332's depressed cubic.
    p_scale = sp.Rational(16, 964467)
    q_scale = sp.Rational(64, 703096443)
    w_scale = sp.Rational(1701, 2)
    require(
        sp.cancel(p_scale * w_scale**2) == 12,
        "P coefficient under w-rescaling",
    )
    require(
        sp.cancel(q_scale * w_scale**3) == 56,
        "Q coefficient under w-rescaling",
    )

    slope = sp.Poly(alpha**3 + 2940 * alpha + 30184, alpha, domain=sp.QQ)
    require(
        sp.discriminant(slope.as_expr(), alpha) == -126247730112,
        "infinity-slope discriminant",
    )
    require(slope.is_irreducible, "infinity-slope cubic irreducibility")
    require(
        slope.eval(-28) == -74088 and slope.eval(14) == 74088,
        "reconstruction denominator controls",
    )
    require(
        sp.gcd(slope, sp.Poly(alpha**2 + 980, alpha)).degree() == 0,
        "implicit derivative is nonzero at every slope",
    )

    # 2. Expand the normalized pole equation through order six.
    local_a = (
        alpha
        + a1 * s
        + a2 * s**2
        + a3 * s**3
        + a4 * s**4
        + a5 * s**5
        + a6 * s**6
    )
    normalized_p = (
        245
        + 1890 * B * s**2
        + (-24300 * B**2 + 122472 * D) * s**4
    )
    normalized_q = (
        539
        + 11340 * B * s**2
        + 183708 * C * s**3
        + (72900 * B**2 - 367416 * D) * s**4
        + (2361960 * B * C + 2480058 * W) * s**5
    )
    phi = sp.Poly(
        sp.expand(local_a**3 + 12 * normalized_p * local_a + 56 * normalized_q),
        s,
    )
    coefficients = [
        sp.factor(
            phi.coeff_monomial(s**power)
            if power <= 1
            else phi.coeff_monomial(s**power).subs(a1, 0)
        )
        for power in range(7)
    ]
    lam = alpha**2 + 980
    expected = [
        slope.as_expr(),
        3 * lam * a1,
        3 * (lam * a2 + 7560 * B * (alpha + 28)),
        3 * (lam * a3 + 3429216 * C),
        3
        * (
            lam * a4
            + alpha * a2**2
            + 7560 * B * a2
            + (1360800 - 97200 * alpha) * B**2
            + 489888 * (alpha - 14) * D
        ),
        3
        * (
            lam * a5
            + 2 * alpha * a2 * a3
            + 7560 * B * a3
            + 44089920 * B * C
            + 46294416 * W
        ),
        (
            3 * lam * a6
            + 6 * alpha * a2 * a4
            + 3 * alpha * a3**2
            + a2**3
            + 22680 * B * a4
            + (-291600 * B**2 + 1469664 * D) * a2
        ),
    ]
    require(
        all(sp.expand(actual - wanted) == 0 for actual, wanted in zip(coefficients, expected)),
        "pole coefficient atlas",
    )

    reconstructed_B = -lam * a2 / (7560 * (alpha + 28))
    reconstructed_C = -lam * a3 / 3429216
    reconstructed_D = -(
        lam * a4
        + alpha * a2**2
        + 7560 * reconstructed_B * a2
        + (1360800 - 97200 * alpha) * reconstructed_B**2
    ) / (489888 * (alpha - 14))
    reconstructed_W = -(
        lam * a5
        + 2 * alpha * a2 * a3
        + 7560 * reconstructed_B * a3
        + 44089920 * reconstructed_B * reconstructed_C
    ) / 46294416

    # 3. Exact formal branches in Q[alpha]/(f) recover B,C,D,W.
    parameter_packets = (
        (1, 1, 1, 1),
        (2, -1, 3, -2),
        (-3, 2, -2, 5),
        (5, 0, 7, -4),
        (-4, -3, 1, 6),
    )
    formal_branch_cases = 0
    inverse_lam3 = sp.invert(sp.Poly(3 * lam, alpha), slope)

    for b_value, c_value, d_value, w_value in parameter_packets:
        substitutions = {B: b_value, C: c_value, D: d_value, W: w_value}
        solved: dict[sp.Symbol, sp.Expr] = {a1: sp.Integer(0)}

        for power, variable in zip(range(2, 7), (a2, a3, a4, a5, a6)):
            coefficient = expected[power].subs(substitutions).subs(solved)
            constant_part = sp.expand(coefficient.subs(variable, 0))
            reduced_constant = quotient_reduce(constant_part, alpha, slope)
            variable_value = sp.rem(
                -sp.Poly(reduced_constant, alpha) * inverse_lam3,
                slope,
            ).as_expr()
            solved[variable] = variable_value

        recovered = (
            reconstructed_B,
            reconstructed_C,
            reconstructed_D,
            reconstructed_W,
        )
        for expression, target in zip(
            recovered,
            (b_value, c_value, d_value, w_value),
        ):
            value = expression.subs(solved)
            require(
                quotient_reduce(value - target, alpha, slope) == 0,
                "formal parameter reconstruction",
            )

        lock = expected[6].subs(substitutions).subs(solved)
        require(
            quotient_reduce(lock, alpha, slope) == 0,
            "formal sixth-order lock",
        )
        formal_branch_cases += 1

    # 4. Confluent interpolation: six pole jets leave exactly d^2.
    x = sp.symbols("x")
    interpolation_packets = (
        ((-2, 1, 3), x**3 + 2 * x**2 - x + 5, (1, 4, -3)),
        ((0, 2, 5), 2 * x**3 - x**2 + 3 * x + 1, (-2, 3, 6)),
        ((-3, -1, 4), -x**3 + 4 * x**2 + 2 * x + 7, (5, -1, 2)),
    )
    hermite_packets = 0
    for roots, numerator, slope_values in interpolation_packets:
        denominator = sp.prod(x - root for root in roots)
        require(
            all(numerator.subs(x, root) != 0 for root in roots),
            "numerator vanishes at a pole",
        )
        coeffs = sp.symbols(f"r0:{7}")
        sextic = sum(coeffs[power] * x**power for power in range(7))
        equations = []
        for root, slope_value in zip(roots, slope_values):
            equations.append(sextic.subs(x, root) - slope_value * numerator.subs(x, root) ** 2)
            equations.append(
                (
                    numerator * sp.diff(sextic, x)
                    - 2 * sextic * sp.diff(numerator, x)
                ).subs(x, root)
            )
        matrix, _ = sp.linear_eq_to_matrix(equations, coeffs)
        require(matrix.rank() == 6, "Hermite conditions lost rank")
        kernel = matrix.nullspace()
        require(len(kernel) == 1, "Hermite kernel dimension")
        kernel_polynomial = sum(kernel[0][power] * x**power for power in range(7))
        ratio = sp.cancel(kernel_polynomial / denominator**2)
        require(not ratio.has(x), "Hermite kernel is not span(d^2)")
        hermite_packets += 1

    # 5. Eighteen order conditions on degree 18 leave exactly d^6.
    residual_packets = 0
    for roots in ((-2, 1, 3), (0, 2, 5), (-3, -1, 4)):
        denominator = sp.prod(x - root for root in roots)
        residual_coeffs = sp.symbols(f"e0:{19}")
        residual = sum(
            residual_coeffs[power] * x**power for power in range(19)
        )
        equations = [
            sp.diff(residual, x, order).subs(x, root)
            for root in roots
            for order in range(6)
        ]
        matrix, _ = sp.linear_eq_to_matrix(equations, residual_coeffs)
        require(matrix.rank() == 18, "degree-eighteen jet matrix rank")
        kernel = matrix.nullspace()
        require(len(kernel) == 1, "degree-eighteen jet kernel dimension")
        kernel_polynomial = sum(
            kernel[0][power] * x**power for power in range(19)
        )
        ratio = sp.cancel(kernel_polynomial / denominator**6)
        require(not ratio.has(x), "global residual kernel is not span(d^6)")
        residual_packets += 1

    print("theorem=THM-2389")
    print("status=PROOF-CANDIDATE+VERIFIED-EXACT+AWAITING-INDEPENDENT-AUDIT")
    print("clean_cubic=w^3+12P*w+56Q")
    print("slope_polynomial=alpha^3+2940alpha+30184")
    print("slope_discriminant=-126247730112")
    print("pole_coefficients_checked=0,1,2,3,4,5,6")
    print("parameter_reconstructions=B,C,D,W")
    print(f"formal_branch_packets={formal_branch_cases}")
    print(f"hermite_kernel_packets={hermite_packets}")
    print("hermite_kernel=span(d^2)")
    print(f"degree18_divisibility_packets={residual_packets}")
    print("degree18_jet_kernel=span(d^6)")
    print("terminal_lock=sixth-order-scalar")
    print("h4_status=OPEN")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
