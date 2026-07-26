#!/usr/bin/env python3
"""Exact companion for THM-2376's degree-eighteen H2 cube-locus closure.

Starting from THM-2357's explicit p3 and q5, this script reconstructs the
quadratic-ring normal form forced by THM-2360 and THM-2371.  It performs
the norm and scalar-trace coefficient elimination, treats every pivot
used in the proof, and certifies the final univariate coprimality both
with SymPy and with a small custom Sylvester determinant modulo 17.

The proof works over characteristic zero.  Modular arithmetic is used
only for the final nonzero-resultant certificate after degree preservation
has been checked.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_numerator(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Poly:
    numerator = sp.together(expression).as_numer_denom()[0]
    content, polynomial = sp.Poly(numerator, *variables).primitive()
    require(content != 0, "primitive numerator unexpectedly vanished")
    return polynomial


def signature(
    polynomial: sp.Poly,
    variables: tuple[sp.Symbol, ...],
) -> tuple[int, ...]:
    expression = polynomial.as_expr()
    return (
        sp.total_degree(expression),
        *(sp.degree(expression, variable) for variable in variables),
        len(polynomial.terms()),
    )


def factor_order(
    polynomial: sp.Poly,
    factor: sp.Poly,
) -> tuple[int, sp.Poly]:
    order = 0
    quotient = polynomial
    while sp.rem(quotient, factor).is_zero:
        quotient = sp.exquo(quotient, factor)
        order += 1
    return order, quotient.primitive()[1]


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    values = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    size = len(values)
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if values[row][column] % prime != 0
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            determinant = -determinant
        pivot_value = values[column][column] % prime
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, size):
            multiplier = values[row][column] * inverse % prime
            if multiplier:
                values[row] = [
                    (
                        values[row][index]
                        - multiplier * values[column][index]
                    )
                    % prime
                    for index in range(size)
                ]
    return determinant % prime


def sylvester_resultant_mod(
    first: sp.Poly,
    second: sp.Poly,
    variable: sp.Symbol,
    prime: int,
) -> tuple[int, int]:
    first_mod = sp.Poly(first, variable, modulus=prime)
    second_mod = sp.Poly(second, variable, modulus=prime)
    first_degree = first_mod.degree()
    second_degree = second_mod.degree()
    first_coefficients = [
        int(coefficient) % prime for coefficient in first_mod.all_coeffs()
    ]
    second_coefficients = [
        int(coefficient) % prime for coefficient in second_mod.all_coeffs()
    ]
    matrix: list[list[int]] = []
    for shift in range(second_degree):
        matrix.append(
            [0] * shift
            + first_coefficients
            + [0] * (second_degree - 1 - shift)
        )
    for shift in range(first_degree):
        matrix.append(
            [0] * shift
            + second_coefficients
            + [0] * (first_degree - 1 - shift)
        )
    custom = determinant_mod(matrix, prime)
    library = int(sp.resultant(first_mod, second_mod, variable)) % prime
    return custom, library


def linear_determinant(
    first: sp.Poly,
    second: sp.Poly,
    variable: sp.Symbol,
    remaining: sp.Symbol,
) -> sp.Poly:
    first_in_variable = sp.Poly(first.as_expr(), variable)
    second_in_variable = sp.Poly(second.as_expr(), variable)
    require(
        first_in_variable.degree() == 1
        and second_in_variable.degree() == 1,
        "linear determinant received a nonlinear polynomial",
    )
    determinant = sp.expand(
        first_in_variable.coeff_monomial(variable)
        * second_in_variable.coeff_monomial(1)
        - second_in_variable.coeff_monomial(variable)
        * first_in_variable.coeff_monomial(1)
    )
    return sp.Poly(determinant, remaining).primitive()[1]


def main() -> None:
    y, bvar, cvar = sp.symbols("y B C")
    avar, root, mvar, nvar, wvar, kvar = sp.symbols(
        "a R m n w k"
    )
    scale_c, scale_lambda = sp.symbols("c lambda")
    ratio, delta = sp.symbols("x delta")

    leading_p = sp.Integer(245)
    leading_q_trace = sp.Integer(3773)
    leading_h2 = sp.Integer(73060029)
    require(
        leading_h2
        == 4 * leading_p**3 + leading_q_trace**2,
        "degree-ten leading identity changed",
    )

    p3 = 35 * (y + 1) * (54 * bvar + 7 * y**2 + 7)
    q5 = 7 * y * (
        1620 * bvar * y**2
        + 1620 * bvar * y
        + 2430 * bvar
        + 26244 * cvar * (y + 1)
        + 77 * (y**4 + y**3 + y**2 + y)
        + 182
    )
    p3_poly = sp.Poly(p3, y)
    q5_poly = sp.Poly(q5, y)
    require(
        p3_poly.degree() == 3
        and p3_poly.LC() == 245
        and p3_poly.coeff_monomial(y**2) == 245
        and q5_poly.degree() == 5
        and q5_poly.LC() == 539
        and q5_poly.coeff_monomial(y**4) == 539
        and q5_poly.coeff_monomial(1) == 0,
        "THM-2357 p3/q5 coefficient atlas changed",
    )

    # Write the monic quadratic h as z^2-Delta.  The chosen lift over y=1
    # is encoded by R=z(1)+t(1), hence Delta=R*(a+2-R).
    z = y + avar / 2
    cover_delta = root * (avar + 2 - root)
    h2_monic = sp.expand(z**2 - cover_delta)

    # ell=(z+t)(z+t-R) is the polynomial-ring image of the selected
    # Laurent factor s(s-r).  Write ell=ell_scalar+ell_t*t.
    ell_scalar = sp.expand(2 * z**2 - cover_delta - root * z)
    ell_t = sp.expand(2 * z - root)
    ell_norm = sp.factor(ell_scalar**2 - h2_monic * ell_t**2)
    require(
        ell_norm == -2 * root * cover_delta * (y - 1),
        "selected root-fibre norm changed",
    )

    # U3 has nonzero constant.  After conjugating if necessary and scaling,
    # s^-2 U3 has this unique leading shape.
    cube_scalar = z**2 + mvar * z + nvar
    cube_t = -z + wvar
    cube_norm_poly = sp.Poly(
        sp.expand(cube_scalar**2 - h2_monic * cube_t**2),
        y,
    )
    require(
        cube_norm_poly.degree() == 3
        and sp.factor(cube_norm_poly.LC() - 2 * (mvar + wvar)) == 0,
        "normalized cubic-factor norm changed",
    )

    w_solution = kvar - mvar
    lambda_solution = 2 * kvar / leading_p
    norm_y2_after_w = sp.factor(
        (
            cube_norm_poly.coeff_monomial(y**2)
            - scale_lambda * p3_poly.coeff_monomial(y**2)
        ).subs(
            {
                wvar: w_solution,
                scale_lambda: lambda_solution,
            }
        )
    )
    n_solution = (
        root**2
        - root * avar
        - 2 * root
        - 3 * avar * kvar
        + kvar**2
        - 2 * kvar * mvar
        + 2 * kvar
    ) / 2
    require(
        sp.factor(norm_y2_after_w - 2 * (nvar - n_solution)) == 0,
        "norm y^2 equation stopped solving for n",
    )

    norm_substitutions = {
        wvar: w_solution,
        scale_lambda: lambda_solution,
        nvar: n_solution,
    }
    norm_y1 = sp.factor(
        (
            cube_norm_poly.coeff_monomial(y)
            - scale_lambda * p3_poly.coeff_monomial(y)
        ).subs(norm_substitutions)
    )
    require(
        sp.factor(sp.diff(norm_y1, bvar) + sp.Rational(108, 7) * kvar)
        == 0,
        "norm y equation stopped having the recorded nonzero B pivot",
    )
    b_solution = sp.factor(sp.solve(norm_y1, bvar)[0])
    expected_b_solution = -sp.Rational(7, 216) * (
        -4 * root**2 * kvar
        + 2 * root**2 * mvar
        + 4 * root * avar * kvar
        - 2 * root * avar * mvar
        + 8 * root * kvar
        - 4 * root * mvar
        + 3 * avar**2 * kvar
        + 6 * avar * kvar * mvar
        - 4 * avar * kvar
        - 2 * kvar**2 * mvar
        + 4 * kvar * mvar**2
        - 4 * kvar * mvar
        + 4 * kvar
    ) / kvar
    require(
        sp.factor(b_solution - expected_b_solution) == 0,
        "norm y equation stopped solving for B with only the k pivot",
    )

    all_norm_substitutions = {
        **norm_substitutions,
        bvar: b_solution,
    }
    norm_constant = (
        cube_norm_poly.coeff_monomial(1)
        - scale_lambda * p3_poly.coeff_monomial(1)
    ).subs(all_norm_substitutions)

    # Scalar and t components of ell*u^3.
    cube_third_scalar = sp.expand(
        cube_scalar**3
        + 3 * cube_scalar * cube_t**2 * h2_monic
    )
    cube_third_t = sp.expand(
        3 * cube_scalar**2 * cube_t
        + cube_t**3 * h2_monic
    )
    product_scalar_poly = sp.Poly(
        sp.expand(
            ell_scalar * cube_third_scalar
            + h2_monic * ell_t * cube_third_t
        ),
        y,
    )
    product_t_poly = sp.Poly(
        sp.expand(
            ell_scalar * cube_third_t
            + ell_t * cube_third_scalar
        ),
        y,
    )
    require(
        product_scalar_poly.degree() == 5
        and product_t_poly.degree() == 4,
        "linear-times-cube bidegrees changed",
    )
    ring_a, ring_b, ring_c, ring_d, ring_h = sp.symbols(
        "ring_a ring_b ring_c ring_d ring_h"
    )
    require(
        sp.expand(
            (ring_a * ring_c + ring_h * ring_b * ring_d) ** 2
            - ring_h * (ring_a * ring_d + ring_b * ring_c) ** 2
            - (ring_a**2 - ring_h * ring_b**2)
            * (ring_c**2 - ring_h * ring_d**2)
        )
        == 0,
        "quadratic-ring norm multiplicativity changed",
    )
    require(
        sp.expand(
            (
                ring_a**3
                + 3 * ring_a * ring_b**2 * ring_h
            )
            ** 2
            - ring_h
            * (
                3 * ring_a**2 * ring_b
                + ring_b**3 * ring_h
            )
            ** 2
            - (ring_a**2 - ring_h * ring_b**2) ** 3
        )
        == 0,
        "quadratic-ring cube norm identity changed",
    )

    leading_product = sp.factor(
        product_scalar_poly.LC().subs(norm_substitutions)
    )
    root_delta = root * cover_delta
    leading_pivot = kvar**3 - root_delta
    require(
        sp.factor(leading_product - 2 * leading_pivot) == 0,
        "scalar trace leading coefficient changed",
    )
    trace_drop = sp.factor(
        (
            product_scalar_poly.coeff_monomial(y**4)
            - product_scalar_poly.coeff_monomial(y**5)
        ).subs(norm_substitutions)
    )
    second_factor = root + 4 * avar - 3 * kvar + 6 * mvar - 4
    require(
        sp.factor(
            trace_drop
            + (kvar**3 + root_delta) * second_factor
        )
        == 0,
        "scalar y^4-y^5 factor split changed",
    )

    c_solution = leading_q_trace / (2 * leading_pivot)
    norm_scale = sp.Poly(
        sp.expand(
            leading_p**3 * leading_pivot**2
            - leading_q_trace**2 * root_delta * kvar**3
        ),
        kvar,
    )
    first_branch = sp.Poly(kvar**3 + root_delta, kvar)
    first_branch_remainder = sp.factor(
        sp.rem(norm_scale, first_branch).as_expr()
    )
    require(
        first_branch_remainder
        == leading_h2 * root_delta**2,
        "norm scale stopped excluding the k^3+R*Delta branch",
    )

    # The first factor is impossible under R*Delta!=0.  On the second
    # factor solve for m and build the four remaining primitive residuals.
    m_solution = (4 - root - 4 * avar + 3 * kvar) / 6
    norm_constant_polynomial = primitive_numerator(
        norm_constant.subs(mvar, m_solution),
        (avar, root, kvar),
    )

    def trace_equation(power: int) -> sp.Expr:
        expression = (
            scale_c * product_scalar_poly.coeff_monomial(y**power)
            - 7 * q5_poly.coeff_monomial(y**power)
        )
        expression = expression.subs(all_norm_substitutions)
        expression = expression.subs(scale_c, c_solution)
        return expression.subs(mvar, m_solution)

    trace_y3_polynomial = primitive_numerator(
        trace_equation(3),
        (avar, root, kvar),
    )
    trace_y0_polynomial = primitive_numerator(
        trace_equation(0),
        (avar, root, kvar),
    )
    trace_y2 = sp.factor(trace_equation(2))
    require(
        sp.factor(sp.diff(trace_y2, cvar) + 49 * 26244) == 0,
        "C pivot in the scalar y^2 equation changed",
    )
    cvar_solution = sp.solve(trace_y2, cvar)[0]
    trace_y1_raw = primitive_numerator(
        trace_equation(1).subs(cvar, cvar_solution),
        (avar, root, kvar),
    )
    safe_leading_factor = sp.Poly(-leading_pivot, avar, root, kvar)
    trace_y1_order, trace_y1_polynomial = factor_order(
        trace_y1_raw,
        safe_leading_factor,
    )
    require(
        trace_y1_order == 1,
        "scalar y^1 equation changed its safe leading-pivot order",
    )

    residuals = [
        norm_constant_polynomial,
        trace_y3_polynomial,
        trace_y0_polynomial,
        trace_y1_polynomial,
    ]
    require(
        [signature(item, (avar, root, kvar)) for item in residuals]
        == [
            (4, 3, 4, 2, 24),
            (6, 3, 6, 5, 39),
            (8, 5, 8, 3, 84),
            (8, 5, 7, 5, 78),
        ],
        "pre-ratio residual signatures changed",
    )

    # Ratio coordinates use only the inherited nonzero divisors R and k:
    # x=k/R and delta=(a+2-R)/R.
    ratio_substitutions = {
        kvar: ratio * root,
        avar: root * (delta + 1) - 2,
    }
    ratio_residuals: list[sp.Poly] = []
    root_orders: list[int] = []
    for polynomial in residuals:
        transformed = sp.Poly(
            sp.expand(polynomial.as_expr().subs(ratio_substitutions)),
            delta,
            root,
            ratio,
        )
        root_order = min(
            monomial[1] for monomial, _ in transformed.terms()
        )
        quotient = sp.Poly(
            sp.cancel(transformed.as_expr() / root**root_order),
            delta,
            root,
            ratio,
        ).primitive()[1]
        root_orders.append(root_order)
        ratio_residuals.append(quotient)
    require(
        root_orders == [3, 4, 3, 4]
        and [
            signature(item, (delta, root, ratio))
            for item in ratio_residuals
        ]
        == [
            (5, 3, 1, 2, 17),
            (8, 3, 2, 5, 24),
            (13, 5, 5, 3, 70),
            (12, 5, 4, 5, 56),
        ],
        "ratio-coordinate residual atlas changed",
    )

    scale_quadratic = sp.Poly(
        125 * delta**2
        - 371 * delta * ratio**3
        + 125 * ratio**6,
        delta,
        ratio,
    )
    require(
        sp.factor(
            norm_scale.as_expr().subs(ratio_substitutions)
            - 7**6 * root**6 * scale_quadratic.as_expr()
        )
        == 0,
        "norm scale stopped inducing the ratio quadratic",
    )
    ratio_norm = sp.Poly(ratio_residuals[0].as_expr(), root)
    pivot_a = sp.factor(ratio_norm.coeff_monomial(root))
    pivot_b = sp.factor(ratio_norm.coeff_monomial(1))
    expected_pivot_a = (
        delta**2 + 5 * delta * ratio + 5 * delta + ratio
    ) * (
        delta * ratio + 2 * delta + 2 * ratio + 1
    )
    expected_pivot_b = -2 * (
        delta + 3 * ratio + 2
    ) * (
        2 * delta * ratio + 3 * delta + ratio
    )
    require(
        ratio_norm.degree() == 1
        and sp.factor(pivot_a - expected_pivot_a) == 0
        and sp.factor(pivot_b - expected_pivot_b) == 0,
        "R-pivot factorization changed",
    )

    # Complementary pivot chart A=0: equation R*A+B=0 also forces B=0.
    pivot_resultant = sp.Poly(
        sp.resultant(pivot_a, pivot_b, delta),
        ratio,
    ).primitive()[1]
    require(
        sp.factor(
            pivot_resultant.as_expr() - ratio * (ratio + 1) ** 8
        )
        == 0,
        "A=B=0 pivot resultant changed",
    )
    require(
        sp.factor(pivot_a.subs(ratio, -1))
        == (delta - 1) ** 2 * (delta + 1)
        and sp.factor(pivot_b.subs(ratio, -1))
        == -2 * (delta - 1) **2
        and scale_quadratic.as_expr().subs(
            {ratio: -1, delta: 1}
        )
        == 621,
        "pivot hostile point or its scale obstruction changed",
    )

    # Main chart A!=0: solve R, then reduce the three trace residuals
    # modulo the monic-up-to-a-constant scale quadratic.
    root_solution = sp.cancel(-pivot_b / pivot_a)
    raw_main_signatures = [
        (12, 7, 9, 52),
        (23, 16, 13, 192),
        (20, 13, 13, 146),
    ]
    remainder_signatures = [
        (24, 1, 24, 38),
        (53, 1, 53, 92),
        (44, 1, 44, 74),
    ]
    expected_ratio_orders = [3, 5, 5]
    main_remainders: list[sp.Poly] = []
    stripped_remainders: list[sp.Poly] = []
    scale_over_fraction_field = sp.Poly(
        scale_quadratic.as_expr(),
        delta,
        domain=sp.QQ.frac_field(ratio),
    )
    for index, residual in enumerate(ratio_residuals[1:]):
        raw_main = primitive_numerator(
            residual.as_expr().subs(root, root_solution),
            (delta, ratio),
        )
        require(
            signature(raw_main, (delta, ratio))
            == raw_main_signatures[index],
            "main-chart substituted numerator signature changed",
        )
        remainder_over_fraction_field = sp.Poly(
            raw_main.as_expr(),
            delta,
            domain=sp.QQ.frac_field(ratio),
        ).rem(scale_over_fraction_field)
        coefficient_denominator = sp.lcm(
            [
                sp.denom(coefficient)
                for coefficient
                in remainder_over_fraction_field.all_coeffs()
            ]
        )
        remainder = primitive_numerator(
            remainder_over_fraction_field.as_expr()
            * coefficient_denominator,
            (delta, ratio),
        )
        require(
            signature(remainder, (delta, ratio))
            == remainder_signatures[index],
            "main-chart scale remainder signature changed",
        )
        ratio_order, stripped = factor_order(
            remainder,
            sp.Poly(ratio, delta, ratio),
        )
        require(
            ratio_order == expected_ratio_orders[index]
            and sp.degree(stripped.as_expr(), delta) == 1,
            "main-chart ratio order or linearity changed",
        )
        main_remainders.append(remainder)
        stripped_remainders.append(stripped)

    determinant_30 = linear_determinant(
        main_remainders[0],
        main_remainders[1],
        delta,
        ratio,
    )
    determinant_31 = linear_determinant(
        main_remainders[0],
        main_remainders[2],
        delta,
        ratio,
    )
    order_30, obstruction_62 = factor_order(
        determinant_30,
        sp.Poly(ratio, ratio),
    )
    order_31, obstruction_53 = factor_order(
        determinant_31,
        sp.Poly(ratio, ratio),
    )
    require(
        order_30 == 12
        and order_31 == 12
        and obstruction_62.degree() == 62
        and obstruction_53.degree() == 53
        and sp.gcd(obstruction_62, obstruction_53).degree() == 0,
        "final characteristic-zero obstruction changed",
    )

    certificate_prime = 17
    obstruction_62_mod = sp.Poly(
        obstruction_62,
        ratio,
        modulus=certificate_prime,
    )
    obstruction_53_mod = sp.Poly(
        obstruction_53,
        ratio,
        modulus=certificate_prime,
    )
    require(
        obstruction_62_mod.degree() == 62
        and obstruction_53_mod.degree() == 53,
        "final modular certificate lost degree",
    )
    custom_resultant, library_resultant = sylvester_resultant_mod(
        obstruction_62,
        obstruction_53,
        ratio,
        certificate_prime,
    )
    require(
        custom_resultant == library_resultant == 11,
        "custom and library final resultants changed",
    )

    # Hostile prime: degree drops and a spurious modular gcd appears.  This
    # check guards against silently replacing the chosen good prime by an
    # invalid one.  Positive controls separately test the resultant code.
    hostile_prime = 11
    hostile_62 = sp.Poly(obstruction_62, ratio, modulus=hostile_prime)
    hostile_53 = sp.Poly(obstruction_53, ratio, modulus=hostile_prime)
    require(
        hostile_62.degree() == 61
        and hostile_53.degree() == 52
        and sp.gcd(hostile_62, hostile_53).degree() == 11,
        "hostile-prime control changed",
    )
    common_control = sp.Poly(
        (ratio**2 + 1) * (ratio + 1),
        ratio,
    )
    common_control_2 = sp.Poly(
        (ratio**2 + 1) * (ratio + 2),
        ratio,
    )
    coprime_control = sp.Poly(ratio**2 + 1, ratio)
    coprime_control_2 = sp.Poly(ratio + 1, ratio)
    require(
        sylvester_resultant_mod(
            common_control,
            common_control_2,
            ratio,
            certificate_prime,
        )
        == (0, 0)
        and sylvester_resultant_mod(
            coprime_control,
            coprime_control_2,
            ratio,
            certificate_prime,
        )
        == (2, 2),
        "modular resultant controls changed",
    )

    print("THM-2376 degree-18 H2 coprime cube-locus exact companion")
    print("p3/q5 degrees and leading coefficients: (3,245), (5,539)")
    print("quadratic cover: h=z^2-Delta; Delta=R*(a+2-R)")
    print("selected lift: Norm(ell)=-2*R*Delta*(y-1)")
    print("normalized cube: u=z^2+m*z+n+(-z+w)*t; k=m+w")
    print("norm pivots: lambda=2*k/245; n and B solved; k!=0")
    print("trace pivot: c=3773/(2*(k^3-R*Delta)); denominator nonzero")
    print(
        "factor split: "
        "(k^3+R*Delta)*(R+4*a-3*k+6*m-4)=0"
    )
    print(
        "norm scale kills first factor by "
        "4*245^3+3773^2=73060029"
    )
    print(
        "ratio scale: "
        "125*delta^2-371*delta*x^3+125*x^6=0"
    )
    print(
        "R pivot: A=(delta^2+5*delta*x+5*delta+x)"
        "*(delta*x+2*delta+2*x+1)"
    )
    print(
        "R constant: "
        "B0=-2*(delta+3*x+2)*(2*delta*x+3*delta+x)"
    )
    print("pivot resultant: x*(x+1)^8; hostile scale value: 621")
    print(
        "main raw signatures: "
        "(12,7,9,52), (23,16,13,192), (20,13,13,146)"
    )
    print(
        "main remainder signatures: "
        "(24,1,24,38), (53,1,53,92), (44,1,44,74)"
    )
    print("ratio orders: 3,5,5; every stripped remainder is delta-linear")
    print("final determinant x-orders: 12,12")
    print("final obstruction degrees: 62,53; gcd over Q: 1")
    print(
        "mod 17 degree-preserving resultant: "
        f"custom={custom_resultant}, library={library_resultant}"
    )
    print("hostile prime 11: degrees 61,52; spurious gcd degree 11")
    print(
        "saturation inventory: "
        "R!=0, Delta!=0, k!=0, k^3-R*Delta!=0; A=0 audited"
    )
    print("VERDICT: the degree-18 H2 coprime cube locus is empty")
    print("SCOPE: H4, degree eighteen as a whole, JC(2), and DC(2) remain open")


if __name__ == "__main__":
    main()
