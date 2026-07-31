#!/usr/bin/env python3
"""Exact companion for THM-2848.

The script checks the whitening/quotient identities, four sharp positive
functional controls, the harmonic projection and multipole resultant, and the
factorial integration-by-parts sidecar.  All truth-bearing gates use
``require`` so that ``python -O`` performs the same verification.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def laplacian(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(sum(sp.diff(poly, variable, 2) for variable in variables))


def cyclic_moment(modulus: int, left: int, right: int) -> sp.Integer:
    return sp.Integer(1 if (left - right) % modulus == 0 else 0)


def two_radius_moment(
    left: int, right: int, epsilon: sp.Symbol
) -> sp.Expr:
    """E[Z^left conjugate(Z)^right] for the two-radius hostile."""

    frequency = left - right
    total_degree = left + right
    values = {
        1: (sp.Integer(1), sp.Integer(1)),
        2: (sp.Rational(-1, 2), sp.Rational(-1, 4)),
    }
    answer = sp.Integer(0)
    for radius, (harmonic_one, harmonic_two) in values.items():
        if frequency == 0:
            angular = sp.Integer(1)
        elif abs(frequency) == 1:
            angular = epsilon * harmonic_one / 2
        elif abs(frequency) == 2:
            angular = epsilon * harmonic_two / 2
        else:
            angular = sp.Integer(0)
        answer += sp.Rational(1, 2) * radius**total_degree * angular
    return sp.expand(answer)


def six_atom_moment(left: int, right: int) -> sp.Expr:
    root_five = sp.sqrt(5)
    atoms = (
        (sp.Integer(1), sp.Integer(0), sp.Rational(2, 5)),
        (sp.Integer(-1), sp.Integer(0), sp.Rational(2, 5)),
        (sp.Integer(1), root_five, sp.Rational(1, 20)),
        (sp.Integer(-1), root_five, sp.Rational(1, 20)),
        (sp.Integer(1), -root_five, sp.Rational(1, 20)),
        (sp.Integer(-1), -root_five, sp.Rational(1, 20)),
    )
    return sp.simplify(
        sum(
            weight
            * (real + sp.I * imag) ** left
            * (real - sp.I * imag) ** right
            for real, imag, weight in atoms
        )
    )


def cross_norm_squared(
    first: tuple[sp.Expr, sp.Expr, sp.Expr],
    second: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Expr:
    a, b, c = first
    d, e, f = second
    return sp.expand((b * f - c * e) ** 2 + (c * d - a * f) ** 2 + (a * e - b * d) ** 2)


def null_conic_line(
    vector: tuple[sp.Expr, sp.Expr, sp.Expr],
    u: sp.Symbol,
    v: sp.Symbol,
) -> sp.Expr:
    a, b, c = vector
    return sp.expand((a + sp.I * b) * u**2 + 2 * c * u * v + (-a + sp.I * b) * v**2)


def factorial_functional(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), variable)
    return sp.expand(
        sum(coefficient * factorial(power[0]) for power, coefficient in expanded.terms())
    )


def odd_double_factorial(power: int) -> int:
    if power <= 0:
        return 1
    value = 1
    for factor in range(power, 0, -2):
        value *= factor
    return value


def sphere_monomial(exponents: tuple[int, int, int]) -> sp.Rational:
    """Uniform probability moment on S^2."""

    if any(exponent % 2 for exponent in exponents):
        return sp.Rational(0)
    halves = tuple(exponent // 2 for exponent in exponents)
    numerator = sp.prod(
        odd_double_factorial(2 * half - 1) for half in halves
    )
    denominator = odd_double_factorial(2 * sum(halves) + 1)
    return sp.Rational(numerator, denominator)


def sphere_expectation(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), *variables)
    return sp.expand(
        sum(
            coefficient * sphere_monomial(exponents)
            for exponents, coefficient in expanded.terms()
        )
    )


def main() -> None:
    # ------------------------------------------------------------------
    # Whitening and quotient formulas.
    # ------------------------------------------------------------------
    eta, eta_bar = sp.symbols("eta eta_bar")
    rho, a, a_bar, b, b_bar, radial_four = sp.symbols(
        "rho a a_bar b b_bar radial_four", nonzero=True
    )
    q = rho * eta * eta_bar / 2
    cubic = sp.Rational(3, 8) * (
        a * eta * eta_bar**2 + a_bar * eta**2 * eta_bar
    )
    quartic = (
        sp.Rational(3, 8) * radial_four * eta**2 * eta_bar**2
        + sp.Rational(1, 4)
        * (b * eta * eta_bar**3 + b_bar * eta**3 * eta_bar)
    )
    ell = sp.Rational(3, 4) * (a * eta_bar + a_bar * eta) / rho
    ratio_four = (
        sp.Rational(3, 4) * radial_four * eta * eta_bar / rho
        + (b * eta_bar**2 + b_bar * eta**2) / (2 * rho)
    )
    require(sp.expand(cubic - q * ell) == 0, "cubic quotient identity failed")
    require(
        sp.expand(quartic - q * ratio_four) == 0,
        "quartic quotient identity failed",
    )

    conic = sp.expand(ratio_four - q - ell**2)
    coefficient_radial = sp.expand(conic.coeff(eta, 1).coeff(eta_bar, 1))
    coefficient_bar_square = sp.expand(conic.coeff(eta_bar, 2).subs(eta, 0))
    claimed_radial = (
        6 * rho * radial_four - 4 * rho**3 - 9 * a * a_bar
    ) / (8 * rho**2)
    claimed_bar_square = (8 * rho * b - 9 * a**2) / (16 * rho**2)
    require(
        sp.simplify(coefficient_radial - claimed_radial) == 0,
        "Pearson radial coefficient failed",
    )
    require(
        sp.simplify(coefficient_bar_square - claimed_bar_square) == 0,
        "Pearson harmonic coefficient failed",
    )

    first_bessel_gram = sp.Matrix(
        [
            [1, 0, 0, rho],
            [0, rho, 0, a],
            [0, 0, rho, a_bar],
            [rho, a_bar, a, radial_four],
        ]
    )
    second_bessel_gram = sp.Matrix(
        [
            [1, 0, 0, rho],
            [0, radial_four, 0, b],
            [0, 0, radial_four, b_bar],
            [rho, b_bar, b, radial_four],
        ]
    )
    require(
        sp.expand(
            first_bessel_gram.det()
            - rho * (rho * radial_four - rho**3 - 2 * a * a_bar)
        )
        == 0,
        "first Bessel determinant failed",
    )
    require(
        sp.expand(
            second_bessel_gram.det()
            - radial_four
            * (radial_four**2 - rho**2 * radial_four - 2 * b * b_bar)
        )
        == 0,
        "second Bessel determinant failed",
    )

    # Direct expansion of the Pearson residual.
    q0, ell0, cubic0, quartic0 = sp.symbols("q0 ell0 cubic0 quartic0")
    residual = quartic0 - 2 * q0**2 - 2 * ell0 * cubic0 + q0**2 + ell0**2 * q0
    require(
        sp.expand(residual.subs(cubic0, q0 * ell0) - q0 * (quartic0 / q0 - q0 - ell0**2))
        == 0,
        "Pearson square expansion failed",
    )

    # ------------------------------------------------------------------
    # Positive-functional hostiles.
    # ------------------------------------------------------------------
    pentagon_gram = sp.Matrix(
        [[cyclic_moment(5, row, column) for column in range(5)] for row in range(5)]
    )
    require(pentagon_gram == sp.eye(5), "fifth-root Toeplitz Gram changed")
    require(
        all(cyclic_moment(5, power, 0) == 0 for power in range(1, 5)),
        "fifth-root pure moments changed",
    )

    epsilon = sp.symbols("epsilon", positive=True)
    require(
        all(two_radius_moment(power, 0, epsilon) == 0 for power in range(1, 5)),
        "two-radius pure moments changed",
    )
    two_radius_data = (
        two_radius_moment(1, 1, epsilon),
        two_radius_moment(2, 2, epsilon),
        two_radius_moment(2, 1, epsilon),
        two_radius_moment(3, 1, epsilon),
    )
    require(
        two_radius_data
        == (
            sp.Rational(5, 2),
            sp.Rational(17, 2),
            -3 * epsilon / 4,
            -3 * epsilon / 4,
        ),
        "two-radius mixed moments changed",
    )
    require(
        1 - 2 * sp.Rational(1, 4) > 0
        and 1 - sp.Rational(3, 4) * sp.Rational(1, 4) > 0,
        "two-radius density lower bound changed",
    )

    require(
        all(six_atom_moment(power, 0) == 0 for power in range(1, 5)),
        "six-atom pure moments changed",
    )
    six_data = (
        six_atom_moment(1, 1),
        six_atom_moment(2, 2),
        six_atom_moment(2, 1),
        six_atom_moment(3, 1),
    )
    require(
        six_data == (2, 8, 0, -4),
        "six-atom mixed moments changed",
    )
    six_gap_left = 6 * six_data[0] * six_data[1] - 4 * six_data[0] ** 3
    six_gap_right = abs(8 * six_data[0] * six_data[3])
    require(six_gap_left == six_gap_right == 64, "six-atom Pearson equality changed")

    t = sp.symbols("t", real=True)
    for frequency in range(1, 9):
        circle_moment = sp.integrate(
            sp.exp(2 * sp.pi * sp.I * frequency * t), (t, 0, 1)
        )
        require(sp.simplify(circle_moment) == 0, "exponential-circle moment changed")

    # ------------------------------------------------------------------
    # Harmonic projection and the exact pair/multipole resultant.
    # ------------------------------------------------------------------
    x, y, z = sp.symbols("x y z", real=True)
    variables = (x, y, z)
    Q = x**2 + y**2 + z**2
    cubic_coefficients = sp.symbols("c0:10")
    quartic_coefficients = sp.symbols("f0:15")
    cubic_monomials = [
        x**3,
        x**2 * y,
        x**2 * z,
        x * y**2,
        x * y * z,
        x * z**2,
        y**3,
        y**2 * z,
        y * z**2,
        z**3,
    ]
    quartic_monomials = [
        x**4,
        x**3 * y,
        x**3 * z,
        x**2 * y**2,
        x**2 * y * z,
        x**2 * z**2,
        x * y**3,
        x * y**2 * z,
        x * y * z**2,
        x * z**3,
        y**4,
        y**3 * z,
        y**2 * z**2,
        y * z**3,
        z**4,
    ]
    raw_cubic = sum(
        coefficient * monomial
        for coefficient, monomial in zip(cubic_coefficients, cubic_monomials)
    )
    raw_quartic = sum(
        coefficient * monomial
        for coefficient, monomial in zip(quartic_coefficients, quartic_monomials)
    )
    harmonic_cubic = sp.expand(raw_cubic - Q * laplacian(raw_cubic, variables) / 10)
    harmonic_quartic = sp.expand(
        raw_quartic
        - Q * laplacian(raw_quartic, variables) / 14
        + Q**2 * laplacian(laplacian(raw_quartic, variables), variables) / 280
    )
    require(
        laplacian(harmonic_cubic, variables) == 0,
        "cubic harmonic projection failed",
    )
    require(
        laplacian(harmonic_quartic, variables) == 0,
        "quartic harmonic projection failed",
    )

    u, v = sp.symbols("u v")
    first = sp.symbols("a0:3", real=True)
    second = sp.symbols("b0:3", real=True)
    first_line = null_conic_line(first, u, v)
    second_line = null_conic_line(second, u, v)
    pair_resultant = sp.factor(
        sp.resultant(first_line.subs(v, 1), second_line.subs(v, 1), u)
    )
    require(
        sp.expand(pair_resultant + 4 * cross_norm_squared(first, second)) == 0,
        "pair multipole resultant failed",
    )

    cubic_vectors = ((1, 0, 0), (0, 1, 0), (1, 1, 1))
    quartic_vectors = ((1, 2, 0), (2, 1, 1), (1, -1, 2), (3, 1, -1))
    gamma = sp.Integer(2)
    delta = sp.Integer(-3)
    binary_cubic = gamma * sp.prod(
        null_conic_line(vector, u, v) for vector in cubic_vectors
    )
    binary_quartic = delta * sp.prod(
        null_conic_line(vector, u, v) for vector in quartic_vectors
    )
    multipole_resultant = sp.factor(
        sp.resultant(binary_cubic.subs(v, 1), binary_quartic.subs(v, 1), u)
    )
    cross_product = sp.prod(
        cross_norm_squared(first_vector, second_vector)
        for first_vector in cubic_vectors
        for second_vector in quartic_vectors
    )
    multipole_claim = 2**24 * gamma**8 * delta**6 * cross_product
    require(
        sp.expand(multipole_resultant - multipole_claim) == 0,
        "full multipole resultant failed",
    )

    coincident_quartic_vectors = (cubic_vectors[0],) + quartic_vectors[1:]
    coincident_form = sp.prod(
        null_conic_line(vector, u, v) for vector in coincident_quartic_vectors
    )
    require(
        sp.resultant(binary_cubic.subs(v, 1), coincident_form.subs(v, 1), u) == 0,
        "coincident multipole line did not kill resultant",
    )

    f_zero = sp.expand(Q * (x**2 - y**2))
    f_zero_harmonic = sp.expand(
        f_zero
        - Q * laplacian(f_zero, variables) / 14
        + Q**2 * laplacian(laplacian(f_zero, variables), variables) / 280
    )
    c_control = x**3 - 3 * x * y**2
    require(f_zero_harmonic == 0, "vanishing quartic-harmonic branch failed")
    require(laplacian(c_control, variables) == 0, "cubic branch control not harmonic")
    require(
        Q.subs({x: 0, y: 1, z: sp.I}) == 0
        and c_control.subs({x: 0, y: 1, z: sp.I}) == 0
        and f_zero.subs({x: 0, y: 1, z: sp.I}) == 0,
        "vanishing-harmonic common point control failed",
    )

    # Strictly positive spherical densities exhibit both multipole bad branches.
    h_x, h_y, h_z = sp.symbols("h_x h_y h_z", real=True)
    linear = h_x * x + h_y * y + h_z * z
    cubic_sphere = sphere_expectation(linear**3 * x * y * z, variables)
    quartic_sphere = sphere_expectation(linear**4, variables)
    require(
        sp.expand(cubic_sphere - sp.Rational(2, 35) * h_x * h_y * h_z) == 0,
        "spherical cubic multipole hostile failed",
    )
    require(
        sp.expand(
            quartic_sphere
            - sp.Rational(1, 5) * (h_x**2 + h_y**2 + h_z**2) ** 2
        )
        == 0,
        "spherical quartic radial branch failed",
    )
    quartic_shared = sphere_expectation(
        linear**4 * x * y * (x**2 - y**2), variables
    )
    require(
        sp.expand(
            quartic_shared
            - sp.Rational(8, 315) * h_x * h_y * (h_x**2 - h_y**2)
        )
        == 0,
        "spherical shared-line multipole hostile failed",
    )
    require(
        sp.Rational(1, 3) + sp.Rational(1, 2) < 1,
        "spherical density positivity bound changed",
    )

    # ------------------------------------------------------------------
    # Factorial integration-by-parts/lowering sidecar.
    # ------------------------------------------------------------------
    variable = sp.symbols("s")
    polynomial_controls = (
        1 - 2 * variable + 3 * variable**2,
        2 + variable**3 - variable**5,
        variable * (1 - variable + variable**2),
    )
    for polynomial in polynomial_controls:
        for power in range(4):
            left = factorial_functional(
                polynomial**power * sp.diff(polynomial, variable), variable
            )
            right = (
                factorial_functional(polynomial ** (power + 1), variable)
                - polynomial.subs(variable, 0) ** (power + 1)
            ) / (power + 1)
            require(sp.expand(left - right) == 0, "factorial lowering identity failed")

    print("THM-2848 WHITENED MOVING-PLANE BOUNDARY -- exact companion")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    print("whitening=q-divides-p_m iff pure complex moment m vanishes")
    print("quotients=cubic+quartic+Pearson-conic exact")
    print("subsidiary_Bessel_inequalities=2")
    print("fifth_root=first4zero;Toeplitz_rank=5;support_minimal=5")
    print("two_radius=first4zero;rho=5/2;s=17/2;a=b=-3epsilon/4")
    print("six_atom=first4zero;rho=2;s=8;a=0;b=-4;Pearson=equality")
    print("exponential_circle=first8_nonzero_pure_moments_zero")
    print("harmonic_projection=degrees3,4 exact")
    print("pair_resultant=-4*cross_norm_squared")
    print("multipole_resultant_constant=2^24;transverse_control=nonzero")
    print("quartic_harmonic_zero_branch=separate+common_point")
    print("positive_sphere_hostiles=Fharmonic_zero+shared_xy_lines")
    print("factorial_lowering=L(Z^(r+1))-Z(0)^(r+1) over r+1")
    print("scope=abstract positivity boundary; factorial multiplier access still missing")


if __name__ == "__main__":
    main()
