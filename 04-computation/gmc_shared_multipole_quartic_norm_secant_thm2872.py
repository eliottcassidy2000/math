#!/usr/bin/env python3
"""Exact companion for THM-2872's shared-multipole reduction.

The script checks the universal binary remainder identities, the positive
quartic norm and its gauge degree, the adjacent-coordinate Euler conic, the
response-secant endpoint and midpoint formulas, and the chart/isotropy/
stability hostiles.  It does not assert the still-open TP3/secant inequality
or exclude a shared cubic--quartic line.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factorial_readout(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), variable)
    return sp.expand(
        sum(
            coefficient * factorial(exponents[0])
            for exponents, coefficient in expanded.terms()
        )
    )


def main() -> None:
    # ------------------------------------------------------------------
    # Binary cubic and quartic remainders.
    # ------------------------------------------------------------------
    z = sp.symbols("z")
    g0, g1, g2 = sp.symbols("g0 g1 g2", nonzero=True)
    t0, t1, t2, t3 = sp.symbols("t0 t1 t2 t3")
    A0, A1, A2, A3, A4 = sp.symbols("A0 A1 A2 A3 A4")
    q_binary = g0 + 2 * g1 * z + g2 * z**2
    cubic = t0 + 3 * t1 * z + 3 * t2 * z**2 + t3 * z**3
    quartic = A0 + 4 * A1 * z + 6 * A2 * z**2 + 4 * A3 * z**3 + A4 * z**4
    field = sp.QQ.frac_field(
        g0,
        g1,
        g2,
        t0,
        t1,
        t2,
        t3,
        A0,
        A1,
        A2,
        A3,
        A4,
    )

    def remainder(poly: sp.Expr) -> sp.Expr:
        return sp.Poly(poly, z, domain=field).rem(
            sp.Poly(q_binary, z, domain=field)
        ).as_expr()

    C0 = 2 * g0 * g1 * t3 - 3 * g0 * g2 * t2 + g2**2 * t0
    C1 = (
        -g0 * g2 * t3
        + 4 * g1**2 * t3
        - 6 * g1 * g2 * t2
        + 3 * g2**2 * t1
    )
    R0 = (
        A0 * g2**3
        - 6 * A2 * g0 * g2**2
        + 8 * A3 * g0 * g1 * g2
        + A4 * g0**2 * g2
        - 4 * A4 * g0 * g1**2
    )
    R1 = 4 * (
        A1 * g2**3
        - 3 * A2 * g1 * g2**2
        - A3 * g0 * g2**2
        + 4 * A3 * g1**2 * g2
        + A4 * g0 * g1 * g2
        - 2 * A4 * g1**3
    )
    require(
        sp.cancel(remainder(cubic) - (C0 + C1 * z) / g2**2) == 0,
        "cubic remainder changed",
    )
    require(
        sp.cancel(remainder(quartic) - (R0 + R1 * z) / g2**3) == 0,
        "quartic remainder changed",
    )

    gram_determinant = g0 * g2 - g1**2
    quartic_norm = g2 * R0**2 - 2 * g1 * R0 * R1 + g0 * R1**2
    square_certificate = (
        (g2 * R0 - g1 * R1) ** 2 + gram_determinant * R1**2
    ) / g2
    require(
        sp.cancel(quartic_norm - square_certificate) == 0,
        "quartic norm square certificate changed",
    )
    require(
        sp.cancel(sp.resultant(q_binary, quartic, z) - quartic_norm / g2**3)
        == 0,
        "quartic resultant normalization changed",
    )

    # At the upper root zeta=(-g1+i sqrt(Dg))/g2, the formal squared
    # remainder is N4/g2^7.  The Hermitian norm of U+zeta V is 2Dg/g2.
    root_modulus = (
        R0**2
        - 2 * g1 * R0 * R1 / g2
        + g0 * R1**2 / g2
    ) / g2**6
    require(
        sp.cancel(root_modulus - quartic_norm / g2**7) == 0,
        "quartic root norm changed",
    )
    hermitian_norm = 2 * gram_determinant / g2
    omega = sp.cancel(root_modulus / hermitian_norm**4)
    require(
        sp.cancel(
            omega
            - quartic_norm / (16 * g2**3 * gram_determinant**4)
        )
        == 0,
        "intrinsic Omega normalization changed",
    )

    scale = sp.symbols("scale", nonzero=True)
    scaled_norm = quartic_norm.subs(
        {
            g0: scale**2 * g0,
            g1: scale**2 * g1,
            g2: scale**2 * g2,
            A0: scale**4 * A0,
            A1: scale**4 * A1,
            A2: scale**4 * A2,
            A3: scale**4 * A3,
            A4: scale**4 * A4,
        },
        simultaneous=True,
    )
    require(
        sp.expand(scaled_norm - scale**22 * quartic_norm) == 0,
        "quartic norm gauge degree changed",
    )

    # ------------------------------------------------------------------
    # Adjacent-coordinate Euler conic.
    # ------------------------------------------------------------------
    a0, p, gap_q, r = sp.symbols("a0 p gap_q r")
    n1, n2, n3 = sp.symbols("n1 n2 n3")
    exponents = (
        a0,
        a0 + p,
        a0 + p + gap_q,
        a0 + p + gap_q + r,
    )
    B1 = sp.Matrix((-1, 1, 0, 0))
    B2 = sp.Matrix((0, -1, 1, 0))
    B3 = sp.Matrix((0, 0, -1, 1))
    U = n3 * B1 - n1 * B3
    V = n3 * B2 - n2 * B3
    euler = sp.diag(*exponents)
    real_determinant = sp.det(
        sp.Matrix.hstack(U, V, euler * U, euler * V)
    )
    euler_conic = (
        p * r * n2 * (n1 + n2 + n3)
        - gap_q * (p + gap_q + r) * n1 * n3
    )
    require(
        sp.expand(real_determinant + n3**2 * euler_conic) == 0,
        "Euler adjacent-normal determinant changed",
    )
    require(
        sp.expand(euler_conic.subs({p: 1, gap_q: 1, r: 1, n1: 1, n2: 1, n3: 1}))
        == 0,
        "Euler-normalization hostile changed",
    )

    delta_h = 4 * gram_determinant * n3**2 * euler_conic / g2**2
    euler_normalized = sp.cancel(root_modulus / delta_h**2)
    require(
        sp.cancel(
            euler_normalized
            - quartic_norm
            / (
                16
                * g2**3
                * gram_determinant**2
                * n3**4
                * euler_conic**2
            )
        )
        == 0,
        "Euler-transverse quartic normalization changed",
    )

    # ------------------------------------------------------------------
    # Cubic transport reversals and quartic endpoint/midpoint holonomy.
    # ------------------------------------------------------------------
    cubic_t1 = (
        2 * g1 * t0 / g0 + g0 * t3 / g2
    ) / 3
    cubic_t2 = (
        g2 * t0 / g0 + 2 * g1 * t3 / g2
    ) / 3
    delta2_U = sp.cancel(cubic_t1 / g1 - 2 * t0 / (3 * g0))
    delta2_V = sp.cancel(cubic_t2 / g1 - 2 * t3 / (3 * g2))
    require(
        sp.cancel(delta2_U - t3 * g0 / (3 * g1 * g2)) == 0,
        "first cubic transport reversal changed",
    )
    require(
        sp.cancel(delta2_V - t0 * g2 / (3 * g1 * g0)) == 0,
        "second cubic transport reversal changed",
    )

    delta3_U = A1 / g1 - A0 / (2 * g0)
    delta3_V = A3 / g1 - A4 / (2 * g2)
    kappa_U = sp.cancel(delta3_U / delta2_U)
    kappa_V = sp.cancel(delta3_V / delta2_V)
    alpha = t0 / g0
    beta = t3 / g2
    left_r1 = (2 * A1 * g0 - A0 * g1) / g0**2
    right_r1 = (2 * A3 * g2 - A4 * g1) / g2**2
    require(
        sp.cancel(beta * kappa_U - sp.Rational(3, 2) * left_r1) == 0,
        "left endpoint secant changed",
    )
    require(
        sp.cancel(alpha * kappa_V - sp.Rational(3, 2) * right_r1) == 0,
        "right endpoint secant changed",
    )

    r0_symbol, r1_symbol, r2_symbol = sp.symbols("r0_symbol r1_symbol r2_symbol")
    endpoint_substitution = {
        A0: g0 * r0_symbol,
        A1: (g1 * r0_symbol + g0 * r1_symbol) / 2,
        A3: (g2 * r1_symbol + g1 * r2_symbol) / 2,
        A4: g2 * r2_symbol,
    }
    midpoint_defect = sp.cancel(
        quartic.subs(z, 1) / q_binary.subs(z, 1)
        - A0 / g0
        - A4 / g2
        - 2 * left_r1
    )
    missing_middle = sp.expand(
        quartic
        - q_binary
        * (r0_symbol + 2 * r1_symbol * z + r2_symbol * z**2)
    ).subs(endpoint_substitution, simultaneous=True)
    require(
        sp.Poly(sp.cancel(missing_middle), z).terms()
        == [((2,), sp.Poly(sp.cancel(missing_middle), z).coeff_monomial(z**2))],
        "endpoint substitution left more than the middle coefficient",
    )
    require(
        sp.cancel(
            midpoint_defect.subs(endpoint_substitution, simultaneous=True)
            - missing_middle.subs(z, 1) / q_binary.subs(z, 1)
        )
        == 0,
        "midpoint/Jensen defect changed",
    )

    # ------------------------------------------------------------------
    # Chart and complex-isotropy hostiles.
    # ------------------------------------------------------------------
    variable = sp.symbols("s")
    f_slots = tuple(
        variable**index / sp.factorial(index) for index in range(4)
    )
    blocks = tuple(
        f_slots[index] - f_slots[index - 1] for index in range(1, 4)
    )
    complex_isotropic = blocks[0] + (-1 + sp.I) * blocks[1] / 2
    require(
        factorial_readout(complex_isotropic**2, variable) == 0,
        "complex g2-isotropy hostile changed",
    )
    chart_boundary_U = -blocks[2]
    chart_boundary_V = -blocks[2]
    require(
        sp.expand(chart_boundary_U - chart_boundary_V) == 0,
        "n3=0 false-chart hostile changed",
    )
    valid_chart_U = blocks[1] - blocks[0]
    valid_chart_V = blocks[2]
    require(
        sp.Poly(valid_chart_U, variable) != sp.Poly(valid_chart_V, variable),
        "cyclic repair chart became dependent",
    )
    degenerate_R0, degenerate_R1 = sp.symbols(
        "degenerate_R0 degenerate_R1"
    )
    degenerate_norm = (
        degenerate_R0**2
        - 2 * degenerate_R0 * degenerate_R1
        + degenerate_R1**2
    )
    require(
        degenerate_norm.subs(
            {degenerate_R0: 1, degenerate_R1: 1}
        )
        == 0,
        "degenerate Gram hostile did not kill the norm",
    )
    require(
        sp.rem(1 + z, (1 + z) ** 2, domain=sp.QQ) != 0,
        "degenerate Gram hostile lost its nonzero remainder",
    )
    require(
        factorial_readout(blocks[0] * sp.diff(blocks[0], variable), variable)
        == 0,
        "d0 boundary response hostile changed",
    )

    # ------------------------------------------------------------------
    # The natural moment-Jacobian is neither stable nor Lorentzian.
    # ------------------------------------------------------------------
    u, v, w = sp.symbols("u v w")
    Z = u * blocks[0] + v * blocks[1] + w * blocks[2]
    moment_jacobian = sp.expand(
        sp.det(
            sp.Matrix(
                [
                    [
                        factorial_readout(blocks[column] * Z**power, variable)
                        for column in range(3)
                    ]
                    for power in range(1, 4)
                ]
            )
        )
    )
    content, primitive_jacobian = sp.primitive(moment_jacobian)
    primitive_poly = sp.Poly(primitive_jacobian, u, v, w)
    require(content == 12, "moment-Jacobian content changed")
    require(
        len(primitive_poly.terms()) == 28
        and all(coefficient > 0 for _, coefficient in primitive_poly.terms()),
        "moment-Jacobian positive support changed",
    )
    specialization = sp.expand(primitive_jacobian.subs({u: 1, w: 1}))
    expected_specialization = (
        217 * v**6
        + 4802 * v**5
        + 43835 * v**4
        + 210800 * v**3
        + 560815 * v**2
        + 779518 * v
        + 441213
    )
    require(
        specialization == expected_specialization,
        "stability hostile specialization changed",
    )
    hostile_discriminant = sp.discriminant(specialization, v)
    expected_discriminant = (
        -2**21
        * 3**6
        * 5**10
        * 7**2
        * 2281
        * 7508739397
    )
    require(
        hostile_discriminant == expected_discriminant < 0,
        "stability hostile discriminant changed",
    )
    fourth_derivative = sp.diff(primitive_jacobian, u, 1, w, 3)
    hostile_hessian = sp.hessian(fourth_derivative, (u, v, w))
    require(
        hostile_hessian.det() == -43366523385600000
        and sp.trace(hostile_hessian) == 38634480 > 0,
        "Lorentzian hostile Hessian changed",
    )

    print("THM-2872 SHARED MULTIPOLE QUARTIC NORM AND SECANT REDUCTION")
    print("status=PROVED-ALGEBRAIC-REDUCTION+VERIFIED-EXACT")
    print("binary_cubic_remainder=(C0+C1*z)/g2^2")
    print("binary_quartic_remainder=(R0+R1*z)/g2^3")
    print("N4=g2*R0^2-2*g1*R0*R1+g0*R1^2;positive_Gram_norm=True")
    print("Omega=N4/(16*g2^3*Dg^4);normal_degree=22/22")
    print("Euler_conic=p*r*n2*(n1+n2+n3)-q*(p+q+r)*n1*n3")
    print("real_Euler_determinant=-n3^2*Euler_conic")
    print("cubic_transport_reversals=2;strict_on_positive_boundary_free_cones")
    print("quartic_endpoint_holonomy=beta*kappa_U-alpha*kappa_V")
    print("remaining_quartic_condition=one_midpoint_Jensen_defect")
    print("cyclic_charts=3;chart_boundary+complex_isotropy+Euler_hostiles=PASS")
    print("degenerate_Gram+zero_cross+bottom_d0_hostiles=PASS")
    print("moment_Jacobian_terms=28;all_positive=True;stable=False;Lorentzian=False")
    print(f"stability_hostile_discriminant={hostile_discriminant}")
    print(f"Lorentzian_hostile_Hessian_det={hostile_hessian.det()}")
    print("scope=exact reduction only; shared-line branch and TP3/secant sign remain open")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
