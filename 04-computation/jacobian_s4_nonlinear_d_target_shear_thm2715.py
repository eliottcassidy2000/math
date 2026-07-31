#!/usr/bin/env python3
"""Exact companion for THM-2715.

Classify polynomial graph solutions for the nonlinear target family

    A=x^2-2y,  B=y^2-2xz,  d=z,  (U,V)=(A,B+H(d)).

The proof is symbolic plus an all-degree leading-term argument recorded in
the theorem.  This companion checks the integrated PDE, the linear family,
the complete quadratic coefficient comparison, both quadratic sections and
their common triangular target map, and the quartic contradiction.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    x, y, t, c, kappa = sp.symbols("x y t c kappa")
    z = sp.symbols("z")
    alpha, beta, gamma = sp.symbols("alpha beta gamma", nonzero=True)

    # The coordinate change y=t+x^2/2 makes f_x+x f_y the x-derivative
    # at fixed t.  Verify the coordinate Jacobian identities directly for
    # an unrestricted symbolic function F(x,t).
    F = sp.Function("F")(x, t)
    fx_old = sp.diff(F, x) - x * sp.diff(F, t)
    fy_old = sp.diff(F, t)
    L = sp.factor(fx_old + x * fy_old)
    G = sp.factor(x * fx_old + x**2 * fy_old + F - x * (t + x**2 / 2))
    require(L == sp.diff(F, x), "fixed-t derivative identity changed")
    require(
        sp.factor(G - (x * sp.diff(F, x) + F - x * t - x**3 / 2)) == 0,
        "graph Euler row changed",
    )

    # Linear H=beta*z+gamma.  With b=-beta/2, this is the shifted cubic
    # from THM-2709, now checked in the nonlinear-family normalization.
    b = -beta / 2
    f_linear = sp.expand(
        (x - b) * t / 2 + (x - b) * (x**2 + b**2) / 8 + c
    )
    h_linear = beta * z + gamma
    pde_linear = sp.factor(
        x * sp.diff(f_linear, x)
        + f_linear
        - x * t
        - x**3 / 2
        - sp.diff(h_linear, z).subs(z, f_linear) * sp.diff(f_linear, x) / 2
        - c
    )
    require(pde_linear == 0, "linear target-shear graph changed")
    v_linear = sp.factor(
        (t + x**2 / 2) ** 2
        - 2 * x * f_linear
        + h_linear.subs(z, f_linear)
    )
    v_linear_expected = sp.expand(
        (t + b**2 / 2) ** 2 - 2 * c * (x + b) + gamma
    )
    require(
        sp.factor(v_linear - v_linear_expected) == 0,
        "linear triangular target changed",
    )

    # A general quadratic graph ansatz against H=alpha*z^2+beta*z+gamma.
    # Coefficient comparison is complete once degree theory forces deg_x f=2.
    P, Q = sp.symbols("P Q", nonzero=True)
    R = sp.symbols("R")
    K = sp.symbols("K")
    f_quad_ansatz = P * x**2 + Q * x + R
    h_quad = alpha * z**2 + beta * z + gamma
    integrated_residual = sp.Poly(
        sp.expand(
            x * f_quad_ansatz
            - h_quad.subs(z, f_quad_ansatz) / 2
            - (x**2 * t / 2 + x**4 / 8 + c * x + K)
        ),
        x,
    )
    quad_coefficients = tuple(
        sp.factor(integrated_residual.coeff_monomial(x**j))
        for j in range(4, -1, -1)
    )
    expected_quad_coefficients = (
        -alpha * P**2 / 2 - sp.Rational(1, 8),
        P * (1 - alpha * Q),
        Q - alpha * (Q**2 + 2 * P * R) / 2 - beta * P / 2 - t / 2,
        R - alpha * Q * R - beta * Q / 2 - c,
        -alpha * R**2 / 2 - beta * R / 2 - gamma / 2 - K,
    )
    require(
        all(
            sp.factor(left - right) == 0
            for left, right in zip(quad_coefficients, expected_quad_coefficients)
        ),
        "quadratic coefficient system changed",
    )
    r_forced = sp.factor((1 / alpha - beta * P - t) / (2 * alpha * P))
    require(
        sp.factor(
            expected_quad_coefficients[2].subs(
                {Q: 1 / alpha, R: r_forced}
            )
        )
        == 0,
        "quadratic x^2 row did not force R",
    )
    require(
        sp.factor(
            expected_quad_coefficients[3].subs(
                {Q: 1 / alpha, R: r_forced}
            )
            + beta / (2 * alpha)
            + c
        )
        == 0,
        "quadratic coefficient wall changed",
    )

    # Parametrize the two roots by p: alpha=-1/(4p^2), beta=-2alpha*c.
    p = sp.symbols("p", nonzero=True)
    alpha_p = -1 / (4 * p**2)
    beta_p = -2 * alpha_p * c
    f_plus = sp.expand(p * x**2 - 4 * p**2 * x + 2 * p * t + 8 * p**3 + c)
    f_minus = sp.expand(f_plus.subs(p, -p))
    h_p = sp.expand(alpha_p * z**2 + beta_p * z + gamma)

    target_expected = sp.expand(
        -2 * c * x
        - 8 * p**2 * t
        + gamma
        + c**2 / (4 * p**2)
        - 16 * p**4
    )
    for label, section in (("plus", f_plus), ("minus", f_minus)):
        pde = sp.factor(
            x * sp.diff(section, x)
            + section
            - x * t
            - x**3 / 2
            - sp.diff(h_p, z).subs(z, section) * sp.diff(section, x) / 2
            - c
        )
        require(pde == 0, f"quadratic {label} section failed PDE")
        target = sp.factor(
            (t + x**2 / 2) ** 2
            - 2 * x * section
            + h_p.subs(z, section)
        )
        require(
            sp.factor(target - target_expected) == 0,
            f"quadratic {label} target failed triangularization",
        )

    U = -2 * t
    V = target_expected
    require(
        sp.factor(
            sp.det(
                sp.Matrix(
                    [
                        [sp.diff(U, x), sp.diff(U, t)],
                        [sp.diff(V, x), sp.diff(V, t)],
                    ]
                )
            )
            + 4 * c
        )
        == 0,
        "quadratic target Jacobian changed",
    )
    x_inverse = sp.factor(
        (-8 * p**2 * (-U / 2) + gamma + c**2 / (4 * p**2) - 16 * p**4 - V)
        / (2 * c)
    )
    require(sp.factor(x_inverse - x) == 0, "quadratic target inverse changed")

    # Degree four is the only other numerical possibility.  The x^4 row
    # makes P constant and nonzero; the x^3 row then makes R constant.  The
    # x^2 row cannot equal t/2.  Verify the two displayed coefficient rows.
    a4, a3, a2, a1, a0 = sp.symbols("a4 a3 a2 a1 a0", nonzero=True)
    P4, R4 = sp.symbols("P4 R4", nonzero=True)
    f_quartic_ansatz = P4 * x + R4
    h_four = a4 * z**4 + a3 * z**3 + a2 * z**2 + a1 * z + a0
    residual_four = sp.Poly(
        sp.expand(
            x * f_quartic_ansatz
            - h_four.subs(z, f_quartic_ansatz) / 2
            - (x**2 * t / 2 + x**4 / 8 + c * x + K)
        ),
        x,
    )
    require(
        sp.factor(residual_four.coeff_monomial(x**4)
                  + a4 * P4**4 / 2 + sp.Rational(1, 8)) == 0,
        "quartic x^4 row changed",
    )
    require(
        sp.factor(residual_four.coeff_monomial(x**3)
                  + P4**3 * (4 * a4 * R4 + a3) / 2) == 0,
        "quartic x^3 row changed",
    )
    r4_forced = -a3 / (4 * a4)
    x2_after = sp.factor(
        residual_four.coeff_monomial(x**2).subs(R4, r4_forced)
    )
    require(
        sp.diff(x2_after, t) == -sp.Rational(1, 2),
        "quartic x^2 coefficient lost the free t obstruction",
    )

    print("THM-2715 NONLINEAR d-TARGET SHEAR GRAPH CLASSIFICATION")
    print("target=(A,B+H(d)); A=x^2-2y; B=y^2-2xd; d=f(x,y)")
    print("fixed_t_coordinate=t=y-x^2/2")
    print("integrated_equation=x*f-H(f)/2=x^2*t/2+x^4/8+c*x+K(t)")
    print("kappa=-4c_nonzero")
    print("deg_H_le_1=unique_shifted_cubic_triangular")
    print("deg_H_2=possible_iff_beta+2*alpha*c=0")
    print("quadratic_sections=exactly_two:p^2=-1/(4alpha)")
    print("quadratic_f_p=p*x^2-4p^2*x+2p*t+8p^3+c")
    print("quadratic_target_V=-2c*x-8p^2*t+gamma+c^2/(4p^2)-16p^4")
    print("quadratic_target_inverse=polynomial")
    print("deg_H_ge_3=no_polynomial_graph_solution")
    print("all_surviving_maps=triangular_polynomial_automorphisms")
    print("scope=THIS_GRAPH_AND_TARGET_FAMILY_NOT_ARBITRARY_NONLINEAR_TARGET_NOT_JC2")
    print("PASS")


if __name__ == "__main__":
    main()
