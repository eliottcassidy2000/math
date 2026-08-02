#!/usr/bin/env python3
"""Exact symbolic companion for THM-3276."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def t_valuation(expression, t):
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, t)
    denominator_poly = sp.Poly(denominator, t)
    require(numerator_poly.as_expr() != 0 and denominator_poly.as_expr() != 0,
            "valuation requested for zero")
    numerator_order = min(monomial[0] for monomial, _ in numerator_poly.terms())
    denominator_order = min(monomial[0] for monomial, _ in denominator_poly.terms())
    return numerator_order - denominator_order


def main():
    x, t, lam = sp.symbols("X t lambda", nonzero=True)

    # Exact 1+3 tame model.  The derivative star D=f_T and the fixed
    # idempotent e0 type the diagonal inverse-different line K*D^(-1).
    g = x**3 - t
    f = sp.expand((x - 1) * g)
    derivative = sp.diff(f, x)
    d0 = sp.factor(g.subs(x, 1))
    e0 = sp.factor(g / d0)
    require(sp.factor(e0.subs(x, 1) - 1) == 0,
            "fixed idempotent does not evaluate to one")
    require(sp.rem(sp.Poly(e0.as_numer_denom()[0], x), sp.Poly(g, x)).as_expr() == 0,
            "fixed idempotent does not vanish on the cubic factor")

    # The determinant-one gauge on inverse-spectral numerators.  It has
    # q0=lambda^3 D0 and qC=lambda^(-1) DC; its inverse cofactor therefore
    # has physical Jacobian values lambda^(-3) and lambda.
    q_lam = sp.factor(
        lam**-1 * derivative + (lam**3 - lam**-1) * d0 * e0
    )
    q0 = sp.factor(q_lam.subs(x, 1))
    require(sp.factor(q0 - lam**3 * d0) == 0,
            "fixed inverse-spectral numerator scale failed")
    cubic_q_difference = sp.rem(
        sp.Poly(sp.together(q_lam - lam**-1 * derivative).as_numer_denom()[0], x),
        sp.Poly(g, x),
    ).as_expr()
    require(sp.expand(cubic_q_difference) == 0,
            "moving inverse-spectral numerator scale failed")

    fixed_jacobian = lam**-3
    moving_jacobian = lam
    pointed_ratio = sp.factor(moving_jacobian / fixed_jacobian)
    determinant_gauge = sp.factor(fixed_jacobian * moving_jacobian**3)
    require(pointed_ratio == lam**4, "pointed ratio is not lambda^4")
    require(determinant_gauge == 1, "four-sheet determinant was not preserved")
    line_equation = sp.factor(lam**3 - lam**-1)
    require(sp.factor(line_equation - (lam**4 - 1) / lam) == 0,
            "diagonal line equation failed")

    # THM-3064's shifted resultant is exactly the scalar equation cutting
    # out the diagonal line on the cubic field.
    numerator = sp.factor(
        sp.resultant(
            g,
            q0 * (x - 1) * sp.diff(g, x) - d0 * q_lam,
            x,
        )
    )
    denominator = sp.factor(d0**3 * sp.resultant(g, q_lam, x))
    shifted_resultant = sp.factor(numerator / denominator)
    require(sp.factor(shifted_resultant - (lam**4 - 1) ** 3) == 0,
            "shifted resultant does not cut out the diagonal line")

    # Arbitrarily deep finite-jet hostiles.  lambda_N=1+t^N makes q_N and
    # D, as well as their inverse cofactors relative to D^(-1), agree modulo
    # t^N.  Yet the exact pointed ratio is not one and Delta has order 3N.
    jet_checks = 0
    for n in range(1, 17):
        lam_n = 1 + t**n
        fixed_scale_n = sp.factor(lam_n**3)
        moving_scale_n = sp.factor(lam_n**-1)
        cofactor_fixed_scale_n = sp.factor(lam_n**-3)
        cofactor_moving_scale_n = lam_n
        rho_n = sp.factor(lam_n**4)
        delta_n = sp.factor((rho_n - 1) ** 3)

        require(t_valuation(fixed_scale_n - 1, t) == n,
                f"fixed q jet order failed at N={n}")
        require(t_valuation(moving_scale_n - 1, t) == n,
                f"moving q jet order failed at N={n}")
        require(t_valuation(cofactor_fixed_scale_n - 1, t) == n,
                f"fixed cofactor jet order failed at N={n}")
        require(t_valuation(cofactor_moving_scale_n - 1, t) == n,
                f"moving cofactor jet order failed at N={n}")
        require(t_valuation(rho_n - 1, t) == n,
                f"pointed ratio order failed at N={n}")
        require(t_valuation(delta_n, t) == 3 * n,
                f"shifted norm order failed at N={n}")
        require(sp.factor(cofactor_fixed_scale_n * cofactor_moving_scale_n**3) == 1,
                f"cofactor determinant gauge failed at N={n}")

        q_n = sp.factor(q_lam.subs(lam, lam_n))
        q_difference = sp.cancel(q_n - derivative)
        for coefficient in sp.Poly(q_difference, x).all_coeffs():
            if coefficient != 0:
                require(t_valuation(coefficient, t) >= n,
                        f"polynomial q jet failed at N={n}")

        specialized_resultant = sp.factor(shifted_resultant.subs(lam, lam_n))
        require(sp.factor(specialized_resultant - delta_n) == 0,
                f"specialized resultant failed at N={n}")
        require(delta_n != 0, f"finite-jet hostile collapsed at N={n}")
        jet_checks += 1

    # One explicit fifth-jet control, reported without series ambiguity.
    lam5 = 1 + t**5
    rho5_minus_one = sp.factor(lam5**4 - 1)
    delta5 = sp.factor(rho5_minus_one**3)
    require(t_valuation(rho5_minus_one, t) == 5, "N=5 ratio control failed")
    require(t_valuation(delta5, t) == 15, "N=5 norm control failed")

    print("THM3276 inverse-different diagonal line and finite-jet hostile exact companion")
    print("typed_line=L_(f,xi)=K*f_T(xi)^(-1) inside inverse_different")
    print("line_incidence: c in L_(f,xi) iff c*f_T(xi) is diagonal")
    print("gauge=(lambda^-3 on fixed,lambda on cubic) determinant=1")
    print("gauge_meets_diagonal_line iff lambda^4=1")
    print("shifted_resultant=(lambda^4-1)^3")
    print(f"finite_jet_hostile_checks={jet_checks} N_range=1..16")
    print("lambda_N=1+t^N: q_N/f_T and c_N/(f_T^-1) equal 1 mod t^N")
    print("lambda_N=1+t^N: v(rho_N-1)=N and v(Delta_N)=3N")
    print("N=5_control: same_through_order_4=True rho_defect_order=5 Delta_order=15")
    print(f"generic_q_lambda={q_lam}")
    print(f"generic_resultant_numerator={numerator}")
    print(f"generic_resultant_denominator={denominator}")


if __name__ == "__main__":
    main()
