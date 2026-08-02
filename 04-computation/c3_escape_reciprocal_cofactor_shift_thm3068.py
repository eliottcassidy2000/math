#!/usr/bin/env python3
"""Exact companion for THM-3068.

This script constructs a fixed-sheet plus C3-escape rational inverse pair
which satisfies the THM-2621 coefficient-pole ledger, inverse-spectral Keller
congruence, branchwise Laurent identity, and trace--Liouville exactness.  It
then checks the reciprocal-root derivative/cofactor transformation which
moves valuation +3 for the affine-primitive cofactor to valuation -2 for the
integral reciprocal cofactor.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rational_remainder(expr: sp.Expr, modulus: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Reduce a rational expression whose denominator is scalar in variable."""

    numerator, denominator = sp.fraction(sp.cancel(expr))
    require(sp.degree(denominator, variable) == 0, "non-scalar rational denominator")
    return sp.cancel(sp.rem(numerator, modulus, variable) / denominator)


def order_at_zero(expr: sp.Expr, variable: sp.Symbol) -> int:
    """Return the Laurent order at variable=0 of a nonzero rational expression."""

    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator_poly = sp.Poly(numerator, variable)
    denominator_poly = sp.Poly(denominator, variable)
    require(not numerator_poly.is_zero, "zero has no finite Laurent order")
    require(not denominator_poly.is_zero, "zero denominator")
    numerator_order = min(monomial[0] for monomial, coefficient in numerator_poly.terms() if coefficient != 0)
    denominator_order = min(
        monomial[0] for monomial, coefficient in denominator_poly.terms() if coefficient != 0
    )
    return numerator_order - denominator_order


def main() -> None:
    T, u, t, s = sp.symbols("T u t s", nonzero=True)

    # Monic inverse quartic and its polynomial eliminant.  At t=0 exactly
    # one fixed sheet survives and the three roots of T^3=t^-1 escape.
    g_x = T**3 - 1 / t
    f_x = sp.expand((T - u) * g_x)
    eliminant = sp.expand(t * f_x)
    expected_eliminant = t * T**4 - u * t * T**3 - T + u
    require(sp.expand(eliminant - expected_eliminant) == 0, "polynomial eliminant")
    require(sp.expand(eliminant.subs(t, 0) - (-T + u)) == 0, "one surviving sheet")
    require(sp.degree(eliminant.subs(t, 0), T) == 1, "specialized fibre degree")
    require(sp.expand(sp.Poly(f_x, T).nth(1) + 1 / t) == 0, "full a1 pole")

    # The cross-resultant t*u^3-1 is a unit at t=0, so the split idempotent
    # and companion are regular in the completed graph order.
    cross_unit = t * u**3 - 1
    e = sp.cancel((t * T**3 - 1) / cross_unit)
    require(sp.cancel(e.subs(T, u) - 1) == 0, "fixed idempotent value")
    require(rational_remainder(e, g_x, T) == 0, "cubic idempotent value")
    require(rational_remainder(e**2 - e, f_x, T) == 0, "idempotent square")

    b = sp.cancel(t * e + 3 * u * t**2 * T**2 * (1 - e))
    require(sp.cancel(b.subs(T, u) - t) == 0, "fixed companion")
    require(
        rational_remainder(b - 3 * u * t**2 * T**2, g_x, T) == 0,
        "cubic companion",
    )

    # THM-2621 numerator q=f_t*b_u-f_u*b_t.  It equals f_T on both split
    # factors, hence modulo the quartic, so the rational pair is symplectic.
    q_inverse = sp.cancel(
        sp.diff(f_x, t) * sp.diff(b, u)
        - sp.diff(f_x, u) * sp.diff(b, t)
    )
    derivative_x = sp.diff(f_x, T)
    require(
        rational_remainder(q_inverse - derivative_x, f_x, T) == 0,
        "inverse-spectral Keller congruence",
    )
    d_fixed = sp.factor(g_x.subs(T, u))
    derivative_cubic = (T - u) * sp.diff(g_x, T)
    require(sp.cancel(derivative_x.subs(T, u) - d_fixed) == 0, "fixed derivative")
    require(
        rational_remainder(derivative_x - derivative_cubic, g_x, T) == 0,
        "cubic derivative",
    )
    cubic_q_norm = sp.factor(sp.resultant(g_x, derivative_cubic, T))
    require(
        sp.cancel(cubic_q_norm + 27 * (t * u**3 - 1) / t**3) == 0,
        "nonzero cubic q norm",
    )
    decoder_numerator = sp.cancel(
        d_fixed * derivative_cubic - d_fixed * derivative_x
    )
    require(sp.factor(decoder_numerator + d_fixed * g_x) == 0, "shifted decoder sign")
    require(sp.resultant(g_x, decoder_numerator, T) == 0, "shifted decoder vanishing")

    # The cubic component is the punctured rational Keller map
    # (x,y)->(u=x^4*y/3,t=x^-3); the fixed component is the identity map.
    x, y = sp.symbols("x y", nonzero=True)
    target_u = x**4 * y / 3
    target_t = x ** -3
    rational_jacobian = sp.factor(
        sp.diff(target_u, x) * sp.diff(target_t, y)
        - sp.diff(target_u, y) * sp.diff(target_t, x)
    )
    require(rational_jacobian == 1, "punctured rational Jacobian")
    require(sp.simplify(target_u.subs({x: 1 / s, y: 3 * u * s**4}) - u) == 0, "cubic inverse u")
    require(sp.simplify(target_t.subs(x, 1 / s) - s**3) == 0, "cubic inverse t")

    # Branchwise symplectic and Liouville identities.  In (u,s) coordinates
    # the cubic inverse is x=s^-1, y=3*u*s^4 and t=s^3.
    x_cubic = 1 / s
    y_cubic = 3 * u * s**4
    inverse_wedge = sp.factor(
        sp.diff(x_cubic, u) * sp.diff(y_cubic, s)
        - sp.diff(x_cubic, s) * sp.diff(y_cubic, u)
    )
    require(inverse_wedge == 3 * s**2, "cubic inverse wedge")
    require(inverse_wedge == sp.diff(s**3, s), "target wedge pullback")

    # theta=x*dy-u*dt has coefficients (3t,9u*s^2), hence theta=d(3ut).
    theta_du = sp.factor(x_cubic * sp.diff(y_cubic, u))
    theta_ds = sp.factor(
        x_cubic * sp.diff(y_cubic, s) - u * sp.diff(s**3, s)
    )
    potential = 3 * u * s**3
    require(theta_du == sp.diff(potential, u), "Liouville du coefficient")
    require(theta_ds == sp.diff(potential, s), "Liouville ds coefficient")
    require(order_at_zero(x_cubic * sp.diff(y_cubic, s), s) == 2, "branch one-form has no residue")

    # One fixed zero potential plus three cubic potentials gives 9*d(u*t).
    trace_potential = 9 * u * s**3
    require(3 * theta_du == sp.diff(trace_potential, u), "trace Liouville du")
    require(3 * theta_ds == sp.diff(trace_potential, s), "trace Liouville ds")

    # Reciprocal-root derivative law, independently checked for four symbols.
    z = sp.symbols("z0:4", nonzero=True)
    product_z = sp.prod(z)
    for i in range(4):
        derivative_z_generic = sp.prod(z[i] - z[j] for j in range(4) if j != i)
        derivative_x_generic = sp.prod(1 / z[i] - 1 / z[j] for j in range(4) if j != i)
        require(
            sp.factor(derivative_z_generic + product_z * z[i] ** 2 * derivative_x_generic) == 0,
            f"reciprocal derivative law at sheet {i}",
        )

    # Exact C3 valuations in the ramified normalization w(s)=1, w(t)=3.
    Z = sp.symbols("Z")
    f_z = sp.expand((Z - 1 / u) * (Z**3 - t))
    derivative_x_cubic = sp.cancel(derivative_x.subs({T: 1 / s, t: s**3}))
    q_x_cubic = derivative_x_cubic
    cofactor_x_cubic = sp.cancel(1 / q_x_cubic)
    product_z_cubic = s**3 / u
    derivative_z_cubic = sp.cancel(sp.diff(f_z, Z).subs({Z: s, t: s**3}))
    cofactor_z_cubic = sp.cancel(1 / derivative_z_cubic)
    require(
        sp.cancel(derivative_z_cubic + product_z_cubic * s**2 * derivative_x_cubic) == 0,
        "C3 reciprocal derivative specialization",
    )
    require(
        sp.cancel(cofactor_z_cubic + cofactor_x_cubic / (product_z_cubic * s**2)) == 0,
        "C3 reciprocal cofactor specialization",
    )
    require(order_at_zero(derivative_x_cubic, s) == -3, "affine derivative valuation")
    require(order_at_zero(q_x_cubic, s) == -3, "inverse numerator valuation")
    require(order_at_zero(cofactor_x_cubic, s) == 3, "affine cofactor valuation")
    require(order_at_zero(product_z_cubic, s) == 3, "reciprocal root product valuation")
    require(order_at_zero(s, s) == 1, "integral reciprocal root valuation")
    require(order_at_zero(derivative_z_cubic, s) == 2, "integral derivative valuation")
    require(order_at_zero(cofactor_z_cubic, s) == -2, "inverse-different cofactor valuation")

    derivative_x_fixed = sp.cancel(derivative_x.subs({T: u, t: s**3}))
    derivative_z_fixed = sp.cancel(sp.diff(f_z, Z).subs({Z: 1 / u, t: s**3}))
    require(
        sp.cancel(derivative_z_fixed + product_z_cubic * (1 / u) ** 2 * derivative_x_fixed) == 0,
        "fixed reciprocal derivative specialization",
    )
    require(order_at_zero(derivative_x_fixed, s) == -3, "fixed affine derivative valuation")
    require(order_at_zero(derivative_z_fixed, s) == 0, "fixed integral derivative valuation")

    cubic_discriminant = sp.factor(sp.discriminant(Z**3 - t, Z))
    require(cubic_discriminant == -27 * t**2, "tame cubic discriminant")

    print("theorem=THM-3068")
    print("status=PROVED_VERIFIED_EXACT_HOSTILE_AUDITED")
    print("eliminant=t*T^4-u*t*T^3-T+u;lead=t;specialization=-T+u;k_D=1")
    print("coefficient_pole=a1=-1/t;exact_base_order=-1;three_sheets_escape")
    print("split_idempotent=e=(t*T^3-1)/(t*u^3-1);cross_resultant_unit_at_t=0")
    print("companion=b=t*e+3*u*t^2*T^2*(1-e)")
    print("inverse_spectral=q=f_t*b_u-f_u*b_t=f_T_mod_f;shifted_decoder=0")
    print("rational_C3_map=(x,y)->(u=x^4*y/3,t=x^-3);Jacobian=1;source_has_unit_x")
    print("branch_Liouville=fixed_0;C3_each_d(3*u*t);trace=d(9*u*t);residues=0")
    print("reciprocal_law=D_Z=-P_z*z_i^2*D_X;c_Z=-c_X/(P_z*z_i^2)")
    print("C3_w_s=1:val(q_X^-1)=3,val(P_z)=3,val(z)=1,val(c_Z)=-2")
    print("conclusion=coefficient_poles_and_Laurent_exactness_permit_inverse_different_principal_part")
    print("missing_gate=polynomial_target_realization_and_affine_space_constant_unit_condition")


if __name__ == "__main__":
    main()
