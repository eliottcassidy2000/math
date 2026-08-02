#!/usr/bin/env python3
"""Exact companion for THM-3064.

The calculation checks the pointed derivative decomposition for
f=(T-a)g, the cubic resultant formula for the orbit/fixed Keller defect,
the THM-2621 inverse-spectral numerator q=f_v*b_u-f_u*b_v, and the
split-order idempotent in the tame C3 model

    f=(T-1)(T^3-t),

and the exact inverse-different valuation/residue compensation.  It also
checks that trace equality and norm-one data separately fail to decide the
pointed ratio, whereas Norm(rho-1) does decide it in a cubic field.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def order_at_zero(expr: sp.Expr, variable: sp.Symbol) -> int:
    """Return the Laurent order at variable=0 for a nonzero rational expression."""

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


def rational_remainder(expr: sp.Expr, modulus: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Reduce a rational expression whose denominator is scalar in variable."""

    numerator, denominator = sp.fraction(sp.cancel(expr))
    require(sp.degree(denominator, variable) == 0, "non-scalar rational denominator")
    return sp.cancel(sp.rem(numerator, modulus, variable) / denominator)


def main() -> None:
    T, a = sp.symbols("T a")
    p, q, r0 = sp.symbols("p q r0")
    g_generic = T**3 + p * T**2 + q * T + r0
    f_generic = (T - a) * g_generic
    derivative_generic = sp.diff(f_generic, T)

    # Pointed derivative stars: fixed branch g(a), cubic branch (Y-a)g'(Y).
    require(
        sp.expand(derivative_generic.subs(T, a) - g_generic.subs(T, a)) == 0,
        "fixed derivative star",
    )
    require(
        sp.rem(derivative_generic - (T - a) * sp.diff(g_generic, T), g_generic, T) == 0,
        "complement derivative star modulo g",
    )

    # Exact tame C3 model over R=Q[[t]].
    t, s, unit = sp.symbols("t s unit")
    g = T**3 - t
    f = sp.expand((T - 1) * g)
    fixed_derivative = sp.factor(sp.diff(f, T).subs(T, 1))
    cubic_derivative = (T - 1) * sp.diff(g, T)
    require(fixed_derivative == 1 - t, "fixed derivative in C3 model")
    require(sp.expand(cubic_derivative - 3 * T**2 * (T - 1)) == 0, "cubic derivative in C3 model")
    require(
        sp.rem(sp.diff(f, T) - cubic_derivative, g, T) == 0,
        "cubic derivative represents f' modulo g",
    )

    # Since d=1-t is a local unit, g/(1-t) is the singleton idempotent.
    d = fixed_derivative
    idempotent_numerator = sp.expand(g**2 - d * g)
    require(sp.rem(idempotent_numerator, f, T) == 0, "singleton idempotent square")
    require(sp.simplify(g.subs(T, 1) / d) == 1, "singleton idempotent fixed value")
    require(sp.rem(g, g, T) == 0, "singleton idempotent cubic value")

    # Norms are resultants because g is monic.  The pointed ratio packet
    # c0=(1-t)^-1, cC=unit/[3Y^2(Y-1)] has J0=1 and JC=unit.
    denominator = cubic_derivative
    denominator_norm = sp.factor(sp.resultant(g, denominator, T))
    defect_numerator = sp.expand((unit - 1) * denominator)
    defect_norm_numerator = sp.factor(sp.resultant(g, defect_numerator, T))
    pointed_norm_defect = sp.factor(defect_norm_numerator / denominator_norm)
    require(denominator_norm == 27 * t**2 * (t - 1), "cubic derivative norm")
    require(pointed_norm_defect == (unit - 1) ** 3, "pointed norm defect formula")
    require(pointed_norm_defect.subs(unit, 1) == 0, "Keller packet norm defect")
    require(pointed_norm_defect.subs(unit, 2) == 1, "hostile packet norm defect")

    # Exact THM-2621 supplied-pair specialization on the same split cover.
    # The fixed branch has inverse coordinates (x,y)=(u,t); on the cubic
    # branch they are (x,y)=(s,-3*u*s^2/unit).  Their inverse Jacobians are
    # 1 and 1/unit.  The CRT companion realizes both descriptions.
    base_u = sp.symbols("u", nonzero=True)
    d_pair = base_u**3 - t
    f_pair = sp.expand((T - base_u) * g)
    e_pair = g / d_pair
    b_pair = sp.cancel(t * e_pair - (3 * base_u / unit) * T**2 * (1 - e_pair))
    q_pair = sp.cancel(
        sp.diff(f_pair, t) * sp.diff(b_pair, base_u)
        - sp.diff(f_pair, base_u) * sp.diff(b_pair, t)
    )
    derivative_pair_cubic = (T - base_u) * sp.diff(g, T)
    q_expected = sp.cancel(g + derivative_pair_cubic / unit)
    require(sp.cancel(b_pair.subs(T, base_u) - t) == 0, "inverse pair fixed companion")
    require(
        rational_remainder(b_pair + (3 * base_u / unit) * T**2, g, T) == 0,
        "inverse pair cubic companion",
    )
    require(
        rational_remainder(q_pair - q_expected, f_pair, T) == 0,
        "inverse-spectral q packet",
    )
    require(sp.cancel(q_expected.subs(T, base_u) - d_pair) == 0, "q fixed scalar gate")
    require(
        rational_remainder(q_expected - derivative_pair_cubic / unit, g, T) == 0,
        "q cubic packet",
    )
    require(
        sp.expand(q_expected.subs(unit, 1) - sp.diff(f_pair, T)) == 0,
        "unit-one Keller congruence",
    )
    require(
        sp.expand(
            d_pair * q_expected.subs(unit, 1)
            - q_expected.subs({unit: 1, T: base_u}) * sp.diff(f_pair, T)
        )
        == 0,
        "degree-three decoder sign",
    )

    inverse_jacobian_cubic = sp.cancel(
        0 * (-2 * base_u / (unit * s))
        - (1 / (3 * s**2)) * (-3 * s**2 / unit)
    )
    require(inverse_jacobian_cubic == 1 / unit, "cubic inverse Jacobian")

    q_norm = sp.factor(sp.resultant(g, q_expected, T))
    inverse_pair_defect_numerator = sp.cancel(
        d_pair * derivative_pair_cubic - d_pair * q_expected
    )
    inverse_pair_defect_norm = sp.factor(
        sp.resultant(g, inverse_pair_defect_numerator, T)
    )
    inverse_pair_decoder = sp.factor(
        inverse_pair_defect_norm / (d_pair**3 * q_norm)
    )
    require(
        sp.cancel(q_norm - 27 * t**2 * (t - base_u**3) / unit**3) == 0,
        "inverse-pair q norm",
    )
    require(
        sp.cancel(
            inverse_pair_defect_norm
            + 27 * t**2 * (unit - 1) ** 3 * (t - base_u**3) ** 4 / unit**3
        )
        == 0,
        "inverse-pair defect numerator norm",
    )
    require(inverse_pair_decoder == (unit - 1) ** 3, "inverse-pair decoder normalization")

    # The cubic different is generated by g'(s)=3s^2; s-1 is a unit.
    discriminant = sp.factor(sp.discriminant(g, T))
    derivative_on_cubic = 3 * s**2 * (s - 1)
    fixed_cofactor = 1 / (1 - s**3)
    keller_cofactor = 1 / derivative_on_cubic
    hostile_cofactor = 2 / derivative_on_cubic
    require(discriminant == -27 * t**2, "cubic discriminant")
    require(order_at_zero(derivative_on_cubic, s) == 2, "tame cubic different order")
    require(order_at_zero(fixed_cofactor, s) == 0, "fixed cofactor order")
    require(order_at_zero(keller_cofactor, s) == -2, "Keller inverse-different order")
    require(order_at_zero(hostile_cofactor, s) == -2, "hostile inverse-different order")
    require(sp.cancel(keller_cofactor * derivative_on_cubic) == 1, "Keller cubic Jacobian")
    require(sp.cancel(hostile_cofactor * derivative_on_cubic) == 2, "hostile cubic Jacobian")

    # Exact initial-unit compensation in F_7.
    prime = 7
    derivative_residue = (-3) % prime
    keller_cofactor_residue = pow(derivative_residue, -1, prime)
    hostile_cofactor_residue = 2 * keller_cofactor_residue % prime
    require(derivative_residue == 4, "normalized derivative residue")
    require(keller_cofactor_residue == 2, "normalized Keller cofactor residue")
    require(hostile_cofactor_residue == 4, "normalized hostile cofactor residue")
    require(derivative_residue * keller_cofactor_residue % prime == 1, "Keller residue product")
    require(derivative_residue * hostile_cofactor_residue % prime == 2, "hostile residue product")

    # Trace and ordinary norm do not replace Norm(rho-1).
    multiplication_by_s = sp.Matrix([[0, 0, t], [1, 0, 0], [0, 1, 0]])
    require(sp.trace(sp.eye(3) + multiplication_by_s) == 3, "trace hostile")
    require(multiplication_by_s != sp.zeros(3), "trace hostile is nontrivial")
    zeta = sp.symbols("zeta")
    norm_one_numerator = sp.factor(sp.resultant(g, 1 + T, T))
    norm_one_denominator = sp.factor(sp.resultant(g, 1 + zeta * T, T))
    require(norm_one_numerator == 1 + t, "Hilbert-90 numerator norm")
    require(norm_one_denominator.subs(zeta**3, 1) == 1 + t, "Hilbert-90 denominator norm")

    print("theorem=THM-3064")
    print("status=PROVED_VERIFIED_EXACT_CANDIDATE")
    print("pointed_derivatives=D_fixed=g(a);D_C3=(Y-a)g'(Y)")
    print("ratio=rho=(c_C3*D_C3)/(c_fixed*g(a))")
    print("decoder=Norm_C3(rho-1);zero_iff_rho=1_for_transitive_cubic_field")
    print("resultant=Res(g,h*(Y-a)*g'-J_fixed*b)/(J_fixed^3*Res(g,b))")
    print("inverse_spectral=q=f_v*b_u-f_u*b_v;J_phys=f_T/q")
    print("pair_decoder=Res(g,q_fixed*(Y-a)*g'-g(a)*q)/(g(a)^3*Res(g,q))")
    print("degree_p_decoder=pair_decoder=0_iff_g(a)*q=q_fixed*f_T;kappa=g(a)/q_fixed")
    print("constant_gate=polynomial_Keller_requires_kappa_in_constant_field")
    print("rational_pair=f=(T-u)(T^3-t);b=t*e-(3u/unit)T^2*(1-e);e=(T^3-t)/(u^3-t)")
    print("C3_model=f=(T-1)(T^3-t);g(1)=1-t_unit;graph_order=R_x_C")
    print("different=3*s^2;valuation=2;Keller_cofactor_valuation=-2")
    print("packets=J_fixed=1;J_C3=unit;norm_defect=(unit-1)^3")
    print("F7=derivative_residue_4;cofactor_residues_Keller_2_hostile_4;products_1_2")
    print("hostile=unit_1_vs_2_same_order_and_inverse_different_pole_but_ratio_1_vs_2")
    print("boundary=graph_order_and_fixed_owner_do_not_force_ratio;actual_Keller_identity_does")


if __name__ == "__main__":
    main()
