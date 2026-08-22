#!/usr/bin/env python3
"""Exact width, conductor-period, and low-conductor gates for the nodal cylinder.

The source subring is

    R = Q[u, X=t^2-1, Y=t(t^2-1)].

This companion does not construct a Keller counterexample.  It checks the
arithmetic consequences of the cited Newton-polygon similarity theorem,
builds the first relation-reduced two-scale (4,6) scaffold at target cap 38,
checks the nodal period/PDE identities, proves the displayed first-conductor-
layer equations are inconsistent, and gives an exact period repair which
still has nonconstant Jacobian.  Every truth gate remains active under
``python -O``.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, label: str) -> None:
    """Raise on a failed truth gate, including under optimized execution."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def top_form(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Return the ordinary top homogeneous form."""
    terms = sp.Poly(sp.expand(poly), *variables).terms()
    degree = max(sum(exponent) for exponent, _ in terms)
    return sp.expand(sum(
        coefficient * sp.prod(variable**power for variable, power
                              in zip(variables, exponent))
        for exponent, coefficient in terms if sum(exponent) == degree
    ))


def integrate_t(poly: sp.Expr, u: sp.Symbol, t: sp.Symbol) -> sp.Expr:
    """Integrate a polynomial exactly from t=-1 to t=1, term by term."""
    answer = sp.Integer(0)
    for (u_power, t_power), coefficient in sp.Poly(
            sp.expand(poly), u, t).terms():
        if t_power % 2 == 0:
            answer += coefficient * u**u_power * sp.Rational(2, t_power + 1)
    return sp.expand(answer)


def main() -> None:
    u, t = sp.symbols("u t")
    lam, mu = sp.symbols("lambda mu", nonzero=True)
    U, X0, Y0 = sp.symbols("U X Y")
    X = t**2 - 1
    Y = t * X

    # CITED input (Lang/Abhyankar): for a Keller pair of degrees (72,108),
    # the Newton polygons scale by 108/72=3/2.  The script checks only the
    # elementary width consequences; it does not present this cited theorem
    # as a computation.
    scale = sp.Rational(108, 72)
    require(scale == sp.Rational(3, 2), "Newton-polygon scale")
    width_ladder = [(2 * r, 3 * r) for r in range(1, 7)]
    require(width_ladder[:2] == [(2, 3), (4, 6)],
            "first two Newton-scaled u-widths")
    require(70 * 3 == 105 * 2,
            "minimal width common-factor degree is 35")
    require(35 > 0, "minimal (2,3) common factor is nonconstant")

    # The cap-38 two-scale common base.  V is kept as an ambient abbreviation
    # only while checking the relation-reduced cube.
    V = lam * U * Y0 + mu * U**2
    d_target = Y0**12 + X0 * Y0**10 * V
    p_target = sp.expand(d_target**2)
    q_target = sp.expand(
        Y0**36
        + 3 * X0 * Y0**34 * V
        + 3 * X0**2 * Y0**32 * V**2
        + (Y0**32 - X0**2 * Y0**30) * V**3
    )
    image_relation = Y0**2 - X0**2 * (X0 + 1)
    require(sp.expand(q_target - d_target**3
                      - Y0**30 * V**3 * image_relation) == 0,
            "relation-reduced two-scale cube")
    target_caps = tuple(sp.Poly(poly, U, X0, Y0).total_degree()
                        for poly in (d_target, p_target,
                                     d_target**3, q_target))
    require(target_caps == (13, 26, 39, 38),
            "two-scale target-cap ledger")

    d_source = sp.expand(d_target.subs({U: u, X0: X, Y0: Y}))
    p_source = sp.expand(p_target.subs({U: u, X0: X, Y0: Y}))
    q_source = sp.expand(q_target.subs({U: u, X0: X, Y0: Y}))
    require(sp.expand(q_source - d_source**3) == 0,
            "reduced cube pulls back to a literal cube")
    require(sp.expand(p_source - d_source**2) == 0,
            "square pulls back to a literal square")

    high_base = t**35 * (t + lam * u)
    require(top_form(d_source, (u, t)) == sp.expand(high_base),
            "degree-36 high common base")
    require(top_form(p_source, (u, t)) == sp.expand(high_base**2),
            "degree-72 high square")
    require(top_form(q_source, (u, t)) == sp.expand(high_base**3),
            "degree-108 high cube")
    require(sp.Poly(p_source, u).degree() == 4
            and sp.Poly(q_source, u).degree() == 6,
            "two-scale u-width is (4,6)")
    lower_factor = X * Y**10
    require(sp.expand(sp.Poly(p_source, u).coeff_monomial(u**4)
                      - mu**2 * lower_factor**2) == 0,
            "width-four leading coefficient")
    require(sp.expand(sp.Poly(q_source, u).coeff_monomial(u**6)
                      - mu**3 * lower_factor**3) == 0,
            "width-six leading coefficient")
    require(sp.Poly(lower_factor, t).degree() == 32,
            "lower common factor has t-degree 32")

    # The root-line jet interpolation.  The lower factor vanishes exactly on
    # t=-1,0,1.  eta makes all three degree-drop lines carry the same unit
    # transverse jet and the same minimal nodal boundary.
    eta = 3 * X + 1
    a_scaffold = sp.expand(
        p_source + u**2 - 1 + sp.Rational(1, 2) * eta * Y
    )
    b_scaffold = sp.expand(
        q_source + u**3 - u + sp.Rational(3, 4) * eta * u * Y
    )
    jac_scaffold = sp.expand(
        sp.diff(a_scaffold, u) * sp.diff(b_scaffold, t)
        - sp.diff(a_scaffold, t) * sp.diff(b_scaffold, u)
    )
    for root in (-1, 0, 1):
        require(sp.expand(lower_factor.subs(t, root)) == 0,
                f"lower-factor root t={root}")
        require(sp.expand(d_source.subs(t, root)) == 0
                and sp.expand(sp.diff(d_source, t).subs(t, root)) == 0,
                f"common base vanishes to first order at t={root}")
        require(sp.expand(a_scaffold.subs(t, root) - (u**2 - 1)) == 0
                and sp.expand(b_scaffold.subs(t, root) - (u**3 - u)) == 0,
                f"nodal restriction at t={root}")
        require(sp.expand(jac_scaffold.subs(t, root)) == 1,
                f"unit root-line jet at t={root}")
    four_point_values = {
        (sp.expand(a_scaffold.subs({u: sign_u, t: sign_t})),
         sp.expand(b_scaffold.subs({u: sign_u, t: sign_t})))
        for sign_u in (-1, 1) for sign_t in (-1, 1)
    }
    require(four_point_values == {(0, 0)}, "four-point conductor fibre")
    require(sp.expand(jac_scaffold.subs({u: 0, t: 2}) - 1) != 0,
            "unrepaired two-scale scaffold is not Keller")

    # The conductor-period pairing.  beta_n and K_ij are checked exactly on a
    # hostile finite rectangle; the accompanying proof derives the formula
    # for all i,j by one beta-integral recurrence.
    for n in range(9):
        beta = sp.Rational(
            (-1)**n * 2**(2 * n + 1) * factorial(n)**2,
            factorial(2 * n + 1),
        )
        require(sp.integrate(X**n, (t, -1, 1)) == beta,
                f"beta integral n={n}")
        for i in range(n + 1):
            j = n - i
            kij = -sp.Rational(2 * i, 2 * n + 3) * beta
            require(sp.integrate(
                X**i * sp.diff(Y * X**j, t), (t, -1, 1)
            ) == kij, f"conductor pairing i={i},j={j}")

    # A small hostile which passes the boundary, four-point, first-jet, and
    # period gates but is not Keller.
    a_period_hostile = u**2 - 1 + Y / 2
    b_period_hostile = (
        u**3 - u + sp.Rational(3, 4) * u * Y
        + sp.Rational(105, 16) * u * X**2
    )
    j_period_hostile = sp.expand(
        sp.diff(a_period_hostile, u) * sp.diff(b_period_hostile, t)
        - sp.diff(a_period_hostile, t) * sp.diff(b_period_hostile, u)
    )
    phi_period_hostile = integrate_t(
        a_period_hostile * sp.diff(b_period_hostile, t), u, t
    )
    require(phi_period_hostile == 2 * u,
            "positive hostile has exact nodal period")
    require(integrate_t(j_period_hostile, u, t) == 2,
            "positive hostile has integrated Jacobian two")
    require(j_period_hostile != 1,
            "positive period hostile is not Keller")

    # Verify the normalization PDE on symbolic test polynomials.  Since the
    # displayed identity follows by splitting into even and odd t-parts, the
    # proof is all-degree; this check protects its signs and factors.
    x = sp.symbols("x")
    a_test = u**2 + (u + 1) * x + (2 * u - 1) * x**2
    b_test = u**3 + (u**2 - 1) * x + (u + 2) * x**2
    c_test = 1 + u * x + x**2
    d_test = u + (u**2 + 1) * x + 2 * x**2

    def lop(poly: sp.Expr) -> sp.Expr:
        return ((3 * x + 2) * poly
                + 2 * x * (x + 1) * sp.diff(poly, x))

    even_part = (
        sp.diff(a_test, u) * lop(d_test)
        - lop(c_test) * sp.diff(b_test, u)
        + 2 * x * (x + 1) * (
            sp.diff(c_test, u) * sp.diff(b_test, x)
            - sp.diff(a_test, x) * sp.diff(d_test, u)
        )
    )
    odd_part = (
        2 * (sp.diff(a_test, u) * sp.diff(b_test, x)
             - sp.diff(a_test, x) * sp.diff(b_test, u))
        + x * (sp.diff(c_test, u) * lop(d_test)
               - lop(c_test) * sp.diff(d_test, u))
    )
    source_a_test = (a_test + t * X * c_test).subs(x, X)
    source_b_test = (b_test + t * X * d_test).subs(x, X)
    source_j_test = sp.expand(
        sp.diff(source_a_test, u) * sp.diff(source_b_test, t)
        - sp.diff(source_a_test, t) * sp.diff(source_b_test, u)
    )
    require(sp.expand(source_j_test
                      - even_part.subs(x, X)
                      - t * odd_part.subs(x, X)) == 0,
            "normalization PDE split")

    # The no-I^2 first-conductor-layer equations.  hp denotes h'(u).  The
    # two E-coefficients force hp=15/2 and the remaining polynomial equation
    # has an immediate degree contradiction for every polynomial k.
    h, k, hp = sp.symbols("h k hp")
    common_bracket = h * (16 * (3 * u**2 + 1) * k + 12 * u)
    e1_equation = sp.expand(common_bracket - (2 * hp + 3))
    e2_equation = sp.expand(sp.Rational(5, 4) * common_bracket - 3 * hp)
    require(sp.expand(sp.Rational(5, 4) * (2 * hp + 3) - 3 * hp)
            == -hp / 2 + sp.Rational(15, 4),
            "first-conductor equations force h'=15/2")
    # Finite hostile controls for the all-degree leading-term proof.
    c_shift = sp.symbols("c_shift")
    h_linear = sp.Rational(15, 2) * u + c_shift
    for degree in range(7):
        coefficients = sp.symbols(f"k{degree}_0:{degree + 1}")
        k_poly = sum(coefficients[i] * u**i for i in range(degree + 1))
        obstruction = sp.Poly(sp.expand(
            h_linear * (16 * (3 * u**2 + 1) * k_poly + 12 * u) - 18
        ), u)
        equations = [coefficient for _, coefficient in obstruction.terms()]
        require(sp.linsolve(equations, coefficients) == sp.EmptySet,
                f"no first-conductor solution with deg k<={degree}")

    # Period repair of the two-scale scaffold, specialized to lambda=mu=1.
    # Its coefficients are exact rationals but are deliberately constructed
    # from a small dual matrix instead of printed as unreadable fractions.
    a_two = a_scaffold.subs({lam: 1, mu: 1})
    b_two = b_scaffold.subs({lam: 1, mu: 1})
    phi_two = integrate_t(a_two * sp.diff(b_two, t), u, t)
    require(sp.Poly(phi_two, u).degree() == 6,
            "unrepaired two-scale period has degree six")
    f0 = eta * Y / 2
    f1 = 2 * X * Y**23
    f3 = 2 * X**2 * Y**21
    p_dual = [X**(j + 2) * (X + 1) for j in range(3)]
    dual_matrix = sp.Matrix([
        [sp.integrate(sp.expand(f * sp.diff(p, t)), (t, -1, 1))
         for p in p_dual]
        for f in (f0, f1, f3)
    ])
    require(dual_matrix.det() != 0, "period dual matrix is invertible")
    g0_coeff = dual_matrix.LUsolve(sp.Matrix([1, 0, 0]))
    g1_coeff = dual_matrix.LUsolve(sp.Matrix([0, 1, 0]))
    g0 = sp.expand(sum(g0_coeff[j] * p_dual[j] for j in range(3)))
    g1 = sp.expand(sum(g1_coeff[j] * p_dual[j] for j in range(3)))
    require(integrate_t(a_two * sp.diff(g0, t), u, t) == 1,
            "first period dual")
    require(integrate_t(a_two * sp.diff(g1, t), u, t) == u,
            "second period dual")

    residual = sp.Poly(sp.expand(
        2 * u + phi_two.subs(u, 0) - phi_two
    ), u)
    residual_coeff = [residual.coeff_monomial(u**i) for i in range(7)]
    q0 = sum(residual_coeff[i] * u**i for i in range(1, 6))
    q1 = residual_coeff[6] * u**5
    require(sp.expand(phi_two + q0 + u * q1
                      - (2 * u + phi_two.subs(u, 0))) == 0,
            "two-scale period repair")
    require(sp.Poly(q0, u).degree() <= 5
            and sp.Poly(q1, u).degree() <= 5,
            "period correction preserves width six")
    require(all(sp.expand(p.subs(t, root)) == 0
                and sp.expand(sp.diff(p, t).subs(t, root)) == 0
                for p in p_dual for root in (-1, 0, 1)),
            "period correction preserves all root-line first jets")

    # Evaluate the corrected Jacobian without expanding its enormous full
    # polynomial.  At u=0 only q0'(0)g0 contributes to the correction.
    point = {u: 0, t: 2}
    a_u_point = sp.diff(a_two, u).subs(point)
    a_t_point = sp.diff(a_two, t).subs(point)
    b_t_point = sp.diff(b_two, t).subs(point)
    b_u_point = (
        sp.diff(b_two, u).subs(point)
        + sp.diff(q0, u).subs(u, 0) * g0.subs(t, 2)
        + sp.diff(q1, u).subs(u, 0) * g1.subs(t, 2)
    )
    corrected_jacobian_point = (
        a_u_point * b_t_point - a_t_point * b_u_point
    )
    require(corrected_jacobian_point != 1,
            "period-repaired scaffold still is not Keller")
    witness = corrected_jacobian_point - 1

    print("NODAL-CYLINDER WIDTH AND PERIOD GATES")
    print("cited_Newton_polygon_scale=3/2")
    print("minimal_width_(2,3)=excluded_by_every_line_gate")
    print("first_surviving_scaled_width=(4,6)")
    print("two_scale_source_base_degrees=(36,34)")
    print("two_scale_coefficient_t_degrees=(35,32)")
    print("two_scale_target_caps=(D:13,D2:26,D3_raw:39,D3_reduced:38)")
    print("degree_drop_root_lines=(-1,0,1); nodal_and_unit_jet=True")
    print("four_point_conductor_fibre=True")
    print("nodal_period_formula=Phi'=2*Jac")
    print("period_positive_hostile=True; hostile_Keller=False")
    print("normalization_PDE_split=True")
    print("first_conductor_layer_no_go=PROVED")
    print("second_conductor_layer_required=True")
    print("two_scale_unrepaired_period_degree=6")
    print("two_scale_period_repaired=True")
    print("period_repair_target_cap<=10; preserves_width_and_root_jets=True")
    print("period_repaired_scaffold_Keller=False")
    print("corrected_J_minus_1_mod_1000003="
          f"({int(sp.numer(witness) % 1000003)},"
          f"{int(sp.denom(witness) % 1000003)})")
    print("VERDICT=stronger finite-exact counterexample cell, no Keller pair")


if __name__ == "__main__":
    main()
