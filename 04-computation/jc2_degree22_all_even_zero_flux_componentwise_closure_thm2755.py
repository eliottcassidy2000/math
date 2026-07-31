#!/usr/bin/env python3
"""Exact algebra companion for the all-even, zero-first-flux closure theorem.

This scratch companion verifies the finite algebra used by the proof draft:

* independent recurrence and multinomial constructions of the Faber rows;
* q-divisibility of every even Phi and response row;
* coprimality and complete support of the residual h=0 intersection;
* the local-order table behind min(v(q),v(s)) >= 4 v(h);
* the universal smooth-infinity exact-prefix resultant; and
* the cancellation R_ex = -(Q/2) A20_ex which refutes the tempting
  "simple response pole" shortcut.

It does not mechanize normalization, properness, or the theorem that a
regular function on a complete integral curve is constant.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int, p: sp.Expr, q: sp.Expr, r: sp.Expr, extra: int = 3
) -> list[sp.Expr]:
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def recurrence_observables(
    degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s)
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    response = sp.factor(
        4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1]
    )
    return phi, psi, response


def multinomial_coefficient(
    degree: int, offset: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> sp.Expr:
    """Coefficient of w^(-offset) in P^(degree/4), independently."""
    target = degree + offset
    exponent = sp.Rational(degree, 4)
    total = sp.Integer(0)
    for a in range(target // 2 + 1):
        for b in range(target // 3 + 1):
            remainder = target - 2 * a - 3 * b
            if remainder < 0 or remainder % 4:
                continue
            c = remainder // 4
            count = a + b + c
            falling = sp.prod(exponent - index for index in range(count))
            scalar = falling / (
                factorial(a) * factorial(b) * factorial(c)
            )
            total += scalar * (2 * d) ** a * q**b * (d**2 - s) ** c
    return sp.factor(total)


def multinomial_observables(
    degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    c1, c2, c3 = (
        multinomial_coefficient(degree, offset, d, q, s)
        for offset in (1, 2, 3)
    )
    return tuple(sp.factor(value) for value in (
        4 * c1,
        4 * c2,
        4 * c3 + 2 * d * c1,
    ))


def local_face(
    expression: sp.Expr, q: sp.Symbol, s: sp.Symbol
) -> tuple[int, sp.Expr]:
    terms = sp.Poly(sp.expand(expression), q, s).terms()
    order = min(sum(monomial) for monomial, _ in terms)
    face = sum(
        coefficient * q ** monomial[0] * s ** monomial[1]
        for monomial, coefficient in terms
        if sum(monomial) == order
    )
    return order, sp.factor(face)


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    basis = sp.groebner(expressions, *variables, order="grevlex")
    return basis.contains(sp.Integer(1))


def main() -> None:
    d, q, s = sp.symbols("d q s")
    degrees = (2, 6, 10, 14, 22)
    rows: dict[int, tuple[sp.Expr, sp.Expr, sp.Expr]] = {}
    for degree in degrees:
        recurrence = recurrence_observables(degree, d, q, s)
        independent = multinomial_observables(degree, d, q, s)
        require(
            all(sp.expand(left - right) == 0
                for left, right in zip(recurrence, independent)),
            f"recurrence/multinomial disagreement in degree {degree}",
        )
        rows[degree] = recurrence
        phi, psi, response = recurrence
        require(
            sp.rem(sp.Poly(phi, q), sp.Poly(q, q)) == 0
            and sp.rem(sp.Poly(response, q), sp.Poly(q, q)) == 0,
            f"even row {degree} lost q-divisibility",
        )
        require(
            sp.expand(phi.subs(q, -q) + phi) == 0
            and sp.expand(psi.subs(q, -q) - psi) == 0
            and sp.expand(response.subs(q, -q) + response) == 0,
            f"even row {degree} lost exact parity",
        )

    phi22, psi22, response22 = rows[22]
    residual_top = sp.factor(phi22 / q)
    require(
        sp.gcd(
            sp.Poly(residual_top, d, q, s),
            sp.Poly(psi22, d, q, s),
        ) == 1,
        "residual top forms acquired a common factor",
    )
    require(
        sp.gcd(sp.Poly(q, d, q, s), sp.Poly(psi22, d, q, s)) == 1,
        "vertical q=0 top component became common with G",
    )
    require(
        sp.factor(residual_top.subs(q, 0)) == -sp.Rational(693, 128) * s**5
        and sp.factor(psi22.subs(q, 0)) == sp.Rational(231, 256) * s**6,
        "q=0 infinity fibre changed",
    )

    # Local order table at P_infty on the d=1 index cover.  The first
    # equation is A20=F23/q; the second remains G24.
    gaps = {22: 0, 14: 8, 10: 12, 6: 16, 2: 20}
    order_rows = []
    faces = {}
    for degree in (22, 14, 10, 6, 2):
        phi, psi, response = (
            expression.subs(d, 1) for expression in rows[degree]
        )
        a_order, a_face = local_face(sp.factor(phi / q), q, s)
        g_order, g_face = local_face(psi, q, s)
        r_order, r_face = local_face(response, q, s)
        order_rows.append((degree, gaps[degree], a_order, g_order, r_order))
        faces[degree] = (a_face, g_face, r_face)

    require(
        tuple(order_rows) == (
            (22, 0, 5, 6, 6),
            (14, 8, 3, 4, 4),
            (10, 12, 2, 3, 3),
            (6, 16, 1, 2, 2),
            (2, 20, 0, 1, 1),
        ),
        "all-even local-order table changed",
    )
    a5, g6, _ = faces[22]
    require(
        sp.gcd(sp.Poly(a5, q, s), sp.Poly(g6, q, s)) == 1,
        "residual P_infty faces acquired a common tangent",
    )

    # If b=min(v(q),v(s))<4a with a=v(h), the degree-22 faces strictly
    # precede every lower column.  At b=4a all rows tie at weights 20/24.
    margin_rows = []
    for degree, gap, a_order, g_order, _ in order_rows[1:]:
        a_margin_at_zero = gap
        a_margin_slope = a_order - 5
        g_margin_at_zero = gap
        g_margin_slope = g_order - 6
        require(
            a_margin_at_zero + 4 * a_margin_slope == 0
            and g_margin_at_zero + 4 * g_margin_slope == 0
            and a_margin_at_zero > 0
            and g_margin_at_zero > 0,
            f"degree {degree} no longer has the sharp slope-four tie",
        )
        margin_rows.append((degree, gap, a_order, g_order))
    require(24 - 6 * 4 == 0, "W h^24 lost the slope-four G tie")

    # The five q!=0 coarse infinity points.  Set q=1 and use the free mu_3
    # invariants rho=d/s^2 and z=s^3.
    rho, z = sp.symbols("rho z")
    smooth_a = 3 * rho - 21 + z * (-84 * rho**2 + 280 * rho - 84)
    smooth_g = (
        3
        + z * (-336 * rho + 560)
        + z**2 * (-224 * rho**3 + 3360 * rho**2 - 3360 * rho + 224)
    )
    smooth_r = (
        3
        + z * (9 * rho**2 - 105 * rho + 70)
        + z**2 * (-84 * rho**3 + 280 * rho**2 - 84 * rho)
    )
    p5 = (
        20141047808 * z**5
        - 14386462720 * z**4
        + 1089822720 * z**3
        - 21288960 * z**2
        - 35910 * z
        + 81
    )
    rho_relation = (
        200070 * rho
        + 930804137984 * z**4
        - 669330178048 * z**3
        + 53434945536 * z**2
        - 1142572032 * z
        - 355401
    )
    lex = sp.groebner([smooth_a, smooth_g], rho, z, order="lex")
    require(
        [sp.factor(poly.as_expr()) for poly in lex.polys]
        == [rho_relation, p5]
        and sp.gcd(p5, sp.diff(p5, z)) == 1,
        "five-point residual infinity quotient changed",
    )
    smooth_jacobian = sp.factor(
        sp.diff(smooth_a, rho) * sp.diff(smooth_g, z)
        - sp.diff(smooth_a, z) * sp.diff(smooth_g, rho)
    )
    require(
        unit_ideal([smooth_a, smooth_g, z], rho, z)
        and unit_ideal([smooth_a, smooth_g, rho], rho, z)
        and unit_ideal([smooth_a, smooth_g, smooth_jacobian], rho, z)
        and unit_ideal([smooth_a, smooth_g, smooth_r], rho, z),
        "smooth infinity nonzero/transverse/response certificate changed",
    )

    # Universal exact-prefix obstruction at any of the five smooth points.
    # These are scaled versions of A20 and G24, with q nonzero.
    phi_inner = sp.factor(phi22 * sp.Rational(512, 33) / q)
    psi_inner = sp.factor(psi22 * sp.Rational(8192, 33))
    omega = sp.symbols("omega")
    smooth_prefix_a = sp.factor(
        phi_inner.subs({d: -omega**2, s: q * omega})
    )
    smooth_prefix_g = sp.factor(
        psi_inner.subs({d: -omega**2, s: q * omega})
    )
    require(
        smooth_prefix_a
        == -8 * omega**2 * q**5 * (56 * omega**3 + 3 * q)
        and smooth_prefix_g
        == q**6 * (
            7168 * omega**6 + 896 * omega**3 * q + 3 * q**2
        ),
        "smooth exact-prefix factors changed",
    )
    smooth_resultant = sp.factor(sp.resultant(
        56 * omega**3 + 3 * q,
        7168 * omega**6 + 896 * omega**3 * q + 3 * q**2,
        q,
    ))
    require(
        smooth_resultant == -76608 * omega**6,
        "smooth exact-prefix resultant changed",
    )

    # Exceptional slope-four face, with all lower even columns and W.
    B, C, D, E, W = sp.symbols("B C D E W")
    coefficients = {22: 1, 14: B, 10: C, 6: D, 2: E}
    a_ex = sp.factor(sum(
        coefficients[degree] * faces[degree][0] for degree in coefficients
    ))
    g_ex = sp.factor(sum(
        coefficients[degree] * faces[degree][1] for degree in coefficients
    ) - W)
    r_ex = sp.factor(sum(
        coefficients[degree] * faces[degree][2] for degree in coefficients
    ))
    require(
        sp.expand(r_ex + q * a_ex / 2) == 0,
        "exceptional response cancellation identity changed",
    )

    # A simple rational-primitive pole does not force U constant.
    x = sp.symbols("x")
    source_u = x**2
    source_response = -1 / x
    require(
        sp.factor(source_u * sp.diff(source_response, x)) == 1,
        "simple finite rational-primitive hostile changed",
    )

    print("ALL-EVEN LAMBDA-ZERO EXACT COMPANION")
    print("faber_rows=recurrence_equals_independent_multinomial degrees=2,6,10,14,22")
    print("exact_factorization=F23(lambda=0)=q*A20; every even R_j divisible_by_q")
    print("regular_sequence_faces=gcd(A20_infty,G24_infty)=1;gcd(q,G24_infty)=1")
    print("h0_support=P_infty_plus_five_smooth_coarse_points")
    print("smooth_points=q_nonzero,s_nonzero,d_nonzero,transverse,R22_nonzero")
    print(f"local_order_rows={tuple(order_rows)}")
    print("valuation_gate=b=min(vq,vs)<4a forces coprime A5/G6 faces;therefore b>=4a")
    print("consequence_at_P_infty=v(q/h^3)>=v(h)>0 on every nonvertical branch")
    print(f"smooth_exact_prefix_resultant={smooth_resultant}")
    print("exceptional_face_identity=R_ex=-(Q/2)*A20_ex;simple_pole_claim=REFUTED")
    print("simple_pole_source_hostile=U=x^2,R=-1/x,U*Rprime=1")
    print("proof_sidecar=global_regular_q_on_complete_component_then_q=0_then_third_flux_contradiction")
    print("PASS")


if __name__ == "__main__":
    main()
