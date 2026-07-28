#!/usr/bin/env python3
"""Exact local referee for the all-even zero-first-flux degree-22 edge.

The script reconstructs the Faber observables from their Laurent recurrence.
It verifies the two facts used in the componentwise closure argument:

* below slope four, the two degree-six forms at P_infty are coprime; and
* on the all-even lambda=0 family, the homogeneous response combination
  R_25+(d/2)F_23 has local weight at least 28 once q and s have weight four.

The source-normalization, rational-primitive, and projective-regularity steps
are mathematical arguments recorded in the companion reflection/theorem, not
finite claims inferred by this program.
"""

from __future__ import annotations

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


def observables(
    degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s)
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    response = sp.factor(
        4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1]
    )
    return phi, psi, response


def local_order(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), q, s)
    return min(sum(monomial) for monomial, _ in polynomial.terms())


def local_face(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), q, s)
    order = local_order(expression, q, s)
    return sp.factor(
        sum(
            coefficient * q ** monomial[0] * s ** monomial[1]
            for monomial, coefficient in polynomial.terms()
            if sum(monomial) == order
        )
    )


def h_order(expression: sp.Expr, h: sp.Symbol, *others: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), h, *others)
    return min(monomial[0] for monomial, _ in polynomial.terms())


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    return sp.groebner(expressions, *variables, order="grevlex").contains(sp.Integer(1))


def main() -> None:
    h, d, q, s = sp.symbols("h d q s")
    Q, S = sp.symbols("Q S")
    a2, a6, a10, a14, W = sp.symbols("a2 a6 a10 a14 W")

    even_degrees = (2, 6, 10, 14, 22)
    rows = {degree: observables(degree, d, q, s) for degree in even_degrees}
    local_orders = {2: 1, 6: 2, 10: 3, 14: 4, 22: 6}
    gaps = {degree: 22 - degree for degree in even_degrees}

    # Every even row lies on the same slope-four face.  Its response differs
    # from -d Phi/2 by one additional local (q,s)-order.
    expected_differences = {
        2: sp.Integer(0),
        6: -q**3 / 4,
        10: 5 * q**3 * s / 8,
        14: 7 * q**3 * (d * q**2 - 5 * s**2) / 32,
        22: -33
        * q**3
        * (6 * d**2 * q**4 - 84 * d * q**2 * s**2 + 3 * q**4 * s + 70 * s**4)
        / 1024,
    }
    for degree in even_degrees:
        phi, psi, response = rows[degree]
        k = local_orders[degree]
        require(
            local_order(phi.subs(d, 1), q, s) == k,
            f"Phi_{degree} local order changed",
        )
        require(
            local_order(psi.subs(d, 1), q, s) == k,
            f"Psi_{degree} local order changed",
        )
        difference = sp.factor(response + d * phi / 2)
        require(
            sp.expand(difference - expected_differences[degree]) == 0,
            f"D_{degree}=R_{degree}+d Phi_{degree}/2 changed",
        )
        if difference != 0:
            require(
                local_order(difference.subs(d, 1), q, s) >= k + 1,
                f"D_{degree} lost its extra local order",
            )
        require(
            gaps[degree] + 4 * k == 24,
            f"degree {degree} left the common slope-four face",
        )

    # Finite hostile control for the all-degree lemma proved from the formal
    # coefficient recurrence in THM-2752.  This is evidence through m=38,
    # not the source of the quantified proof.
    general_even_control = tuple(range(2, 39, 4))
    for degree in general_even_control:
        phi, psi, response = observables(degree, d, q, s)
        k = (degree + 2) // 4
        phi_local = phi.subs(d, 1)
        psi_local = psi.subs(d, 1)
        require(local_order(phi_local, q, s) == k, f"Phi_{degree} all-m order changed")
        require(local_order(psi_local, q, s) == k, f"Psi_{degree} all-m order changed")
        phi_general_face = local_face(phi_local, q, s)
        psi_general_face = local_face(psi_local, q, s)
        require(
            sp.gcd(
                sp.Poly(phi_general_face, q, s),
                sp.Poly(psi_general_face, q, s),
            ).total_degree()
            == 0,
            f"degree {degree} all-m faces acquired a common direction",
        )
        require(
            sp.Poly(psi_general_face, q, s).coeff_monomial(s**k) != 0,
            f"Psi_{degree} lost pure s^{k}",
        )
        require(
            sp.Poly(phi_general_face, q, s).coeff_monomial(q**k) != 0
            or sp.Poly(psi_general_face, q, s).coeff_monomial(q**k) != 0,
            f"degree {degree} lost pure q^{k}",
        )
        difference = sp.factor(response + d * phi / 2)
        if difference == 0:
            continue
        require(
            local_order(difference.subs(d, 1), q, s) >= k + 1,
            f"degree {degree} lost all-m response order gain",
        )
        require(
            all(
                monomial[1] >= 3
                for monomial, _ in sp.Poly(difference, d, q, s).terms()
            ),
            f"degree {degree} lost q^3 divisibility",
        )

    # If min(v(q),v(s))<4v(h), all lower even rows and W h^24 are
    # strictly higher than the degree-six top faces.  Those top faces have no
    # common projective direction.  Unequal valuations are already killed by
    # the pure q^6 and s^6 terms in Psi; equal valuations reduce to this gcd.
    phi22, psi22, _ = rows[22]
    require(
        sp.gcd(sp.Poly(phi22, d, q, s), sp.Poly(psi22, d, q, s)) == 1,
        "top fluxes acquired a common h=0 component",
    )
    phi_face = local_face(phi22.subs(d, 1), q, s)
    psi_face = local_face(psi22.subs(d, 1), q, s)
    require(local_order(phi22.subs(d, 1), q, s) == 6, "Phi22 face order changed")
    require(local_order(psi22.subs(d, 1), q, s) == 6, "Psi22 face order changed")
    require(
        sp.gcd(sp.Poly(phi_face, q, s), sp.Poly(psi_face, q, s)) == 1,
        "P_infty degree-six faces acquired a common direction",
    )
    psi_poly = sp.Poly(psi_face, q, s)
    require(psi_poly.coeff_monomial(q**6) != 0, "Psi22 lost pure q^6")
    require(psi_poly.coeff_monomial(s**6) != 0, "Psi22 lost pure s^6")

    # Recheck the complete coefficient-independent top boundary rather than
    # importing it through an odd-seed theorem.  Besides P_infty (q=0), the
    # q=1 coarse quotient has exactly five transverse points, and R22 is
    # nonzero at all five.
    _, _, response22 = rows[22]
    require(
        sp.factor(psi22.subs(q, 0)) == sp.Rational(231, 256) * s**6,
        "the q=0 top fibre changed",
    )
    require(
        unit_ideal([phi22.subs(q, 1), psi22.subs(q, 1), s], d, s),
        "a q-nonzero top point with s=0 appeared",
    )
    rho, z, test_s = sp.symbols("rho z test_s")
    phi_inner = sp.factor(phi22 * sp.Rational(512, 33) / q)
    psi_inner = sp.factor(psi22 * sp.Rational(8192, 33))
    response_inner = sp.factor(-response22 * sp.Rational(1024, 33) / q)
    smooth_phi = 3 * rho - 21 + z * (-84 * rho**2 + 280 * rho - 84)
    smooth_psi = (
        3
        + z * (-336 * rho + 560)
        + z**2 * (-224 * rho**3 + 3360 * rho**2 - 3360 * rho + 224)
    )
    smooth_response = (
        3
        + z * (9 * rho**2 - 105 * rho + 70)
        + z**2 * (-84 * rho**3 + 280 * rho**2 - 84 * rho)
    )
    require(
        sp.expand(
            sp.factor(
                phi_inner.subs({q: 1, d: rho * test_s**2, s: test_s})
                / test_s**2
            ).subs(test_s**3, z)
            - smooth_phi
        )
        == 0,
        "five-point Phi quotient changed",
    )
    require(
        sp.expand(
            sp.factor(
                psi_inner.subs({q: 1, d: rho * test_s**2, s: test_s})
            ).subs(test_s**3, z)
            - smooth_psi
        )
        == 0,
        "five-point Psi quotient changed",
    )
    require(
        sp.expand(
            sp.factor(
                response_inner.subs({q: 1, d: rho * test_s**2, s: test_s})
                / test_s
            ).subs(test_s**3, z)
            - smooth_response
        )
        == 0,
        "five-point response quotient changed",
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
    lex_basis = sp.groebner([smooth_phi, smooth_psi], rho, z, order="lex")
    require(
        [sp.factor(poly.as_expr()) for poly in lex_basis.polys]
        == [rho_relation, p5],
        "five-point lexicographic basis changed",
    )
    require(sp.gcd(p5, sp.diff(p5, z)) == 1, "five-point polynomial not squarefree")
    smooth_jacobian = sp.factor(
        sp.diff(smooth_phi, rho) * sp.diff(smooth_psi, z)
        - sp.diff(smooth_phi, z) * sp.diff(smooth_psi, rho)
    )
    require(unit_ideal([smooth_phi, smooth_psi, z], rho, z), "smooth point has s=0")
    require(unit_ideal([smooth_phi, smooth_psi, rho], rho, z), "smooth point has d=0")
    require(
        unit_ideal([smooth_phi, smooth_psi, smooth_jacobian], rho, z),
        "five-point intersection lost transversality",
    )
    require(
        unit_ideal([smooth_phi, smooth_psi, smooth_response], rho, z),
        "R22 vanished at a smooth top point",
    )

    coefficients = {22: 1, 14: a14, 10: a10, 6: a6, 2: a2}
    F23 = sp.factor(
        sum(
            coefficients[degree] * h ** gaps[degree] * rows[degree][0]
            for degree in coefficients
        )
    )
    G24 = sp.factor(
        sum(
            coefficients[degree] * h ** gaps[degree] * rows[degree][1]
            for degree in coefficients
        )
        - W * h**24
    )
    R25 = sp.factor(
        sum(
            coefficients[degree] * h ** gaps[degree] * rows[degree][2]
            for degree in coefficients
        )
    )

    # Work on the d=1 local index cover and encode the valuation conclusion by
    # q=h^4 Q, s=h^4 S.  Both fluxes begin in h^24, while the exact homogeneous
    # cancellation begins in h^28.  On F23=0 it is exactly R25.
    slope_four = {d: 1, q: h**4 * Q, s: h**4 * S}
    scaled_F = sp.expand(F23.subs(slope_four))
    scaled_G = sp.expand(G24.subs(slope_four))
    cancellation = sp.factor(R25 + d * F23 / 2)
    scaled_cancellation = sp.expand(cancellation.subs(slope_four))
    require(h_order(scaled_F, h, Q, S) == 24, "F23 slope-four order changed")
    require(h_order(scaled_G, h, Q, S) == 24, "G24 slope-four order changed")
    require(
        h_order(scaled_cancellation, h, Q, S) >= 28,
        "R25+dF23/2 lost slope-four regularization",
    )
    require(sp.rem(scaled_F, h**24, h) == 0, "F23 is not divisible by h^24")
    require(sp.rem(scaled_G, h**24, h) == 0, "G24 is not divisible by h^24")
    require(
        sp.rem(scaled_cancellation, h**28, h) == 0,
        "response cancellation is not divisible by h^28",
    )

    # The five other h=0 points are the already audited smooth points.  The
    # exact-prefix residues d=-omega^2 and s=q omega cannot satisfy both top
    # fluxes there.
    omega = sp.symbols("omega")
    smooth_phi = sp.factor(phi_inner.subs({d: -omega**2, s: q * omega}))
    smooth_psi = sp.factor(psi_inner.subs({d: -omega**2, s: q * omega}))
    require(
        smooth_phi == -8 * omega**2 * q**5 * (56 * omega**3 + 3 * q),
        "smooth exact-prefix Phi factor changed",
    )
    require(
        smooth_psi
        == q**6 * (7168 * omega**6 + 896 * omega**3 * q + 3 * q**2),
        "smooth exact-prefix Psi factor changed",
    )
    smooth_resultant = sp.factor(
        sp.resultant(
            56 * omega**3 + 3 * q,
            7168 * omega**6 + 896 * omega**3 * q + 3 * q**2,
            q,
        )
    )
    require(
        smooth_resultant == -76608 * omega**6,
        "smooth exact-prefix resultant changed",
    )

    print("DEGREE22 ALL-EVEN ZERO-FIRST-FLUX RESPONSE REGULARIZATION")
    print("even_degrees=2,6,10,14,22")
    print("local_orders=1,2,3,4,6")
    print("common_face_invoice=(22-m)+4*k_m=24")
    print("all_m_4ell_plus_2_control=m_2_through_38:q3,orders,gcd,pure_terms")
    print("P_infty_top_degree6_gcd=1")
    print("universal_h0_boundary=P_infty_plus_five_transverse_coarse_points")
    print("five_smooth_points=R22_nonzero")
    print("unequal_subfour_slopes_killed_by=Psi22_pure_q6_and_s6")
    for degree in even_degrees:
        print(f"D{degree}={sp.sstr(expected_differences[degree])}")
    print("homogeneous_identity=R25+(d/2)F23=sum_even_a_m*h^(22-m)*D_m")
    print("branch_valuation=min(vq,vs)>=4*v(h)")
    print("on_F23_zero:v(R25)>=28*v(h)>25*v(h)")
    print("P_infty_response=regular_and_vanishing")
    print(f"smooth_prefix_resultant={smooth_resultant}")
    print("all_even_zero_flux_edge=componentwise_empty_subject_to_projective_argument")
    print("PASS")


if __name__ == "__main__":
    main()
