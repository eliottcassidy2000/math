#!/usr/bin/env python3
"""Exact referee for the degree-22 highest-odd componentwise closure.

This companion reconstructs the degree-22 Faber flux/response faces, computes
the complete h=0 support, and checks the two exact-prefix incompatibilities
used on a hypothetical one-infinity physical image component.  It is a
finite exact sidecar to the proof; it does not infer the source-normalization
or rational-primitive steps computationally.
"""

from __future__ import annotations

from math import gcd

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


def local_face(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> tuple[int, sp.Expr]:
    terms = sp.Poly(sp.expand(expression), q, s).terms()
    order = min(sum(monomial) for monomial, _ in terms)
    face = sum(
        coefficient * q**monomial[0] * s**monomial[1]
        for monomial, coefficient in terms
        if sum(monomial) == order
    )
    return order, sp.factor(face)


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    basis = sp.groebner(expressions, *variables, order="grevlex")
    return basis.contains(sp.Integer(1))


def main() -> None:
    d, q, s = sp.symbols("d q s")
    phi22, psi22, response22 = observables(22, d, q, s)

    # Remove only the declared nonzero scalar factors from the top forms.
    phi_inner = sp.factor(phi22 * sp.Rational(512, 33) / q)
    psi_inner = sp.factor(psi22 * sp.Rational(8192, 33))
    response_inner = sp.factor(-response22 * sp.Rational(1024, 33) / q)
    require(
        phi_inner
        == -84 * d**2 * q**4 * s
        + 3 * d * q**6
        + 280 * d * q**2 * s**3
        - 21 * q**4 * s**2
        - 84 * s**5,
        "Phi22 top form changed",
    )
    require(
        psi_inner
        == -224 * d**3 * q**6
        + 3360 * d**2 * q**4 * s**2
        - 336 * d * q**6 * s
        - 3360 * d * q**2 * s**4
        + 3 * q**8
        + 560 * q**4 * s**3
        + 224 * s**6,
        "Psi22 top form changed",
    )
    require(
        response_inner
        == -84 * d**3 * q**4 * s
        + 9 * d**2 * q**6
        + 280 * d**2 * q**2 * s**3
        - 105 * d * q**4 * s**2
        - 84 * d * s**5
        + 3 * q**6 * s
        + 70 * q**2 * s**4,
        "R22 top form changed",
    )
    require(
        sp.gcd(sp.Poly(q * phi_inner, d, q, s), sp.Poly(psi_inner, d, q, s))
        == 1,
        "top fluxes acquired a common h=0 component",
    )

    # q=0 leaves only P_infty=[d:q:s]=[1:0:0].
    require(sp.factor(psi_inner.subs(q, 0)) == 224 * s**6, "q=0 fibre changed")

    # On q=1, divide the nonzero s^2 from Phi and pass to the coarse
    # mu_3-invariants rho=d/s^2 and z=s^3.
    rho, z = sp.symbols("rho z")
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
    test_s = sp.symbols("test_s")
    require(
        sp.expand(
            sp.factor(
                phi_inner.subs({q: 1, d: rho * test_s**2, s: test_s})
                / test_s**2
            ).subs(test_s**3, z)
            - smooth_phi
        )
        == 0,
        "smooth Phi quotient changed",
    )
    require(
        sp.expand(
            sp.factor(
                psi_inner.subs({q: 1, d: rho * test_s**2, s: test_s})
            ).subs(test_s**3, z)
            - smooth_psi
        )
        == 0,
        "smooth Psi quotient changed",
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
        "smooth response quotient changed",
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
        "smooth infinity intersection lost transversality",
    )
    require(
        unit_ideal([smooth_phi, smooth_psi, smooth_response], rho, z),
        "R22 vanished at a smooth infinity point",
    )

    # The P_infty faces and all highest-odd response poles.
    d1_phi_order, d1_phi_face = local_face(phi22.subs(d, 1), q, s)
    d1_psi_order, d1_psi_face = local_face(psi22.subs(d, 1), q, s)
    d1_response_order, d1_response_face = local_face(response22.subs(d, 1), q, s)
    cross_face = q * s * (q**2 - 3 * s**2) * (3 * q**2 - s**2)
    g2_face = (q**2 - s**2) * (q**4 - 14 * q**2 * s**2 + s**4)
    require((d1_phi_order, d1_psi_order, d1_response_order) == (6, 6, 6),
            "P_infty local orders changed")
    require(sp.expand(d1_phi_face + sp.Rational(231, 128) * cross_face) == 0,
            "P_infty Phi face changed")
    require(sp.expand(d1_psi_face + sp.Rational(231, 256) * g2_face) == 0,
            "P_infty Psi face changed")
    require(sp.expand(d1_response_face - sp.Rational(231, 256) * cross_face) == 0,
            "P_infty response face changed")

    pole_rows: list[tuple[int, int, int, int, int]] = []
    for selected in range(21, 0, -2):
        gap = 22 - selected
        branch_gcd = gcd(gap, 6)
        selected_phi, _, selected_response = observables(selected, d, q, s)
        phi_zero = sp.factor(selected_phi.subs({d: 1, q: 0, s: 0}))
        response_zero = sp.factor(selected_response.subs({d: 1, q: 0, s: 0}))
        require(phi_zero != 0 and response_zero != 0, f"odd row {selected} lost transversality")
        restricted_response = sp.factor(
            d1_response_face - (response_zero / phi_zero) * d1_phi_face
        )
        require(
            sp.gcd(
                sp.Poly(restricted_response, q, s), sp.Poly(d1_psi_face, q, s)
            )
            == 1,
            f"odd row {selected} response meets a G2 tangent",
        )
        coarse_branches = 3 * branch_gcd
        h_order = 6 // branch_gcd
        response_pole = (150 - 6 * gap) // branch_gcd
        require(response_pole >= 8, f"odd row {selected} response pole too small")
        require(gap < 24, f"odd row {selected} no longer beats the h^4 L remainder")
        pole_rows.append((selected, gap, branch_gcd, coarse_branches, response_pole))
        require(coarse_branches * h_order == 18, "P_infty h-divisor degree changed")

    # Exact-prefix leading identities at a finite source pole.  At a smooth
    # infinity point, W^2+d=0 and qW-s=0 are incompatible with both top
    # fluxes when q,d are nonzero.
    W = sp.symbols("W")
    smooth_prefix_phi = sp.factor(phi_inner.subs({d: -W**2, s: q * W}))
    smooth_prefix_psi = sp.factor(psi_inner.subs({d: -W**2, s: q * W}))
    require(
        smooth_prefix_phi == -8 * W**2 * q**5 * (56 * W**3 + 3 * q),
        "smooth exact-prefix Phi factor changed",
    )
    require(
        smooth_prefix_psi
        == q**6 * (7168 * W**6 + 896 * W**3 * q + 3 * q**2),
        "smooth exact-prefix Psi factor changed",
    )
    smooth_resultant = sp.factor(
        sp.resultant(
            56 * W**3 + 3 * q,
            7168 * W**6 + 896 * W**3 * q + 3 * q**2,
            q,
        )
    )
    require(smooth_resultant == -76608 * W**6, "smooth prefix resultant changed")

    # At P_infty, W represents the residue omega and q,s represent the first
    # nonzero tangent coefficients q_*,s_*.  The exact-prefix relations
    # omega^2=-1 and s_*=q_*omega miss G2: its face reduces to 32 q_*^6.
    pinfty_prefix = sp.expand(g2_face.subs(s, q * W))
    pinfty_remainder = sp.rem(
        sp.Poly(pinfty_prefix / q**6, W), sp.Poly(W**2 + 1, W)
    ).as_expr()
    require(pinfty_remainder == 32, "P_infty exact-prefix/G2 obstruction changed")

    print("DEGREE22 HIGHEST-ODD COMPONENTWISE CLOSURE")
    print("h0_support=P_infty_plus_five_smooth_coarse_points")
    print("smooth_quotient_p5_degree=5 squarefree=YES unique_rho=YES")
    print("smooth_points=q_nonzero,d_nonzero,transverse,R22_nonzero")
    print("h_divisor_degree=3*g*(6/g)+5=23")
    print("smooth_poles=q_aff:3,R_aff:25")
    print("P_infty_response_pole_min=8")
    print("source_one_pole=>physical_component_has_exactly_one_infinity_point")
    print("pullback_response_order_gt_1=>U_nonconstant=>unique_source_pole_finite")
    print("finite_pole_exact_prefix=h^2*H=W^2+d;h^4*L=qW-s")
    print(f"smooth_prefix_resultant={smooth_resultant}")
    print("P_infty_tangent_after_omega2=-1_and_sstar=qstar*omega=32*qstar^6")
    print("all_odd_response_members=componentwise_empty_including_reducible_nonreduced")
    print("residual=all_even_zero_first_flux_only_not_JC2")
    print("PASS")


if __name__ == "__main__":
    main()
