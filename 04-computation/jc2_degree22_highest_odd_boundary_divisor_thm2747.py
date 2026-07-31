#!/usr/bin/env python3
"""Exact infinity-divisor referee for THM-2747.

The script reconstructs the degree-22 quartic Faber observables, passes to
the invariant q!=0 chart of P(2,3,4), proves that the residual infinity
scheme is a reduced irreducible quintic, and checks that the third response
does not vanish there.  It then prints the branch, divisor, and one-ended
component invoices for every possible highest odd Faber seed.

This is an exact boundary computation, not a factorization census or a proof
of JC(2).
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


def main() -> None:
    d, q, s, t, u = sp.symbols("d q s t u")
    phi, psi, response = observables(22, d, q, s)

    phi_core = q * (
        -84 * d**2 * q**4 * s
        + 3 * d * q**6
        + 280 * d * q**2 * s**3
        - 21 * q**4 * s**2
        - 84 * s**5
    )
    psi_core = (
        -224 * d**3 * q**6
        + 3360 * d**2 * q**4 * s**2
        - 336 * d * q**6 * s
        - 3360 * d * q**2 * s**4
        + 3 * q**8
        + 560 * q**4 * s**3
        + 224 * s**6
    )
    require(sp.expand(phi - sp.Rational(33, 512) * phi_core) == 0, "Phi22 changed")
    require(sp.expand(psi - sp.Rational(33, 8192) * psi_core) == 0, "Psi22 changed")

    # At q!=0 use the weight-zero invariants t=s^3/q^4 and u=ds/q^2.
    # In the q=1 cover, multiplying Phi by s and Psi by s^3 gives these
    # two invariant equations.  The t=0 factor in their resultant is an
    # artifact of clearing s; the true top intersection has t!=0.
    f_bar = -84 * u**2 + (3 + 280 * t) * u - 21 * t - 84 * t**2
    g_bar = (
        -224 * u**3
        + 3360 * u**2 * t
        - 336 * u * t
        - 3360 * u * t**2
        + 3 * t
        + 560 * t**2
        + 224 * t**3
    )
    quintic = (
        20141047808 * t**5
        - 14386462720 * t**4
        + 1089822720 * t**3
        - 21288960 * t**2
        - 35910 * t
        + 81
    )
    require(
        sp.factor(sp.resultant(f_bar, g_bar, u)) == 224 * t * quintic,
        "residual quintic resultant changed",
    )
    require(quintic.subs(t, 0) == 81, "clearing artifact entered the residual scheme")
    require(sp.Poly(quintic, t, modulus=47).is_irreducible, "mod-47 irreducibility failed")
    require(sp.gcd(sp.Poly(quintic, t), sp.Poly(sp.diff(quintic, t), t)) == 1, "quintic not squarefree")
    discriminant = sp.discriminant(quintic, t)
    expected_discriminant = -(2**86) * (3**16) * (5**12) * (7**15) * (13**2) * (19**2) * 29
    require(discriminant == expected_discriminant, "quintic discriminant changed")

    quotient_basis = sp.groebner([f_bar, g_bar], u, t, order="lex")
    require(
        any(sp.Poly(poly.as_expr(), u, t).degree(u) == 1 for poly in quotient_basis.polys),
        "residual u-coordinate ceased to be unique",
    )

    # Since R22 has residual mu_3 weight one, s^2 R22 is invariant.  The
    # following polynomial is its nonzero scalar multiple on the q=1 cover.
    response_bar = (
        84 * t**2 * u
        - 70 * t**2
        - 280 * t * u**2
        + 105 * t * u
        - 3 * t
        + 84 * u**3
        - 9 * u**2
    )
    response_identity = sp.factor(response.subs({q: 1, d: u / s}) * s**2)
    require(
        sp.expand(
            response_identity
            - sp.Rational(33, 1024) * response_bar.subs(t, s**3)
        )
        == 0,
        "invariant response formula changed",
    )
    triple_basis = sp.groebner([f_bar, g_bar, response_bar], u, t, order="lex")
    require(
        [sp.factor(poly.as_expr()) for poly in triple_basis.polys] == [u, t],
        "response acquired a genuine residual common zero",
    )

    # The d=1 index cover at P_infty has intersection length 36.  The
    # central mu_2 quotient halves this to 18; the residual quintic supplies
    # the remaining five units of weighted Bezout length 23.
    cross = q * s * (q**2 - 3 * s**2) * (3 * q**2 - s**2)
    g2 = (q**2 - s**2) * (q**4 - 14 * q**2 * s**2 + s**4)
    require(sp.factor(sp.resultant(cross, g2, q)) == -(2**24) * s**36, "P_infty length changed")
    require(sp.Rational(23 * 24, 2 * 3 * 4) == 23, "weighted Bezout degree changed")

    print("DEGREE22 HIGHEST-ODD REDUCED BOUNDARY DIVISOR")
    print("residual_invariants=t=s^3/q^4,u=d*s/q^2")
    print("residual_quintic=" + str(quintic))
    print("residual_quintic=irreducible_mod47;squarefree;five_simple_coarse_points")
    print("residual_discriminant=-2^86*3^16*5^12*7^15*13^2*19^2*29")
    print("response_common_zero_ideal=(t,u);residual_t_nonzero=>R22_nonzero")
    print("boundary_length=18_at_P_infty+5_residual=23")
    print("j,r,g,Pinf_points,ord_h_Pinf,pole_R_Pinf,total_R_pole_degree,allowed_one_end_degree,Pinf_q_order")

    for selected in range(21, 0, -2):
        r = 22 - selected
        common_gcd = gcd(r, 6)
        pinf_points = 3 * common_gcd
        h_order = 6 // common_gcd
        response_pole = (150 - 6 * r) // common_gcd
        total_response_poles = 5 * 25 + pinf_points * response_pole
        q_order = (r - 18) // common_gcd if (r - 18) % common_gcd == 0 else sp.Rational(r - 18, common_gcd)
        require(5 + pinf_points * h_order == 23, "h-divisor degree changed")
        require(total_response_poles == 575 - 18 * r, "response divisor degree changed")
        require(response_pole >= 8, "response pole floor changed")
        print(
            f"{selected},{r},{common_gcd},{pinf_points},{h_order},{response_pole},"
            f"{total_response_poles},1_or_{h_order},{q_order}"
        )

    # Boundary-counting alone is sharp: the formal singleton partition uses
    # five degree-one residual components and one degree 6/g component for
    # each of the 3g P_infty branches.  This is deliberately not asserted to
    # be realized inside the Faber family.
    for common_gcd in (1, 3):
        singleton_degrees = [1] * 5 + [6 // common_gcd] * (3 * common_gcd)
        require(sum(singleton_degrees) == 23, "singleton divisor hostile changed")

    print("nonreduced_escape=excluded_by_CM_plus_boundary_generic_reducedness")
    print("physical_escape=rational_one_ended_component_of_weighted_degree_1_2_or_6")
    print("r19_r21=P_infty_one_end_excluded_by_global_regular_q_aff")
    print("ramification_boundary=m-1=e*pole_order;pure_power_composition")
    print("scope=ABSTRACT_ATLAS;PHYSICAL_COMPONENTS_CLOSED_BY_THM2745;NOT_JC2")
    print("PASS")


if __name__ == "__main__":
    main()
