#!/usr/bin/env python3
"""Exact G2/third-response sidecar to THM-2719.

This is a finite exact scout, not a JC(2) closure.  It reconstructs the first,
second, and response Laurent observables in degrees 21 and 22 and computes
their initial forms at the unique infinity point of THM-2719.
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
    response = sp.factor(4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1])
    return phi, psi, response


def initial_form(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> tuple[int, sp.Expr]:
    terms = sp.Poly(sp.expand(expression), q, s).terms()
    order = min(sum(monomial) for monomial, _ in terms)
    face = sum(
        coefficient * q**monomial[0] * s**monomial[1]
        for monomial, coefficient in terms
        if sum(monomial) == order
    )
    return order, sp.factor(face)


def main() -> None:
    d, q, s = sp.symbols("d q s")
    phi21, _, response21 = observables(21, d, q, s)
    phi22, psi22, response22 = observables(22, d, q, s)

    phi21_zero = sp.factor(phi21.subs({d: 1, q: 0, s: 0}))
    response21_zero = sp.factor(response21.subs({d: 1, q: 0, s: 0}))
    require(phi21_zero == sp.Rational(88179, 131072), "Phi21 transverse row changed")
    require(response21_zero == sp.Rational(323323, 1048576), "R21 row changed")
    ratio = sp.factor(response21_zero / phi21_zero)
    require(ratio == sp.Rational(11, 24), "R21/Phi21 ratio changed")

    phi_order, phi_face = initial_form(phi22.subs(d, 1), q, s)
    psi_order, psi_face = initial_form(psi22.subs(d, 1), q, s)
    response_order, response_face = initial_form(response22.subs(d, 1), q, s)
    cross_face = q * s * (q**2 - 3 * s**2) * (3 * q**2 - s**2)
    g2_face = (q**2 - s**2) * (q**4 - 14 * q**2 * s**2 + s**4)
    require(phi_order == psi_order == response_order == 6, "local order-six face changed")
    require(
        sp.expand(phi_face + sp.Rational(231, 128) * cross_face) == 0,
        "Phi22 face changed",
    )
    require(
        sp.expand(psi_face + sp.Rational(231, 256) * g2_face) == 0,
        "Psi22 face changed",
    )
    require(
        sp.expand(response_face - sp.Rational(231, 256) * cross_face) == 0,
        "R22 face changed",
    )

    # Solving F23=0 gives h_6=-Phi22_6/(a21 Phi21(0)).  The a21 factor
    # cancels in the degree-25 response, leaving this exact initial form.
    solved_response_face = sp.factor(response_face - ratio * phi_face)
    require(
        sp.expand(solved_response_face - sp.Rational(1771, 1024) * cross_face)
        == 0,
        "solved third-response face changed",
    )
    require(
        sp.gcd(sp.Poly(solved_response_face, q, s), sp.Poly(psi_face, q, s)) == 1,
        "third response acquired a G2 tangent branch",
    )

    imaginary_unit = sp.I
    g2_complex = sp.expand(
        ((q + imaginary_unit * s) ** 6 + (q - imaginary_unit * s) ** 6) / 2
    )
    require(sp.expand(g2_complex - g2_face) == 0, "dihedral/G2 identity changed")

    h_order = 6
    physical_q_pole = 3 * h_order - 1
    response_pole = 25 * h_order - response_order
    require((physical_q_pole, response_pole) == (17, 144), "pole invoice changed")

    print("DEGREE22 ODD G2 THIRD-RESPONSE FACE")
    print("cover_tangent=Re((q+i*s)^6)=I2(6)=G2_reflection_arrangement")
    print("central_C2_pairs_six_lines_to_three_coarse_branches")
    print("Phi21_at_infinity=88179/131072")
    print("R21_at_infinity=323323/1048576")
    print("R21_over_Phi21=11/24")
    print("top_orders=Phi22:6,Psi22:6,R22:6")
    print("solved_R25_face=1771/1024*q*s*(q^2-3s^2)*(3q^2-s^2)")
    print("gcd(solved_R25_face,Psi22_face)=1")
    print("each_coarse_branch_poles=physical_q:17,affine_response:144")
    print("scope=FINITE_EXACT_SIDECAR_NOT_THIRD_FLUX_CLOSURE_NOT_JC2")
    print("PASS")


if __name__ == "__main__":
    main()
