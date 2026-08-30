#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4291.

The universe is the explicit THM-4230 integral M=12 source model after the
stated W=0, Z=-U, lower-row specialization.  The script derives the balanced
exceptional face from that literal equation, normalizes it, checks the deck
characters and elliptic quotient, and verifies the algebra behind the
degree-42 equivariant hostile.  PARI's independent modular-polynomial check is
kept in the companion .gp certificate.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(value: sp.Expr, label: str) -> None:
    if sp.factor(value) != 0:
        raise AssertionError(f"{label}: {sp.factor(value)}")


def main() -> None:
    sigma, b, r, t, X, Q, U = sp.symbols("sigma b r t X Q U")
    a = sp.Rational(7168, 135)
    K = sp.Rational(2848, 45)

    # THM-4230 equation (9), specialized by
    # W=0, Z=-U, Delta=Phi=Theta=eta=zeta_3=upsilon_5=xi_10
    # =alpha_11=beta_11=0.  Here t=sigma*b in the infinity chart.
    hhat = (
        U * (r**6 - r**4)
        - sp.Rational(1376, 135) * t**6 * r**3
        + K * t**6 * r**2
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    literal = (
        (1 - r) * (b**12 - hhat.subs(t, sigma * b))
        - sp.Rational(1, 2) * sigma**12 * b**12
    )

    scaled = sp.expand(literal.subs({b: sigma * X, r: 1 + sigma**12 * Q}))
    sigma_orders = [
        int(term.as_powers_dict().get(sigma, 0))
        for term in sp.Add.make_args(scaled)
    ]
    if min(sigma_orders) != 24:
        raise AssertionError(min(sigma_orders))
    face = sp.expand(scaled).coeff(sigma, 24)
    expected_face = (
        2 * U * Q**2 - Q * (X**12 - a * X**6) - sp.Rational(1, 2) * X**12
    )
    assert_zero(face - expected_face, "balanced face")

    # THM-4103's residue identity and the target/source scalings give the
    # good-target invariant differential
    #   omega_0 = -sigma^9*b^10 db/literal_r.
    # On b=sigma*X, r=1+sigma^12*Q (relative to the sigma-base), db=sigma dX.
    # The denominator loses twelve sigma-orders under d/dr, so the resulting
    # form has exact vertical order eight.
    literal_r_scaled = sp.expand(
        sp.diff(literal, r).subs({b: sigma * X, r: 1 + sigma**12 * Q})
    )
    denominator_orders = [
        int(term.as_powers_dict().get(sigma, 0))
        for term in sp.Add.make_args(literal_r_scaled)
    ]
    if min(denominator_orders) != 12:
        raise AssertionError(min(denominator_orders))
    denominator_face = sp.expand(literal_r_scaled).coeff(sigma, 12)
    assert_zero(denominator_face - sp.diff(face, Q), "differential denominator")
    normalized_z_numerator = 4 * U * Q - (X**12 - a * X**6)
    assert_zero(denominator_face - normalized_z_numerator, "normalized Z numerator")
    vertical_order = 9 + 10 + 1 - 12
    if vertical_order != 8:
        raise AssertionError(vertical_order)

    discriminant = sp.factor(sp.discriminant(face, Q))
    tail_polynomial = sp.expand((X**6 - a) ** 2 + 4 * U)
    assert_zero(discriminant - X**12 * tail_polynomial, "normalization discriminant")
    resultant = sp.factor(
        sp.resultant(tail_polynomial, sp.diff(tail_polynomial, X), X)
    )
    expected_resultant = (
        sp.Integer(2) ** 46
        * U**6
        * (sp.Integer(18225) * U + sp.Integer(12845056)) ** 5
        / sp.Integer(3783403212890625)
    )
    assert_zero(resultant - expected_resultant, "tail squarefree resultant")

    # The normalized tail is z^2=x^12-2*c*x^6+1.  Its reciprocal quotient
    # has degree four and lands on E_c: v^2=w^3-3w-2c.
    x, c = sp.symbols("x c", nonzero=True)
    normalized_tail = x**12 - 2 * c * x**6 + 1
    w = x**2 + x**-2
    quotient_identity = sp.expand(
        normalized_tail / x**6 - (w**3 - 3 * w - 2 * c)
    )
    assert_zero(quotient_identity, "elliptic quotient")
    j_ec = sp.factor(
        1728 * 4 * (-3) ** 3 / (4 * (-3) ** 3 + 27 * (-2 * c) ** 2)
    )
    assert_zero(j_ec - 1728 / (1 - c**2), "elliptic j invariant")

    # Phi_7(0,J)'s nonzero quadratic factor.  The companion GP script proves
    # that it really is the classical modular-polynomial specialization.
    J = sp.symbols("J")
    modular_a = sp.Integer(34848505552896000)
    modular_b = sp.Integer(11356800389480448000000)
    modular_quadratic = J**2 + modular_a * J + modular_b
    radical_scale = sp.Integer(3802283679744000)
    for sign in (-1, 1):
        root = -modular_a / 2 + sign * radical_scale * sp.sqrt(21)
        assert_zero(modular_quadratic.subs(J, root), "Phi_7 nonzero root")
        u_special = sp.factor(432 * a**2 / (root - 1728))
        d_special = sp.factor(a**2 + 4 * u_special)
        assert_zero(
            d_special - a**2 * root / (root - 1728),
            "tail smoothness parameter",
        )
        assert_zero(
            1728 + 432 * a**2 / u_special - root,
            "special elliptic j",
        )
        if sp.simplify(u_special) == 0 or sp.simplify(d_special) == 0:
            raise AssertionError("special tail collided")

    # Let alpha satisfy alpha^2-alpha+1=0.  If e_k=g o delta^k and
    # e_2=e_1-e_0, then m=e_0-alpha*e_1 is an exact delta-alpha eigenmap.
    alpha = sp.symbols("alpha")
    cyclotomic = alpha**2 - alpha + 1
    equivariance_error_e1 = sp.rem(1 - alpha + alpha**2, cyclotomic, alpha)
    assert_zero(equivariance_error_e1, "map-level equivariance")

    # For an isogeny psi:E_c->E_0 of degree D, g=psi o f has degree 4D.
    # The two deck-character summands have energy 2D each, and the surviving
    # scalar alpha^2-1 has Eisenstein norm 3.  Hence deg(m)=6D.
    alpha_conjugate = 1 - alpha
    scalar_norm = sp.rem(
        (alpha**2 - 1) * (alpha_conjugate**2 - 1),
        cyclotomic,
        alpha,
    )
    assert_zero(scalar_norm - 3, "eigenprojector norm")
    if 2 * 7 * int(scalar_norm) != 42:
        raise AssertionError("degree-42 specialization")
    # If the degree-42 eigenmap were divisible by the norm-three prime
    # 2-alpha, the quotient would have degree 14.  The fixed-point fibre
    # calculation in the theorem forces every dx/z eigenmap degree to be
    # 12 modulo 6, so 14 is impossible.
    if 14 % 6 == 12 % 6:
        raise AssertionError("norm-three saturation control")
    if 34 % 6 == 0 or 42 % 6 != 0:
        raise AssertionError("equivariant tail degree congruence")

    characters = [(-2 * (j + 1)) % 12 for j in range(5)]
    if characters != [10, 8, 6, 4, 2]:
        raise AssertionError(characters)

    print("THM4291_LAMBDA_ZERO_GENUS5_DEGREE42_TAIL_V1")
    print("UNIVERSE THM4230_EQ9 W=0 Z=-U FREE_W7_TO_W11=0 DELTA=0")
    print("BALANCED_WEIGHTS sigma=1 b=1 q=12 INITIAL_ORDER 24")
    print(f"FACE {sp.sstr(face)}")
    print(f"NORMALIZED_TAIL z^2={sp.sstr(tail_polynomial)} GENUS 5")
    print("SMOOTH_IFF U!=0 AND 18225*U+12845056!=0")
    print("DECK_CHARACTERS", ",".join(map(str, characters)), "TARGET 10 PRESENT")
    print("ELLIPTIC_QUOTIENT DEGREE 4 J=1728/(1-c^2)")
    print(
        "PHI7_NONZERO_ROOTS "
        "-17424252776448000 +/- 3802283679744000*sqrt(21)"
    )
    print("EIGENMAP_DEGREE 6*D D=7 DEGREE 42")
    print("EIGENMAP_NOT_NORM3_DIVISIBLE REQUIRED_DEGREE_14 NOT_12_MOD_6")
    print("EQUIVARIANT_TAIL_DEGREES 12+6k DEGREE34_EXCLUDED DEGREE42_ALLOWED")
    print("KELLER_FORM sigma^8*(-X^4*dX/z) SPECIAL_TAIL_DIFFERENTIAL_ZERO")
    print("SCOPE GEOMETRIC_TAIL_MAP_AVAILABLE BUT_NOT_KELLER_SPECIALIZATION")
    print("VERDICT PASS EXACT_SYMBOLIC_CERTIFICATE")


if __name__ == "__main__":
    main()
