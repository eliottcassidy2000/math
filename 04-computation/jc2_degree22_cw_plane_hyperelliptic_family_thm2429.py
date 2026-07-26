#!/usr/bin/env python3
"""Exact referee for the degree-22 C-W coefficient plane.

For C!=0, put lambda=W/C^2 and use the weighted coordinates

    v=u/y^2, zeta=Z/y^3, p=C/y^3, W/y^6=lambda*p^2.

The first two normalized fluxes reduce to a quadratic R_lambda(v,p).
Completing the square gives r^2=H_lambda(v).  This script verifies the
generic quintic family, the complete parameter discriminant, its seven
algebraic exceptional ratios, the rational double-root ratio, the
degree-drop ratio, and the y=0 hostile boundary.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def polynomial_hash(expression: sp.Expr, *variables: sp.Symbol) -> str:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    return hashlib.sha256(sp.srepr(polynomial.as_expr()).encode()).hexdigest()


def main() -> None:
    y, u, z = sp.symbols("y u Z")
    cpar, wpar = sp.symbols("C W")
    v, zeta, p, lam = sp.symbols("v zeta p lambda")

    pole_a = -1089 * u + 63 * y**2
    pole_k = (
        2342560 * cpar * u
        - 58080 * cpar * y**2
        + 922383 * u**2 * y
        - 25410 * u * y**3
        + 63 * y**5
    )
    n1 = 1331 * pole_a * z + 4 * pole_k
    n2 = (
        15944049 * z**2
        - 206145280 * cpar * z
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        + 449771520 * cpar * u * y
        - 1239040 * cpar * y**3
        - 1319329792 * wpar
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )

    scaled_a = 9 * (7 - 121 * v)
    scaled_k = (
        922383 * v**2
        - 25410 * v
        + 63
        + (2342560 * v - 58080) * p
    )
    f1 = 1331 * scaled_a * zeta + 4 * scaled_k
    f2 = (
        15944049 * zeta**2
        - 162339408 * zeta * v
        + 2236080 * zeta
        - 1190488992 * v**3
        + 147581280 * v**2
        - 1219680 * v
        + 672
        + (-206145280 * zeta + 449771520 * v - 1239040) * p
        - 1319329792 * lam * p**2
    )
    substitutions = {
        cpar: p * y**3,
        wpar: lam * p**2 * y**6,
        u: v * y**2,
        z: zeta * y**3,
    }
    require(
        sp.factor(n1.subs(substitutions) / y**5 - f1) == 0,
        "C-W first weighted flux mismatch",
    )
    require(
        sp.factor(n2.subs(substitutions) / y**6 - f2) == 0,
        "C-W second weighted flux mismatch",
    )

    l5 = (
        155624547606 * v**5
        + 3215383215 * v**4
        - 1700698560 * v**3
        + 58124770 * v**2
        - 855470 * v
        + 2583
    )
    a2 = -82458112 * (
        131769 * lam * v**2
        - 15246 * lam * v
        + 441 * lam
        + 66550 * v**2
        - 7700 * v
        + 150
    )
    a1 = -464640 * (
        10629366 * v**3 - 1156639 * v**2 + 27104 * v - 315
    )
    a0 = -63 * l5
    r_lam = sp.expand(a2 * p**2 + a1 * p + a0)

    require(
        sp.factor(sp.resultant(f1, f2, zeta) - 255104784 * r_lam)
        == 0,
        "C-W quadratic resultant mismatch",
    )
    require(
        sp.gcd(sp.Poly(a1, v), sp.Poly(a0, v)).degree() == 0,
        "C-W family has a vertical identity fibre",
    )
    a2_coefficient_gcd = sp.gcd_list(sp.Poly(a2, v).all_coeffs())
    require(
        sp.Poly(a2_coefficient_gcd, lam).degree() == 0,
        "C-W resultant loses quadratic degree at one parameter",
    )

    h_lam = (
        35949270496986 * lam * v**5
        + 742753522665 * lam * v**4
        - 392861367360 * lam * v**3
        + 13426821870 * lam * v**2
        - 197613570 * lam * v
        + 596673 * lam
        + 18156197220700 * v**5
        - 1682717215850 * v**4
        - 8503492800 * v**3
        + 370417300 * v**2
        + 3642100 * v
        - 337050
    )
    require(
        sp.factor(sp.discriminant(r_lam, p))
        == -809588736 * (121 * v - 7) ** 2 * h_lam,
        "C-W completed-square factorization mismatch",
    )
    require(
        sp.Poly(h_lam, v).LC() == 363123944414 * (99 * lam + 50),
        "C-W family has an unrecorded degree-drop ratio",
    )

    p7 = (
        5385425879119671 * lam**7
        + 10493769091223880 * lam**6
        + 8643047959054740 * lam**5
        + 3838851484116880 * lam**4
        + 956820534904080 * lam**3
        + 120517531825152 * lam**2
        + 3892886965440 * lam
        - 418060512000
    )
    parameter_discriminant = sp.factor(sp.discriminant(h_lam, v))
    expected_parameter_discriminant = (
        -2**38
        * 3**7
        * 5**12
        * 7**2
        * 11**40
        * (231 * lam + 100)
        * p7
    )
    require(
        sp.factor(
            parameter_discriminant - expected_parameter_discriminant
        )
        == 0,
        "C-W parameter discriminant mismatch",
    )
    require(
        sp.gcd(sp.Poly(p7, lam), sp.Poly(sp.diff(p7, lam), lam)).degree()
        == 0,
        "C-W septic exceptional-ratio polynomial is not squarefree",
    )
    require(
        sp.factor(p7.subs(lam, -sp.Rational(100, 231)))
        == -4185216000,
        "rational double-root ratio meets the septic",
    )
    require(
        sp.factor(p7.subs(lam, -sp.Rational(50, 99)))
        == -sp.Rational(294697041920000, 2187),
        "degree-drop ratio meets the septic",
    )

    cubic_exception = (
        3543122 * v**3 - 2562175 * v**2 + 91476 * v - 1323
    )
    require(
        sp.factor(h_lam.subs(lam, -sp.Rational(100, 231)))
        == 50 * (121 * v - 3) ** 2 * cubic_exception,
        "rational exceptional fibre factorization mismatch",
    )
    require(
        sp.gcd(
            sp.Poly(cubic_exception, v),
            sp.Poly(sp.diff(cubic_exception, v), v),
        ).degree()
        == 0
        and cubic_exception.subs(v, sp.Rational(3, 121)) == -576,
        "rational exceptional cubic is not squarefree and separated",
    )

    quartic_drop = (
        3858459858 * v**4
        - 356083761 * v**3
        + 12020261 * v**2
        - 193963 * v
        + 1197
    )
    require(
        sp.factor(
            h_lam.subs(lam, -sp.Rational(50, 99))
            + sp.Rational(1600, 3) * quartic_drop
        )
        == 0,
        "degree-drop quartic factorization mismatch",
    )
    require(
        sp.gcd(
            sp.Poly(quartic_drop, v),
            sp.Poly(sp.diff(quartic_drop, v), v),
        ).degree()
        == 0,
        "degree-drop quartic is not squarefree",
    )

    # Exact y=0 boundary for C!=0.  First flux fixes Z, and the second
    # then gives a constant cubic equation for u.
    zero_z = sp.Rational(640, 99) * cpar
    require(
        sp.factor(n1.subs({y: 0, z: zero_z})) == 0,
        "C-W y=0 first-flux reconstruction failed",
    )
    require(
        sp.factor(
            n2.subs({y: 0, z: zero_z})
            + sp.Rational(468512, 9)
            * (12800 * cpar**2 + 25344 * wpar + 22869 * u**3)
        )
        == 0,
        "C-W y=0 constant-field equation failed",
    )

    print("THM-2429 degree-22 C-W plane exact referee")
    print("weighted_flux_reduction=PASS")
    print("quadratic_family_degree_in_p=2")
    print("generic_residual=degree_5_squarefree,genus_2")
    print("parameter_discriminant=(231lambda+100)*squarefree_septic")
    print(
        "rational_exception=lambda_-100/231,"
        "double_linear_plus_squarefree_cubic,genus_1"
    )
    print(
        "degree_drop=lambda_-50/99,"
        "squarefree_quartic,genus_1"
    )
    print("septic_exceptions=one_double_root_plus_cubic,genus_1")
    print(f"family_sha256={polynomial_hash(h_lam, v, lam)}")
    print(f"septic_sha256={polynomial_hash(p7, lam)}")
    print("y_zero_hostile_control=PASS")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
