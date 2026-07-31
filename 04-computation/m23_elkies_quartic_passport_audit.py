#!/usr/bin/env python3
"""Exact branch-passport audit for Elkies's degree-23 M23 polynomial.

This does *not* by itself distinguish M23 from A23.  It independently checks
the exact algebraic input to Elkies's certified finite-field/Weil argument:

    P = P2^2 P3 P4^4 = P7 P8^2 + tau

over F = Q(g), g^4 + g^3 + 9g^2 - 10g + 8 = 0, with all displayed factors
coprime.  Thus the finite branch cycles have types 1^3 2^2 4^4 and 1^7 2^8;
infinity has type 23.  It also checks sqrt(-23) lies in F, which makes the
discriminant of P(x)-t a square and hence puts geometric monodromy in A23.

Source:
  N. D. Elkies, "The complex polynomials P(x) with
  Gal(P(x)-t) ~= M23", Open Book Series 1 (2013), 359--367,
  DOI 10.2140/obs.2013.1.359.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    x, g = sp.symbols("x g")
    minpoly = sp.Poly(g**4 + g**3 + 9 * g**2 - 10 * g + 8, g, domain=QQ)
    require(minpoly.is_irreducible, "quartic field polynomial must be irreducible")

    field = QQ.alg_field_from_poly(minpoly)
    a = field.ext

    eta = (2 * a**3 + 4 * a**2 + 16 * a - 7) / 3
    require(field.from_sympy(eta**2 + 23) == field.zero, "eta^2 must equal -23")

    p2 = sp.Poly(
        (8 * a**3 + 16 * a**2 - 20 * a + 20) * x**2
        - (7 * a**3 + 17 * a**2 - 7 * a + 76) * x
        - 13 * a**3
        + 25 * a**2
        - 107 * a
        + 596,
        x,
        domain=field,
    )
    p3 = sp.Poly(
        8 * (31 * a**3 + 405 * a**2 - 459 * a + 333) * x**3
        + (941 * a**3 + 1303 * a**2 - 1853 * a + 1772) * x
        + 85 * a**3
        - 385 * a**2
        + 395 * a
        - 220,
        x,
        domain=field,
    )
    p4 = sp.Poly(
        32 * (4 * a**3 - 69 * a**2 + 74 * a - 49) * x**4
        + 32 * (21 * a**3 + 53 * a**2 - 68 * a + 58) * x**3
        - 8 * (97 * a**3 + 95 * a**2 - 145 * a + 148) * x**2
        + 8 * (41 * a**3 - 89 * a**2 - a + 140) * x
        - 123 * a**3
        + 391 * a**2
        - 93 * a
        + 3228,
        x,
        domain=field,
    )

    p = p2**2 * p3 * p4**4
    require(p.degree() == 23, "P must have degree 23")

    p8, derivative_remainder = sp.div(p.diff(), 23 * p2 * p4**3, domain=field)
    require(derivative_remainder.is_zero, "P'/(23 P2 P4^3) must be exact")
    require(p8.degree() == 8, "P8 must have degree 8")

    tau = (
        sp.Rational(2**38 * 3**17, 23**3)
        * (47323 * a**3 - 1084897 * a**2 + 7751 * a - 711002)
    )
    require(field.from_sympy(tau) != field.zero, "tau must be nonzero")

    p7, second_remainder = sp.div(
        p - sp.Poly(tau, x, domain=field), p8**2, domain=field
    )
    require(second_remainder.is_zero, "(P-tau)/P8^2 must be exact")
    require(p7.degree() == 7, "P7 must have degree 7")

    factors = (p2, p3, p4, p7, p8)
    for i, left in enumerate(factors):
        require(sp.gcd(left, left.diff()).degree() == 0, f"P{i + 2} not squarefree")
        for j, right in enumerate(factors[i + 1 :], i + 1):
            require(
                sp.gcd(left, right).degree() == 0,
                f"factors {i} and {j} not coprime",
            )

    finite_defect_zero = 2 * (2 - 1) + 4 * (4 - 1)
    finite_defect_tau = 8 * (2 - 1)
    infinity_defect = 23 - 1
    require(
        finite_defect_zero + finite_defect_tau + infinity_defect == 2 * 23 - 2,
        "Riemann--Hurwitz defect sum must be 44",
    )

    # Since n=23 is odd, the leading t^22 coefficient of disc_x(P-t) is
    # -23^23*lc(P)^22.  The factorization above gives exactly the two finite
    # discriminant zeros, with exponents 14 and 8.  As eta^2=-23,
    #
    #   -23^23*lc(P)^22 = (eta * 23^11 * lc(P)^11)^2.
    #
    # Hence the whole discriminant is a square in F(t).
    leading_coefficient = field.from_sympy(p.LC())
    square_root_lead = (
        field.from_sympy(eta * 23**11) * leading_coefficient**11
    )
    discriminant_lead = (
        -field.from_sympy(sp.Integer(23) ** 23) * leading_coefficient**22
    )
    require(
        square_root_lead**2 == discriminant_lead,
        "discriminant leading coefficient must be a square",
    )

    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print("status=ELKIES_M23_PASSPORT_EXACT")
    print("field_minpoly=g^4+g^3+9*g^2-10*g+8")
    print("field_contains=eta^2=-23")
    print("degree_P=23")
    print("identity=P2^2*P3*P4^4=P7*P8^2+tau")
    print("branch_types=1^3*2^2*4^4 ; 1^7*2^8 ; 23")
    print("defects=14+8+22=44")
    print("discriminant_shape=square_constant*t^14*(t-tau)^8")
    print("group_boundary=passport_and_square_discriminant_leave_M23_vs_A23")
    print("group_separator=Elkies_effective_Chebotarev_5_set_Weil_bound")
    print(f"script_sha256={source_hash}")


if __name__ == "__main__":
    main()
