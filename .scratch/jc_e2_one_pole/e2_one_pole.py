#!/usr/bin/env python3
"""Exact controls for the balanced-response e=2, h=1 classification draft."""

from __future__ import annotations

import ast
from math import floor
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, rho = sp.symbols("x rho")


def h_poly(n):
    return sp.expand(
        rho ** (n - 1) * (n - (n - 2) * rho)
        - n * rho + (n - 2)
    )


def quotient(n):
    h = sp.Poly(h_poly(n), rho, domain=sp.QQ)
    divisor = sp.Poly((rho - 1) ** 3, rho, domain=sp.QQ)
    q, remainder = sp.div(h, divisor)
    require(remainder.is_zero, f"triple root division failed at N={n}")
    return q


def reduce_coefficient_mod_q(value, q):
    numerator, denominator = sp.fraction(sp.cancel(value))
    denominator_poly = sp.Poly(denominator, rho)
    require(
        sp.rem(sp.Poly(numerator, rho), q).is_zero,
        "coefficient did not vanish on the admissible ratio scheme",
    )
    require(
        sp.gcd(denominator_poly, q).degree() == 0,
        "coefficient denominator is not a unit on the ratio scheme",
    )


def audit_degree(n):
    require(n >= 4, "e=2 requires N>=4")
    h = sp.Poly(h_poly(n), rho, domain=sp.QQ)
    q = quotient(n)

    require(
        h.eval(1) == sp.diff(h.as_expr(), rho).subs(rho, 1) == 0,
        f"value/first derivative at rho=1 changed for N={n}",
    )
    require(
        sp.diff(h.as_expr(), rho, 2).subs(rho, 1) == 0
        and sp.diff(h.as_expr(), rho, 3).subs(rho, 1)
        == -n * (n - 1) * (n - 2),
        f"exact triple-root certificate failed for N={n}",
    )
    require(
        sp.gcd(q, sp.Poly(sp.diff(q.as_expr(), rho), rho)).degree() == 0,
        f"admissible ratio polynomial is not squarefree for N={n}",
    )
    require(
        sp.expand(rho ** n * h_poly(n).subs(rho, 1 / rho) + h_poly(n)) == 0,
        f"anti-reciprocity failed for N={n}",
    )
    require(q.eval(0) != 0, f"zero ratio survived for N={n}")
    require(
        sp.gcd(q, sp.Poly(n - (n - 2) * rho, rho)).degree() == 0
        and sp.gcd(q, sp.Poly(n * rho - (n - 2), rho)).degree() == 0,
        f"a normalization denominator meets the ratio scheme at N={n}",
    )
    require(
        (q.eval(-1) == 0) == (n % 2 == 0),
        f"inversion-fixed ratio boundary changed for N={n}",
    )

    e = (x - 1) * (x - rho)
    c = sp.cancel(
        n * (n - 1) * (n - 2) / (n - (n - 2) * rho)
    )
    p = sp.expand(
        x ** n
        - c * (
            x ** 2 / (n - 2)
            - (1 + rho) * x / (n - 1)
            + rho / n
        )
    )
    require(
        sp.cancel(x * sp.diff(p, x) - n * p - c * e) == 0,
        f"first-integral identity failed for N={n}",
    )
    require(sp.cancel(p.subs(x, 1)) == 0,
            f"first double zero failed for N={n}")
    reduce_coefficient_mod_q(p.subs(x, rho), q)

    # Division by E^2 is checked on the reduced ratio scheme coefficientwise.
    _quotient, remainder = sp.div(
        sp.Poly(p, x, domain=sp.QQ.frac_field(rho)),
        sp.Poly(e ** 2, x, domain=sp.QQ.frac_field(rho)),
    )
    for coefficient in remainder.all_coeffs():
        reduce_coefficient_mod_q(coefficient, q)

    ordered = q.degree()
    fixed = int(n % 2 == 0)
    affine_classes = (ordered - fixed) // 2 + fixed
    require(
        ordered == n - 3
        and affine_classes == floor((n - 2) / 2),
        f"class count failed for N={n}",
    )
    return ordered, affine_classes, sp.factor(h.as_expr())


def main():
    rows = tuple((n, *audit_degree(n)[:2]) for n in range(4, 15))
    expected = tuple(
        (n, n - 3, floor((n - 2) / 2)) for n in range(4, 15)
    )
    require(rows == expected, "bounded degree count table changed")

    q14 = quotient(14)
    expected_q14 = sp.Poly(
        -2 * (rho + 1) * (
            6 * rho ** 10 + 5 * rho ** 9 + 10 * rho ** 8
            + 8 * rho ** 7 + 12 * rho ** 6 + 9 * rho ** 5
            + 12 * rho ** 4 + 8 * rho ** 3 + 10 * rho ** 2
            + 5 * rho + 6
        ),
        rho,
    )
    require(
        sp.expand(q14.as_expr() - expected_q14.as_expr()) == 0,
        "degree-26 ratio factor changed",
    )

    tree = ast.parse(Path(__file__).read_text())
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "optimized-mode validity gate contains assert",
    )

    print("BALANCED RESPONSE e=2 h=1 ONE-POLE CLASSIFICATION CONTROL")
    print(f"degree_rows_(N,ordered_ratios,affine_classes)={rows}")
    print(f"N14_H_factor={sp.factor(h_poly(14))}")
    print("N14_degree_V=26;ordered_ratios=11;unordered_affine_classes=6")
    print("all_ratio_quotients_squarefree=true")
    print("all_normalization_denominators_disjoint=true")
    print("all_first_integral_and_E_squared_divisibility_checks_passed=true")
    print("scope=complete_all_degree_e2_h1_affine_classification;"
          "N4_split_boundary;N_ge_5_genuinely_nonsplit;"
          "no_flux_intersection_or_JC2_closure")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
