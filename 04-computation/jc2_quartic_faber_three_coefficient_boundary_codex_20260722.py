#!/usr/bin/env python3
"""Exact referee for THM-2129's quartic boundary classification.

The paper proof is all-degree.  This companion independently checks the
Laurent-tail flux identity, compares a differential recurrence with the direct
multinomial coefficient formula, and computes exact Groebner certificates for
the first odd and twice-odd degrees.

Tournament Analysis is not natural here: the faithful carrier is an ordered
triple of consecutive coefficients plus its recurrence.  Turning the three
coordinates into pairwise wins would forget adjacency, the resonant index, and
the proper-square exceptional point.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


DEGREE_LIMIT = 14


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def falling(alpha: sp.Rational, length: int) -> sp.Rational:
    answer = sp.Rational(1)
    for j in range(length):
        answer *= alpha - j
    return answer


def direct_coefficient(
    n: int, degree: int, beta: sp.Symbol, gamma: sp.Symbol
) -> sp.Expr:
    """[u^degree](1+4u+beta*u^2+gamma*u^3)^(n/4)."""

    alpha = sp.Rational(n, 4)
    answer = sp.Rational(0)
    for ell in range(degree // 3 + 1):
        for j in range((degree - 3 * ell) // 2 + 1):
            i = degree - 2 * j - 3 * ell
            used = i + j + ell
            answer += (
                falling(alpha, used)
                * 4**i
                * beta**j
                * gamma**ell
                / (factorial(i) * factorial(j) * factorial(ell))
            )
    return sp.expand(answer)


def original_boundary_tail_coefficient(
    n: int, degree: int, beta: sp.Symbol, gamma: sp.Symbol
) -> sp.Expr:
    """Coefficient in the original depressed boundary quartic branch."""

    alpha = sp.Rational(n, 4)
    p = beta - 6
    q = gamma - 2 * beta + 8
    r = beta - gamma - 3
    answer = sp.Rational(0)
    for ell in range(degree // 4 + 1):
        for j in range((degree - 4 * ell) // 3 + 1):
            remainder = degree - 4 * ell - 3 * j
            if remainder % 2:
                continue
            i = remainder // 2
            used = i + j + ell
            answer += (
                falling(alpha, used)
                * p**i
                * q**j
                * r**ell
                / (factorial(i) * factorial(j) * factorial(ell))
            )
    return sp.expand(answer)


def recurrence_coefficients(
    n: int, limit: int, beta: sp.Symbol, gamma: sp.Symbol
) -> list[sp.Expr]:
    """Coefficients from T F'=(n/4)T'F, through the requested degree."""

    alpha = sp.Rational(n, 4)
    coefficients = [sp.Integer(1)]
    for k in range(limit):
        c_k = coefficients[k]
        c_km1 = coefficients[k - 1] if k >= 1 else sp.Integer(0)
        c_km2 = coefficients[k - 2] if k >= 2 else sp.Integer(0)
        next_coefficient = -(
            4 * (k - alpha) * c_k
            + beta * (k - 1 - 2 * alpha) * c_km1
            + gamma * (k - 2 - 3 * alpha) * c_km2
        ) / (k + 1)
        coefficients.append(sp.expand(next_coefficient))
    return coefficients


def polynomial_part_in_z(expression: sp.Expr, z: sp.Symbol) -> sp.Expr:
    answer = sp.Integer(0)
    for term in sp.expand(expression).as_ordered_terms():
        exponent = term.as_powers_dict().get(z, sp.Integer(0))
        if exponent >= 0:
            answer += term
    return sp.expand(answer)


def check_tail_flux_identity() -> None:
    z, p, q = sp.symbols("z p q")
    px, qx, rx = sp.symbols("px qx rx")
    c = sp.symbols("c1:7")
    cx = sp.symbols("c1x:7x")
    p_z = 4 * z**3 + 2 * p * z + q
    p_x = px * z**2 + qx * z + rx
    tail = sum(c[j - 1] * z ** (-j) for j in range(1, 7))
    tail_x = sum(cx[j - 1] * z ** (-j) for j in range(1, 7))
    hamiltonian_tail = p_z * tail_x - p_x * sp.diff(tail, z)
    actual = polynomial_part_in_z(hamiltonian_tail, z)
    expected = (
        4 * cx[0] * z**2
        + 4 * cx[1] * z
        + 4 * cx[2]
        + 2 * p * cx[0]
        + px * c[0]
    )
    require(sp.expand(actual - expected) == 0, "quartic tail-flux identity")


def first_power_in_ideal(
    basis: sp.GroebnerBasis, polynomial: sp.Expr, maximum: int
) -> int | None:
    for exponent in range(1, maximum + 1):
        _, remainder = basis.reduce(sp.expand(polynomial**exponent))
        if sp.expand(remainder) == 0:
            return exponent
    return None


def check_boundary_ideals() -> list[tuple[int, str]]:
    beta, gamma = sp.symbols("beta gamma")
    summaries: list[tuple[int, str]] = []
    for n in range(1, DEGREE_LIMIT + 1):
        if n % 4 == 0:
            continue
        coefficients = recurrence_coefficients(n, n + 2, beta, gamma)
        for degree, coefficient in enumerate(coefficients):
            direct = direct_coefficient(n, degree, beta, gamma)
            require(
                sp.expand(coefficient - direct) == 0,
                f"recurrence/direct mismatch n={n}, degree={degree}",
            )

        triple = coefficients[n : n + 3]
        original_c1 = original_boundary_tail_coefficient(n, n + 1, beta, gamma)
        original_c2 = original_boundary_tail_coefficient(n, n + 2, beta, gamma)
        require(
            sp.expand(original_c1 - coefficients[n + 1]) == 0,
            f"translated residue mismatch n={n}",
        )
        require(
            sp.expand(original_c2 - coefficients[n + 1] - coefficients[n + 2]) == 0,
            f"translated second-tail mismatch n={n}",
        )
        basis = sp.groebner(triple, beta, gamma, order="lex", domain=sp.QQ)
        if n % 2 == 1:
            require(
                any(poly.as_expr() == 1 for poly in basis.polys),
                f"odd boundary ideal is not unit at n={n}",
            )
            summaries.append((n, "unit"))
            continue

        require(n % 4 == 2, f"unexpected degree class n={n}")
        require(
            all(sp.expand(poly.subs({beta: 4, gamma: 0})) == 0 for poly in triple),
            f"proper-square point missing at n={n}",
        )
        beta_power = first_power_in_ideal(basis, beta - 4, 2 * n + 6)
        gamma_power = first_power_in_ideal(basis, gamma, 2 * n + 6)
        require(beta_power is not None, f"no beta radical certificate at n={n}")
        require(gamma_power is not None, f"no gamma radical certificate at n={n}")
        summaries.append((n, f"radical=(beta-4,gamma); powers={beta_power},{gamma_power}"))
    return summaries


def main() -> None:
    check_tail_flux_identity()
    summaries = check_boundary_ideals()
    print("THM-2129 exact quartic three-coefficient boundary referee")
    print(f"sympy={sp.__version__}")
    print("quartic_tail_flux_identity: PASS")
    print(f"recurrence_translation_and_Groebner: n<= {DEGREE_LIMIT}: PASS")
    print("tournament_analysis=not_applicable_ordered_triple_and_resonance_sidecars")
    for n, summary in summaries:
        print(f"n={n}: {summary}")
    print("VERDICT: PASS")


if __name__ == "__main__":
    main()
