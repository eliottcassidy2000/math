#!/usr/bin/env python3
"""Exact hostile controls for the proved THM-2778 all-degree closure.

The all-degree quantifiers in the audit report are proved symbolically there.
This script supplies finite exact controls for every algebraic interface used
by that proof, without importing the proposed theorem's discovery script.
"""

from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    require(condition, message)
    GATES += 1


def coefficient_recurrence(
    exponent: sp.Expr,
    quartic: dict[int, sp.Expr],
    last_index: int,
) -> list[sp.Expr]:
    """Coefficients of (1+sum quartic[r] z^r)^exponent."""

    values = [sp.Integer(1)]
    for index in range(1, last_index + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * values[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index
        values.append(sp.factor(value))
    return values


def top_observables(
    degree: int,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, list[sp.Expr]]:
    coefficients = coefficient_recurrence(
        sp.Rational(degree, 4),
        {2: 2 * d, 3: q, 4: d**2 - s},
        degree + 3,
    )
    return (
        sp.factor(4 * coefficients[degree + 1]),
        sp.factor(4 * coefficients[degree + 2]),
        sp.factor(
            4 * coefficients[degree + 3]
            + 2 * d * coefficients[degree + 1]
        ),
        coefficients,
    )


def perturbation_coefficient(k: int, index: int, q: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    """Coefficient at z^index of (q z^3-s z^4)^k/(1+z^2)."""

    value = sp.Integer(0)
    for q_count in range(k + 1):
        numerator_degree = 4 * k - q_count
        difference = index - numerator_degree
        if difference >= 0 and difference % 2 == 0:
            value += (
                sp.binomial(k, q_count)
                * q**q_count
                * (-s) ** (k - q_count)
                * (-1) ** (difference // 2)
            )
    return sp.expand(value)


def exact_prefix_residuals(k: int, r: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Independent recurrence construction of THM-2760's P_k,Q_k."""

    exponent = sp.Rational(2 * k - 1, 2)
    coefficients = coefficient_recurrence(
        exponent,
        {2: -2, 3: r, 4: 1 - r},
        4 * k,
    )
    divisor = sp.Poly(r**k, r)
    first, first_remainder = sp.div(
        sp.Poly(sp.expand(coefficients[4 * k - 1]), r),
        divisor,
    )
    second, second_remainder = sp.div(
        sp.Poly(sp.expand(coefficients[4 * k]), r),
        divisor,
    )
    require(first_remainder.is_zero and second_remainder.is_zero, f"k={k}: r^k")
    return sp.factor(first.as_expr()), sp.factor(second.as_expr())


def faber_at_original_section(
    degree: int,
    beta: sp.Symbol,
    capital_c: sp.Symbol,
    capital_e: sp.Symbol,
    q: sp.Symbol,
) -> sp.Expr:
    d, s, w = sp.symbols("d s w")
    coefficients = coefficient_recurrence(
        sp.Rational(degree, 4),
        {2: 2 * d, 3: q, 4: d**2 - s},
        degree,
    )
    faber = sum(
        coefficients[index] * w ** (degree - index)
        for index in range(degree + 1)
    )
    return sp.factor(
        faber.subs(
            {
                w: beta,
                d: capital_c - beta**2,
                s: q * beta - capital_e,
            }
        )
    )


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    return sp.groebner(expressions, *variables, order="grevlex").contains(
        sp.Integer(1)
    )


def local_order_or_infinity(
    expression: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
) -> int | None:
    """Return None for the zero polynomial, interpreted as infinite order."""

    polynomial = sp.Poly(sp.expand(expression), first, second)
    if polynomial.is_zero:
        return None
    return min(sum(monomial) for monomial, _ in polynomial.terms())


def main() -> None:
    d, q, s, r, y = sp.symbols("d q s r y")
    beta, capital_c, capital_e = sp.symbols("beta C E")

    checked_k = tuple(range(1, 31))
    checked_degrees = tuple(4 * k - 2 for k in checked_k)

    # Uniform recurrence-tail interface.  The only potentially non-tail
    # predecessor at n=M+4 is c_M in the quartic step, whose multiplier is 0.
    for k, degree in zip(checked_k, checked_degrees):
        exponent = sp.Rational(degree, 4)
        terminal_multiplier = sp.expand((exponent + 1) * 4 - (degree + 4))
        gate(terminal_multiplier == 0, f"M={degree}: terminal multiplier")
        for index in range(degree + 5, degree + 13):
            predecessors = tuple(index - step for step in (2, 3, 4))
            gate(
                min(predecessors) >= degree + 1,
                f"M={degree}: tail induction predecessor escaped",
            )

    # Exact P_infty leading faces, their coprimality, syzygy, pure monomials,
    # exact-prefix residues, and the q=0 top classification.
    for k, degree in zip(checked_k, checked_degrees):
        phi_face = perturbation_coefficient(k, 4 * k - 1, q, s)
        psi_face = perturbation_coefficient(k, 4 * k, q, s)
        next_face = perturbation_coefficient(k, 4 * k + 1, q, s)
        phi_model = sp.expand(
            sum(
                sp.binomial(k, a)
                * (-1) ** ((a - 1) // 2)
                * q**a
                * (-s) ** (k - a)
                for a in range(1, k + 1, 2)
            )
        )
        psi_model = sp.expand(
            sum(
                sp.binomial(k, a)
                * (-1) ** (a // 2)
                * q**a
                * (-s) ** (k - a)
                for a in range(0, k + 1, 2)
            )
        )
        gate(sp.expand(phi_face - phi_model) == 0, f"M={degree}: Phi face")
        gate(sp.expand(psi_face - psi_model) == 0, f"M={degree}: Psi face")
        gate(sp.expand(next_face + phi_face) == 0, f"M={degree}: response syzygy")
        gate(
            sp.gcd(sp.Poly(phi_face, q, s), sp.Poly(psi_face, q, s)).total_degree()
            == 0,
            f"M={degree}: face gcd",
        )
        gate(
            sp.Poly(psi_face, q, s).coeff_monomial(s**k) != 0,
            f"M={degree}: pure s^k",
        )
        pure_q_flux = phi_face if k % 2 else psi_face
        gate(
            sp.Poly(pure_q_flux, q, s).coeff_monomial(q**k) != 0,
            f"M={degree}: pure q^k",
        )
        for beta_residue in (sp.I, -sp.I):
            gate(
                sp.expand(phi_face.subs({q: 1, s: beta_residue})) != 0,
                f"M={degree}: exact-prefix Phi at beta={beta_residue}",
            )
            gate(
                sp.expand(psi_face.subs({q: 1, s: beta_residue})) != 0,
                f"M={degree}: exact-prefix Psi at beta={beta_residue}",
            )
        q_zero_coefficient = 4 * sp.binomial(sp.Rational(degree, 4), k) * (-1) ** k
        gate(q_zero_coefficient != 0, f"M={degree}: q=0 Psi coefficient")

    # Independent exact-prefix gcd replays and root-zero lane classification.
    for k, degree in zip(checked_k, checked_degrees):
        first, second = exact_prefix_residuals(k, r)
        gate(sp.Poly(first, r).eval(0) != 0, f"M={degree}: P_k(0)")
        gate(sp.Poly(second, r).eval(0) != 0, f"M={degree}: Q_k(0)")
        gate(
            sp.gcd(sp.Poly(first, r), sp.Poly(second, r)).degree() == 0,
            f"M={degree}: exact-prefix residual gcd",
        )

        first_zero = (4 * k - 1) % 3 != 0
        second_zero = (4 * k) % 3 != 0
        exceptional = k % 3 == 2
        gate(
            (first_zero and second_zero) == exceptional,
            f"M={degree}: root-zero congruence",
        )
        if exceptional:
            response_index = sp.Rational(4 * k + 1, 3)
            response_coefficient = 4 * sp.binomial(
                sp.Rational(2 * k - 1, 2),
                response_index,
            )
            section_coefficient = sp.binomial(
                sp.Rational(degree, 4),
                sp.Rational(degree, 3),
            )
            gate(response_coefficient != 0, f"M={degree}: root-zero response")
            gate(section_coefficient != 0, f"M={degree}: root-zero section")

    # Every odd row constant and its response ratio; every lower even row has
    # slope-four alignment.  This includes all rows in the largest test bank.
    largest_degree = checked_degrees[-1]
    for odd_degree in range(1, largest_degree, 2):
        c_value = 4 * sp.binomial(
            sp.Rational(odd_degree, 2),
            sp.Rational(odd_degree + 1, 2),
        )
        next_value = (
            4
            * sp.binomial(
                sp.Rational(odd_degree, 2),
                sp.Rational(odd_degree + 3, 2),
            )
            + 2
            * sp.binomial(
                sp.Rational(odd_degree, 2),
                sp.Rational(odd_degree + 1, 2),
            )
        )
        expected_ratio = sp.Rational(odd_degree + 1, 2 * (odd_degree + 3))
        gate(c_value != 0, f"j={odd_degree}: odd Phi constant")
        gate(
            sp.factor(next_value / c_value - expected_ratio) == 0,
            f"j={odd_degree}: odd response ratio",
        )
        gate(
            sp.Rational(1, 2) + expected_ratio != 0,
            f"j={odd_degree}: cancellation multiplier",
        )
        vertical_lead = 4 * sp.binomial(
            sp.Rational(odd_degree, 2),
            sp.Rational(odd_degree + 1, 2),
        )
        gate(vertical_lead != 0, f"j={odd_degree}: vertical triangular lead")

    for top_degree in checked_degrees:
        top_k = (top_degree + 2) // 4
        for lower_degree in range(2, top_degree, 4):
            lower_k = (lower_degree + 2) // 4
            gate(
                (top_degree - lower_degree) + 4 * lower_k == top_degree + 2,
                f"M={top_degree},j={lower_degree}: slope-four alignment",
            )
            gate(lower_k < top_k, f"M={top_degree},j={lower_degree}: lower row")

    # Full section reconstruction through degree 22 checks the weight grading
    # and the exact pure-q coefficient independently of the proposed scout.
    section_degrees = tuple(j for j in range(1, 23) if j % 4 != 0)
    for degree in section_degrees:
        section = faber_at_original_section(degree, beta, capital_c, capital_e, q)
        polynomial = sp.Poly(sp.expand(section), beta, capital_c, capital_e, q)
        for monomial, coefficient in polynomial.terms():
            gate(
                monomial[0] + 2 * monomial[1] + 4 * monomial[2] + 3 * monomial[3]
                == degree,
                f"E_{degree}(0): weight grading",
            )
        pure_q = polynomial.coeff_monomial(q ** (degree // 3)) if degree % 3 == 0 else 0
        expected_q = (
            sp.binomial(sp.Rational(degree, 4), sp.Rational(degree, 3))
            if degree % 3 == 0
            else 0
        )
        gate(pure_q == expected_q, f"E_{degree}(0): pure-q coefficient")

    # Direct low-degree hostile saturation: the common top triple has no
    # projective point with q or s nonzero.  This is a finite check of the
    # uniform square-tail classification, not its proof.
    saturated_degrees = (2, 6, 10, 14, 18)
    for degree in saturated_degrees:
        phi, psi, response, _ = top_observables(degree, d, q, s)
        k = (degree + 2) // 4
        phi_order = local_order_or_infinity(phi.subs(d, 1), q, s)
        psi_order = local_order_or_infinity(psi.subs(d, 1), q, s)
        difference_order = local_order_or_infinity(
            (response + phi / 2).subs(d, 1),
            q,
            s,
        )
        gate(phi_order == k, f"M={degree}: full Phi local order")
        gate(psi_order == k, f"M={degree}: full Psi local order")
        gate(
            difference_order is None or difference_order >= k + 1,
            f"M={degree}: full response syzygy local order",
        )
        gate(
            unit_ideal([phi, psi, response, 1 - y * q], d, q, s, y),
            f"M={degree}: common triple acquired q-nonzero support",
        )
        gate(
            unit_ideal([phi, psi, response, 1 - y * s], d, q, s, y),
            f"M={degree}: common triple acquired s-nonzero support",
        )

    # Minimal M=2 hostile: write the entire bank explicitly.
    phi_2, psi_2, response_2, _ = top_observables(2, d, q, s)
    gate(sp.expand(phi_2 - 2 * q) == 0, "M=2 Phi")
    gate(sp.expand(psi_2 + 2 * s) == 0, "M=2 Psi")
    gate(sp.expand(response_2 + d * q) == 0, "M=2 response")
    gate(
        local_order_or_infinity((response_2 + phi_2 / 2).subs(d, 1), q, s)
        is None,
        "M=2 zero response syzygy was not treated as infinite order",
    )

    source = Path(__file__).read_text(encoding="utf-8")
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "optimized mode could erase a validity gate",
    )

    print("THM-2778 ALL-DEGREE SPLIT CLOSURE CERTIFICATE")
    print(f"checked_k=1..{checked_k[-1]}")
    print(f"checked_M={checked_degrees[0]}..{checked_degrees[-1]} step 4")
    print(f"section_reconstruction=1..22 excluding multiples of 4")
    print(f"triple_saturation={','.join(map(str, saturated_degrees))}")
    print(f"exact_gates={GATES}")
    print("tail_square_interface=PASS")
    print("nonzero_prefix_gcd_interface=PASS")
    print("rootzero_section_interface=PASS")
    print("full_bank_slope_four_interface=PASS")
    print("vertical_and_M2_interfaces=PASS")
    print("AUDIT_EXACT_CONTROLS=PASS")


if __name__ == "__main__":
    main()
