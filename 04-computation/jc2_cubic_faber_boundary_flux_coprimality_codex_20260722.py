#!/usr/bin/env python3
"""Exact referee for THM-2118.

The theorem is all-degree and proved on paper by Lagrange inversion and the
hypergeometric ODE.  This companion independently compares the original
THM-2084 Faber sums with the quadratic coefficient model, checks both
differential identities exactly, and runs a second Fraction-only polynomial
gcd census.

Tournament Analysis is not natural here: the objects are two adjacent scalar
coefficient polynomials and an ODE, with no intrinsic pairwise binary relation
or tie Hamiltonian path.  Forcing a tournament would discard the derivative
and singular-point data that the proof uses.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial

import sympy as sp


SYMBOLIC_LIMIT = 80
FRACTION_LIMIT = 200


def generalized_binomial_fraction(x: Fraction, k: int) -> Fraction:
    out = Fraction(1)
    for j in range(k):
        out *= x - j
    return out / factorial(k)


def trim(poly: list[Fraction]) -> list[Fraction]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_divmod(
    numerator: list[Fraction], denominator: list[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    num = trim(numerator[:])
    den = trim(denominator[:])
    if den == [0]:
        raise ZeroDivisionError("zero polynomial")
    if len(num) < len(den):
        return [Fraction(0)], num
    quotient = [Fraction(0)] * (len(num) - len(den) + 1)
    while num != [0] and len(num) >= len(den):
        shift = len(num) - len(den)
        coeff = num[-1] / den[-1]
        quotient[shift] = coeff
        for j, value in enumerate(den):
            num[j + shift] -= coeff * value
        trim(num)
    return trim(quotient), trim(num)


def monic_gcd(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    a = trim(left[:])
    b = trim(right[:])
    while b != [0]:
        _, remainder = poly_divmod(a, b)
        a, b = b, remainder
    lead = a[-1]
    return trim([value / lead for value in a])


def quadratic_coefficient_fraction(n: int, degree: int) -> list[Fraction]:
    """Coefficients in s of [x^degree](1+3x+s*x^2)^(n/3)."""

    alpha = Fraction(n, 3)
    out = [Fraction(0)] * (degree // 2 + 1)
    for j in range(degree // 2 + 1):
        i = degree - 2 * j
        total = i + j
        out[j] = (
            generalized_binomial_fraction(alpha, total)
            * Fraction(factorial(total), factorial(i) * factorial(j))
            * 3**i
        )
    return trim(out)


def exact_fraction_census() -> None:
    for n in range(1, FRACTION_LIMIT + 1):
        if n % 3 == 0:
            continue
        e = quadratic_coefficient_fraction(n, n)
        phi = [3 * value for value in quadratic_coefficient_fraction(n, n + 1)]
        gcd = monic_gcd(e, phi)
        assert gcd == [Fraction(1)], (n, gcd)


def original_faber_boundary(n: int, s: sp.Symbol) -> sp.Expr:
    alpha = sp.Rational(n, 3)
    p = s - 3
    q = 2 - s
    total = sp.Integer(0)
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            if 2 * i + 3 * j <= n:
                total += (
                    sp.binomial(alpha, i + j)
                    * sp.binomial(i + j, i)
                    * p**i
                    * q**j
                )
    return sp.cancel(sp.expand(total))


def original_faber_flux(n: int, s: sp.Symbol) -> sp.Expr:
    alpha = sp.Rational(n, 3)
    p = s - 3
    q = 2 - s
    total = sp.Integer(0)
    for i in range((n + 1) // 2 + 1):
        for j in range((n + 1) // 3 + 1):
            if 2 * i + 3 * j == n + 1:
                total += (
                    3
                    * sp.binomial(alpha, i + j)
                    * sp.binomial(i + j, i)
                    * p**i
                    * q**j
                )
    return sp.cancel(sp.expand(total))


def quadratic_coefficient_sympy(n: int, degree: int, s: sp.Symbol) -> sp.Expr:
    alpha = sp.Rational(n, 3)
    total = sp.Integer(0)
    for j in range(degree // 2 + 1):
        i = degree - 2 * j
        count = i + j
        total += (
            sp.binomial(alpha, count)
            * sp.binomial(count, j)
            * 3**i
            * s**j
        )
    return sp.cancel(sp.expand(total))


def exact_symbolic_census() -> dict[int, tuple[sp.Expr, sp.Expr]]:
    s = sp.symbols("s")
    samples: dict[int, tuple[sp.Expr, sp.Expr]] = {}
    for n in range(1, SYMBOLIC_LIMIT + 1):
        if n % 3 == 0:
            continue
        e_faber = original_faber_boundary(n, s)
        phi_faber = original_faber_flux(n, s)
        e = quadratic_coefficient_sympy(n, n, s)
        phi = 3 * quadratic_coefficient_sympy(n, n + 1, s)

        assert sp.expand(e_faber - e) == 0, n
        assert sp.expand(phi_faber - phi) == 0, n

        first_order = (
            (n + 1) * phi
            - s * (9 - 4 * s) * sp.diff(e, s)
            - 2 * n * (s - 3) * e
        )
        assert sp.expand(first_order) == 0, n

        ode = (
            s * (9 - 4 * s) * sp.diff(e, s, 2)
            + (9 - 6 * n + (4 * n - 6) * s) * sp.diff(e, s)
            - n * (n - 1) * e
        )
        assert sp.expand(ode) == 0, n

        at_zero = 3**n * sp.binomial(sp.Rational(n, 3), n)
        at_quarter = sp.Rational(3, 2) ** n * sp.binomial(
            sp.Rational(2 * n, 3), n
        )
        assert sp.simplify(e.subs(s, 0) - at_zero) == 0, n
        assert sp.simplify(e.subs(s, sp.Rational(9, 4)) - at_quarter) == 0, n
        assert at_zero != 0 and at_quarter != 0, n
        assert sp.Poly(e, s, domain=sp.QQ).gcd(
            sp.Poly(phi, s, domain=sp.QQ)
        ).degree() == 0, n

        if n in {1, 2, 4, 13, 40, 80}:
            samples[n] = (sp.factor(e), sp.factor(phi))
    return samples


def main() -> None:
    samples = exact_symbolic_census()
    exact_fraction_census()
    print("THM-2118 exact cubic Faber boundary-flux referee")
    print(f"sympy={sp.__version__}")
    print(f"symbolic_original-vs-quadratic_and_ODE: n<= {SYMBOLIC_LIMIT}: PASS")
    print(f"independent_Fraction_gcd_census: n<= {FRACTION_LIMIT}: PASS")
    print("tournament_analysis=not_applicable_no_intrinsic_binary_pair_relation")
    for n in sorted(samples):
        e, phi = samples[n]
        print(f"sample n={n}: deg(e)={sp.degree(e)} deg(phi)={sp.degree(phi)}")
    print("VERDICT: PASS")


if __name__ == "__main__":
    main()
