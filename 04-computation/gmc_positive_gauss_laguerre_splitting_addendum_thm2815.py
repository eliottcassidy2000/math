#!/usr/bin/env python3
"""Exact companion for the positive-splitting addendum to THM-2815.

The script never approximates a Laguerre root.  It works in the etale
quotient Q[x]/(ell_N), represents the Gauss weights by one residue-class
polynomial omega_N, and recovers the atomic sums as exact quotient traces.
Sturm root counts certify that all N roots are positive and squarefree.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


x = sp.symbols("x")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def monic_laguerre(n: int) -> sp.Expr:
    """ell_n=(-1)^n n! L_n^(0), monic over Z."""
    return sp.expand(
        sum(
            (-1) ** (n + k)
            * factorial(n)
            * sp.binomial(n, k)
            * x**k
            / factorial(k)
            for k in range(n + 1)
        )
    )


def standard_laguerre(n: int) -> sp.Expr:
    """The standard L_n^(0), normalized by integral norm one."""
    return sp.expand((-1) ** n * monic_laguerre(n) / factorial(n))


def factorial_functional(poly: sp.Expr) -> sp.Expr:
    result = sp.Integer(0)
    for (degree,), coefficient in sp.Poly(sp.expand(poly), x).terms():
        result += coefficient * factorial(degree)
    return sp.simplify(result)


def quotient_remainder(poly: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    return sp.rem(
        sp.Poly(sp.expand(poly), x, domain=sp.QQ),
        sp.Poly(modulus, x, domain=sp.QQ),
    ).as_expr()


def quotient_trace(poly: sp.Expr, modulus: sp.Expr, dimension: int) -> sp.Expr:
    """Trace of multiplication by poly on Q[x]/(modulus)."""
    result = sp.Integer(0)
    for column in range(dimension):
        remainder = sp.Poly(
            quotient_remainder(poly * x**column, modulus),
            x,
            domain=sp.QQ,
        )
        result += remainder.coeff_monomial(x**column)
    return sp.simplify(result)


def weight_residue(n: int, ell: sp.Expr) -> sp.Expr:
    """omega(x_i)=(n!)^2/[x_i ell'(x_i)^2] at every root x_i."""
    derivative = sp.diff(ell, x)
    denominator = sp.Poly(x * derivative**2, x, domain=sp.QQ)
    modulus = sp.Poly(ell, x, domain=sp.QQ)
    inverse = sp.invert(denominator, modulus).as_expr()
    return quotient_remainder(factorial(n) ** 2 * inverse, ell)


max_n = 12
positive_root_counts = []
weight_degrees = []
moment_checks = 0
first_failure_checks = 0
coprime_successors = 0
christoffel_darboux_checks = 0
alternate_weight_checks = 0

for n in range(1, max_n + 1):
    ell = monic_laguerre(n)
    ell_poly = sp.Poly(ell, x, domain=sp.QQ)
    require(ell_poly.LC() == 1, f"monicity failed at N={n}")
    require(
        sp.gcd(ell_poly, ell_poly.diff()).degree() == 0,
        f"squarefreeness failed at N={n}",
    )
    positive_count = ell_poly.count_roots(0, sp.oo)
    require(positive_count == n, f"positive-root count failed at N={n}")
    positive_root_counts.append(positive_count)

    require(
        tuple(factorial_functional(x**r * ell) for r in range(n))
        == (0,) * n,
        f"finite-difference orthogonality failed at N={n}",
    )
    require(
        factorial_functional(ell**2) == factorial(n) ** 2,
        f"Laguerre norm failed at N={n}",
    )

    standard = standard_laguerre(n)
    kernel_diagonal = sp.expand(
        sum(standard_laguerre(k) ** 2 for k in range(n))
    )
    require(
        quotient_remainder(
            kernel_diagonal - x * sp.diff(standard, x) ** 2,
            ell,
        )
        == 0,
        f"Christoffel-Darboux diagonal failed at N={n}",
    )
    christoffel_darboux_checks += 1

    successor = standard_laguerre(n + 1)
    require(
        quotient_remainder(
            (n + 1) ** 2 * successor**2
            - x**2 * sp.diff(standard, x) ** 2,
            ell,
        )
        == 0,
        f"alternate Gauss-weight formula failed at N={n}",
    )
    alternate_weight_checks += 1

    omega = weight_residue(n, ell)
    weight_degrees.append(sp.Poly(omega, x).degree())
    require(
        quotient_trace(omega, ell, n) == 1,
        f"weight mass failed at N={n}",
    )
    for degree in range(2 * n):
        actual = quotient_trace(omega * x**degree, ell, n)
        require(
            actual == factorial(degree),
            f"Gauss moment failed at N={n}, degree={degree}",
        )
        moment_checks += 1

    degree_2n = quotient_trace(omega * x ** (2 * n), ell, n)
    expected_2n = factorial(2 * n) - factorial(n) ** 2
    require(
        degree_2n == expected_2n,
        f"first failed moment changed at N={n}",
    )
    require(
        quotient_trace(omega * ell**2, ell, n) == 0,
        f"quadrature did not kill ell_N^2 at N={n}",
    )
    first_failure_checks += 1

    hankel = sp.Matrix(
        [[factorial(i + j) for j in range(n)] for i in range(n)]
    )
    require(
        hankel.det() == sp.prod(factorial(j) ** 2 for j in range(n)),
        f"Hankel determinant failed at N={n}",
    )

    if n >= 2:
        ell_previous = monic_laguerre(n - 2)
        ell_current = monic_laguerre(n - 1)
        recurrence = sp.expand(
            ell
            - (
                (x - (2 * (n - 1) + 1)) * ell_current
                - (n - 1) ** 2 * ell_previous
            )
        )
        require(recurrence == 0, f"monic recurrence failed at N={n}")
    ell_next = monic_laguerre(n + 1)
    require(
        sp.gcd(
            sp.Poly(ell, x, domain=sp.QQ),
            sp.Poly(ell_next, x, domain=sp.QQ),
        ).degree()
        == 0,
        f"successive Laguerre polynomials share a root at N={n}",
    )
    coprime_successors += 1

# Exact N=2 algebraic positive split and scalar-nullity hostile.
ell_2 = monic_laguerre(2)
omega_2 = weight_residue(2, ell_2)
require(ell_2 == x**2 - 4 * x + 2, "N=2 node polynomial changed")
require(omega_2 == 1 - x / 4, "N=2 weight residue changed")
require(
    quotient_trace(omega_2 * (x - 1), ell_2, 2) == 0,
    "N=2 scalar-nullity hostile changed",
)
require(
    sp.gcd(
        sp.Poly(x - 1, x, domain=sp.QQ),
        sp.Poly(ell_2, x, domain=sp.QQ),
    ).degree()
    == 0,
    "N=2 hostile unexpectedly vanished at a node",
)

print("GMC POSITIVE GAUSS-LAGUERRE FINITE CARRIER")
print(
    "status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED; "
    "addendum_to=THM-2815"
)
print(f"N_range=1..{max_n}")
print(f"positive_root_counts={tuple(positive_root_counts)}")
print(f"weight_residue_degrees={tuple(weight_degrees)}")
print(f"exact_moment_trace_checks={moment_checks}")
print(f"first_failure_checks={first_failure_checks}")
print(f"successive_coprime_checks={coprime_successors}")
print(f"christoffel_darboux_checks={christoffel_darboux_checks}")
print(f"alternate_weight_checks={alternate_weight_checks}")
print(f"N2_node_polynomial={ell_2}")
print(f"N2_weight_residue={omega_2}")
print("N2_scalar_null_hostile=x-1; atomic_values_nonzero=1")
print(
    f"N12_degree24_quadrature_moment="
    f"{factorial(24)-factorial(12)**2}"
)
print(
    "scope=unique positive N-atom real splitting through degree 2N-1; "
    "cutoff-dependent algebraic nodes; no tower quotient or nullity separation"
)
print("ALL EXACT CHECKS PASSED")
