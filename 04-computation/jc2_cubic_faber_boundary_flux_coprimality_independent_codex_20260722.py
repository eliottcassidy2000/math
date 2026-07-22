#!/usr/bin/env python3
"""Exact referee for THM-2118 (all-degree cubic Faber boundary-flux law).

Universe: exact rational polynomial arithmetic in Q[z,p,q,a,u,x].  The
all-degree quantifiers are proved in THM-2118.  This referee independently
checks the translated coefficient formula, boundary-flux identity,
hypergeometric ODE, gcd conclusion, and valuation coefficients through n=80.
It also checks the n=13 resultant, excluded 3|n boundary, and representative
power-free cubic faces.

Reproduce:
    python3 04-computation/jc2_cubic_faber_boundary_flux_coprimality_codex_20260722.py
    python3 -O 04-computation/jc2_cubic_faber_boundary_flux_coprimality_codex_20260722.py

The check function raises RuntimeError explicitly, so optimized mode executes
the same tests as normal mode.
"""

from __future__ import annotations

from functools import lru_cache
from math import gcd

import sympy as sp


N_MAX = 80
z, p, q, a, u, x = sp.symbols("z p q a u x")


def require(condition: bool, label: str) -> None:
    """Fail loudly in both normal and optimized Python modes."""
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")


@lru_cache(maxsize=None)
def faber(n: int) -> sp.Expr:
    """Polynomial part at z=infinity of (z^3+pz+q)^(n/3)."""
    alpha = sp.Rational(n, 3)
    return sp.expand(
        sum(
            sp.binomial(alpha, i + j)
            * sp.binomial(i + j, i)
            * p**i
            * q**j
            * z ** (n - 2 * i - 3 * j)
            for i in range(n // 2 + 1)
            for j in range(n // 3 + 1)
            if 2 * i + 3 * j <= n
        )
    )


@lru_cache(maxsize=None)
def tail(n: int, offset: int = 1) -> sp.Expr:
    """Three times [z^(-offset)] (z^3+pz+q)^(n/3)."""
    alpha = sp.Rational(n, 3)
    weight = n + offset
    return sp.expand(
        sum(
            3
            * sp.binomial(alpha, i + j)
            * sp.binomial(i + j, i)
            * p**i
            * q**j
            for i in range(weight // 2 + 1)
            for j in range(weight // 3 + 1)
            if 2 * i + 3 * j == weight
        )
    )


def shifted_coefficient(n: int, k: int) -> sp.Expr:
    """[u^k](1+3u+(a+3)u^2)^(n/3), by an independent sum."""
    alpha = sp.Rational(n, 3)
    return sp.expand(
        sum(
            sp.binomial(alpha, k - j)
            * sp.binomial(k - j, j)
            * 3 ** (k - 2 * j)
            * (a + 3) ** j
            for j in range(k // 2 + 1)
        )
    )


def boundary_pair(n: int) -> tuple[sp.Expr, sp.Expr]:
    """The original Faber definitions e_n and phi_n."""
    e_n = sp.expand(faber(n).subs({z: 1, p: a, q: -1 - a}))
    phi_n = sp.expand(tail(n, 1).subs({p: a, q: -1 - a}))
    return e_n, phi_n


def primitive_integer_poly(expr: sp.Expr) -> sp.Poly:
    """Primitive integer polynomial with positive leading coefficient."""
    _, cleared = sp.Poly(expr, a, domain=sp.QQ).clear_denoms()
    _, primitive = cleared.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def multiplicity_gcd(expr: sp.Expr) -> int:
    """Gcd of irreducible multiplicities over Q; 1 means power-free."""
    _, factors = sp.factor_list(expr, gens=(x, z))
    require(bool(factors), f"nonconstant factorization for {expr}")
    value = 0
    for _, exponent in factors:
        value = gcd(value, int(exponent))
    return value


checked = 0
for n in range(1, N_MAX + 1):
    if n % 3 == 0:
        continue

    e_n, phi_n = boundary_pair(n)
    translated_e = shifted_coefficient(n, n)
    translated_phi = 3 * shifted_coefficient(n, n + 1)
    require(sp.expand(e_n - translated_e) == 0, f"translated e at n={n}")
    require(sp.expand(phi_n - translated_phi) == 0, f"translated phi at n={n}")

    delta = (a + 3) * (4 * a + 3)
    first_order = sp.expand(
        (n + 1) * phi_n - 2 * n * a * e_n + delta * sp.diff(e_n, a)
    )
    require(first_order == 0, f"boundary-flux identity at n={n}")

    ode = sp.expand(
        delta * sp.diff(e_n, a, 2)
        - (2 * n - 3) * (2 * a + 3) * sp.diff(e_n, a)
        + n * (n - 1) * e_n
    )
    require(ode == 0, f"hypergeometric ODE at n={n}")
    require(sp.Poly(sp.gcd(e_n, phi_n), a).degree() == 0, f"gcd at n={n}")
    require(e_n.subs(a, -3) != 0, f"singular point -3 at n={n}")
    require(e_n.subs(a, sp.Rational(-3, 4)) != 0, f"singular point -3/4 at n={n}")

    # Exact nonzero coefficients used in the rho>2H valuation regime.
    phi_poly = sp.Poly(tail(n, 1), p, q)
    if n % 2:
        minimal_j_monomial = p ** ((n + 1) // 2)
    else:
        minimal_j_monomial = p ** ((n - 2) // 2) * q
    require(
        phi_poly.coeff_monomial(minimal_j_monomial) != 0,
        f"large-rho extremal coefficient at n={n}",
    )

    # Exact nonzero coefficient used when rho<2H.
    r = n // 3
    low_rho_boundary = sp.expand(faber(n).subs({z: 1, p: 0, q: -1}))
    expected_low_rho = (-1) ** r * sp.binomial(sp.Rational(n, 3) - 1, r)
    require(low_rho_boundary == expected_low_rho, f"small-rho boundary at n={n}")
    require(low_rho_boundary != 0, f"small-rho nonzero at n={n}")

    # Exact parity split used to remove a remaining finite pole of p.
    if n % 2:
        require(
            phi_poly.coeff_monomial(p ** ((n + 1) // 2)) != 0,
            f"odd p-pole coefficient at n={n}",
        )
    else:
        e_poly = sp.Poly(faber(n), z, p, q)
        require(
            e_poly.coeff_monomial(p ** (n // 2))
            == sp.binomial(sp.Rational(n, 3), n // 2),
            f"even p-pole coefficient at n={n}",
        )

    checked += 1

print(f"translated Faber boundary/flux formulas: {checked} degrees through n={N_MAX} PASS")
print(f"first-order identity and hypergeometric ODE: {checked} degrees PASS")
print(f"exact coprimality and singular-point controls: {checked} degrees PASS")
print(f"uniform valuation and parity coefficients: {checked} degrees PASS")


# The degree-thirteen exact resultant agrees with THM-2110.
e_13, phi_13 = boundary_pair(13)
e_13_primitive = primitive_integer_poly(e_13)
phi_13_primitive = primitive_integer_poly(phi_13)
resultant_13 = int(sp.resultant(e_13_primitive.as_expr(), phi_13_primitive.as_expr(), a))
expected_resultant_13 = -(3**21) * (17**5) * (23**2)
require(resultant_13 == expected_resultant_13, "degree-thirteen resultant")
print(f"degree-thirteen primitive resultant: {resultant_13} PASS")


# The excluded 3|n boundary is real: the fractional power is then polynomial.
for n in (3, 6, 9, 12, 15):
    e_n, phi_n = boundary_pair(n)
    require(e_n == 0 and phi_n == 0, f"excluded divisible degree n={n}")
print("excluded 3-divisible common-zero boundary: n=3,6,9,12,15 PASS")


# Representative instances of each possible depressed-cubic leading face.
faces = (
    z**3 + x**2 * z,
    z**3 + x**3,
    z**3 + x**2 * z + x**3,
)
for face in faces:
    require(multiplicity_gcd(face) == 1, f"power-free cubic face {face}")
require(multiplicity_gcd(z**3) == 3, "proper-cube hostile control")
print("depressed-cubic power-free faces and proper-cube hostile control: PASS")

print("RESULT: PASS")
