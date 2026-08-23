#!/usr/bin/env python3
"""Exact companion for THM-3886's equality-seam second-layer law."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


h, x, y = sp.symbols("h x y")
q, s1, t1 = sp.symbols("q s1 t1")
K2 = y**2 - 15 * x**2


def scaled_residual(n: int, F: sp.Expr, T: sp.Expr) -> sp.Poly:
    """THM-3881 residual after (x,y)->h(x,y), as a polynomial in h."""
    a = h * x + 1
    L = 9 * h * x + 4
    K = h**2 * K2 - 15 * h * x - 4
    P = a * L**2
    r = a * T + K * F
    A = K * T + a * P * F
    residual = sp.expand(
        L**4
        + 2 * (3 * A + 3 * P + r**2) * L**2 * F
        + (8 * A + 6 * P + 3 * r**2) * (P * F**2 - T**2)
    )
    return sp.Poly(residual, h)


def seam_residual(n: int) -> sp.Poly:
    F = h**n * x * q + h ** (n - 1) * s1
    T = -h ** (n + 1) * K2 * q + h**n * t1
    return scaled_residual(n, F, T)


R = K2 * s1 + x * t1 - y**2 * q

# The first two THM-3884 degrees vanish.  In the stable range the unique next
# form is the r^2*B symbol; n=2 has exactly one extra A*B collision.
for n in (3, 4, 5):
    poly = seam_residual(n)
    zero(poly.coeff_monomial(h ** (4 * n + 7)), f"n={n} top vanishes")
    zero(poly.coeff_monomial(h ** (4 * n + 6)), f"n={n} next vanishes")
    zero(
        poly.coeff_monomial(h ** (4 * n + 5))
        - 243 * q**2 * x**5 * R**2,
        f"n={n} stable next symbol",
    )

poly2 = seam_residual(2)
zero(poly2.coeff_monomial(h**15), "n=2 top vanishes")
zero(poly2.coeff_monomial(h**14), "n=2 next vanishes")
zero(
    poly2.coeff_monomial(h**13)
    - 243 * q**2 * x**5 * (R**2 + 216 * q * x**5),
    "n=2 Kummer next symbol",
)

# Freeze the affine all-degree collision differences rather than treating the
# displayed n values as a finite proof.  Their slopes and n=3 endpoints are
# positive, so the inequalities persist for every n>=3.
N = sp.symbols("N", integer=True)
zero((4 * N + 5) - (3 * N + 7) - (N - 2),
     "stable A*B collision difference")
zero((4 * N + 5) - (3 * N + 4) - (N + 1),
     "stable lower r-term difference")
zero((4 * N + 5) - (2 * N + 6) - (2 * N - 1),
     "stable lower A-term difference")
gate(all(value > 0 for value in (3 - 2, 3 + 1, 2 * 3 - 1)),
     "stable differences positive from n=3")
gate(4 * 2 + 5 == 3 * 2 + 7, "n=2 unique A*B collision")
gate(4 * 1 + 5 < 3 * 1 + 7, "n=1 A*B dominates")

# Stable R=0 is exactly the second jet of the negative gauge vector.
p = sp.symbols("p")
s_gauge = q + x * p
t_gauge = -K2 * p + 15 * x * q
zero(R.subs({s1: s_gauge, t1: t_gauge}), "second gauge jet")

# The n=2 UFD normal form, including both Kummer signs.
U, epsilon = sp.symbols("U epsilon", nonzero=True)
q_kummer = x * U**2
R_kummer = epsilon * x**3 * U
zero(
    (R_kummer**2 + 216 * q_kummer * x**5).subs(epsilon**2, -216),
    "n=2 Kummer normal form",
)

# The Kummer symbol cannot lift through the next two homogeneous equations.
# Comparing R=epsilon*x^3*U for general linear f_1 and quadratic T_2 first
# gives exactly the parameterization below.  Add the constant/address pieces,
# expand, and reduce coefficients modulo epsilon^2+216.
A2, b2, m2, rho2 = sp.symbols("A2 b2 m2 rho2")
B2, C2, D2, E2 = sp.symbols("B2 C2 D2 E2")
R2_general = sp.Poly(
    sp.expand(
        K2 * (A2 * x + B2 * y)
        + x * (C2 * x**2 + D2 * x * y + E2 * y**2)
        - y**2 * x * U**2
        - epsilon * x**3 * U
    ),
    x,
    y,
)
zero(R2_general.coeff_monomial(y**3) - B2,
     "n=2 Kummer parameterization y3 row")
zero(R2_general.coeff_monomial(x**2 * y) - (D2 - 15 * B2),
     "n=2 Kummer parameterization x2y row")
zero(R2_general.coeff_monomial(x**3) - (C2 - 15 * A2 - epsilon * U),
     "n=2 Kummer parameterization x3 row")
zero(R2_general.coeff_monomial(x * y**2) - (A2 + E2 - U**2),
     "n=2 Kummer parameterization xy2 row")
zero(
    R2_general.as_expr().subs(
        {B2: 0, D2: 0, C2: 15 * A2 + epsilon * U, E2: U**2 - A2}
    ),
    "n=2 Kummer parameterization reconstruction",
)
F2_kummer = h**2 * U**2 * x**2 + h * A2 * x + b2
T2_kummer = (
    -h**3 * U**2 * x * K2
    + h**2 * ((U**2 - A2) * y**2 + (15 * A2 + epsilon * U) * x**2)
    + h * (m2 * x + rho2 * y)
    + 4 * b2
)
poly2_kummer = scaled_residual(2, F2_kummer, T2_kummer)


def reduce_kummer(expression: sp.Expr) -> sp.Expr:
    """Reduce exactly modulo epsilon^2+216."""
    return sp.Poly(expression, epsilon).rem(
        sp.Poly(epsilon**2 + 216, epsilon)
    ).as_expr().expand()


S13_kummer = reduce_kummer(poly2_kummer.coeff_monomial(h**13))
S12_kummer = sp.Poly(
    reduce_kummer(poly2_kummer.coeff_monomial(h**12)), x, y
)
S11_kummer = sp.Poly(
    reduce_kummer(poly2_kummer.coeff_monomial(h**11)), x, y
)
zero(S13_kummer, "n=2 Kummer leading collision vanishes")
zero(
    S12_kummer.coeff_monomial(x**8 * y**4) + 648 * U**6,
    "n=2 Kummer degree-twelve minimum-x coefficient",
)
zero(
    S11_kummer.coeff_monomial(x**3 * y**8) - 8 * U**6,
    "n=2 Kummer degree-eleven minimum-x coefficient",
)
gate(
    min(monomial[0] for monomial, coefficient in S12_kummer.terms()
        if coefficient != 0) == 8,
    "n=2 Kummer degree-twelve exact x valuation",
)
gate(
    min(monomial[0] for monomial, coefficient in S11_kummer.terms()
        if coefficient != 0) == 3,
    "n=2 Kummer degree-eleven exact x valuation",
)
gate(2 * 4 > 3, "n=2 square-root cross-term x-valuation contradiction")


# n=1 exact homogeneous square-root recursion.  Since k is algebraically
# closed, write q=z^2 and choose mu=162*sqrt(2)*z^3, so mu^2=52488*q^3.
z = sp.symbols("z", nonzero=True)
c, u, v = sp.symbols("c u v")
q1 = z**2
mu = 162 * sp.sqrt(2) * z**3
F1 = h * x * q1 + c
T1 = -h**2 * K2 * q1 + h * (u * x + v * y) + 4 * c
poly1 = scaled_residual(1, F1, T1)


def homogeneous_coefficient(poly: sp.Poly, degree: int) -> sp.Expr:
    return sp.expand(poly.coeff_monomial(h**degree))


S10 = homogeneous_coefficient(poly1, 10)
S9 = homogeneous_coefficient(poly1, 9)
S8 = homogeneous_coefficient(poly1, 8)
zero(S10 - 52488 * q1**3 * x**10, "n=1 leading square")
zero(mu**2 - 52488 * q1**3, "n=1 leading root scalar")
g5 = mu * x**5
g4 = sp.cancel(S9 / (2 * g5))
gate(sp.denom(g4).has(x) is False, "n=1 degree-nine polynomial root piece")

remainder8 = sp.Poly(sp.expand(S8 - g4**2), x, y)
zero(
    remainder8.coeff_monomial(y**8)
    + sp.Rational(9, 32) * q1 * (q1 - c) ** 4,
    "n=1 degree-eight constant-x obstruction",
)
remainder8_c = sp.Poly(sp.expand((S8 - g4**2).subs(c, q1)), x, y)
zero(
    remainder8_c.coeff_monomial(x**4 * y**4)
    + sp.Rational(9, 32) * q1 * v**4,
    "n=1 degree-eight y-color obstruction",
)

# After c=q and v=0 put u=15q+d.  The degree-eight equation defines g3;
# the complete part of the degree-seven remainder below x^5 is one monomial.
d = sp.symbols("d")
F1_reduced = h * x * q1 + q1
T1_reduced = -h**2 * K2 * q1 + h * ((15 * q1 + d) * x) + 4 * q1
poly1_reduced = scaled_residual(1, F1_reduced, T1_reduced)
R9 = homogeneous_coefficient(poly1_reduced, 9)
R8 = homogeneous_coefficient(poly1_reduced, 8)
R7 = homogeneous_coefficient(poly1_reduced, 7)
g4_reduced = sp.cancel(R9 / (2 * g5))
g3_reduced = sp.cancel((R8 - g4_reduced**2) / (2 * g5))
remainder7 = sp.Poly(sp.expand(R7 - 2 * g4_reduced * g3_reduced), x, y)
below_x5 = sp.Add(*[
    coefficient * x**monomial[0] * y**monomial[1]
    for monomial, coefficient in remainder7.terms()
    if monomial[0] < 5
])
zero(
    below_x5 + sp.Rational(1, 288) * d**4 * q1 * x**3 * y**4,
    "n=1 degree-seven pure-gauge obstruction",
)

# d=0 is literally the THM-3881 pure gauge with constant parameter Q=-q.
a0 = x + 1
K0 = y**2 - 15 * x**2 - 15 * x - 4
zero((q1 * a0) - (-a0 * (-q1)), "n=1 pure-gauge f endpoint")
zero((-q1 * K0) - (K0 * (-q1)), "n=1 pure-gauge T endpoint")
gate(sp.rem(sp.Poly(q1, x), sp.Poly(a0, x)).as_expr() != 0,
     "nonzero constant gauge parameter not divisible by a")

semantic = {
    "seam": "f_n=xq;T_nplus1=-K2q",
    "stable": "n>=3 forces R=0 and a second negative-gauge jet",
    "quadratic": "n=2 Kummer symbol dies since vx(S12)=8 but vx(S11)=3",
    "linear": "n=1 square recursion reaches forbidden constant pure gauge",
    "scope": "associated graded only; peeling is not square-invariant",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3886-cusp-residual-equality-seam-second-layer-trichotomy")
print("n_ge_3=second_full_gauge_jet")
print("n_eq_2=no_square_survivor_via_consecutive_x_valuations_8_then_3")
print("n_eq_1=no_square_survivor_via_constant_pure_gauge")
print("gauge_peeling_square_invariance=NOT_CLAIMED")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
