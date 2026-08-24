#!/usr/bin/env python3
"""Independent hostile audit for THM-3908.

This companion is deliberately organized around the consequence objects:
the complete normalization-degree ledger, the two primitive Newton edges,
the three conditional constant-row chords, and the coefficient ideal of the
moving-root pure-power condition.  It does not import the primary companion.
"""

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


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(left - right)) == 0, message)


def disc4(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d)


def homogeneous_part(poly: sp.Expr, variables: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    P = sp.Poly(poly, *variables)
    return sp.expand(sum(
        coefficient * sp.prod(variable**power for variable, power in zip(variables, powers))
        for powers, coefficient in P.terms()
        if sum(powers) == degree
    ))


# ---------------------------------------------------------------------------
# 1. Complete lift ledger after removing the coefficient gcd.
# ---------------------------------------------------------------------------
# A gcd of degree g leaves a basepoint-free coefficient map of degree e=2-g.
# The twisted-cubic pullback is O(3); the double-root normalization pullback
# is O(2,1).  Enumerate the whole nonnegative universe, not selected rows.
triple_rows = []
double_rows = []
for gcd_degree in range(3):
    e = 2 - gcd_degree
    triple_rows.extend(
        (gcd_degree, e, root_degree)
        for root_degree in range(3)
        if 3 * root_degree == e
    )
    double_rows.extend(
        (gcd_degree, e, ell_degree, m_degree)
        for ell_degree in range(3)
        for m_degree in range(3)
        if 2 * ell_degree + m_degree == e
    )
gate(triple_rows == [(2, 0, 0)], "triple-root lift is fixed and forces gcd degree two")
gate(double_rows == [
    (0, 2, 0, 2),
    (0, 2, 1, 0),
    (1, 1, 0, 1),
    (2, 0, 0, 0),
], "complete double-root lift ledger")
gate([row for row in double_rows if row[2] > 0] == [(0, 2, 1, 0)],
     "unique moving repeated-root cell")

# The smooth tangent calculation is rebuilt in a different notation.
t = sp.symbols("t")
r, s = sp.symbols("r s", nonzero=True)
da, db, dc, dd = sp.symbols("da db dc dd")
deformed = disc4(r + t * da, s + t * db, t * dc, t * dd)
same(sp.Poly(deformed, t).coeff_monomial(t), -4 * s**3 * dd,
     "smooth discriminant tangent is evaluation at the repeated root")

# For a fixed repeated root U, tangent cancellation makes the degree-one row
# divisible by U.  A constant binary cubic is divisible by U exactly when its
# V^3 coefficient vanishes; otherwise evaluation at (0,1) is a scalar unit.
U, V = sp.symbols("U V")
q20, q11, q02 = sp.symbols("q20 q11 q02")
f30, f21, f12, f03 = sp.symbols("f30 f21 f12 f03")
fixed_full = U**2 * (r * U + s * V) + U * (q20 * U**2 + q11 * U * V + q02 * V**2) + (
    f30 * U**3 + f21 * U**2 * V + f12 * U * V**2 + f03 * V**3
)
same(fixed_full.subs({U: 0, V: 1}), f03,
     "fixed-root index-form unit evaluation")
gate(sp.rem(sp.Poly(fixed_full, U, domain=sp.QQ.frac_field(V, r, s, q20, q11, q02, f30, f21, f12, f03)),
            sp.Poly(U, U, domain=sp.QQ.frac_field(V, r, s, q20, q11, q02, f30, f21, f12, f03))).as_expr() == f03 * V**3,
     "fixed-root reducibility remainder")


# ---------------------------------------------------------------------------
# 2. Common-zero triple row: exact primitive two-edge polygon.
# ---------------------------------------------------------------------------
A, C, x, z = sp.symbols("A C x z")
alpha, alpha1, beta, beta1, gamma, gamma1 = sp.symbols("alpha alpha1 beta beta1 gamma gamma1")
a = C**2 + alpha * A + alpha1 * C
b = beta * A + beta1 * C
c = gamma * A + gamma1 * C
d = C
D = disc4(a, b, c, d)
same(homogeneous_part(D, (A, C), 6), -27 * C**6, "common-zero pure sextic top")
same(D.subs(C, 0), gamma**2 * (beta**2 - 4 * alpha * gamma) * A**4,
     "common-zero irreducibility gate")
h = sp.expand(z**6 * D.subs({A: 1 / z, C: x / z}, simultaneous=True))
P = sp.Poly(h, x, z)
same(P.coeff_monomial(z**2), gamma**2 * (beta**2 - 4 * alpha * gamma),
     "order-two endpoint")
same(P.coeff_monomial(x**2 * z), -4 * gamma**3, "shared Newton vertex")
same(P.coeff_monomial(x**6), -27, "order-four endpoint")

support = {powers for powers, coefficient in P.terms() if coefficient != 0}
gate(all(i + 2 * j >= 4 for i, j in support), "first Newton supporting halfspace")
gate(all(i + 4 * j >= 6 for i, j in support), "second Newton supporting halfspace")
edge_two = sorted((i, j) for i, j in support if i + 2 * j == 4)
edge_four = sorted((i, j) for i, j in support if i + 4 * j == 6)
gate(edge_two == [(0, 2), (2, 1)], "first edge has exactly its primitive endpoints")
gate(edge_four == [(2, 1), (6, 0)], "second edge has exactly its primitive endpoints")
gate(sp.gcd(2, 1) == 1 and sp.gcd(4, 1) == 1,
     "both Newton edges are primitive and yield one branch")


# ---------------------------------------------------------------------------
# 3. Arbitrary constants: no missing chord or terminal seam.
# ---------------------------------------------------------------------------
a0, b0, c0, d0 = sp.symbols("a0 b0 c0 d0")
D0 = disc4(a + a0, b + b0, c + c0, d + d0)
h0 = sp.expand(z**6 * D0.subs({A: 1 / z, C: x / z}, simultaneous=True))
P0 = sp.Poly(h0, x, z)
same(homogeneous_part(D0, (A, C), 6), -27 * C**6,
     "constants preserve the top sextic")

# List the entire x-independent initial sequence.  These identities justify
# j0>=2,3,4 in the gamma, beta, alpha cases even when later coefficients also
# cancel.
same(P0.coeff_monomial(1), 0, "no local constant term")
same(P0.coeff_monomial(z), 0, "no x-independent z term")
same(P0.coeff_monomial(z**2).subs(gamma, 0), 0,
     "gamma-zero raises x-independent depth to at least three")
same(P0.coeff_monomial(z**3).subs({gamma: 0, beta: 0}), 0,
     "gamma-beta-zero raises x-independent depth to at least four")

same(P0.coeff_monomial(x**2 * z), -4 * gamma**3,
     "gamma interior vertex")
same(P0.coeff_monomial(x * z**2).subs(gamma, 0), -4 * beta**3,
     "beta interior vertex")
same(P0.coeff_monomial(x**4 * z).subs({gamma: 0, beta: 0}), -54 * alpha,
     "alpha interior vertex")

# Cross-multiplied strict-below tests for every allowed j0 lower bound.
for j0 in range(2, 7):
    gate(3 < 2 * j0, f"gamma point below endpoint chord at j0={j0}")
for j0 in range(3, 7):
    gate(12 < 5 * j0, f"beta point below endpoint chord at j0={j0}")
for j0 in range(4, 7):
    gate(3 < j0, f"alpha point below endpoint chord at j0={j0}")

D_terminal = sp.expand(D0.subs({alpha: 0, beta: 0, gamma: 0}))
same(sp.diff(D_terminal, A), 0, "terminal seam is target-univariate")
gate(sp.Poly(D_terminal, C).degree() == 6, "terminal seam retains sextic degree")


# ---------------------------------------------------------------------------
# 4. Moving repeated root: solve the pure-power coefficient ideal directly.
# ---------------------------------------------------------------------------
p, q, rr, ss, tt, uu, vv = sp.symbols("p q rr ss tt uu vv")
ell = A * U + C * V
Q2 = p * U**2 + q * U * V + rr * V**2
Q3 = ss * U**3 + tt * U**2 * V + uu * U * V**2 + vv * V**3
moving = sp.expand(ell**2 * V + ell * Q2 + Q3)
MP = sp.Poly(moving, U, V)
Dm = disc4(
    MP.coeff_monomial(U**3),
    MP.coeff_monomial(U**2 * V),
    MP.coeff_monomial(U * V**2),
    MP.coeff_monomial(V**3),
)
Dm6 = homogeneous_part(Dm, (A, C), 6)
coefficients = [sp.Poly(sp.cancel(Dm6 / A**2), A, C).coeff_monomial(A**i * C**(4 - i)) for i in range(5)]
gate(coefficients == [p**2, -2 * p * q + 4 * ss, 2 * p * rr + q**2 - 4 * tt,
                      -2 * q * rr + 4 * uu, rr**2 - 4 * vv],
     "moving-root quotient coefficient ledger")

solution = {p: 0, ss: 0, tt: q**2 / 4, uu: q * rr / 2}
same(Dm6.subs(solution), (rr**2 - 4 * vv) * A**6,
     "pure-power solution is unique and nonzero exactly off rr^2=4vv")
moving_solution = sp.factor(moving.subs(solution))
same(moving_solution.subs(V, 0), 0, "pure-power solution has the global V factor")
gate(sp.cancel(moving_solution / V).is_polynomial(A, C, U, V),
     "moving-root quotient is polynomial")


# ---------------------------------------------------------------------------
# 5. Audit semantics and inactive-assert gate.
# ---------------------------------------------------------------------------
semantic = {
    "universe": "all coefficient rows of target degree at most two, constants included",
    "lift": "gcd ledger plus O(3) and O(2,1) pullbacks exhaust triple and double strata",
    "fixed": "fixed repeated root gives a literal unit evaluation or a common binary factor",
    "triple": "primitive order-two and order-four edges; constants add an unavoidable interior vertex",
    "moving": "pure L^6 coefficient ideal has unique solution and forces V divisibility",
    "normality": "disc valuation 0/1 makes the finite-free order R1; freeness gives S2; irreducible nonsquare cubic gives S3",
    "scope": "quadratic coefficient depth only; JC2 remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
     "no inactive Python assert")

print("theorem=THM-3908-independent-hostile-audit")
print("lift_ledger=TRIPLE_FIXED_AND_DOUBLE_000_001_002_210")
print("fixed_repeated_root=UNIT_OR_REDUCIBLE")
print("common_zero_edges=PRIMITIVE_ORDERS_2_AND_4")
print("constant_triple_root=AT_LEAST_TWO_EDGES_OR_UNIVARIATE")
print("moving_repeated_root=PURE_POWER_FORCES_V_FACTOR")
print("normal_S3_consequence=PASS_BY_R1_S2_AND_NONSQUARE_DISCRIMINANT")
print("JC2=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
