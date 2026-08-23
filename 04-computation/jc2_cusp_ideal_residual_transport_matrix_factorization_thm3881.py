#!/usr/bin/env python3
"""Exact companion for THM-3881's cusp-ideal residual transport packet."""

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


# Universal residual transport.  No quotient or branch-ring assumption is
# used in this identity.
P0, Q0, A0, r0 = sp.symbols("P0 Q0 A0 r0")
P1 = P0 + 2 * A0 + r0**2
Q1 = Q0 + 3 * r0 * P0 + 3 * r0 * A0 + r0**3
transport = (
    P0**3
    - Q0**2
    + 2 * (3 * A0 + 3 * P0 + r0**2) * (A0 * P0 - r0 * Q0)
    + (8 * A0 + 6 * P0 + 3 * r0**2) * (A0**2 - r0**2 * P0)
)
zero(sp.expand(P1**3 - Q1**2 - transport), "universal residual transport")


# The THM-3869/3872 positive representative and its norm-form branch law.
x, y = sp.symbols("x y")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
P = a * L**2
Q = K * L**2
delta = sp.expand(
    81 * x**5
    + 90 * x**4
    + 25 * x**3
    + 30 * x**2 * y**2
    + 30 * x * y**2
    - y**4
    + 8 * y**2
)
zero(a**3 * L**2 - K**2 - delta, "branch norm form")
zero(P**3 - Q**2 - delta * L**4, "base square residual")


# Exact cusp-value-zero ideal and the three minimum mixed lifts from THM-3872.
j1 = a * x
j2 = a * y
j3 = y**2 + 4 * x
A1 = x * K
A2 = y * K
A3 = (15 * x + 4) * K + a * P
zero(j3 - (a * (15 * x + 4) + K), "third ideal generator split")
zero(
    A3
    - (
        81 * x**4
        + 9 * x**3
        - 44 * x**2
        + 15 * x * y**2
        - 16 * x
        + 4 * y**2
    ),
    "third mixed lift compact form",
)
j_groebner = sp.groebner([j1, j2, j3], y, x, order="lex")
expected_j_groebner = sp.groebner(
    [4 * x + y**2, x * y + y, x**2 + x], y, x, order="lex"
)
gate(j_groebner == expected_j_groebner, "cusp ideal Groebner presentation")
zero(j_groebner.reduce(delta)[1], "branch equation belongs to cusp ideal")


# The two Hilbert--Burch generators.  One is a genuine kernel; the other is
# exactly a delta-valued change of mixed lift.
zero(-y * j1 + x * j2, "first ideal syzygy")
zero(-y * A1 + x * A2, "first mixed syzygy")
zero(4 * j1 + y * j2 - a * j3, "second ideal syzygy")
zero(4 * A1 + y * A2 - a * A3 + delta, "delta-lift syzygy")
hilbert_burch = sp.Matrix([[-y, 4], [x, y], [0, -a]])
minors = [
    hilbert_burch.extract(rows, [0, 1]).det()
    for rows in ((1, 2), (0, 2), (0, 1))
]
gate(
    [sp.expand(minor) for minor in minors]
    == [sp.expand(-j1), sp.expand(j2), sp.expand(-j3)],
    "Hilbert--Burch signed minors",
)


# Every coefficient presentation contracts to two polynomial coordinates.
f1, f2, f3 = sp.symbols("f1 f2 f3")
f = f3
T = x * f1 + y * f2 + (15 * x + 4) * f3
r = f1 * j1 + f2 * j2 + f3 * j3
A = f1 * A1 + f2 * A2 + f3 * A3
zero(r - (a * T + K * f), "rank-two r coordinate")
zero(A - (K * T + a * P * f), "rank-two A coordinate")
zero(T.subs({x: 0, y: 0}) - 4 * f, "address condition")

matrix = sp.Matrix([[a, K], [K, a * P]])
adjugate = sp.Matrix([[a * P, -K], [-K, a]])
zero(matrix.det() - delta, "matrix determinant")
gate(
    (matrix * adjugate - delta * sp.eye(2)).applyfunc(sp.expand) == sp.zeros(2),
    "left matrix factorization",
)
gate(
    (adjugate * matrix - delta * sp.eye(2)).applyfunc(sp.expand) == sp.zeros(2),
    "right matrix factorization",
)


# The two defect sidecars become a linear coordinate and a binary norm.
C_numerator = sp.expand(A * P - r * Q)
B_numerator = sp.expand(A**2 - r**2 * P)
zero(C_numerator - delta * L**2 * f, "linear sidecar")
zero(B_numerator - delta * (P * f**2 - T**2), "quadratic norm sidecar")


# Substitute the sidecars into the universal transport law.  This is the
# exact quartic square equation left by arbitrary polynomial coefficients.
residual = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * f
    + (8 * A + 6 * P + 3 * r**2) * (P * f**2 - T**2)
)
zero(
    sp.expand(
        (P + 2 * A + r**2) ** 3
        - (Q + 3 * r * P + 3 * r * A + r**3) ** 2
        - delta * residual
    ),
    "specialized transported residual",
)
zero(residual.subs({f1: 0, f2: 0, f3: 0}) - L**4,
     "base-point positive control")


# The THM-3872 constant span embeds as f=w and
# T=ux+vy+w(15x+4), without losing any parameter.
u, v, w = sp.symbols("u v w")
constant_T = u * x + v * y + w * (15 * x + 4)
zero(T.subs({f1: u, f2: v, f3: w}) - constant_T,
     "constant-span T slice")
zero(f.subs(f3, w) - w, "constant-span f slice")


# Gauge action: add q(4,y,-a) to a coefficient presentation.  It fixes r,
# sends A to A-delta*q, and acts on the pair by (T,f)+(K,-a)q.
q = sp.symbols("q")
T_gauge = x * (f1 + 4 * q) + y * (f2 + y * q) + (15 * x + 4) * (f3 - a * q)
f_gauge = f3 - a * q
zero(T_gauge - (T + K * q), "pair gauge T")
zero(f_gauge - (f - a * q), "pair gauge f")
zero(a * T_gauge + K * f_gauge - r, "gauge fixes r")
zero(K * T_gauge + a * P * f_gauge - (A - delta * q),
     "gauge shifts mixed lift")


# Pure lift gauge through the positive base point.  These identities feed the
# theorem's Mason arm obstruction and the independent L=0 square-class gate.
pure_gauge_residual = sp.expand(
    L**4 - 6 * P**2 * q + 12 * delta * P * q**2 - 8 * delta**2 * q**3
)
zero(
    residual.subs({f1: 4 * q, f2: y * q, f3: -a * q})
    - pure_gauge_residual,
    "pure gauge residual",
)
q0 = sp.symbols("q0")
zero(
    pure_gauge_residual.subs({x: -1, q: q0})
    - (625 - 8 * (y**2 - 4) ** 4 * q0**3),
    "a-zero Mason boundary",
)
qL = sp.symbols("qL")
zero(
    pure_gauge_residual.subs({x: -sp.Rational(4, 9), q: qL})
    + 8 * (y**2 - sp.Rational(8, 27)) ** 4 * qL**3,
    "L-zero square-class boundary",
)


# The entire T=0 lane has a branch-factor form.  The theorem combines this
# identity with the normalization UFD and a constant Pell factorization.
F = sp.symbols("F")
r_F = K * F
A_F = a * P * F
T_zero_residual = sp.expand(
    L**4
    + 2 * (3 * A_F + 3 * P + r_F**2) * L**2 * F
    + (8 * A_F + 6 * P + 3 * r_F**2) * P * F**2
)
T_zero_inner = sp.expand(
    L**2 * (1 + a * F) ** 3 * (1 + 3 * a * F)
    - delta * F**3 * (2 + 3 * a * F)
)
zero(T_zero_residual - L**2 * T_zero_inner,
     "T-zero branch-factor residual")
zero(
    (P + 2 * A_F + r_F**2) ** 3
    - (Q + 3 * r_F * P + 3 * r_F * A_F + r_F**3) ** 2
    - delta * T_zero_residual,
    "T-zero direct residual reconstruction",
)
z, U, V = sp.symbols("z U V")
zero(3 * (1 + z) - (1 + 3 * z) - 2,
     "T-zero coprime-factor Bezout constant")
zero((V - sp.sqrt(3) * U) * (V + sp.sqrt(3) * U) - (V**2 - 3 * U**2),
     "T-zero Pell factorization")
zero(delta.subs(x, -1) + (y**2 - 4) ** 2,
     "delta is not divisible by a")
t_norm = sp.symbols("t_norm")
x_norm = t_norm**4 - 2 * t_norm**2
y_norm = 3 * t_norm**5 - 5 * t_norm**3
zero(delta.subs({x: x_norm, y: y_norm}), "normalization lies on delta")
gate(
    (x_norm.subs(t_norm, 0), y_norm.subs(t_norm, 0), a.subs(x, x_norm).subs(t_norm, 0))
    == (0, 0, 1),
    "T-zero address fixes z at zero",
)


semantic = {
    "universal": "two sidecars transport every translated cubic residual",
    "branch": "delta=(x+1)^3(9x+4)^2-(y2-15x2-15x-4)^2",
    "module": "r=aT+Kf;A=KT+aPf;T(0,0)=4f(0,0)",
    "matrix": "[[a,K],[K,aP]] has determinant delta",
    "sidecars": "C=L2f;B=Pf2-T2",
    "gauge": "(T,f)->(T+Kq,f-aq);r fixed;A->A-delta*q",
    "T_zero": "square survivor forces delta divides f by branch Pell factorization",
    "scope": "arbitrary polynomial T,f quartic square equation remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization")
print("universal_transport_sidecars=2")
print("cusp_ideal_presentation_rank=3;syzygy_rank=2")
print("contracted_pair_rank=2;matrix_determinant=delta")
print("sidecars=C=L^2*f;B=P*f^2-T^2")
print("delta_lift_gauge=(T,f)+(K,-a)*q")
print("pure_gauge_survivor_requires=a_divides_q")
print("T_zero_square_survivor_requires=delta_divides_f")
print("arbitrary_polynomial_coefficient_square_equation=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
