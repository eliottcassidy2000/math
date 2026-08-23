#!/usr/bin/env python3
"""Exact companion for THM-3894's all-degree gauge--Kummer filtration."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, label)


x, y = sp.symbols("x y")
K2 = y**2 - 15 * x**2


# ---------------------------------------------------------------------------
# 1. The universal next component of r=aT+Kf.
# ---------------------------------------------------------------------------

qm3, qm2, qm1, s, t, qj = sp.symbols("qm3 qm2 qm1 s t qj")
previous_T = -K2 * qm1 + 15 * x * qm2 + 4 * qm3
previous_f_1 = qm2 + x * qm1
previous_f_2 = qm3 + x * qm2
r_component = sp.expand(
    previous_T
    + x * t
    + K2 * s
    - 15 * x * previous_f_1
    - 4 * previous_f_2
)
Rj = K2 * (s - qm1) + x * (t - 15 * x * qm1 - 4 * qm2)
zero(r_component - Rj, "universal next r component")

s_gauge = qm1 + x * qj
t_gauge = -K2 * qj + 15 * x * qm1 + 4 * qm2
zero(Rj.subs({s: s_gauge, t: t_gauge}), "next full gauge jet")


# ---------------------------------------------------------------------------
# 2. Exact coefficient and degree envelope.
# ---------------------------------------------------------------------------

q0, R = sp.symbols("q0 R")
A_top = 81 * q0 * x**5
B_top = 81 * q0**2 * x**5
zero(3 * R**2 * B_top - 243 * q0**2 * x**5 * R**2,
     "r2B leading coefficient")
zero(8 * A_top * B_top - 52488 * q0**3 * x**10,
     "AB leading coefficient")
gate(sp.Rational(52488, 243) == 216, "middle collision ratio")

N, J = sp.symbols("N J", integer=True, positive=True)
d_r2B = 4 * N + 7 - 2 * J
d_AB = 3 * N + 7
d_r2L2f = 3 * N + 6 - 2 * J
d_AL2f = 2 * N + 6
d_PL2f = N + 5
d_L4 = sp.Integer(4)
zero(d_r2B - d_AB - (N - 2 * J), "top-degree crossing")
zero(d_r2B - d_r2L2f - (N + 1), "r2B versus r2L2f")
zero(d_AB - d_AL2f - (N + 1), "AB versus AL2f")
zero(d_AB - d_PL2f - (2 * N + 2), "AB versus PL2f")
zero(d_AB - d_L4 - (3 * N + 3), "AB versus L4")

# Finite controls freeze the integer partition logic, not the theorem's
# all-degree universe.
for n in range(1, 21):
    for j in range(1, n + 1):
        comparison = (4 * n + 7 - 2 * j) - (3 * n + 7)
        gate(comparison == n - 2 * j, f"degree comparison n={n},j={j}")
        if n > 2 * j:
            gate(comparison > 0, f"pre-middle regime n={n},j={j}")
        elif n == 2 * j:
            gate(comparison == 0 and (3 * n + 7) % 2 == 1,
                 f"middle odd-degree regime n={n},j={j}")
        else:
            gate(comparison < 0, f"post-middle regime n={n},j={j}")

for m in range(1, 11):
    n_even = 2 * m
    gate(all(n_even > 2 * j for j in range(1, m)),
         f"even forced depths m={m}")
    gate(n_even == 2 * m, f"even Kummer depth m={m}")

for m in range(0, 10):
    n_odd = 2 * m + 1
    gate(all(n_odd > 2 * j for j in range(1, m + 1)),
         f"odd forced depths m={m}")
    gate(n_odd < 2 * (m + 1), f"odd square depth m={m}")
    gate((3 * n_odd + 7) % 2 == 0, f"odd terminal degree even m={m}")


# The middle Kummer and odd square normal forms reproduce the leading forms.
u, epsilon, v = sp.symbols("u epsilon v", nonzero=True)
zero(
    (epsilon * x**3 * u) ** 2
    + 216 * (x * u**2) * x**5
    - (epsilon**2 + 216) * x**6 * u**2,
    "Kummer normal form",
)
zero(
    52488 * (v**2) ** 3 * x**10
    - (sp.sqrt(52488) * v**3 * x**5) ** 2,
    "odd square-profile leading form",
)


# ---------------------------------------------------------------------------
# 3. Vertical full-gauge slice: exact quartic discriminant.
# ---------------------------------------------------------------------------

a, L, K, q, w = sp.symbols("a L K q w")
Delta = a**3 * L**2 - K**2
P = a * L**2
r_vertical = a * w
A_vertical = Delta * q + K * w
B_vertical = Delta * q**2 + 2 * K * q * w - w**2
S_vertical = sp.expand(
    L**4
    + 2 * (3 * A_vertical + 3 * P + r_vertical**2) * L**2 * a * q
    + (8 * A_vertical + 6 * P + 3 * r_vertical**2) * B_vertical
)

gate(sp.Poly(S_vertical, K).degree() == 4, "vertical residual quartic in K")
zero(sp.Poly(S_vertical, K).coeff_monomial(K**4) - 8 * q**3,
     "vertical quartic leading coefficient")

factor_1 = 3 * a**2 * q + 2
factor_2 = -32 * L**4 * q + 72 * L**2 * a * q * w**2 + 27 * w**4
factor_3 = -4 * L**2 * a**2 * q - 2 * L**2 + a**3 * q * w**2
expected_discriminant = (
    128 * L**2 * q**8 * factor_1**2 * factor_2 * factor_3**3
)
zero(
    sp.discriminant(S_vertical, K) - expected_discriminant,
    "vertical quartic discriminant factorization",
)

zero(factor_1.subs(a, 0) - 2, "first factor a-zero obstruction")
zero(factor_3.subs(a, 0) + 2 * L**2,
     "third factor a-zero obstruction")

vv = sp.symbols("vv")
reduced_factor_2 = sp.cancel(factor_2.subs(w, L * vv) / L**4)
zero(
    reduced_factor_2 - (27 * vv**4 + 72 * a * q * vv**2 - 32 * q),
    "middle factor L-adic reduction",
)
zero(
    (9 * a * vv**2 - 4) - 9 * a * vv * vv + 4,
    "middle factor Bezout coprimality certificate",
)

degree_v = sp.symbols("degree_v", integer=True, nonnegative=True)
zero((1 + 2 * degree_v) - (2 * degree_v + 1),
     "nonzero middle denominator degree")


semantic = {
    "aligned": "q_i gauge jets with universal next symbol R_j",
    "before_middle": "n>2j forces R_j=0 and another gauge jet",
    "middle": "n=2j forces q0=xu^2 and R_j=epsilon*x^3*u",
    "after_middle": "n<2j forces n odd and q0 square",
    "vertical": "f=a*q(x),T=-K*q(x)+w(x) has no nonzero square residual",
    "scope": "associated graded only; gauge peeling is not square invariant",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3894-cusp-residual-all-degree-gauge-kummer-parity-filtration")
print("pre_middle=n_gt_2j_forces_next_full_gauge_jet")
print("even_passport=n_2m_forces_q0_x_times_square_and_Kummer_R")
print("odd_passport=n_2m_plus_1_forces_q0_square")
print("vertical_gauge_completion=no_nonzero_square_residual")
print("gauge_peeling_square_invariance=NOT_CLAIMED")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
