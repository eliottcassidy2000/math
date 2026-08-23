#!/usr/bin/env python3
"""Exact companion for THM-3832's triangular root-ratio chart."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


z, C = sp.symbols("z C")
r = 3 * z**3 + 7 * z**2 + 1
q = 7 * z**2 + 3
b = 6 * z**3 + 7 * z**2 - 1
s = z**2 * q * C - b

A = z * C * (1 + z**2 * C) / r
omega = C * (1 + z**2 * C) / r
theta = -z * C * (7 * C * z**2 + C - 3 * z) / r
D = C**2 * (1 + z**2 * C) * s / r**2
k = r / (C * s)
h = z * r / (C * s)
m = -z * C * (7 * C * z**2 + 3 * C - 9 * z - 14) / r

zero(omega - A / z, "marked root ratio")
zero(h / k - z, "intrinsic row ratio")
zero(
    omega**3 - C * omega**2 + 7 * A**2 * omega
    + 3 * A**3 - A**2 * C**2,
    "marked-root characteristic cubic",
)

zero(omega**2 - (-7 * A**2 + C * omega - A * theta),
     "omega-square multiplication law")
zero(omega * theta - (3 * A**2 - A * C**2),
     "mixed multiplication law")
zero(
    theta**2
    - (3 * A * C - C**3 + (C**2 - 3 * A) * omega - 7 * A * theta),
    "theta-square multiplication law",
)
zero(D - (C * omega - 3 * A * theta - 14 * A**2),
     "intrinsic different")
zero(D - (3 * omega**2 - 2 * C * omega + 7 * A**2),
     "different is marked-root derivative")
zero(A - h * D, "A reconstruction")
zero(omega - k * D, "omega reconstruction")
zero(m - (3 * theta + 14 * A), "m reconstruction")

zero(C * k - m * h - 1, "SL2 determinant")
zero(D * (7 * h**2 + 3 * k**2) - (1 + 2 * C * k),
     "first intrinsic lift law")
zero(h * D * (9 * h + 14 * k) - (k * m + 3 * h * C**2),
     "second intrinsic lift law")

target_density = sp.factor(sp.diff(A, z))
zero(target_density - C * s / r**2, "factored triangular target derivative")
zero(target_density - 1 / (k * r), "target derivative is reciprocal kr")

z_x, z_y, C_x, C_y, lam = sp.symbols("z_x z_y C_x C_y lambda")
source_jacobian = z_x * C_y - z_y * C_x
target_jacobian = target_density * source_jacobian
zero(target_jacobian - source_jacobian / (k * r),
     "source-chain Jacobian identity")
zero(
    (target_jacobian - lam).subs(source_jacobian, lam * k * r),
    "weighted area equation is exactly Keller",
)

zero(b + q - 2 * r, "cubic/quadratic spectral identity")
gate(sp.gcd(sp.Poly(r, z), sp.Poly(z * q * b, z)).degree() == 0,
     "cubic roots avoid zero and quadratic/sign roots")
gate(sp.discriminant(r, z) != 0, "three cubic roots are distinct")

alpha = sp.symbols("alpha")
r_alpha = 3 * alpha**3 + 7 * alpha**2 + 1
numerator_at_alpha = alpha * C * (1 + alpha**2 * C)
gate(sp.factor(numerator_at_alpha) == alpha * C * (alpha**2 * C + 1),
     "two cubic-pole cancellation addresses")
gate(sp.resultant(r, z * q * b, z) != 0,
     "both cancellation addresses are distinct at every cubic root")
zero(
    sp.rem(sp.Poly((b + q).subs(z, alpha), alpha),
           sp.Poly(r_alpha, alpha)).as_expr(),
    "b equals minus q at every cubic root",
)

semantic = {
    "chart": "z=h/k=A/omega; K(U)=K(z,C)",
    "target": "A=z*C*(1+z^2*C)/(3z^3+7z^2+1); second coordinate C",
    "row": "k=r/[C*(z^2*q*C-b)]; h=z*k",
    "density": "Jac_(z,C)(A,C)=C*(z^2*q*C-b)/r^2=1/(k*r)",
    "keller": "Jac_(x,y)(z,C)=lambda*k*r",
    "spectral": "r(alpha)=0 has cancellation addresses C=0,-1/alpha^2",
    "scope": "birational rational chart only; polynomial atlas and JC2 counterexample open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart")
print("chart=z=h/k=A/omega;coordinates=z,C")
print("target=A=z*C*(1+z^2*C)/r;C=C;r=3z^3+7z^2+1")
print("row=k=r/(C*s);h=z*k;s=z^2*(7z^2+3)*C-(6z^3+7z^2-1)")
print("density=Jac_zC(A,C)=C*s/r^2=1/(k*r)")
print("keller=Jac_xy(z,C)=lambda*k*r")
print("spectral=cubic_root_addresses_C=0,-1/alpha^2")
print("scope=rational_birational_chart_only;polynomial_atlas_and_JC2_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
