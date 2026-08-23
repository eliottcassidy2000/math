#!/usr/bin/env python3
"""Exact companion for THM-3884's residual total-degree filtration."""

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
    gate(sp.expand(expression) == 0, message)


x, y, T, F = sp.symbols("x y T F")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
P = a * L**2
delta = a**3 * L**2 - K**2
r = a * T + K * F
A = K * T + a * P * F

# THM-3881's exact residual, now regarded as a polynomial in the two module
# coordinates T,F with coefficients in k[x,y].
residual = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * F
    + (8 * A + 6 * P + 3 * r**2) * (P * F**2 - T**2)
)
poly = sp.Poly(residual, T, F)


def coefficient_degree(coefficient: sp.Expr) -> int:
    return int(sp.Poly(coefficient, x, y).total_degree())


def leading_form(coefficient: sp.Expr) -> sp.Expr:
    p = sp.Poly(coefficient, x, y)
    degree = int(p.total_degree())
    return sp.Add(*[
        coeff * x**monomial[0] * y**monomial[1]
        for monomial, coeff in p.terms()
        if sum(monomial) == degree
    ])


# Freeze every nonzero (T-degree,F-degree,coefficient-degree) cell.  This is
# the exact universe used in the all-degree inequalities in the proof.
degree_cells = {
    monomial: coefficient_degree(coefficient)
    for monomial, coefficient in poly.terms()
}
expected_cells = {
    (4, 0): 2,
    (3, 1): 3,
    (3, 0): 2,
    (2, 2): 5,
    (2, 1): 4,
    (2, 0): 3,
    (1, 3): 6,
    (1, 2): 5,
    (1, 1): 4,
    (0, 4): 7,
    (0, 3): 7,
    (0, 2): 6,
    (0, 1): 5,
    (0, 0): 4,
}
gate(degree_cells == expected_cells, "complete coefficient-degree table")
zero(poly.coeff_monomial(F**4) - 3 * a * K**2 * L**2,
     "exact F4 coefficient")
zero(poly.coeff_monomial(T * F**3) - 6 * K * a**2 * L**2,
     "exact TF3 coefficient")
zero(poly.coeff_monomial(T**2 * F**2) - 3 * delta,
     "exact T2F2 coefficient")
zero(poly.coeff_monomial(T**4) + 3 * a**2,
     "exact T4 coefficient")

# If m<=n and n>=1, (0,4) is uniquely maximal.  Replacing m by n gives the
# largest possible degree for every competing cell; the difference is affine
# in n with nonnegative slope, so checking n=1 proves the all-n inequality.
for (i, j), coefficient_degree_value in degree_cells.items():
    if (i, j) == (0, 4):
        continue
    slope = 4 - i - j
    intercept = 7 - coefficient_degree_value
    gate(slope >= 0, f"nonnegative all-n slope {(i, j)}")
    gate(slope + intercept >= 1, f"strict m<=n gap {(i, j)}")

top_f4 = poly.coeff_monomial(F**4)
zero(
    leading_form(top_f4) - 243 * x**3 * (y**2 - 15 * x**2) ** 2,
    "m<=n unique leading coefficient",
)

# On m=n+1, exactly the F^4, T F^3, and T^2 F^2 cells reach degree 4n+7.
equality_top_cells = set()
for (i, j), coefficient_degree_value in degree_cells.items():
    slope = 4 - i - j
    intercept = 7 - (i + coefficient_degree_value)
    # Difference (4n+7)-[i(n+1)+jn+c] = slope*n+intercept.
    if slope == 0 and intercept == 0:
        equality_top_cells.add((i, j))
    else:
        gate(slope >= 0, f"equality nonnegative slope {(i, j)}")
        gate(slope + intercept >= 1, f"strict equality gap {(i, j)}")
gate(
    equality_top_cells == {(0, 4), (1, 3), (2, 2)},
    "complete equality top-cell packet",
)

Tm, Fn = sp.symbols("Tm Fn")
equality_top = sp.expand(
    leading_form(poly.coeff_monomial(F**4)) * Fn**4
    + leading_form(poly.coeff_monomial(T * F**3)) * Tm * Fn**3
    + leading_form(poly.coeff_monomial(T**2 * F**2)) * Tm**2 * Fn**2
)
K2 = y**2 - 15 * x**2
zero(
    equality_top - 243 * x**3 * Fn**2 * (K2 * Fn + x * Tm) ** 2,
    "equality leading-form factorization",
)
gate(sp.gcd(x, K2) == 1, "x coprime to K2")

# The seam equation K2*Fn+x*Tm=0 is exactly the leading part of the
# THM-3881 gauge vector (K,-a).  The full gauge also preserves the address.
q = sp.symbols("q")
Fn_seam = x * q
Tm_seam = -K2 * q
zero(K2 * Fn_seam + x * Tm_seam, "leading seam relation")
zero(Tm_seam + K2 * q, "gauge cancels leading T form")
zero(Fn_seam - x * q, "gauge cancels leading F form")
zero(K.subs({x: 0, y: 0}) * q - 4 * (-a.subs({x: 0, y: 0}) * q),
     "gauge vector preserves address")

semantic = {
    "universe": "THM3881 residual in k[x,y][T,f]",
    "degree_floor": "nonconstant f square survivor has degT>=degf+1",
    "equality_top": "243*x^3*f_n^2*(K2*f_n+x*T_m)^2",
    "equality_seam": "f_n=x*q and T_m=-K2*q",
    "interpretation": "leading gauge cancellation; square invariance and termination not claimed",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3884-cusp-residual-total-degree-leading-gauge-filtration")
print("coefficient_degree_cells=14")
print("nonconstant_f_square_survivor_requires=degT_at_least_degf_plus_1")
print("equality_top=243*x^3*f_n^2*(K2*f_n+x*T_m)^2")
print("equality_survivor_requires=f_n=x*q;T_m=-K2*q")
print("leading_direction=THM3881_gauge_(K,-a)")
print("gauge_square_invariance_or_termination=NOT_CLAIMED")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
