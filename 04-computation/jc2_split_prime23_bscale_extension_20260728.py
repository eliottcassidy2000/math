#!/usr/bin/env python3
"""Exact B-scale extension of the split even-Faber prime-23 curve.

The original prime-23 scout normalizes the weight-two constant coefficient
``B_0`` to ``rho^2``.  This companion retains

    b = B_0/rho^2

as an arbitrary constant while keeping ``rho`` nonzero.  It verifies that
the generalized scaled fluxes differ from the normalized ``b=1`` fluxes by
explicit ``t^2`` terms only.  Therefore the complete ``t=0`` fibre, its old
and new local controls, the weighted corner, and the component divisor budget
are unchanged, including at ``b=0``.

This is an exact structural extension of the response-curve input.  The
mathematical lci/component/reducedness argument is recorded separately; this
file does not by itself close a Keller trajectory or the odd Faber bank.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_split_lambda_prime23_scout.py"
BASE_SHA256 = "1bbadb900e27112f57f600b0b5b73a3b85fc5b02f23e3f8f687dac4fa1c41fc3"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited split-lambda scout changed",
)

spec = importlib.util.spec_from_file_location("split_lambda_prime23_bscale", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load lambda scout")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)


# Retain the weight-two scale instead of setting it to one.
b = s.symbols("b")
scale = {
    base.B: b * base.rho**2,
    base.C: base.c * base.rho**3,
    base.D: base.dd * base.rho**4,
    base.E: base.e * base.rho**5,
    base.W: base.w * base.rho**6,
    base.y: base.rho / base.tau,
    base.u: base.rho**2 * base.v / base.tau**2,
    base.Z: base.rho**3 * base.zeta / base.tau**3,
}

f1_b = s.factor(s.cancel(base.N1.subs(scale) * base.tau**5 / base.rho**5))
f2_b = s.factor(s.cancel(base.N2.subs(scale) * base.tau**6 / base.rho**6))
require(
    base.rho not in f1_b.free_symbols and base.rho not in f2_b.free_symbols,
    "arbitrary nonzero rho did not cancel",
)

f1_delta = -616 * base.tau**2 * (b - 1) * (
    4840 * base.v - 1331 * base.zeta - 40
)
f2_delta = 49280 * base.tau**2 * (b - 1) * (
    29282 * base.v**2 - 1452 * base.v + 1331 * base.zeta + 2
)
require(s.expand(f1_b - base.f1 - f1_delta) == 0, "general f1 delta changed")
require(s.expand(f2_b - base.f2 - f2_delta) == 0, "general f2 delta changed")


# The prime-23 weighted relation is independent of the normalization b=1.
scaled_left = s.factor(s.cancel((base.Z * base.N1**4).subs(scale)))
require(
    s.factor(
        scaled_left
        - base.rho**23 * base.zeta * f1_b**4 / base.tau**23
    )
    == 0,
    "general b prime-23 relation changed",
)


# Every b-dependent term vanishes at t=0.  Hence all fixed local controls
# checked by the imported audited scout apply literally.
require(
    s.expand(f1_b.subs(base.tau, 0) - base.f1.subs(base.tau, 0)) == 0,
    "general b changed the fixed first flux",
)
require(
    s.expand(f2_b.subs(base.tau, 0) - base.f2.subs(base.tau, 0)) == 0,
    "general b changed the fixed second flux",
)
require(s.Poly(f1_b, base.zeta).degree() == 1, "general b f1 is not zeta-linear")
require(s.Poly(f2_b, base.zeta).degree() == 2, "general b f2 is not zeta-quadratic")


# Check the weighted corner directly.  Weighted homogenization with h=0
# retains only terms of weights five and six; then t=0 leaves these forms.
def weighted_part(expression: s.Expr, weight: int) -> s.Expr:
    polynomial = s.Poly(
        s.expand(expression),
        base.tau,
        base.v,
        base.zeta,
        b,
        base.c,
        base.dd,
        base.e,
        base.w,
    )
    result = s.Integer(0)
    for monomial, coefficient in polynomial.terms():
        t_power, v_power, zeta_power = monomial[:3]
        if t_power + 2 * v_power + 3 * zeta_power == weight:
            result += coefficient * s.prod(
                variable**power
                for variable, power in zip(polynomial.gens, monomial)
            )
    return s.expand(result)


corner_f1 = s.expand(weighted_part(f1_b, 5).subs(base.tau, 0))
corner_f2 = s.expand(weighted_part(f2_b, 6).subs(base.tau, 0))
require(corner_f1 == -1449459 * base.v * base.zeta, "corner f1 changed")
require(
    corner_f2 == 15944049 * base.zeta**2 - 1190488992 * base.v**3,
    "corner f2 changed",
)
require(
    s.expand((base.zeta * corner_f1**4).subs(base.v, 0)) == 0
    and s.expand((base.zeta * corner_f1**4).subs(base.zeta, 0)) == 0,
    "corner product support changed",
)
require(corner_f2.subs(base.v, 0) != 0, "zeta-coordinate corner entered f2")
require(corner_f2.subs(base.zeta, 0) != 0, "v-coordinate corner entered f2")


# The opposing section budgets are parameter-free and retain one survivor.
survivors = []
for old_count in range(6):
    for new_count in range(4):
        degree = 4 * old_count + new_count
        if degree and 23 * old_count <= 5 * degree and 23 * new_count <= 3 * degree:
            survivors.append((old_count, new_count, degree))
require(survivors == [(5, 3, 23)], "general b component budget changed")


print("split prime-23 arbitrary B-scale extension")
print(f"base_split_lambda_sha256={BASE_SHA256}")
print("scale=B0=b*rho^2:C0=c*rho^3:D0=d*rho^4:E0=e*rho^5:W0=w*rho^6")
print("rho=arbitrary_nonzero_constant:cancelled=True")
print("f1_b-f1_1=-616*t^2*(b-1)*(4840*v-1331*zeta-40)")
print("f2_b-f2_1=49280*t^2*(b-1)*(29282*v^2-1452*v+1331*zeta+2)")
print("prime23_relation=zeta*f1_b^4=(7496192^4*lambda^4/rho^23)*t^23")
print("fixed_fibre_independent_of_b=True:G3*L5^4")
print("fixed_local_controls_inherited_exactly=True")
print("weighted_corner_independent_of_b=True:empty=True")
print("ambient_coordinate_values_independent_of_b=-1190488992,15944049")
print("component_budget_survivor=(5,3,23)_only")
print("b_zero_included=True")
print("scope=EXACT_RESPONSE_CURVE_EXTENSION_NOT_ODD_BANK_NOT_JC2")
print("ALL CHECKS PASSED")
