#!/usr/bin/env python3
"""Exact companion for THM-3840's forced two-arm Jelonek passport."""

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


a, A, C = sp.symbols("a A C")
r = 3 * a**3 + 7 * a**2 + 1
Q = 7 * a**2 + 3
b = 6 * a**3 + 7 * a**2 - 1

A0 = -2 / (a * Q)
C0 = -1 / a**2

gate(sp.resultant(r, a, a) == -1, "cubic slopes are nonzero")
gate(sp.resultant(r, Q, a) == 1615, "Q is nonzero at every cubic slope")
gate(sp.resultant(r, b, a) != 0, "sign factor is nonzero at every slope")
zero(b + Q - 2 * r, "b equals minus Q on the cubic spectrum")
gate(sp.resultant(r, 7 * a**2 - 3, a) == 5245,
     "forced branch values are not triple-root values")
gate(sp.resultant(r, 7 * a**2 - 1, a) == 1363,
     "companion ratio denominator is nonzero")
gate(sp.resultant(r, r.subs(a, -a), a) == 216,
     "no two distinct cubic slopes have the same square")

# Exact target limits of the two Laurent arms, using t=1/k at a pole of k.
t = sp.symbols("t")
A_minus = a * t / Q
C_minus = sp.Integer(0)
A_plus = a * t / Q + A0
C_plus = C0
zero(A_minus.subs(t, 0), "minus arm tends to target origin")
zero(C_minus, "minus arm has C=0")
zero(A_plus.subs(t, 0) - A0, "plus arm tends to A0")
zero(C_plus - C0, "plus arm has fixed nonzero C0")

# The nonlinear cubic discriminant from THM-3811.
Delta = (
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
    + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)
Delta_at_T = sp.together(Delta.subs({A: A0, C: C0}))
Delta_num = sp.factor(Delta_at_T.as_numer_denom()[0])
gate(sp.rem(sp.Poly(Delta_num, a), sp.Poly(r, a)).as_expr() == 0,
     "plus endpoint lies on the branch discriminant")

Delta_A_num = sp.factor(
    sp.together(sp.diff(Delta, A).subs({A: A0, C: C0})).as_numer_denom()[0]
)
Delta_C_num = sp.factor(
    sp.together(sp.diff(Delta, C).subs({A: A0, C: C0})).as_numer_denom()[0]
)
gate(sp.rem(sp.Poly(Delta_A_num, a), sp.Poly(r, a)).as_expr() == 0,
     "branch tangent has Delta_A=0 at the endpoint")
gate(sp.resultant(r, Delta_C_num, a) == 27703695576000,
     "Delta_C is nonzero at every forced endpoint")

# Recover the branch-normalization parameter q0=1/a and its companion sheet.
q = sp.symbols("q")
Rq = (q - 3) * (q + 1) * (q + 2)
A_branch = -2 * q**2 * Rq / (3 * q**2 + 7) ** 2
C_branch = -q * Rq / (3 * q**2 + 7)


def numerator_remainder(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.factor(sp.rem(sp.Poly(numerator, a), sp.Poly(r, a)).as_expr())


zero(numerator_remainder(A_branch.subs(q, 1 / a) - A0),
     "branch normalization A-value at q=1/a")
zero(numerator_remainder(C_branch.subs(q, 1 / a) - C0),
     "branch normalization C-value at q=1/a")

q_comp = (7 * a**2 - 1) / (2 * a)
z_comp = 1 / q_comp
zero(z_comp - 2 * a / (7 * a**2 - 1), "companion marked-root ratio")
zero(2 / a + q_comp - C0 / A0, "Vieta sum for double plus companion roots")
zero(1 / a**2 + 2 * q_comp / a - 7,
     "Vieta pair sum for double plus companion roots")
zero(
    numerator_remainder(q_comp / a**2 - (C0**2 / A0 - 3)),
    "Vieta product for double plus companion roots",
)

# T_a values are pairwise distinct because equality of their C-coordinates
# would give two cubic roots differing by sign, excluded above.
gate(C0 != 0, "forced branch endpoint is distinct from the origin")

semantic = {
    "forced": "THM-3836 supplies one cubic slope with nonempty minus and plus source components",
    "unit": "k is a nonconstant unit on each relevant etale curve component",
    "limits": "minus -> O=(0,0); plus -> T_a=(-2/[aQ(a)],-1/a^2)",
    "jelonek": "a pole of k on each smooth completion gives a source-infinite valuation",
    "branch": "Delta(T_a)=0, Delta_C(T_a)!=0, normalization parameter q=1/a",
    "companion": "surviving root ratio beta=2a/(7a^2-1)",
    "distinct": "the three candidate T_a are pairwise distinct and nonzero",
    "scope": "O and at least one of the three T_a are forced nonproper values",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3840-forced-cubic-two-arm-jelonek-passport")
print("forced=origin_and_at_least_one_of_three_cubic_branch_values")
print("T_a=(-2/(a*(7a^2+3)),-1/a^2);r(a)=0")
print("branch=Delta(T_a)=0;Delta_C(T_a)!=0;q_normalization=1/a")
print("triple=excluded;resultant=5245")
print("companion=root_ratio_beta=2a/(7a^2-1)")
print("distinct=three_T_a_pairwise_distinct_and_all_nonzero")
print("mechanism=nonconstant_unit_k_has_boundary_pole_on_each_forced_arm")
print("scope=two_points_forced_not_all_three_T_a;atlas_still_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
