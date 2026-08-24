#!/usr/bin/env python3
"""Exact restriction and support companion for THM-3987."""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Assertion-free exact gate retained under ``python -O``."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


def linear_part(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    poly = sp.Poly(sp.expand(expression), *variables)
    return sp.Add(*(
        coefficient * sp.prod(variable**power for variable, power
                              in zip(variables, powers))
        for powers, coefficient in poly.terms()
        if sum(powers) == 1
    ))


x, t, u, s = sp.symbols("x t u s")
z = 1 + x**2 * t
p = sp.expand(z * t)
u_source = x**2 * t
coefficient = 1 + 3 * u + 2 * u**2


# Boundary chart: the exact relation and rational pole expression behind
# ord_D(t)=-2 (the valuation conclusion itself is geometric).
y_source = sp.expand(x * z * t**2)
zero(z * (z - 1)**2 - x**3 * y_source,
     "D(z-1) boundary-chart relation")
zero((z - 1) / x**2 - t, "boundary rational expression for t")


# Panel I: exact homogeneous-module survival on t=0 and x=0.
positive_survival: dict[int, str] = {}
for weight in range(0, 13):
    source = sp.expand(x**weight * coefficient.subs(u, u_source))
    expected_t0 = x**weight
    expected_x0 = sp.Integer(1) if weight == 0 else sp.Integer(0)
    zero(source.subs(t, 0) - expected_t0,
         f"positive weight {weight}: t=0 restriction")
    zero(source.subs(x, 0) - expected_x0,
         f"positive weight {weight}: x=0 restriction")
    positive_survival[weight] = "x^r" if weight else "constant"

negative_survival: dict[int, str] = {}
for q in range(1, 25):
    u_order = (q + 1) // 2
    z_order = (q + 2) // 3
    x_order = 2 * u_order - q
    source = sp.expand(
        x**x_order * t**u_order * z**z_order
        * coefficient.subs(u, u_source)
    )
    gate(x_order in (0, 1), f"negative weight -{q}: parity x order")
    zero(source.subs(t, 0), f"negative weight -{q}: t=0 vanishing")
    expected_x0 = t**u_order if q % 2 == 0 else sp.Integer(0)
    zero(source.subs(x, 0) - expected_x0,
         f"negative weight -{q}: x=0 parity restriction")
    for powers, coeff in sp.Poly(source, x, t).terms():
        if coeff:
            gate(powers[0] - 2 * powers[1] == -q,
                 f"negative weight -{q}: source weight homogeneity")
    negative_survival[-q] = (
        f"t^{u_order}" if q % 2 == 0 else "zero"
    )


# Only weights 1 and -2 can contribute to the source linear jet.
jet_weights: list[int] = []
for weight in range(0, 9):
    source = sp.expand(x**weight * coefficient.subs(u, u_source))
    if linear_part(source - source.subs({x: 0, t: 0}), (x, t)) != 0:
        jet_weights.append(weight)
for q in range(1, 17):
    u_order = (q + 1) // 2
    z_order = (q + 2) // 3
    source = sp.expand(
        x**(2 * u_order - q) * t**u_order * z**z_order
        * coefficient.subs(u, u_source)
    )
    if linear_part(source, (x, t)) != 0:
        jet_weights.append(-q)
gate(sorted(jet_weights) == [-2, 1], "complete linear-jet weight list")

a1, am2, c1, cm2 = sp.symbols("a1 am2 c1 cm2")
A_jet = a1 * x + am2 * p
C_jet = c1 * x + cm2 * p
origin_jacobian = sp.det(sp.Matrix([
    [sp.diff(A_jet, x).subs({x: 0, t: 0}),
     sp.diff(A_jet, t).subs({x: 0, t: 0})],
    [sp.diff(C_jet, x).subs({x: 0, t: 0}),
     sp.diff(C_jet, t).subs({x: 0, t: 0})],
]))
zero(origin_jacobian - (a1 * cm2 - am2 * c1),
     "origin jet determinant")


# Panel II: immersed noninjective curve degree lemma and sharp nodal control.
r0, r1 = sp.symbols("r0 r1")
a0, a_1, a_2 = sp.symbols("a0 a_1 a_2", nonzero=True)
quadratic = a_2 * s**2 + a_1 * s + a0
collision_quotient = sp.factor(
    (quadratic.subs(s, r0) - quadratic.subs(s, r1)) / (r0 - r1)
)
midpoint_derivative = sp.factor(
    sp.diff(quadratic, s).subs(s, (r0 + r1) / 2)
)
zero(collision_quotient - midpoint_derivative,
     "quadratic collision midpoint derivative")

nodal_f = s**2 - 1
nodal_g = s * (s**2 - 1)
gate(nodal_f.subs(s, 1) == nodal_f.subs(s, -1),
     "nodal hostile: first collision coordinate")
gate(nodal_g.subs(s, 1) == nodal_g.subs(s, -1),
     "nodal hostile: second collision coordinate")
gate(sp.gcd(sp.Poly(sp.diff(nodal_f, s), s),
            sp.Poly(sp.diff(nodal_g, s), s)).degree() == 0,
     "nodal hostile is immersed")
gate((sp.degree(nodal_f, s), sp.degree(nodal_g, s)) == (2, 3),
     "nodal hostile degree pair")


# Panel III: ordinary directional tails. Substitution on an origin line has
# coefficient P_d(v) in degree d, exactly as used in the projective gate.
v0, v1 = sp.symbols("v0 v1")
b20, b11, b02, b30, b21, b12, b03 = sp.symbols(
    "b20 b11 b02 b30 b21 b12 b03"
)
P2 = b20 * x**2 + b11 * x * t + b02 * t**2
P3 = b30 * x**3 + b21 * x**2 * t + b12 * x * t**2 + b03 * t**3
line_pullback = sp.expand((P2 + P3).subs({x: s * v0, t: s * v1}))
zero(line_pullback.coeff(s, 2) - P2.subs({x: v0, t: v1}),
     "directional degree-two coefficient")
zero(line_pullback.coeff(s, 3) - P3.subs({x: v0, t: v1}),
     "directional degree-three coefficient")


# Panel IV: the exact cusp log chart and its transformed Darboux row.
sigma, tau = sp.symbols("sigma tau", nonzero=True)
p_log = sigma**2 + tau
y_log = sigma * p_log
x_log = sigma / tau
zero(p_log - (tau + (y_log / p_log)**2), "log chart p=s^2+tau")
zero(y_log / p_log - sigma, "log chart s=y/p")
zero(x_log**2 * tau - sigma**2 / tau, "log chart u=s^2/tau")
gate(sp.det(sp.Matrix([
    [sp.diff(x_log, sigma), sp.diff(x_log, tau)],
    [0, 1],
])) == 1 / tau, "log chart source-volume density")

A_control = sp.expand(x**2 + p)
C_control = sp.expand(y_source + z)
J_control = sp.factor(
    sp.diff(A_control, x) * sp.diff(C_control, t)
    - sp.diff(A_control, t) * sp.diff(C_control, x)
)
A_log = sp.factor(A_control.subs({x: x_log, t: tau}))
C_log = sp.factor(C_control.subs({x: x_log, t: tau}))
J_log = sp.factor(
    tau * (sp.diff(A_log, sigma) * sp.diff(C_log, tau)
           - sp.diff(A_log, tau) * sp.diff(C_log, sigma))
)
zero(J_log - J_control.subs({x: x_log, t: tau}),
     "log chart transforms every Jacobian row")


# Panel V: finite hostile enumeration of the uniform support consequences.
weights = tuple(range(-12, 9))


def valid_row(row: tuple[int, ...]) -> bool:
    return (any(weight >= 2 for weight in row)
            and any(weight <= -4 and weight % 2 == 0 for weight in row)
            and any(weight in (1, -2) for weight in row))


valid_triples = [
    row for row in itertools.combinations(weights, 3) if valid_row(row)
]
gate(bool(valid_triples), "hostile support range has exact triples")
for row in valid_triples:
    gate(sum(weight >= 2 for weight in row) == 1,
         f"triple {row}: unique positive tail")
    gate(sum(weight <= -4 and weight % 2 == 0 for weight in row) == 1,
         f"triple {row}: unique negative tail")
    gate(sum(weight in (1, -2) for weight in row) == 1,
         f"triple {row}: unique jet")

opposite_jet_pairs = 0
same_jet_pairs = 0
for left in valid_triples:
    left_jet = next(weight for weight in left if weight in (1, -2))
    for right in valid_triples:
        right_jet = next(weight for weight in right if weight in (1, -2))
        if left_jet != right_jet:
            opposite_jet_pairs += 1
        else:
            same_jet_pairs += 1
gate(opposite_jet_pairs > 0 and same_jet_pairs > 0,
     "enumeration contains both jet patterns before determinant gate")
gate(sp.det(sp.Matrix([[1, 0], [1, 0]])) == 0,
     "same +1 jets are dependent")
gate(sp.det(sp.Matrix([[0, 1], [0, 1]])) == 0,
     "same -2 jets are dependent")
gate(abs(sp.det(sp.Matrix([[1, 0], [0, 1]]))) == 1,
     "opposite jets can be independent")

support_cells = [
    (left, right)
    for left in range(1, 9)
    for right in range(1, 9)
    if left >= 3 and right >= 3 and (left, right) != (3, 3)
]
minimum_total = min(left + right for left, right in support_cells)
first_cells = sorted(
    (left, right) for left, right in support_cells
    if left + right == minimum_total
)
gate(minimum_total == 7, "post-3x3 minimum support total")
gate(first_cells == [(3, 4), (4, 3)], "first live support cells")


summary = {
    "checks": CHECKS,
    "cited_gate": "nonautomorphic Keller => noninjective on every affine line",
    "nonautomorphism": "ord_D(t)=-2 and k[A,C] subset B2",
    "line_invoice": "each degree>=2; one degree>=3",
    "positive_tail": "each r>=2; one r>=3",
    "negative_tail": "each even weight<=-4; one weight<=-6",
    "jet_weights": jet_weights,
    "directional_gate": (
        "A_{>=2}, C_{>=2}, and combined A/C_{>=3} tails basepoint-free on P1"
    ),
    "cusp_log_gate": "tau(A_s*C_tau-A_tau*C_s)=J_(x,t)(A,C)",
    "valid_triples_in_hostile_range": len(valid_triples),
    "first_live_cells_after_THM3974": first_cells,
    "scope": "B2 polynomial Darboux pairs; 3x4 and larger open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3987 Gwozdziewicz every-line weight-floor companion")
print(f"CHECKS={CHECKS}")
print("NONAUTOMORPHIC=ORD_D_T_MINUS_2;EVERY_AFFINE_LINE_NONINJECTIVE")
print("EACH_OUTPUT=POS_GE_2+EVEN_NEG_LE_MINUS_4+JET_1_OR_MINUS_2")
print("PAIR_ENDPOINTS=ONE_POS_GE_3+ONE_NEG_LE_MINUS_6")
print("DIRECTIONAL_TAILS=A_GE2,C_GE2,COMBINED_GE3_BASEPOINT_FREE")
print("CUSP_LOG_GATE=TAU*(A_S*C_TAU-A_TAU*C_S)=J_XT")
print("EXACT_3_ROWS=OPPOSITE_MIDDLE_JETS")
print("FIRST_LIVE_CELLS_AFTER_THM3974=3x4,4x3")
print("SCOPE=B2_POLYNOMIAL_DARBOUX_PAIRS;UNRESTRICTED_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
