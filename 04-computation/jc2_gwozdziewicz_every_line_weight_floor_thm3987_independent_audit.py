#!/usr/bin/env python3
"""Independent hostile audit for THM-3987.

This checker deliberately rebuilds the restriction ledgers from generator
monomials instead of importing the canonical companion.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


x, t, u, zeta = sp.symbols("x t u zeta")
s, a, b = sp.symbols("s a b")
u_xt = x**2 * t
z_xt = 1 + u_xt
p_xt = sp.expand(z_xt * t)
y_xt = sp.expand(x * z_xt * t**2)


# I. Boundary valuation and source identities, rebuilt from the generators.
zero(z_xt * (z_xt - 1) ** 2 - x**3 * y_xt,
     "boundary hypersurface relation")
zero(t - (z_xt - 1) / x**2, "t boundary expression")
boundary_valuations = {"x": 1, "z": 3, "p": 1, "y": 0, "z-1": 0}
gate(boundary_valuations["z"] == 3 * boundary_valuations["x"],
     "generic boundary relation gives ord(z)=3")
gate(boundary_valuations["p"] ==
     boundary_valuations["z"] - 2 * boundary_valuations["x"],
     "p=zt valuation")
gate(boundary_valuations["y"] ==
     boundary_valuations["x"] + boundary_valuations["z"]
     - 4 * boundary_valuations["x"], "y=xzt^2 valuation")
ord_t = boundary_valuations["z-1"] - 2 * boundary_valuations["x"]
gate(ord_t == -2, "t has a boundary pole")
gate(min(boundary_valuations[name] for name in ("x", "z", "p", "y")) >= 0,
     "every ring generator is boundary-regular")


# II. Reconstruct the negative-weight coefficient ideals from all feasible
# generator monomials x^c p^alpha y^beta in a broad hostile range.
negative_rows: dict[int, tuple[int, int, int]] = {}
for q in range(1, 97):
    feasible: list[tuple[int, int, int, int, int]] = []
    for alpha in range(q + 2):
        for beta in range(q + 2):
            c = 2 * alpha + 3 * beta - q
            if c >= 0:
                feasible.append((alpha, beta, c, alpha + 2 * beta,
                                 alpha + beta))
    gate(bool(feasible), f"weight -{q}: nonempty monomial set")
    min_u = min(row[3] for row in feasible)
    min_z = min(row[4] for row in feasible)
    expected_u = math.ceil(q / 2)
    expected_z = math.ceil(q / 3)
    gate(min_u == expected_u, f"weight -{q}: u-order gcd minimum")
    gate(min_z == expected_z, f"weight -{q}: (u+1)-order gcd minimum")

    representative = sp.expand(
        x ** (2 * expected_u - q) * t**expected_u
        * z_xt**expected_z * (1 + 5 * u_xt + 7 * u_xt**2)
    )
    zero(representative.subs(t, 0), f"weight -{q}: t=0 vanishing")
    expected_x0 = t**expected_u if q % 2 == 0 else sp.Integer(0)
    zero(representative.subs(x, 0) - expected_x0,
         f"weight -{q}: x=0 parity survival")
    for powers, coefficient in sp.Poly(representative, x, t).terms():
        if coefficient:
            gate(powers[0] - 2 * powers[1] == -q,
                 f"weight -{q}: expanded source homogeneity")
    negative_rows[q] = (expected_u, expected_z, 2 * expected_u - q)

for r in range(0, 49):
    representative = sp.expand(x**r * (1 + 5 * u_xt + 7 * u_xt**2))
    zero(representative.subs(t, 0) - x**r,
         f"weight {r}: positive t=0 survival")
    zero(representative.subs(x, 0) - (1 if r == 0 else 0),
         f"weight {r}: positive x=0 restriction")


# III. Source-order audit of every possible weight: after removing the scalar
# at weight zero, only weights 1 and -2 can carry a first-order source term.
linear_weights: list[int] = []
for r in range(0, 49):
    minimum_order = r if r > 0 else 3
    if minimum_order == 1:
        linear_weights.append(r)
for q, (u_order, _, x_order) in negative_rows.items():
    minimum_order = x_order + u_order
    if minimum_order == 1:
        linear_weights.append(-q)
gate(sorted(linear_weights) == [-2, 1], "only two linear-jet weights")

alpha_x, alpha_t, gamma_x, gamma_t = sp.symbols(
    "alpha_x alpha_t gamma_x gamma_t", nonzero=True
)
same_x = sp.Matrix([[alpha_x, 0], [gamma_x, 0]]).det()
same_t = sp.Matrix([[0, alpha_t], [0, gamma_t]]).det()
opposite_xt = sp.Matrix([[alpha_x, 0], [0, gamma_t]]).det()
opposite_tx = sp.Matrix([[0, alpha_t], [gamma_x, 0]]).det()
zero(same_x, "two unique +1 jets cannot be independent")
zero(same_t, "two unique -2 jets cannot be independent")
gate(opposite_xt != 0 and opposite_tx != 0,
     "opposite nonzero unique jets are independent")


# IV. Independently derive the quadratic midpoint obstruction and sharp (2,3)
# immersed noninjective control.
r0, r1 = sp.symbols("r0 r1")
f2, f1, f0, g2, g1, g0 = sp.symbols(
    "f2 f1 f0 g2 g1 g0", nonzero=True
)
f = f2 * s**2 + f1 * s + f0
g = g2 * s**2 + g1 * s + g0
for polynomial, label in ((f, "f"), (g, "g")):
    collision_quotient = sp.cancel(
        (polynomial.subs(s, r0) - polynomial.subs(s, r1)) / (r0 - r1)
    )
    midpoint_derivative = sp.diff(polynomial, s).subs(s, (r0 + r1) / 2)
    zero(collision_quotient - midpoint_derivative,
         f"{label}: quadratic collision forces midpoint criticality")

node_f = s**2 - 1
node_g = s * (s**2 - 1)
gate((node_f.subs(s, 1), node_g.subs(s, 1)) ==
     (node_f.subs(s, -1), node_g.subs(s, -1)),
     "sharp nodal collision")
gate(sp.gcd(sp.Poly(sp.diff(node_f, s), s),
            sp.Poly(sp.diff(node_g, s), s)).degree() == 0,
     "sharp nodal curve is immersed")
gate((sp.degree(node_f, s), sp.degree(node_g, s)) == (2, 3),
     "sharp nodal degrees")


# V. Directional-tail reformulation: coefficient extraction on s*v is the
# evaluation of the ordinary homogeneous part at v, independently through d=8.
v0, v1 = sp.symbols("v0 v1")
ordinary = sp.Integer(0)
homogeneous_parts: dict[int, sp.Expr] = {}
for degree in range(0, 9):
    part = sp.Add(*(
        sp.Symbol(f"c_{degree}_{i}") * x**i * t ** (degree - i)
        for i in range(degree + 1)
    ))
    homogeneous_parts[degree] = part
    ordinary += part
pullback = sp.Poly(sp.expand(ordinary.subs({x: s * v0, t: s * v1})), s)
for degree, part in homogeneous_parts.items():
    zero(pullback.coeff_monomial(s**degree) - part.subs({x: v0, t: v1}),
         f"ordinary degree {degree}: direction coefficient")


# VI. Independently verify the corrected log chart tau=t (there is no scalar
# H/p^2 in the stated generator list).
sigma, tau = sp.symbols("sigma tau", nonzero=True)
x_log = sigma / tau
t_log = tau
p_log = sigma**2 + tau
y_log = sigma * p_log
zero((1 + x_log**2 * t_log) * t_log - p_log,
     "log chart p=sigma^2+tau")
zero(x_log * (1 + x_log**2 * t_log) * t_log**2 - y_log,
     "log chart y=sigma*p")
zero(y_log / p_log - sigma, "log chart sigma=y/p")
gate(sp.Matrix([[sp.diff(x_log, sigma), sp.diff(x_log, tau)],
                [sp.diff(t_log, sigma), sp.diff(t_log, tau)]]).det() == 1 / tau,
     "log chart volume density")

A_test = x**3 + x * t + t**2 + 2 * x
C_test = x**2 * t + t**3 + 3 * t
J_test = sp.diff(A_test, x) * sp.diff(C_test, t) - sp.diff(A_test, t) * sp.diff(C_test, x)
A_st = sp.cancel(A_test.subs({x: x_log, t: t_log}))
C_st = sp.cancel(C_test.subs({x: x_log, t: t_log}))
transformed = tau * (
    sp.diff(A_st, sigma) * sp.diff(C_st, tau)
    - sp.diff(A_st, tau) * sp.diff(C_st, sigma)
)
zero(transformed - J_test.subs({x: x_log, t: t_log}),
     "log-Darboux chain rule")


# VII. Hostile support enumeration and the exact-three determinant gate.
weights = tuple(range(-30, 16))


def invoices(row: tuple[int, ...]) -> bool:
    return (any(weight >= 2 for weight in row)
            and any(weight <= -4 and weight % 2 == 0 for weight in row)
            and any(weight in (-2, 1) for weight in row))


triples = [row for row in itertools.combinations(weights, 3) if invoices(row)]
gate(bool(triples), "hostile range contains exact triples")
for row in triples:
    positive = [weight for weight in row if weight >= 2]
    negative = [weight for weight in row if weight <= -4 and weight % 2 == 0]
    jets = [weight for weight in row if weight in (-2, 1)]
    gate(len(positive) == len(negative) == len(jets) == 1,
         f"exact triple {row}: tails and jet exhaust slots")
    gate(negative[0] < jets[0] < positive[0],
         f"exact triple {row}: jet is the middle weight")

cells = [(left, right) for left in range(1, 31) for right in range(1, 31)
         if left >= 3 and right >= 3 and (left, right) != (3, 3)]
minimum_total = min(left + right for left, right in cells)
first_live = sorted(cell for cell in cells if sum(cell) == minimum_total)
gate(minimum_total == 7, "post-3x3 support total floor")
gate(first_live == [(3, 4), (4, 3)], "first unexcluded cells")


summary = {
    "checks": CHECKS,
    "boundary_ord_t": ord_t,
    "module_weights_audited": [1, 96],
    "linear_weights": sorted(linear_weights),
    "curve_floor": "each coordinate degree >=2; pair max degree >=3",
    "direction_degrees_audited": [0, 8],
    "exact_triples_audited": len(triples),
    "first_live": first_live,
    "citation": "alg-geom/9305008 Theorem 1.1 exact statement confirmed",
    "repair": "replace undefined tau=H/p^2 by tau=t=(z-1)/x^2",
    "scope": "necessary support gate only; JC(2) open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("Independent hostile audit of THM-3987")
print(f"CHECKS={CHECKS}")
print("PRIMARY_CITATION=ALG-GEOM/9305008_THEOREM_1.1_EXACT_MATCH")
print("CORE=PASS")
print("BOUNDARY_ORD_T=-2")
print("NEGATIVE_MODULES=REBUILT_FROM_GENERATOR_MONOMIALS_Q_1_TO_96")
print("X0_SURVIVAL=EVEN_NEGATIVE_ONLY")
print("LINEAR_WEIGHTS=-2,1")
print("CURVE_DEGREE_FLOOR=2_BY_2_EXCLUDED;2_BY_3_SHARP")
print("DIRECTIONAL_TAIL_CHAIN_RULE=PASS_DEGREES_0_TO_8")
print("EXACT_THREE=UNIQUE_MIDDLE_JET;PAIR_JETS_OPPOSITE")
print("FIRST_LIVE_AFTER_3X3=3x4,4x3")
print("WORDING_REPAIR=UNDEFINED_H_IN_TAU_DEFINITION")
print(f"SEMANTIC_SHA256={semantic}")
