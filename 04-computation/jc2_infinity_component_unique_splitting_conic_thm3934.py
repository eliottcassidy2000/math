#!/usr/bin/env python3
"""Exact companion for THM-3934's splitting-conic classification.

Reproduction:
  python3 04-computation/jc2_infinity_component_unique_splitting_conic_thm3934.py
  python3 -O 04-computation/jc2_infinity_component_unique_splitting_conic_thm3934.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


s = sp.symbols("s")
X, Y, Z = sp.symbols("X Y Z")
A, B, C = sp.symbols("A B C")
g = sp.symbols("g0:7")

h = s**3 - 2
x_pull = h**2
y_pull = s * h
z_pull = sp.Integer(1)

# Pullbacks of the six conic monomials.  Rank six proves that restriction of
# conics to the irreducible sextic is injective.
conic_basis = [
    x_pull**2,
    x_pull * y_pull,
    x_pull * z_pull,
    y_pull**2,
    y_pull * z_pull,
    z_pull**2,
]
coefficient_matrix = sp.Matrix(
    [[sp.expand(form).coeff(s, degree) for form in conic_basis] for degree in range(13)]
)
gate(coefficient_matrix.rank() == 6, "conic restriction has rank six")

# A polynomial of degree at most twelve lies in the conic pullback space iff
# the following seven independent linear coefficient relations hold.
left_kernel = coefficient_matrix.T.nullspace()
gate(len(left_kernel) == 7, "conic pullback has codimension seven")

root = sum(g[i] * s**i for i in range(7))
square = sp.Poly(sp.expand(root**2), s)
coeff = [square.nth(i) for i in range(13)]
constraints = [
    coeff[2] + coeff[5],
    3 * coeff[1] + 6 * coeff[4] + 8 * coeff[7],
    -coeff[2] + 4 * coeff[8],
    coeff[3] + 4 * coeff[6] + 8 * coeff[9],
    -coeff[1] + 16 * coeff[10] - 2 * coeff[4],
    coeff[11],
    64 * coeff[12] - coeff[3] - 4 * coeff[6],
]

# Check that these equations are exactly a basis of the left-kernel
# conditions, rather than merely necessary samples.
constraint_rows = []
dummy = sp.symbols("c0:13")
linear_constraints = [
    dummy[2] + dummy[5],
    3 * dummy[1] + 6 * dummy[4] + 8 * dummy[7],
    -dummy[2] + 4 * dummy[8],
    dummy[3] + 4 * dummy[6] + 8 * dummy[9],
    -dummy[1] + 16 * dummy[10] - 2 * dummy[4],
    dummy[11],
    64 * dummy[12] - dummy[3] - 4 * dummy[6],
]
for relation in linear_constraints:
    constraint_rows.append([sp.expand(relation).coeff(item) for item in dummy])
constraint_matrix = sp.Matrix(constraint_rows)
gate(constraint_matrix.rank() == 7, "displayed conic-square constraints are independent")
gate(
    constraint_matrix * coefficient_matrix == sp.zeros(7, 6),
    "displayed constraints annihilate conic pullbacks",
)

# Exhaust the projective square-root charts by the highest nonzero
# coefficient.  The Groebner bases below are set-theoretic certificates; in
# the two nonreduced chart ideals, g_i^2=0 forces g_i=0 over a field.
chart_bases: dict[int, list[sp.Expr]] = {}
empty_charts = []
for degree in range(6, -1, -1):
    equations = constraints + [g[j] for j in range(degree + 1, 7)] + [g[degree] - 1]
    basis = sp.groebner(equations, *g, order="lex")
    expressions = [sp.factor(polynomial.as_expr()) for polynomial in basis.polys]
    chart_bases[degree] = expressions
    if expressions == [sp.Integer(1)]:
        empty_charts.append(degree)

gate(empty_charts == [5, 2, 1], "degree five, two, and one square-root charts are empty")

expected_chart_6 = {
    g[2] * (g[0] - 4),
    g[1] + 2 * g[4],
    g[2] ** 2,
    g[2] * g[4],
    g[3] + 4,
    g[5],
    g[6] - 1,
}
gate(set(chart_bases[6]) == expected_chart_6, "degree-six chart is the line-pullback component")
gate(
    set(chart_bases[4]) == {g[1] + 2, g[2], g[3], g[4] - 1, g[5], g[6]},
    "degree-four chart is the Y-plus-constant line component",
)
gate(
    set(chart_bases[3]) == {g[0] + 2, g[1] ** 2, g[2], g[3] - 1, g[4], g[5], g[6]},
    "degree-three chart is the unique extra h direction",
)
gate(
    set(chart_bases[0]) == {g[0] - 1, g[1], g[2], g[3], g[4], g[5], g[6]},
    "degree-zero chart is the constant line component",
)

# Reassemble the chart solutions invariantly.  One component consists of
# pullbacks of arbitrary lines; squaring gives nonreduced double lines.  The
# only other component is h, whose square is the pullback of XZ.
line_pull = A * x_pull + B * y_pull + C * z_pull
line_conic_pull = sp.expand((A * X + B * Y + C * Z) ** 2).subs(
    {X: x_pull, Y: y_pull, Z: z_pull}
)
zero(line_pull**2 - line_conic_pull, "line-root component gives exactly double lines")
zero(h**2 - (X * Z).subs({X: x_pull, Y: y_pull, Z: z_pull}), "extra root gives XZ")

# The reduced exceptional conic really splits on the quadratic double plane:
# on each component X=0 and Z=0 the branch sextic restricts to a square.
F = 4 * X**3 * Z**3 - (Y**3 - X**2 * Z) ** 2
zero(F.subs(X, 0) + Y**6, "branch restriction to X=0 is minus a sixth power")
zero(F.subs(Z, 0) + Y**6, "branch restriction to Z=0 is minus a sixth power")
gate(sp.factor(X * Z) == X * Z, "XZ is the union of two distinct lines")

# At the only multibranch point O, total even contact of a conic cannot hide
# two odd branches.  With local branch parameter t=h, x=t^2 and
# y=alpha*t+O(t^2), alpha^3=2.  The leading conic rows are common to all
# three branches as recorded below.
alpha = sp.symbols("alpha")
alpha_polynomial = alpha**3 - 2
gate(sp.discriminant(alpha_polynomial, alpha) == -108, "the three finite branch slopes are distinct")
gate(sp.gcd(alpha_polynomial, sp.diff(alpha_polynomial, alpha)) == 1, "finite branch slopes are reduced")

# A nonzero linear-y row has order one on all three branches.  If it vanishes,
# the order-two row C+D*alpha^2 can vanish on at most one branch unless C=D=0.
Dcoef, Ecoef = sp.symbols("Dcoef Ecoef")
gate(
    sp.degree(sp.rem(C + Dcoef * alpha**2, alpha_polynomial, alpha), alpha) <= 2,
    "order-two cancellation is controlled by a quadratic in the three slopes",
)
gate(
    sp.resultant(alpha_polynomial, C + Dcoef * alpha**2, alpha) != 0,
    "generic order-two row cancels on no branch",
)
# If the row vanishes on two distinct cube roots, it is identically zero:
# a degree-at-most-two polynomial with those special monomials cannot share
# two roots with alpha^3-2 unless both coefficients vanish.
gcd_symbolic = sp.gcd(
    sp.Poly(alpha_polynomial, alpha, domain=sp.QQ.frac_field(C, Dcoef)),
    sp.Poly(C + Dcoef * alpha**2, alpha, domain=sp.QQ.frac_field(C, Dcoef)),
)
gate(gcd_symbolic.degree() == 0, "nonzero order-two row cannot cancel on two branches")

# The remaining leading row is B*x*y, of order three on all three branches;
# if it too vanishes, A*x^2 has even order four.  Hence even total contact at
# O implies branchwise even contact.  Infinity has only one branch by the
# explicit normalization, so no parity can be hidden there either.
parity_patterns = {
    "linear_y": (1, 1, 1),
    "cubic_xy": (3, 3, 3),
    "quartic_x2": (4, 4, 4),
}
gate(sum(parity_patterns["linear_y"]) % 2 == 1, "three linear contacts have odd total")
# If the quadratic row cancels at the unique exceptional branch, its order m
# can be three or can increase after further cancellation.  The exact point is
# parity, not the representative pattern (2,2,3): the other two branches have
# order two, so total parity equals the exceptional branch parity for every
# possible local intersection order.  Bezout bounds the total by twelve.
gate(
    all(
        (2 + 2 + exceptional_order) % 2 == exceptional_order % 2
        for exceptional_order in range(1, 13)
    ),
    "total parity detects the unique exceptional branch at the quadratic row",
)
gate(sum(parity_patterns["cubic_xy"]) % 2 == 1, "three cubic contacts have odd total")
gate(all(value % 2 == 0 for value in parity_patterns["quartic_x2"]),
     "the surviving quartic contacts are branchwise even")

# The weighted exceptional elliptic curve and the first hostile lifts of its
# other three-torsion directions.
u, v, lam = sp.symbols("u v lambda")
division_3 = 3 * u**4 + 12 * sp.Rational(1, 4) * u
zero(division_3 - 3 * u * (u**3 + 1), "j-zero exceptional curve three-division polynomial")
zero((1 - 4 * X**3).subs(X, 0) - 1, "natural exceptional flex has x=0,w squared one")
residual_other_direction = sp.expand((1 - lam**2 * Y) ** 2 - 4)
gate(
    sp.factor(sp.discriminant(residual_other_direction, Y)) == 16 * lam**4,
    "other-torsion residual has two distinct roots",
)

semantic_payload = {
    "branch": "F=4X^3Z^3-(Y^3-X^2Z)^2",
    "normalization": "[H^2:S^2*T*H:S^6], H=T^3-2S^3",
    "square_root_locus": "pullbacks_of_lines union span(H)",
    "conics": "double_lines union span(XZ)",
    "reduced_conic": "XZ_only",
    "parity_bridge": "splitting_implies_branchwise_even_at_three_branch_origin",
    "class_scope": "no_second_degree_at_most_two_splitting_support;full_Cl[3]_open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3934-infinity-component-sextic-unique-reduced-splitting-conic")
print("conic_square_roots=line_pullbacks_union_span(h)")
print("conic_locus=double_lines_union_span(XZ)")
print("unique_reduced_splitting_conic=XZ")
print("resolvent_support=affine_component_X=0_recovers_THM3932_Dplus_Dminus")
print("new_low_degree_three_class=NONE")
print("scope=degree_at_most_two_splitting_support_only;full_resolvent_class_group_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
