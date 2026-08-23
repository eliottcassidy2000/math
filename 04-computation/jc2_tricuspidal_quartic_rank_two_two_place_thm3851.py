#!/usr/bin/env python3
"""Exact companion for THM-3851's tricuspidal quartic tradeoff."""

from __future__ import annotations

import ast
import hashlib
import json
from itertools import combinations
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


s, t, u = sp.symbols("s t u")
A, B, C = sp.symbols("A B C")

# The tricuspidal quartic and its basepoint-free normalization.
A_st = s**2 * t**2
B_st = t**2 * (s - t) ** 2
C_st = s**2 * (s - t) ** 2
delta = sp.expand((A * B + A * C + B * C) ** 2 - 4 * A * B * C * (A + B + C))
zero(delta.subs({A: A_st, B: B_st, C: C_st}), "normalization lies on quartic")
gate(sp.factor(delta) == delta, "tricuspidal quartic is irreducible")
gate(sp.Poly(delta, A, B, C).total_degree() == 4, "quartic degree")

# No normalization base point, including all three special addresses.
for address, values in {
    "s_zero": {s: 0, t: 1},
    "t_zero": {s: 1, t: 0},
    "s_equals_t": {s: 1, t: 1},
}.items():
    coordinates = [entry.subs(values) for entry in (A_st, B_st, C_st)]
    gate(any(entry != 0 for entry in coordinates), f"no base point at {address}")

# Rational inverse s/t=(AB+AC-BC)/(2AB).
inverse_numerator = A * B + A * C - B * C
inverse_denominator = 2 * A * B
zero(
    t * inverse_numerator.subs({A: A_st, B: B_st, C: C_st})
    - s * inverse_denominator.subs({A: A_st, B: B_st, C: C_st}),
    "birational normalization inverse",
)

# Complete projective singular support: one vertex in each standard chart.
for chart_variable, remaining in (
    (A, (B, C)),
    (B, (A, C)),
    (C, (A, B)),
):
    chart_polynomial = sp.expand(delta.subs(chart_variable, 1))
    chart_groebner = sp.groebner(
        [chart_polynomial, *[sp.diff(chart_polynomial, variable) for variable in remaining]],
        *remaining,
        order="lex",
    )
    first, second = remaining
    expected = sp.groebner([first - second, second**2], first, second, order="lex")
    gate(chart_groebner == expected, f"complete singular chart {chart_variable}")

# Direct A2 leading packets at u=0, u=1, and u=infinity.
v = sp.symbols("v")
zero(
    u**2 / (u - 1) ** 2 - u**2 - u**3 * (2 - u) / (u - 1) ** 2,
    "A2 packet at u zero",
)
zero(
    v**2 - v**2 / (1 + v) ** 2 - v**3 * (2 + v) / (1 + v) ** 2,
    "A2 packet at u one",
)
zero(
    v**2 / (1 - v) ** 2 - v**2 - v**3 * (2 - v) / (1 - v) ** 2,
    "A2 packet at u infinity",
)

# The chosen smooth bitangent and its two distinct normalization addresses.
bitangent_pullback = sp.expand(A_st + B_st + C_st)
bitangent_quadratic = s**2 - s * t + t**2
zero(bitangent_pullback - bitangent_quadratic**2, "split bitangent pullback")
zero(
    delta.subs(C, -A - B) - (A**2 + A * B + B**2) ** 2,
    "target-line branch restriction is a square",
)
gate(sp.discriminant(u**2 - u + 1, u) == -3, "two distinct bitangent addresses")
gate(sp.gcd(u**2 - u + 1, 2 * u - 1) == 1,
     "the two bitangent addresses have distinct projective images")
for values in ({s: 0, t: 1}, {s: 1, t: 0}, {s: 1, t: 1}):
    gate(bitangent_quadratic.subs(values) != 0, "bitangent avoids every cusp")

# Uniform no-fourfold-line proof.  On t=1, a general line pulls back to
# c*u^4-2c*u^3+(a+b+c)u^2-2bu+b.  A finite one-point support would be
# kappa*(u-r)^4.  A root at infinity would give kappa*t^4, constant here.
a, b, c, kappa, r, kappa_inverse = sp.symbols("a b c kappa r kappa_inverse")
line_pullback = sp.expand(
    a * A_st.subs({s: u, t: 1})
    + b * B_st.subs({s: u, t: 1})
    + c * C_st.subs({s: u, t: 1})
)
expected_line_pullback = c * u**4 - 2 * c * u**3 + (a + b + c) * u**2 - 2 * b * u + b
zero(line_pullback - expected_line_pullback, "general line pullback")
finite_fourfold_equations = [
    sp.Poly(line_pullback - kappa * (u - r) ** 4, u).coeff_monomial(u**degree)
    for degree in range(5)
]
finite_fourfold_groebner = sp.groebner(
    [*finite_fourfold_equations, kappa * kappa_inverse - 1],
    a,
    b,
    c,
    kappa,
    r,
    kappa_inverse,
    order="grevlex",
)
gate(finite_fourfold_groebner == sp.groebner(
    [1], a, b, c, kappa, r, kappa_inverse, order="grevlex"),
     "no finite fourfold line address")
infinity_fourfold_equations = [
    sp.Poly(line_pullback - kappa, u).coeff_monomial(u**degree)
    for degree in range(5)
]
infinity_fourfold_groebner = sp.groebner(
    [*infinity_fourfold_equations, kappa * kappa_inverse - 1],
    a,
    b,
    c,
    kappa,
    kappa_inverse,
    order="grevlex",
)
gate(
    infinity_fourfold_groebner
    == sp.groebner([1], a, b, c, kappa, kappa_inverse, order="grevlex"),
    "no infinity fourfold line address",
)


# Weak degree-two del Pezzo lattice.  Coordinates are H,E1,...,E7 with
# intersection diag(1,-1,...,-1).  The exact enumeration is deliberately
# separate from the standard geometric resolution bridge used in the proof.
def lattice_dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(x * y for x, y in zip(left[1:], right[1:]))


def lattice_neg(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-entry for entry in vector)


roots: set[tuple[int, ...]] = set()
for i in range(7):
    for j in range(7):
        if i == j:
            continue
        vector = [0] * 8
        vector[i + 1] = 1
        vector[j + 1] = -1
        roots.add(tuple(vector))
for indices in combinations(range(7), 3):
    vector = [0] * 8
    vector[0] = 1
    for index in indices:
        vector[index + 1] = -1
    roots.add(tuple(vector))
    roots.add(lattice_neg(tuple(vector)))
for omitted in range(7):
    vector = [2] + [-1] * 7
    vector[omitted + 1] = 0
    roots.add(tuple(vector))
    roots.add(lattice_neg(tuple(vector)))
gate(len(roots) == 126, "complete E7 roots")
gate(all(lattice_dot(root, root) == -2 for root in roots), "E7 root squares")

boundary_plus = (0, 0, 0, 0, 0, 0, 0, 1)
anticanonical = (3, -1, -1, -1, -1, -1, -1, -1)
boundary_minus = tuple(
    anticanonical[index] - boundary_plus[index] for index in range(8)
)
gate(lattice_dot(boundary_plus, boundary_plus) == -1, "first boundary square")
gate(lattice_dot(boundary_minus, boundary_minus) == -1, "second boundary square")
gate(lattice_dot(boundary_plus, boundary_minus) == 2, "two bitangent contacts")
gate(sp.Matrix([boundary_plus, boundary_minus]).rank() == 2,
     "boundary classes independent")

orthogonal_roots = sorted(root for root in roots if lattice_dot(root, boundary_plus) == 0)
gate(len(orthogonal_roots) == 72, "boundary-orthogonal E6 roots")
root_lookup = set(orthogonal_roots)
a2_systems: dict[
    tuple[tuple[int, ...], ...], tuple[tuple[int, ...], tuple[int, ...]]
] = {}
for left, right in combinations(orthogonal_roots, 2):
    if lattice_dot(left, right) != 1:
        continue
    root_sum = tuple(left[index] + right[index] for index in range(8))
    subsystem = tuple(sorted({left, right, lattice_neg(left), lattice_neg(right),
                              root_sum, lattice_neg(root_sum)}))
    gate(all(root in root_lookup for root in subsystem), "A2 subsystem closure")
    a2_systems[subsystem] = (left, right)
gate(len(a2_systems) == 120, "complete E6 A2 systems")

triple_configurations: list[
    tuple[
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[tuple[int, ...], tuple[int, ...]],
    ]
] = []
a2_items = list(a2_systems.items())
for first_index, (first_roots, first_pair) in enumerate(a2_items):
    for second_index in range(first_index + 1, len(a2_items)):
        second_roots, second_pair = a2_items[second_index]
        if any(lattice_dot(left, right) != 0 for left in first_roots for right in second_roots):
            continue
        for third_index in range(second_index + 1, len(a2_items)):
            third_roots, third_pair = a2_items[third_index]
            if any(
                lattice_dot(left, right) != 0
                for left in first_roots + second_roots
                for right in third_roots
            ):
                continue
            triple_configurations.append((first_pair, second_pair, third_pair))
gate(len(triple_configurations) == 40, "complete three-A2 configurations")

smith_profiles: dict[tuple[int, ...], int] = {}
for first_pair, second_pair, third_pair in triple_configurations:
    relation_matrix = sp.Matrix(
        [*first_pair, *second_pair, *third_pair, boundary_plus, boundary_minus]
    )
    smith = smith_normal_form(relation_matrix, domain=ZZ)
    profile = tuple(abs(int(smith[index, index])) for index in range(8))
    smith_profiles[profile] = smith_profiles.get(profile, 0) + 1
gate(
    smith_profiles == {(1, 1, 1, 1, 1, 1, 3, 3): 40},
    "every affine quotient is Z3 direct sum Z3",
)

# One marking displays the quotient and its local rank-two address map.
r1 = (1, -1, -1, -1, 0, 0, 0, 0)
r2 = (1, 0, 0, 0, -1, -1, -1, 0)
r3 = (0, 0, 1, -1, 0, 0, 0, 0)
r4 = (0, 1, -1, 0, 0, 0, 0, 0)
r5 = (0, 0, 0, 0, 0, 1, -1, 0)
r6 = (0, 0, 0, 0, 1, -1, 0, 0)
explicit_matrix = sp.Matrix([r1, r2, r3, r4, r5, r6, boundary_plus, boundary_minus])
explicit_smith = smith_normal_form(explicit_matrix, domain=ZZ)
gate(
    tuple(abs(int(explicit_smith[index, index])) for index in range(8))
    == (1, 1, 1, 1, 1, 1, 3, 3),
    "explicit rank-two Smith quotient",
)
x_class = (0, 1, 0, 0, 0, 0, 0, 0)  # E1
y_class = (0, 0, 0, 0, 1, 0, 0, 0)  # E4
three_x_coefficients = (-2, -1, 1, 2, 0, 0, 2, 1)
three_y_coefficients = (-1, -2, 0, 0, 1, 2, 2, 1)
for divisor, coefficients, label in (
    (x_class, three_x_coefficients, "x"),
    (y_class, three_y_coefficients, "y"),
):
    relation = tuple(
        sum(coefficients[row] * explicit_matrix[row, column] for row in range(8))
        for column in range(8)
    )
    gate(relation == tuple(3 * entry for entry in divisor),
         f"explicit order-three relation for {label}")
ordered_a2_pairs = ((r1, r2), (r3, r4), (r5, r6))


def local_a2_character(divisor: tuple[int, ...], pair: tuple[tuple[int, ...], tuple[int, ...]]) -> int:
    return (lattice_dot(divisor, pair[0]) + 2 * lattice_dot(divisor, pair[1])) % 3


x_local = tuple(local_a2_character(x_class, pair) for pair in ordered_a2_pairs)
y_local = tuple(local_a2_character(y_class, pair) for pair in ordered_a2_pairs)
gate(x_local == (1, 1, 0), "first global class local addresses")
gate(y_local == (2, 0, 1), "second global class local addresses")
gate(
    all((vector[0] - vector[1] + vector[2]) % 3 == 0 for vector in (x_local, y_local)),
    "local classes lie in the rank-two anti-sum plane",
)
gate(
    (x_local[1] * y_local[2] - x_local[2] * y_local[1]) % 3 != 0,
    "local address map has rank two",
)
canonical = lattice_neg(anticanonical)
for divisor, label in ((x_class, "x"), (y_class, "y")):
    deck_divisor = tuple(
        lattice_dot(divisor, canonical) * canonical[index] - divisor[index]
        for index in range(8)
    )
    gate(
        tuple(deck_divisor[index] + divisor[index] for index in range(8))
        == anticanonical,
        f"deck inversion of {label} modulo the boundary",
    )


semantic = {
    "branch": "tricuspidal quartic (AB+AC+BC)^2-4ABC(A+B+C)",
    "normalization": "[s2t2:t2(s-t)2:s2(s-t)2];three A2 cusps",
    "bitangent": "A+B+C pulls back to (s2-st+t2)^2;two distinct places",
    "fourfold": "no nonzero line has fourth-power pullback, including infinity root",
    "resolvent": "split-boundary quadratic double plane;Cl=(Z/3)^2;units constants",
    "lattice": "40 compatible 3A2 markings;uniform Smith profile 1^6,3,3",
    "scope": "rank-two/two-place design tradeoff;no cubic Keller atlas claimed",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff")
print("branch=tricuspidal_quartic;finite_singularities=3A2")
print("chosen_bitangent_places=2;fourfold_line=NONE")
print("quadratic_resolvent_affine_class_group=Z3_plus_Z3")
print("root_boundary_configurations=40;smith_profile=1^6_3_3_uniform")
print("local_torsion_rank=2;deck_character=anti_invariant")
print("scope=rank_two_two_place_tradeoff_no_Keller_counterexample")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
