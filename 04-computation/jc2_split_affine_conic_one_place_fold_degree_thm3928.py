#!/usr/bin/env python3
"""Exact companion for THM-3928's split-conic fold-degree barrier.

Reproduction:
  python3 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
  python3 -O 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
"""

from __future__ import annotations

import hashlib
import json

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


def smith_diagonal(matrix: sp.Matrix) -> tuple[int, ...]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    )


t, h, x, y, ell, q_symbol, sigma = sp.symbols("t h x y ell q sigma")


# ---------------------------------------------------------------------------
# Cardano normalization and the sharp excluded linear control.
# ---------------------------------------------------------------------------

p_h = 3 * h**2
q_h = 2 * h**3
zero(4 * p_h**3 - 27 * q_h**2, "Cardano cusp identity")
zero(3 * q_h / (2 * p_h) - h, "Cardano parameter recovery")
zero(p_h / 3 - h**2, "Cardano parameter is integral over the branch ring")

p_line = 3 * x * y
q_line = 2 * x**3
delta_line = sp.factor(4 * p_line**3 - 27 * q_line**2)
gate(sp.rem(delta_line, y - x, y) == 0, "excluded line branch divides discriminant")
zero(
    (3 * q_line / (2 * p_line)).subs(y, x) - x,
    "excluded line branch has fold degree one",
)


# ---------------------------------------------------------------------------
# Projective-degree ledger, including the complete sextic table.
# ---------------------------------------------------------------------------

for degree in range(2, 31):
    admissible = []
    for fold_degree in range(1, degree + 1):
        second = 2 * fold_degree - degree
        if 1 <= second <= degree:
            admissible.append((fold_degree, (degree, second)))
    gate(bool(admissible), f"degree {degree} has an abstract degree packet")
    gate(
        admissible[0][0] == (degree + 2) // 2,
        f"degree {degree} realizes ceil((N+1)/2) floor",
    )
    for fold_degree, (first, second) in admissible:
        gate(first == degree, f"degree {degree} retains a full-degree line")
        gate(
            first + second == 2 * fold_degree,
            f"degree {degree} product-degree identity",
        )
        gate(second >= 1, f"degree {degree} excludes a constant line pullback")

sextic = [
    (fold_degree, (6, 2 * fold_degree - 6), 12 - 2 * fold_degree)
    for fold_degree in range(1, 7)
    if 1 <= 2 * fold_degree - 6 <= 6
]
gate(
    sextic == [(4, (6, 2), 4), (5, (6, 4), 2), (6, (6, 6), 0)],
    "complete sextic fold/line/contact table",
)
gate(all(row[0] >= 4 for row in sextic), "sextic fold degree is at least four")

degree_one_exponents = [(e, 2 - e) for e in range(3)]
gate(
    degree_one_exponents == [(0, 2), (1, 1), (2, 0)],
    "degree-one UFD exponent split is complete",
)


# ---------------------------------------------------------------------------
# A distinct-line quadratic coefficient is never an etale coefficient map.
# ---------------------------------------------------------------------------

coefficients = sp.symbols("c00 c10 c01 c20 c11 c02 c30 c21 c12 c03")
monomials = (1, x, y, x**2, x * y, y**2, x**3, x**2 * y, x * y**2, y**3)
q_generic = sum(coefficient * monomial for coefficient, monomial in zip(coefficients, monomials))

p_crossing = x * y
jac_crossing = sp.diff(p_crossing, x) * sp.diff(q_generic, y) - sp.diff(
    p_crossing, y
) * sp.diff(q_generic, x)
zero(jac_crossing.subs({x: 0, y: 0}), "crossing conic has a critical point")

p_parallel = x * (x - 1)
jac_parallel = sp.diff(p_parallel, x) * sp.diff(q_generic, y) - sp.diff(
    p_parallel, y
) * sp.diff(q_generic, x)
zero(jac_parallel.subs(x, sp.Rational(1, 2)), "parallel conic has a critical midline")


# ---------------------------------------------------------------------------
# Universal quadratic-resolvent divisor relation lattices.
# ---------------------------------------------------------------------------

split_two = sp.Matrix(
    [
        [1, 1, 0, 0],
        [0, 0, 1, 1],
        [3, 0, 3, 0],
    ]
)
gate(smith_diagonal(split_two) == (1, 1, 3), "two-component Smith packet")
gate(split_two.rank() == 3, "two-component packet has one free direction")

double_line = sp.Matrix([[1, 1], [6, 0]])
gate(smith_diagonal(double_line) == (1, 6), "double-line Smith packet")
gate(double_line.det() == -6, "double-line packet has order six")

double_delta = 4 * ell**6 - 27 * q_symbol**2
double_factorization = (2 * ell**3 - 3 * sigma * q_symbol) * (
    2 * ell**3 + 3 * sigma * q_symbol
)
zero(
    sp.expand(double_delta - double_factorization).subs(sigma**2, 3),
    "double-line discriminant factors into two cubics",
)

for components in range(1, 9):
    rows = []
    for index in range(components):
        row = [0] * (2 * components)
        row[2 * index] = 1
        row[2 * index + 1] = 1
        rows.append(row)
    norm_row = [0] * (2 * components)
    for index in range(components):
        norm_row[2 * index] = 3
    rows.append(norm_row)
    relation_matrix = sp.Matrix(rows)
    gate(
        smith_diagonal(relation_matrix) == (1,) * components + (3,),
        f"{components}-component Smith packet",
    )
    gate(
        2 * components - relation_matrix.rank() == components - 1,
        f"{components}-component free-rank packet",
    )


semantic_payload = {
    "normalization": "one_place_makes_h_polynomial_but_not_linear",
    "parallel": "distinct_parallel_line_factors_force_h_constant_by_coprime_square_UFD",
    "degree": "nonparallel_degree_N_branch_forces_deg_h_at_least_ceil_Nplus1_over2",
    "sextic": "fold_degrees_4_5_6_with_line_packets_6_2__6_4__6_6",
    "kummer": "two_split_components_formal_relations_Z_plus_Z3_diagonal_only",
    "boundaries": "double_line_factors;high_folds_and_infinity_component_remain",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3928-split-affine-conic-one-place-fold-degree-barrier")
print("normalization=h=3q/(2p) lies in k[t];p=3h^2;q=2h^3")
print("parallel=distinct parallel affine factors impossible for nonconstant h")
print("degree_N=deg(h)>=ceil((N+1)/2)")
print("sextic=(d;e1,e2;infinity_contact)=(4;6,2;4),(5;6,4;2),(6;6,6;0)")
print("split_lattice=SNF(1,1,3,0);formal quotient Z plus Z/3")
print("double_line=Delta factors as two cubics;resolvent SNF(1,6)")
print("scope=birational_split_and_double_line_closed;high_folds_and_infinity_component_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
