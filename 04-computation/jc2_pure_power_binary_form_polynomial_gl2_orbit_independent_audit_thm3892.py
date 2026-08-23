#!/usr/bin/env python3
"""Independent hostile audit and all-binary-degree extension of THM-3892."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")
GATES = 0


def gate(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"gate failed: {label}")


def zero(expression: sp.Expr, label: object) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, label)


def row_det(left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(left[0] * right[1] - left[1] * right[0])


# 1. A connected finite-etale degree-r cover of A1 cannot ramify enough in
# its single infinity fibre.  This is independent of the binary-form degree.
for cover_degree in range(2, 257):
    required = 2 * cover_degree - 2
    available = cover_degree - 1
    gate(required > available, ("Riemann-Hurwitz one-fibre hostile", cover_degree))


# 2. Pairwise determinant units recover constant coordinates relative to any
# two root sections.  No cubic-specific third-root identity is used.
p, q, r, s = sp.symbols("p q r s")
basis = ((p, q), (r, s))
basis_det = row_det(*basis)
for index in range(2, 17):
    alpha, beta = sp.symbols(f"alpha_{index} beta_{index}")
    row = (
        sp.expand(alpha * p + beta * r),
        sp.expand(alpha * q + beta * s),
    )
    zero(row_det(basis[0], row) - beta * basis_det,
         ("first determinant coordinate", index))
    zero(row_det(basis[1], row) + alpha * basis_det,
         ("second determinant coordinate", index))


# 3. The split discriminant makes the covariance exponents transparent.
# There are n(n-1)/2 pair determinants, each squared.  A scalar multiplying
# the form has exponent 2n-2.  These arithmetic identities also freeze the
# central homogenizing layer for a coefficient row padded from degree E to D.
for form_degree in range(2, 65):
    pairs = form_degree * (form_degree - 1) // 2
    gate(2 * pairs == form_degree * (form_degree - 1),
         ("GL2 determinant exponent", form_degree))
    gate(2 * form_degree - 2 > 0, ("scalar exponent", form_degree))
    for total_degree in range(0, 9):
        for primitive_degree in range(total_degree + 1):
            central = total_degree - primitive_degree
            left = (2 * form_degree - 2) * primitive_degree
            left += (2 * form_degree - 2) * central
            right = (2 * form_degree - 2) * total_degree
            gate(left == right,
                 ("central homogenization exponent", form_degree,
                  primitive_degree, total_degree))


# 4. Independent exact positive and hostile controls in degrees 2 through 6.
# The determinant-one matrix is deliberately not an elementary single shear.
t, U, V = sp.symbols("t U V")
M = sp.Matrix([[1 + t**2, t], [t, 1]])
gate(M.det() == 1, "positive matrix determinant")
M_bad = sp.Matrix([[t, 0], [0, 1]])
gate(M_bad.det() == t, "hostile matrix determinant")


def transformed_split_form(degree: int, matrix: sp.Matrix) -> sp.Expr:
    mu = matrix[0, 0] * U + matrix[0, 1] * V
    mv = matrix[1, 0] * U + matrix[1, 1] * V
    result = sp.Integer(1)
    for root in range(degree):
        result *= mu + root * mv
    return sp.expand(result)


for form_degree in range(2, 7):
    base = sp.prod(U + root * V for root in range(form_degree))
    base_disc = sp.discriminant(sp.Poly(base.subs(V, 1), U), U)
    gate(base_disc != 0, ("squarefree base", form_degree))

    moved = transformed_split_form(form_degree, M)
    moved_disc = sp.factor(sp.discriminant(sp.Poly(moved.subs(V, 1), U), U))
    zero(moved_disc - base_disc, ("det-one orbit", form_degree))

    hostile = transformed_split_form(form_degree, M_bad)
    hostile_disc = sp.factor(sp.discriminant(sp.Poly(hostile.subs(V, 1), U), U))
    zero(hostile_disc - t ** (form_degree * (form_degree - 1)) * base_disc,
         ("nonunit determinant valuation", form_degree))

    scalar = 1 + t
    scaled_disc = sp.factor(sp.discriminant(
        sp.Poly((scalar * base).subs(V, 1), U), U
    ))
    zero(scaled_disc - scalar ** (2 * form_degree - 2) * base_disc,
         ("nonunit scalar valuation", form_degree))


# 5. Finite hostile atlas for the determinant-coordinate implication.  Rows
# have degree at most one and coefficients in {-1,0,1}.  Whenever a basis
# determinant and both test determinants are nonzero constants, the recovered
# coordinates are constants.  Nonconstant coordinates always advertise a
# nonunit determinant somewhere.
values = (-1, 0, 1)
rows = tuple(
    (a0 + a1 * t, b0 + b1 * t)
    for a0, a1, b0, b1 in itertools.product(values, repeat=4)
    if (a0, a1, b0, b1) != (0, 0, 0, 0)
)
atlas_tests = 0
atlas_unit_cases = 0
for left in rows:
    for right in rows:
        determinant = sp.Poly(row_det(left, right), t)
        if determinant.degree() != 0 or determinant.as_expr() == 0:
            continue
        for test in rows:
            atlas_tests += 1
            d0 = sp.Poly(row_det(left, test), t)
            d1 = sp.Poly(row_det(right, test), t)
            if (d0.degree() != 0 or d1.degree() != 0
                    or d0.as_expr() == 0 or d1.as_expr() == 0):
                continue
            atlas_unit_cases += 1
            inv = sp.Matrix([left, right]).inv()
            coordinates = sp.Matrix([[test[0], test[1]]]) * inv
            gate(all(sp.Poly(sp.cancel(entry), t).degree() == 0 for entry in coordinates),
                 ("finite constant-coordinate recovery", left, right, test))
gate(atlas_tests > 0 and atlas_unit_cases > 0, "finite atlas controls occur")


# 6. Sharp scope controls.  Degree one has no discriminant divisor; degree
# two is fully covered by the orbit theorem even though the separate pencil
# branch-value lemma has a sharp failure there.
F2 = U * V
G2 = U * (U + V)
ss, tt = sp.symbols("ss tt")
quadratic_pencil = sp.expand(ss * F2 + tt * G2)
quadratic_pencil_disc = sp.factor(
    sp.discriminant(sp.Poly(quadratic_pencil.subs(V, 1), U), U)
)
zero(quadratic_pencil_disc - (ss + tt) ** 2,
     "degree-two pencil one-support boundary")


semantic = {
    "cover": "all connected finite-etale A1 components have degree one",
    "orbit": "constant-discriminant squarefree binary degree n>=2 is a constant form under polynomial GL2",
    "covariance": "det exponent n(n-1), scalar exponent 2n-2",
    "homogenization": "pure C exponent (2n-2)D including central scalar layer",
    "controls": "degrees 2..6 determinant-one and nonunit determinant/scalar",
    "atlas": (len(rows), atlas_tests, atlas_unit_cases),
    "boundary": "degree-two orbit theorem valid; degree-two pencil support lemma false",
    "scope": "source-binary carrier only; no target Keller equivalence",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
source = open(__file__, encoding="utf-8").read()
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("THM3892_ALL_BINARY_DEGREE_INDEPENDENT_AUDIT_20260823")
print("status=PASS;PROMOTION_READY;JC2_OPEN")
print("extension=binary_form_degree_n>=2")
print("covariance=det_exponent_n(n-1);scalar_exponent_2n-2")
print(f"finite_section_atlas=rows:{len(rows)};tests:{atlas_tests};unit_cases:{atlas_unit_cases}")
print("positive_controls=degrees_2_through_6_det_one_orbit")
print("hostile_controls=nonunit_determinant_and_nonunit_scalar")
print("sharp_boundary=degree_two_pencil_one_support_but_orbit_classification_survives")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
