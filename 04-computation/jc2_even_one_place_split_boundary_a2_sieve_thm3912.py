#!/usr/bin/env python3
"""Assertion-free exact lattice gates for the THM-3912 proof candidate.

The companion freezes the even one-place boundary Gram/SNF calculation,
blowup congruence, A2 discriminant pairing, mod-18 order-three table, and
the sharp degree-six specialization.

Reproduction:
  python3 04-computation/jc2_even_one_place_split_boundary_a2_sieve_thm3912.py
  python3 -O 04-computation/jc2_even_one_place_split_boundary_a2_sieve_thm3912.py
"""

from __future__ import annotations

import ast
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"{label}: failed")


def equal(label: str, left: object, right: object) -> None:
    gate(left == right, f"{label}: {left!r} != {right!r}")


def mod_one(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


# ---------------------------------------------------------------------------
# 1. Universal boundary block K_d for d=2m.
# ---------------------------------------------------------------------------

m = sp.symbols("m", integer=True, positive=True)
K = sp.Matrix([[1 - m, m], [m, 1 - m]])
equal("boundary_determinant", sp.factor(K.det()), 1 - 2 * m)
equal("boundary_sum_square", sum(K), 2)
equal("boundary_sum_eigenvector", K * sp.Matrix([1, 1]), sp.Matrix([1, 1]))
equal("boundary_difference_eigenvector", K * sp.Matrix([1, -1]), (1 - 2 * m) * sp.Matrix([1, -1]))

for degree in range(2, 42, 2):
    half = degree // 2
    Kd = sp.Matrix([[1 - half, half], [half, 1 - half]])
    smith = smith_normal_form(Kd, domain=ZZ)
    diagonal = sorted(abs(int(smith[i, i])) for i in range(2))
    equal(f"boundary_snf_d{degree}", diagonal, [1, degree - 1])


# ---------------------------------------------------------------------------
# 2. Boundary blowups add only a unimodular (-1) summand.
# ---------------------------------------------------------------------------

s, t = sp.symbols("s t", integer=True)
old = sp.Matrix([[s, t], [t, s]])
blown = sp.Matrix(
    [
        [s - 1, t - 1, 1],
        [t - 1, s - 1, 1],
        [1, 1, -1],
    ]
)
# Columns express pullback(B+), pullback(B-), E in the basis
# strict(B+), strict(B-), E.
T = sp.Matrix([[1, 0, 0], [0, 1, 0], [1, 1, 1]])
target = sp.diag(1, 1, -1)
target[:2, :2] = old
equal("intersection_blowup_congruence", sp.simplify(T.T * blown * T), target)
equal("intersection_blowup_basis_unimodular", int(T.det()), 1)

# A blowup at a smooth point of one removed component has the same property.
one_blown = sp.Matrix([[s - 1, 1], [1, -1]])
T_one = sp.Matrix([[1, 0], [1, 1]])
equal("smooth_boundary_blowup_congruence", T_one.T * one_blown * T_one, sp.diag(s, -1))
equal("smooth_boundary_blowup_basis_unimodular", int(T_one.det()), 1)


# ---------------------------------------------------------------------------
# 3. The negative A2 discriminant generator and the intrinsic triple.
# ---------------------------------------------------------------------------

A2 = sp.Matrix([[-2, 1], [1, -2]])
g = sp.Matrix([sp.Rational(1, 3), sp.Rational(2, 3)])
equal("A2_determinant", int(A2.det()), 3)
equal("A2_generator_pairings", A2 * g, sp.Matrix([0, -1]))
equal("A2_generator_square", (g.T * A2 * g)[0], -sp.Rational(2, 3))
equal("three_A2_diagonal_square", 3 * (g.T * A2 * g)[0], -2)

# Support-size k has self-pairing k/3 modulo Z, since -2/3 = 1/3 mod Z.
for support_size in range(7):
    q = mod_one(Fraction(-2 * support_size, 3))
    equal(
        f"A2_support_{support_size}_isotropic",
        q == 0,
        support_size % 3 == 0,
    )


# ---------------------------------------------------------------------------
# 4. Boundary order-three self-pairing and the exact mod-18 permission table.
# ---------------------------------------------------------------------------

def boundary_order_three_q(degree: int) -> Fraction | None:
    n = degree - 1
    if n % 3:
        return None
    half = degree // 2
    return mod_one(Fraction(n * (half - 1), 9))


expected_q = {4: Fraction(1, 3), 10: Fraction(0), 16: Fraction(2, 3)}
for degree in range(2, 110, 2):
    q_boundary = boundary_order_three_q(degree)
    if degree % 18 in expected_q:
        equal(f"boundary_q_mod18_d{degree}", q_boundary, expected_q[degree % 18])
    else:
        equal(f"boundary_no_order_three_d{degree}", q_boundary, None)


def isotropy_permitted(degree: int, a2_blocks: int) -> bool:
    # Any three A2 blocks already contain the intrinsic diagonal class.
    if a2_blocks >= 3:
        return True
    q_boundary = boundary_order_three_q(degree)
    # With no boundary 3-part, support size must be a positive multiple of 3.
    if q_boundary is None:
        return False
    # The boundary coordinate may be zero or nonzero.  We seek a nonzero
    # class, so k=0 with boundary zero is not counted.
    for use_boundary in (False, True):
        for support_size in range(a2_blocks + 1):
            if not use_boundary and support_size == 0:
                continue
            q = Fraction(support_size, 3)
            if use_boundary:
                q += q_boundary
            if mod_one(q) == 0:
                return True
    return False


permission_table = {
    residue: [isotropy_permitted(residue if residue else 18, r) for r in range(4)]
    for residue in range(0, 18, 2)
}
expected_permission_table = {
    0: [False, False, False, True],
    2: [False, False, False, True],
    4: [False, False, True, True],
    6: [False, False, False, True],
    8: [False, False, False, True],
    10: [True, True, True, True],
    12: [False, False, False, True],
    14: [False, False, False, True],
    16: [False, True, True, True],
}
equal("complete_even_mod18_permission_table", permission_table, expected_permission_table)


# ---------------------------------------------------------------------------
# 5. The degree-six specialization behind THM-3911.
# ---------------------------------------------------------------------------

K6 = sp.Matrix([[-2, 3], [3, -2]])
equal("degree_six_boundary_determinant", int(K6.det()), -5)
equal("degree_six_boundary_snf", smith_normal_form(K6, domain=ZZ), sp.diag(1, 5))
equal("degree_six_two_A2_isotropy", isotropy_permitted(6, 2), False)

# One ordinary-quadruple elliptic (-2) block and two A1 (-2) blocks are all
# prime to three.  The complete sharp removed lattice has determinant 360.
sharp_determinant = abs(int(K6.det())) * 2 * (int(A2.det()) ** 2) * (2**2)
equal("sharp_removed_lattice_determinant", sharp_determinant, 360)

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "boundary": "K_d=[[1-m,m],[m,1-m]], SNF=(1,d-1), d=2m",
    "blowups": "removed-boundary blowups add an orthogonal -1 by unimodular congruence",
    "A2": "negative A2 generator has q=-2/3; three-block diagonal is integrally isotropic",
    "mod18": "boundary order-3 q is 1/3,0,2/3 at d=4,10,16 mod18",
    "sieve": "Pic[3]=0 plus no order-3 isotropic discriminant class implies Cl(open)[3]=0",
    "scope": "permission means only lattice-level possibility, not saturation or cover existence",
    "degree6": "boundary determinant -5 and two A2 blocks are 3-anisotropic",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3912-even-one-place-split-boundary-a2-three-torsion-design-sieve")
print("boundary_gram=[[1-m,m],[m,1-m]]")
print("boundary_snf=[1,d-1]")
print("boundary_determinant=1-d")
print("boundary_order3_q_by_d_mod18={4:1/3,10:0,16:2/3}")
print("a2_intrinsic_isotropic_support_min=3")
print("permission_rows_r0_to_r3=d0:[0,0,0,1];d4:[0,0,1,1];d10:[1,1,1,1];d16:[0,1,1,1]")
print("degree6_two_A2_three_torsion=EXCLUDED_BY_LATTICE")
print("permission_is_not_globalization=YES")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
