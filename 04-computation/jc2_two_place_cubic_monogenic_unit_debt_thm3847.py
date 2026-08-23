#!/usr/bin/env python3
"""Assertion-free exact gates for the THM-3847 proof candidate.

Reproduction:
  python3 04-computation/jc2_two_place_cubic_monogenic_unit_debt_thm3847.py
  python3 -O 04-computation/jc2_two_place_cubic_monogenic_unit_debt_thm3847.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, beta, T, Z, X, Y, W, eps = sp.symbols(
    "A C beta T Z X Y W eps"
)
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    if isinstance(value, sp.MatrixBase):
        entries = [sp.factor(entry) for entry in value]
        if any(entry != 0 for entry in entries):
            raise AssertionError(f"{label}: {value}")
        return
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


# The constant-profile deformation of the depressed m=2 tower.
p = sp.Rational(3, 2) + A * C
u = 1 + A * C + A**2 * beta
Delta = (
    -27 * A**2 * beta**2
    + 8 * A * C**3
    - 54 * A * C * beta
    + 9 * C**2
    - 54 * beta
)
F = T**3 - 2 * p * T + 2 * u

equal("tower_discriminant", sp.discriminant(F, T), 4 * A**2 * Delta)
equal("tower_cusp_identity", 8 * p**3 - 27 * u**2, A**2 * Delta)

# The branch is an irreducible quadratic in A.  Its quadratic discriminant is
# a cube of the squarefree profile L, which produces a conic normalization.
L = -2 * C**2 + 9 * beta
disc_A = sp.discriminant(Delta, A)
equal("branch_quadratic_discriminant", disc_A, -8 * L**3)
equal("profile_squarefree_discriminant", sp.discriminant(L, C), 72 * beta)
equal("branch_A_degree", sp.degree(Delta, A), 2)
equal("branch_A_lead", sp.LC(sp.Poly(Delta, A), A), -27 * beta**2)

# If eps^2=-8 and W^2=L, this is the finite birational normalization map.
Bquad = 8 * C**3 - 54 * beta * C
A_norm = (Bquad - eps * L * W) / (54 * beta**2)
norm_num = sp.together(Delta.subs(A, A_norm)).as_numer_denom()[0]
norm_reduced = sp.reduced(
    sp.Poly(norm_num, eps, W, domain=sp.QQ.frac_field(C, beta)),
    [
        sp.Poly(eps**2 + 8, eps, W, domain=sp.QQ.frac_field(C, beta)),
        sp.Poly(W**2 - L, eps, W, domain=sp.QQ.frac_field(C, beta)),
    ],
)[1]
equal("conic_normalization_map", norm_reduced.as_expr(), 0)
equal("conic_projective_infinity_degree", sp.degree(W**2 + 2 * C**2, W), 2)
nonzero(
    "conic_projective_infinity_separable",
    sp.discriminant(W**2 + 2 * C**2, W),
)

# A monic cubic presentation of the finite order, together with the depressed
# tower element t.  Reduction is always performed in k(A,C,beta)[Z]/(G).
G = (
    Z**3
    - 4 * C * Z**2
    + (4 * C**2 + 6 * beta) * Z
    - 4 * A * beta**2
    - 12 * beta * C
)
t = (2 * C * Z - Z**2 - 4 * beta) / (2 * beta)


def reduce_G(value: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(value).as_numer_denom()
    remainder = sp.rem(numerator, G, Z)
    return sp.factor(remainder / denominator)


equal("depressed_tower_relation", reduce_G(F.subs(T, t)), 0)
equal("integral_z_relation", reduce_G(t**2 + t - 2 - A * Z), 0)
equal(
    "mixed_multiplication",
    reduce_G(t * Z),
    Z + 2 * C * t - 2 * C - 2 * A * beta,
)
equal(
    "z_square_multiplication",
    reduce_G(Z**2),
    2 * C * Z - 2 * beta * t - 4 * beta,
)
equal("z_characteristic_polynomial", reduce_G(G), 0)
equal("power_discriminant", sp.discriminant(G, Z), 16 * beta**2 * Delta)

# Verify the full free multiplication table in the basis (1,t,z), including
# all 27 associativity triples.  Vectors are coefficient columns in that basis.
e0 = sp.Matrix([1, 0, 0])
e1 = sp.Matrix([0, 1, 0])
e2 = sp.Matrix([0, 0, 1])
basis = (e0, e1, e2)
table = {
    (0, 0): e0,
    (0, 1): e1,
    (0, 2): e2,
    (1, 1): sp.Matrix([2, -1, A]),
    (1, 2): sp.Matrix([-2 * C - 2 * A * beta, 2 * C, 1]),
    (2, 2): sp.Matrix([-4 * beta, -2 * beta, 2 * C]),
}


def basis_product(i: int, j: int) -> sp.Matrix:
    return table[(min(i, j), max(i, j))]


def multiply(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    ans = sp.zeros(3, 1)
    for i in range(3):
        for j in range(3):
            ans += left[i] * right[j] * basis_product(i, j)
    return ans.applyfunc(sp.factor)


for i in range(3):
    for j in range(3):
        for k in range(3):
            left = multiply(multiply(basis[i], basis[j]), basis[k])
            right = multiply(basis[i], multiply(basis[j], basis[k]))
            equal(f"associativity_{i}{j}{k}", left, right)


def multiplication_matrix(value: sp.Matrix) -> sp.Matrix:
    return sp.Matrix.hstack(*(multiply(value, entry) for entry in basis))


M_t = multiplication_matrix(e1)
M_z = multiplication_matrix(e2)
equal("t_characteristic_polynomial", M_t.charpoly(T).as_expr(), F)
equal("z_characteristic_from_table", M_z.charpoly(Z).as_expr(), G)

# The trace discriminant in (1,t,z), and its binary index form.
trace_matrix = sp.zeros(3)
for i in range(3):
    for j in range(3):
        trace_matrix[i, j] = sp.trace(multiplication_matrix(multiply(basis[i], basis[j])))
equal("trace_discriminant", sp.det(trace_matrix), 4 * Delta)

theta = X * e1 + Y * e2
theta2 = multiply(theta, theta)
index_form = sp.factor(sp.det(sp.Matrix.hstack(e0, theta, theta2)))
index_expected = A * X**3 + 3 * X**2 * Y - 2 * C * X * Y**2 + 2 * beta * Y**3
equal("binary_index_form", index_form, index_expected)
equal("constant_index_value", index_form.subs({X: 0, Y: 1}), 2 * beta)
equal(
    "index_discriminant_identity",
    sp.discriminant(G, Z),
    index_form.subs({X: 0, Y: 1}) ** 2 * sp.det(trace_matrix),
)

# The global generator's different is nonconstant.  On every etale affine
# open it becomes a unit, which is the exact constant-unit debt.
different = sp.diff(G, Z)
equal("different_formula", different, 3 * Z**2 - 8 * C * Z + 4 * C**2 + 6 * beta)
nonzero("different_nonconstant_Z", sp.diff(different, Z))
nonzero("different_nonconstant_C", sp.diff(different, C))

# Hostile control: making the profile nonconstant can lower the quadratic
# radicand, but simultaneously destroys the unit index and adds a branch
# component.  It is a design signal, not a uniform puncture theorem.
b_signal = C * (2 * C + 1) / 9
Delta_signal = sp.factor(Delta.subs(beta, b_signal))
equal("nonconstant_signal_radicand", L.subs(beta, b_signal), C)
equal("nonconstant_signal_extra_component", Delta_signal.subs(C, 0), 0)
equal("nonconstant_signal_index", (2 * beta).subs(beta, b_signal), 2 * b_signal)
nonzero("nonconstant_signal_index_nonunit", sp.diff(2 * b_signal, C))

print("THM3847_TOWER_P", p)
print("THM3847_TOWER_U", u)
print("THM3847_BRANCH", Delta)
print("THM3847_BRANCH_A_DISCRIMINANT", sp.factor(disc_A))
print("THM3847_BRANCH_NORMALIZATION", "W^2=-2*C^2+9*beta")
print("THM3847_BRANCH_PLACES_AT_INFINITY", 2)
print("THM3847_CUBIC_GENERATOR", G)
print("THM3847_BINARY_INDEX", index_form)
print("THM3847_INDEX_REPRESENTS_ONE", "Y=(2*beta)^(-1/3)")
print("THM3847_TRACE_DISCRIMINANT", sp.factor(sp.det(trace_matrix)))
print("THM3847_DIFFERENT", different)
print("THM3847_NONCONSTANT_SIGNAL", sp.factor(Delta_signal))
print(
    "THM3847_SCOPE",
    "constant beta theorem; nonconstant profile is design-only; no JC(2) conclusion",
)
semantic_packet = (
    "constant-profile depressed cubic",
    "irreducible branch with conic normalization and two infinity places",
    "normal finite-flat cubic order",
    "binary index A X^3+3 X^2Y-2CXY^2+2beta Y^3 represents one",
    "global monogenic generator and nonconstant different unit debt",
    "nonconstant profile remains design-only",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
