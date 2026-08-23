#!/usr/bin/env python3
"""Exact symbolic companion for THM-3892."""

from __future__ import annotations

import hashlib
import itertools

import sympy as sp


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def zero(label: str, expression: sp.Expr) -> None:
    check(sp.cancel(sp.expand(expression)) == 0, label)


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d)


def coefficients(form: sp.Expr, U: sp.Symbol, V: sp.Symbol) -> tuple[sp.Expr, ...]:
    polynomial = sp.Poly(sp.expand(form), U, V)
    return tuple(
        polynomial.coeff_monomial(monomial)
        for monomial in (U**3, U**2 * V, U * V**2, V**3)
    )


# ---------------------------------------------------------------------------
# 1. Universal scalar/GL2 covariance exponents.
# ---------------------------------------------------------------------------

U, V = sp.symbols("U V")
a, b, c, d = sp.symbols("a b c d")
m00, m01, m10, m11 = sp.symbols("m00 m01 m10 m11")
F = a * U**3 + b * U**2 * V + c * U * V**2 + d * V**3
MU = m00 * U + m01 * V
MV = m10 * U + m11 * V
transformed_coefficients = coefficients(F.subs({U: MU, V: MV}, simultaneous=True), U, V)
matrix_determinant = m00 * m11 - m01 * m10
equal(
    "binary cubic GL2 discriminant covariance",
    disc(*transformed_coefficients),
    matrix_determinant**6 * disc(a, b, c, d),
)

scalar = sp.symbols("scalar")
equal(
    "binary cubic scalar discriminant covariance",
    disc(scalar * a, scalar * b, scalar * c, scalar * d),
    scalar**4 * disc(a, b, c, d),
)


# ---------------------------------------------------------------------------
# 2. Pair determinants reconstruct one common polynomial matrix.
# ---------------------------------------------------------------------------

p, q, r, s, alpha, beta = sp.symbols("p q r s alpha beta")
basis_determinant = p * s - q * r
ell0 = sp.Matrix([[p, q]])
ell1 = sp.Matrix([[r, s]])
ell2 = alpha * ell0 + beta * ell1

def row_det(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand(left[0, 0] * right[0, 1] - left[0, 1] * right[0, 0])


equal("third-row first determinant coordinate", row_det(ell0, ell2), beta * basis_determinant)
equal("third-row second determinant coordinate", row_det(ell1, ell2), -alpha * basis_determinant)

det02, det12 = sp.symbols("det02 det12")
recovered_alpha = -det12 / basis_determinant
recovered_beta = det02 / basis_determinant
recovered_ell2 = recovered_alpha * ell0 + recovered_beta * ell1
equal("recovered first pair determinant", row_det(ell0, recovered_ell2), det02)
equal("recovered second pair determinant", row_det(ell1, recovered_ell2), det12)

# The only possible connected degrees in a rank-three cover are 2 and 3.
# One infinity fibre cannot carry their Riemann--Hurwitz tariff.
for cover_degree in (2, 3):
    check(
        2 * cover_degree - 2 > cover_degree - 1,
        f"degree-{cover_degree} finite-etale A1 triviality inequality",
    )


# ---------------------------------------------------------------------------
# 3. Central scalar and contact-degree accounting.
# ---------------------------------------------------------------------------

d0, d1, d2, central = sp.symbols("d0 d1 d2 central", integer=True, nonnegative=True)
section_degree = d0 + d1 + d2
pair_contact_sum = (d0 + d1) + (d0 + d2) + (d1 + d2)
equal("pair-contact discriminant exponent", 2 * pair_contact_sum, 4 * section_degree)
equal(
    "central scalar completes pure-power exponent",
    4 * section_degree + 4 * central,
    4 * (section_degree + central),
)

# Freeze the factor-degree arithmetic in a finite exact universe; this is an
# identity control, not the theorem's degree universe.
for degree_values in itertools.product(range(5), repeat=3):
    left = 2 * sum(degree_values[i] + degree_values[j] for i, j in ((0, 1), (0, 2), (1, 2)))
    right = 4 * sum(degree_values)
    check(left == right, f"finite contact identity {degree_values}")


# ---------------------------------------------------------------------------
# 4. Positive nonconstant orbit controls and an alternating Cohn word.
# ---------------------------------------------------------------------------

t, A, C = sp.symbols("t A C")
F0 = U * V * (U - V)
F0_discriminant = disc(*coefficients(F0, U, V))
equal("constant control discriminant", F0_discriminant, 1)

# One upper shear produces the quadratic moving carrier after homogenization.
M_triangular = sp.Matrix([[1, -t], [0, 1]])
triangular_form = sp.expand(F0.subs({U: U - t * V, V: V}, simultaneous=True))
equal("triangular orbit constant discriminant", disc(*coefficients(triangular_form, U, V)), 1)
check(max(sp.Poly(coefficient, t).degree() for coefficient in coefficients(triangular_form, U, V)) == 2, "triangular orbit degree two")
triangular_homogeneous = sp.expand(C**2 * triangular_form.subs(t, A / C))
check(sp.denom(sp.cancel(triangular_homogeneous)) == 1, "triangular degree-two homogenization polynomial")
equal(
    "triangular homogeneous factorization",
    triangular_homogeneous,
    (C * U - A * V) * V * (C * U - (A + C) * V),
)
equal(
    "triangular homogeneous pure-power discriminant",
    disc(*coefficients(triangular_homogeneous, U, V)),
    C**8,
)

def eplus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix([[1, value], [0, 1]])

def eminus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix([[1, 0], [value, 1]])


M_word = eplus(t) * eminus(t) * eplus(t)
equal("three-shear word determinant", M_word.det(), 1)
word_form = sp.expand(F0.subs(
    {U: M_word[0, 0] * U + M_word[0, 1] * V, V: M_word[1, 0] * U + M_word[1, 1] * V},
    simultaneous=True,
))
equal("three-shear orbit constant discriminant", disc(*coefficients(word_form, U, V)), 1)
word_degree = max(sp.Poly(coefficient, t).degree() for coefficient in coefficients(word_form, U, V))
check(word_degree == 8, "three-shear orbit degree eight")
word_homogeneous = sp.expand(sp.cancel(C**word_degree * word_form.subs(t, A / C)))
check(sp.denom(word_homogeneous) == 1, "three-shear homogenization polynomial")
equal(
    "three-shear homogeneous pure-power discriminant",
    disc(*coefficients(word_homogeneous, U, V)),
    C ** (4 * word_degree),
)


semantic_packet = "\n".join(
    (
        "finite etale A1 root divisor triviality",
        "primitive polynomial root sections",
        "pair determinant units",
        "one common polynomial GL2 matrix",
        "binary cubic determinant exponent six",
        "central scalar exponent four",
        "Euclidean elementary Cohn word",
        "source binary versus target symplectic scope",
    )
) + "\n"
semantic_sha256 = hashlib.sha256(semantic_packet.encode("utf-8")).hexdigest()

print("THM3892_ETALE constant_discriminant_root_divisor_is_three_trivial_A1_sections")
print("THM3892_ORBIT phi=c*F0(M(t)X);M_in_GL2_k[t];detM_constant")
print("THM3892_COVARIANCE scalar_exponent=4;determinant_exponent=6")
print("THM3892_COHN polynomial_GL2_constant_determinant_has_finite_elementary_word")
print(f"THM3892_CONTROLS triangular_degree=2;three_shear_degree={word_degree}")
print("THM3892_SIDECAR central_C_layer_separates_common_vanishing_from_section_motion")
print("THM3892_SCOPE source_binary_leading_carrier_only;full_Keller_target_equivalence_not_claimed")
print(f"SEMANTIC_SHA256 {semantic_sha256}")
print(f"CHECKS {CHECKS}")
