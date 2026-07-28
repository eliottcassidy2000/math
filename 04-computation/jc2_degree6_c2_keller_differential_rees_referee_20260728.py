#!/usr/bin/env python3
"""Exact degree-six R_Q/Keller-row and first differential-grade referee.

This is the repaired follow-on to the negative primitive B/Phi/Psi audit.
It derives the fourth Faber diagonal which gives

    R_j = 4 c_(j,3) + p c_(j,1)

in THM-2129, translates it to THM-2194's monic chosen-sheet coordinates,
and takes the canonical primitive Z_(2)-line.  It then kills the free
additive constant by differentiating relative to the constant Faber seed
coefficients, saturates the resulting one-form by its common 2-content, and
computes the first C2 augmentation grade.

The outcome is negative: both the primitive R row and the stronger universal
Kahler differential envelope still see only E5 modulo two.  E3 and E1 remain
in the graded kernel.  The physical Keller equation is weaker than the
universal one-form envelope: it only contracts dR with the actual trajectory
derivation.  Thus the calculation is a stopping certificate for the first
canonical differential grade, not an inherited all-degree Rees lattice.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def v2_integer(value: int) -> int:
    require(value != 0, "v2(0) is not used")
    value = abs(value)
    order = 0
    while value % 2 == 0:
        value //= 2
        order += 1
    return order


def v2_rational(value: sp.Rational) -> int:
    value = sp.Rational(value)
    require(value != 0, "v2(0) is not used")
    return v2_integer(int(value.p)) - v2_integer(int(value.q))


X, a, c, Lambda, Omega = sp.symbols("X a c Lambda Omega")
da, dc, dLambda, dOmega = sp.symbols("da dc dLambda dOmega")
variables = (a, c, Lambda, Omega)
form_variables = variables + (da, dc, dLambda, dOmega)
degrees = (6, 5, 3, 2, 1)

K = 1 + X + c * X**2 / 2
T = sp.expand(K**2 + Lambda * X**3 + Omega * X**4)
t_coefficients = [None] + [T.coeff(X, index) for index in range(1, 5)]


def faber_diagonals(degree: int) -> list[sp.Expr]:
    """Return A_(degree+r), r=0,1,2,3, from the exact recurrence."""
    alpha = sp.Rational(degree, 4)
    coefficient = [sp.Integer(1)]
    for index in range(1, degree + 4):
        value = sum(
            t_coefficients[i]
            * ((alpha + 1) * i - index)
            * coefficient[index - i]
            for i in range(1, min(4, index) + 1)
        ) / index
        coefficient.append(sp.expand(value))
    return [coefficient[degree + offset] for offset in range(4)]


diagonals = {degree: faber_diagonals(degree) for degree in degrees}

# Derive the depressed coefficient rather than assuming it.  Before the
# shift, the relevant leading terms are
#
#   w^4 + 2 a w^3 + a^2(1+c) w^2 + ... .
#
# The depressing shift is w=Z-a/2.
w, Z, q0, r0 = sp.symbols("w Z q0 r0")
predepressed = w**4 + 2 * a * w**3 + a**2 * (1 + c) * w**2 + q0 * w + r0
depressed = sp.Poly(sp.expand(predepressed.subs(w, Z - a / 2)), Z)
require(depressed.coeff_monomial(Z**3) == 0, "depressing shift retained a cubic")
p = sp.expand(depressed.coeff_monomial(Z**2))
require(sp.expand(p - a**2 * (c - sp.Rational(1, 2))) == 0,
        "p/a^2 relation changed")

# If Z=w+h, the old and new Laurent tails obey
# d1=c1, d2=c2-h*c1, d3=c3-2h*c2+h^2*c1.  Solving for c3 gives the displayed
# translation.  Here h/a=1/2.
c1, c2, c3, h = sp.symbols("c1 c2 c3 h")
d1 = c1
d2 = c2 - h * c1
d3 = c3 - 2 * h * c2 + h**2 * c1
require(sp.expand(d3 + 2 * h * d2 + h**2 * d1 - c3) == 0,
        "third Laurent translation changed")

# Normalize by four, as THM-2194 does for Phi and Psi.  If A_r denotes the
# raw diagonal at degree+j, then
#
# R_j/4 = a^(j+3) [A_3+A_2+(1+2c)A_1/8].
H = {
    degree: sp.expand(
        diagonals[degree][3]
        + diagonals[degree][2]
        + (1 + 2 * c) * diagonals[degree][1] / 8
    )
    for degree in degrees
}

expected_H = {
    6: Lambda
    * (-4 * Lambda**2 + 6 * Lambda * c - 3 * Lambda - 12 * Omega * c + 6 * Omega)
    / 64,
    5: 5
    * (
        -512 * Lambda**2 * c
        + 512 * Lambda**2
        - 1024 * Lambda * Omega
        - 128 * Lambda * c**2
        + 128 * Lambda * c
        - 32 * Lambda
        + 1024 * Omega**2
        + 256 * Omega * c**2
        - 256 * Omega * c
        + 64 * Omega
        + 48 * c**4
        - 96 * c**3
        + 72 * c**2
        - 24 * c
        + 3
    )
    / 32768,
    3: (-48 * Lambda**2 + 8 * c**3 - 12 * c**2 + 6 * c - 1) / 512,
    2: -Lambda * (2 * c - 1) / 16,
    1: (-16 * Lambda + 32 * Omega + 4 * c**2 - 4 * c + 1) / 128,
}
require(all(sp.expand(H[j] - expected_H[j]) == 0 for j in degrees),
        "exact R/4 bank changed")

R_over_4_row = [sp.expand(a ** (degree + 3) * H[degree]) for degree in degrees]


def coefficient_list(expressions: list[sp.Expr], polynomial_variables: tuple[sp.Symbol, ...]) -> list[sp.Rational]:
    answer: list[sp.Rational] = []
    for expression in expressions:
        if expression == 0:
            continue
        answer.extend(sp.Poly(expression, *polynomial_variables, domain=sp.QQ).coeffs())
    require(answer, "zero row has no coefficient content")
    return answer


def content_v2(expressions: list[sp.Expr], polynomial_variables: tuple[sp.Symbol, ...]) -> int:
    return min(v2_rational(coefficient) for coefficient in coefficient_list(expressions, polynomial_variables))


def primitive_two_local_row(
    row: list[sp.Expr], polynomial_variables: tuple[sp.Symbol, ...]
) -> tuple[int, list[sp.Expr]]:
    minimum = content_v2(row, polynomial_variables)
    shift = -minimum
    primitive = [sp.expand(2**shift * entry) for entry in row]
    require(content_v2(primitive, polynomial_variables) == 0,
            "row was not made two-locally primitive")
    return shift, primitive


R_shift, R_primitive = primitive_two_local_row(R_over_4_row, variables)
require(R_shift == 15, "primitive R/4 scale changed")


def rational_mod_two(value: sp.Rational) -> int:
    value = sp.Rational(value)
    require(int(value.q) % 2 == 1, "coefficient is not two-integral")
    return int(value.p) % 2


def polynomial_mod_two(expression: sp.Expr, polynomial_variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    polynomial = sp.Poly(expression, *polynomial_variables, domain=sp.QQ)
    answer = sp.Integer(0)
    for monomial, coefficient in polynomial.terms():
        if rational_mod_two(coefficient):
            term = sp.Integer(1)
            for variable, exponent in zip(polynomial_variables, monomial):
                term *= variable**exponent
            answer += term
    return sp.Poly(answer, *polynomial_variables, modulus=2).as_expr()


R_mod_two = [polynomial_mod_two(entry, variables) for entry in R_primitive]
require(
    all(
        sp.Poly(found - wanted, *variables, modulus=2).is_zero
        for found, wanted in zip(R_mod_two, [0, a**8, 0, 0, 0])
    ),
    "primitive R row mod two changed",
)

# The Faber seed coefficients are constants for the base derivation.  Take
# the relative universal differential in a,c,Lambda,Omega only.  This is a
# stronger envelope than the physical Keller equation, which contracts the
# one-form with one actual trajectory derivation.
def relative_differential(expression: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(expression, a) * da
        + sp.diff(expression, c) * dc
        + sp.diff(expression, Lambda) * dLambda
        + sp.diff(expression, Omega) * dOmega
    )


dR_unsaturated = [relative_differential(entry) for entry in R_primitive]
dR_content = content_v2(dR_unsaturated, form_variables)
require(dR_content == 3, "differential R row common content changed")
dR_primitive = [sp.expand(entry / 2**dR_content) for entry in dR_unsaturated]
require(content_v2(dR_primitive, form_variables) == 0,
        "differential R row was not saturated")
dR_mod_two = [polynomial_mod_two(entry, form_variables) for entry in dR_primitive]
require(
    all(
        sp.Poly(found - wanted, *form_variables, modulus=2).is_zero
        for found, wanted in zip(
            dR_mod_two,
            [0, a**7 * da + a**8 * dc, 0, 0, 0],
        )
    ),
    "primitive differential R row mod two changed",
)

# Reconstruct the first three primitive rows to compare their exact column
# depths with the new R and dR rows.
B_row = [sp.expand(a**degree * diagonals[degree][0]) for degree in degrees]
F_row = [sp.expand(a ** (degree + 1) * diagonals[degree][1]) for degree in degrees]
G_row = [
    sp.expand(a ** (degree + 2) * (diagonals[degree][2] + diagonals[degree][1] / 2))
    for degree in degrees
]
B_shift, B_primitive = primitive_two_local_row(B_row, variables)
F_shift, F_primitive = primitive_two_local_row(F_row, variables)
G_shift, G_primitive = primitive_two_local_row(G_row, variables)
require((B_shift, F_shift, G_shift) == (8, 10, 9), "old primitive row scales changed")

odd_indices = (1, 2, 4)  # E5,E3,E1.
rows_for_depth = {
    "B": (B_primitive, variables),
    "F": (F_primitive, variables),
    "G": (G_primitive, variables),
    "R": (R_primitive, variables),
    "dR": (dR_primitive, form_variables),
}
odd_depths = {
    label: tuple(content_v2([row[index]], row_variables) for index in odd_indices)
    for label, (row, row_variables) in rows_for_depth.items()
}
require(
    odd_depths
    == {
        "B": (0, 4, 7),
        "F": (0, 3, 7),
        "G": (0, 4, 7),
        "R": (0, 6, 8),
        "dR": (0, 4, 7),
    },
    "odd-column two-adic tooth profile changed",
)
combined_odd_depth = tuple(min(profile[index] for profile in odd_depths.values()) for index in range(3))
require(combined_odd_depth == (0, 3, 7), "combined odd-column tooth profile changed")

# Exact-square stopping control.  It is the same algebraic square used in
# THM-2194's own referee.  Raw diagonals r=0,1,2 vanish on the primitive odd
# vector, while r=3 is 1/256.  The resulting R value is 1/64, but it is a free
# additive constant and hence has zero derivative/principal part.
square = {a: 1, c: 0, Lambda: 0, Omega: 0}
kernel_weights = {5: 128, 3: 32, 1: 1}
raw_sums = tuple(
    sp.factor(sum(kernel_weights[j] * diagonals[j][offset].subs(square) for j in kernel_weights))
    for offset in range(4)
)
require(raw_sums == (0, 0, 0, sp.Rational(1, 256)),
        "exact-square raw diagonal word changed")
R_over_4_kernel = sp.expand(
    sum(kernel_weights[j] * a ** (j + 3) * H[j] for j in kernel_weights)
)
require(R_over_4_kernel.subs(square) == sp.Rational(1, 256),
        "exact-square R/4 value changed")
R_kernel = sp.expand(4 * R_over_4_kernel)
require(R_kernel.subs(square) == sp.Rational(1, 64),
        "exact-square physical R value changed")
R_gradient = tuple(sp.diff(R_kernel, variable).subs(square) for variable in variables)
require(R_gradient == (sp.Rational(1, 2), sp.Rational(-1, 2), -3, 6),
        "exact-square R tangent changed")

# Literal tangent controls: a nonconstant a-direction sees the R row, whereas
# the constant fibre and the da=dc tangent kill it.  These controls emphasize
# that physical detection depends on a trajectory tangent; it is not seed
# elimination in the augmentation grade.
def contract_gradient(tangent: tuple[int, int, int, int]) -> sp.Rational:
    return sp.factor(sum(left * right for left, right in zip(R_gradient, tangent)))


require(contract_gradient((1, 0, 0, 0)) == sp.Rational(1, 2),
        "positive R tangent control failed")
require(contract_gradient((0, 0, 0, 0)) == 0,
        "constant-fibre R control failed")
require(contract_gradient((1, 1, 0, 0)) == 0,
        "da=dc hostile tangent stopped cancelling")

print("degree-6 C2 R_Q / differential-Rees hostile referee")
print("columns=E6,E5,E3,E2,E1:monic_Faber_normalization")
print("depressed_relation=p/a^2=c-1/2")
print("Rj_over_4=a^(j+3)*(A3+A2+(1+2c)A1/8)")
print("primitive_R_over_4_scale=2^15")
print("mod2_R=[0,a^8,0,0,0]")
print("relative_differential_constants=Faber_seed_coefficients")
print("d_of_primitive_R_common_content=2^3")
print("primitive_dR_effective_scale=2^12*d(R_over_4)")
print("mod2_primitive_dR=[0,a^7*da+a^8*dc,0,0,0]")
print("odd_depths_B=E5:0,E3:4,E1:7")
print("odd_depths_F=E5:0,E3:3,E1:7")
print("odd_depths_G=E5:0,E3:4,E1:7")
print("odd_depths_R=E5:0,E3:6,E1:8")
print("odd_depths_dR=E5:0,E3:4,E1:7")
print("combined_odd_tooth_profile=E5:0,E3:3,E1:7")
print("first_differential_grade_odd_rank=1:E3,E1_survive")
print("exact_square_raw_diagonals_r0_r1_r2_r3=0,0,0,1/256")
print("exact_square_R_over_4=1/256:physical_R=1/64")
print("exact_square_R_is_free_constant=zero_in_derivative_and_finite_place_principal_part")
print("exact_square_R_gradient=(1/2,-1/2,-3,6)")
print("positive_tangent_da=1:R_derivative=1/2")
print("hostile_tangent_da=dc=1:R_derivative=0")
print("physical_constraint=contraction_with_actual_derivation_not_universal_one_form")
print("verdict=FIRST_CANONICAL_DIFFERENTIAL_GRADE_STILL_FAILS_TO_RAISE_E3_E1")
print("later_Rees_weights_E3:3,E1:7=EXTRA_STRUCTURE_NOT_INHERITED")
print("scope=FINITE_EXACT_STOPPING_CERTIFICATE_NOT_ALL_REES_LATTICES_NOT_JC2")
print("ALL CHECKS PASSED")
