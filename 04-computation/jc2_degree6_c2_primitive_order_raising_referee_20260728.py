#!/usr/bin/env python3
"""Exact hostile referee for THM-2206's degree-six C2 order-raising test.

The tested object is the degree-six monic Faber column bank

    E6, E5, E3, E2, E1

with THM-2194's normalized boundary, first-flux/4, and second-flux/4
rows.  For each rational polynomial row we take its canonical primitive
Z_(2)-line: multiply by the unique power of two which makes every
coefficient 2-integral and at least one coefficient a 2-adic unit.  This
choice is unique up to a Z_(2)-unit and therefore gives a well-defined
map on the first C2 augmentation grade.

The exact conclusion is negative.  The three primitive rows reduce modulo
two to multiples of the E5 coordinate, so the E3 and E1 anti-invariant
classes are not raised.  At the exact square T=(1+X)^2 there is also a
primitive characteristic-zero kernel vector

    128 E5 + 32 E3 + E1.

Adding the normalized top seed E6 does not change any of the three
observables there.  Hence row saturation or alternative rowwise clearing
cannot remove this exact kernel.  A different Rees/source lattice or an
additional observable is extra structure, not a consequence of this test.

This refutes only the primitive-row test named in THM-2206.  It does not
refute every possible integral order-raising construction and proves no
Jacobian-conjecture statement.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd

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
variables = (a, c, Lambda, Omega)
degrees = (6, 5, 3, 2, 1)

K = 1 + X + c * X**2 / 2
T = sp.expand(K**2 + Lambda * X**3 + Omega * X**4)
t_coefficients = [None] + [sp.expand(T).coeff(X, index) for index in range(1, 5)]


def faber_row(degree: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Derive B_j,F_j,G_j from the exact coefficient recurrence."""
    alpha = sp.Rational(degree, 4)
    coefficient = [sp.Integer(1)]
    for index in range(1, degree + 3):
        value = sum(
            t_coefficients[i]
            * ((alpha + 1) * i - index)
            * coefficient[index - i]
            for i in range(1, min(4, index) + 1)
        ) / index
        coefficient.append(sp.expand(value))
    return (
        coefficient[degree],
        coefficient[degree + 1],
        sp.expand(coefficient[degree + 2] + coefficient[degree + 1] / 2),
    )


bank = {degree: faber_row(degree) for degree in degrees}

# THM-2194 equation (14): the three chosen-sheet contributions have powers
# a^j, a^(j+1), a^(j+2).  Monic normalization means that the columns are not
# independently rescaled.
rational_rows = {
    "boundary": [sp.expand(a**degree * bank[degree][0]) for degree in degrees],
    "phi_over_4": [sp.expand(a ** (degree + 1) * bank[degree][1]) for degree in degrees],
    "psi_over_4": [sp.expand(a ** (degree + 2) * bank[degree][2]) for degree in degrees],
}


def coefficient_list(row: list[sp.Expr]) -> list[sp.Rational]:
    answer: list[sp.Rational] = []
    for entry in row:
        if entry == 0:
            continue
        answer.extend(sp.Poly(entry, *variables, domain=sp.QQ).coeffs())
    require(answer, "zero row has no primitive line")
    return answer


def primitive_two_local_row(row: list[sp.Expr]) -> tuple[int, list[sp.Expr]]:
    """Return 2^shift*row with coefficientwise minimum v2 equal to zero."""
    minimum = min(v2_rational(coefficient) for coefficient in coefficient_list(row))
    shift = -minimum
    require(shift >= 0, "the audited Faber rows should require denominator clearing")
    primitive = [sp.expand(2**shift * entry) for entry in row]
    valuations = [v2_rational(coefficient) for coefficient in coefficient_list(primitive)]
    require(min(valuations) == 0, "row was not made Z_(2)-primitive")
    return shift, primitive


primitive_rows: dict[str, list[sp.Expr]] = {}
primitive_shifts: dict[str, int] = {}
for label, row in rational_rows.items():
    shift, primitive = primitive_two_local_row(row)
    primitive_shifts[label] = shift
    primitive_rows[label] = primitive

require(
    primitive_shifts == {"boundary": 8, "phi_over_4": 10, "psi_over_4": 9},
    "primitive two-local row scales changed",
)


def rational_mod_two(value: sp.Rational) -> int:
    value = sp.Rational(value)
    require(int(value.q) % 2 == 1, "coefficient is not 2-integral")
    return (int(value.p) % 2) * pow(int(value.q) % 2, -1, 2) % 2


def polynomial_mod_two(expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    answer = sp.Integer(0)
    for monomial, coefficient in polynomial.terms():
        residue = rational_mod_two(coefficient)
        if residue:
            term = sp.Integer(1)
            for variable, exponent in zip(variables, monomial):
                term *= variable**exponent
            answer += term
    return sp.Poly(answer, *variables, modulus=2).as_expr()


reduced_rows = {
    label: [polynomial_mod_two(entry) for entry in row]
    for label, row in primitive_rows.items()
}
expected_reductions = {
    "boundary": [0, a**5, 0, 0, 0],
    "phi_over_4": [0, a**6, 0, 0, 0],
    "psi_over_4": [0, a**7 * Lambda, 0, 0, 0],
}
for label in reduced_rows:
    require(
        all(
            sp.Poly(found - wanted, *variables, modulus=2).is_zero
            for found, wanted in zip(reduced_rows[label], expected_reductions[label])
        ),
        f"mod-two primitive row changed: {label}",
    )

# Exact-square hostile/positive algebraic control from THM-2194's own referee.
square_substitution = {a: 1, c: 0, Lambda: 0, Omega: 0}
square_bank = {
    degree: tuple(sp.Rational(value.subs(square_substitution)) for value in bank[degree])
    for degree in degrees
}
require(
    square_bank
    == {
        6: (0, 0, 0),
        5: (sp.Rational(3, 256), sp.Rational(-5, 1024), 0),
        3: (sp.Rational(-1, 16), sp.Rational(3, 128), 0),
        2: (0, 0, 0),
        1: (sp.Rational(1, 2), sp.Rational(-1, 8), 0),
    },
    "exact-square bank changed",
)

normalized_degree_six_vector = {6: 1, 5: 128, 3: 32, 2: 0, 1: 1}
for channel in range(3):
    value = sum(
        normalized_degree_six_vector[degree] * square_bank[degree][channel]
        for degree in degrees
    )
    require(value == 0, "primitive exact-square kernel stopped being a kernel")

# In the canonical generic primitive row scalings, the two nonzero specialized
# odd-column rows are [3,-16,128] and [-5,24,-128].  Their Smith invariant
# factors are 1 and 8, and their primitive integer kernel is [128,32,1].
odd_indices = (1, 2, 4)  # E5,E3,E1 in the declared column order.
square_integer_rows = [
    [int(entry.subs(square_substitution)) for index, entry in enumerate(primitive_rows[label]) if index in odd_indices]
    for label in ("boundary", "phi_over_4")
]
require(square_integer_rows == [[3, -16, 128], [-5, 24, -128]],
        "specialized primitive integer rows changed")
all_entries = [abs(value) for row in square_integer_rows for value in row if value]
d1 = 0
for value in all_entries:
    d1 = gcd(d1, value)
minor_gcd = 0
for left, right in combinations(range(3), 2):
    minor = (
        square_integer_rows[0][left] * square_integer_rows[1][right]
        - square_integer_rows[0][right] * square_integer_rows[1][left]
    )
    minor_gcd = gcd(minor_gcd, abs(minor))
require(d1 == 1 and minor_gcd == 8, "specialized Smith factors changed")
smith_factors = (d1, minor_gcd // d1)
require(smith_factors == (1, 8), "specialized Smith normal form changed")

kernel = (128, 32, 1)
require(
    all(sum(row[index] * kernel[index] for index in range(3)) == 0 for row in square_integer_rows),
    "primitive exact-square odd kernel changed",
)
require(gcd(gcd(kernel[0], kernel[1]), kernel[2]) == 1,
        "exact-square kernel is not primitive")

# Split-deck typing.  On M=Z_2 x Z_2 with swap, epsilon=(1,-1) generates
# I M, and I^j M=2^(j-1) Z_2 epsilon.  The E1 coordinate of the kernel is a
# unit, so 2^(j-1)*kernel*epsilon is never in I^(j+1), although every audited
# observable is exactly zero.
sigma = sp.Matrix([[0, 1], [1, 0]])
identity = sp.eye(2)
Delta = sigma - identity
epsilon = sp.Matrix([1, -1])
require(Delta * sp.Matrix([0, 1]) == epsilon, "epsilon does not generate I M")
for exponent in range(1, 9):
    require(Delta**exponent == (-2) ** (exponent - 1) * Delta,
            "quadratic augmentation law changed")
    require(2 ** (exponent - 1) % 2**exponent != 0,
            "graded kernel representative accidentally raised")

print("degree-6 C2 primitive-row order-raising hostile referee")
print("columns=E6,E5,E3,E2,E1:monic_Faber_normalization")
print("coefficient_ring=Z_(2)[a,c,Lambda,Omega]:formal_chosen_sheet")
print("primitive_row_scales=boundary:2^8,phi_over_4:2^10,psi_over_4:2^9")
print("mod2_boundary=[0,a^5,0,0,0]")
print("mod2_phi_over_4=[0,a^6,0,0,0]")
print("mod2_psi_over_4=[0,a^7*Lambda,0,0,0]")
print("odd_column_rank_over_F2(a,c,Lambda,Omega)=1")
print("odd_graded_survivors=E3,E1")
print("exact_square=T=(1+X)^2:a=1,c=Lambda=Omega=0")
print("exact_square_normalized_bank=E6+128E5+32E3+E1")
print("exact_square_observables=boundary:0,phi_over_4:0,psi_over_4:0")
print("exact_square_odd_kernel=[128,32,1]:primitive_E1_class=True")
print("exact_square_primitive_matrix=[[3,-16,128],[-5,24,-128],[0,0,0]]")
print("exact_square_nonzero_row_smith_factors=1,8:cokernel_has_Z/8_torsion")
print("augmentation_kernel=2^(j-1)*[128,32,1]*epsilon_nonzero_on_every_gr_j")
print("rowwise_rescaling_cannot_remove_an_exact_kernel=True")
print("verdict=CANONICAL_PRIMITIVE_ROW_TEST_FAILS_AT_DEGREE_6")
print("repair_requires=extra_saturated_Rees_source_lattice_or_additional_observable")
print("scope=FINITE_EXACT_REFUTATION_OF_STATED_TEST_NOT_ALL_INTEGRAL_ORDER_RAISING_NOT_JC2")
print("ALL CHECKS PASSED")
