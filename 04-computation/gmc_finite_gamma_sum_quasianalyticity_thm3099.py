#!/usr/bin/env python3
"""Exact finite companion for THM-3099's Gamma-sum zero dichotomy."""

from fractions import Fraction
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(matrix):
    matrix = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(matrix)
    answer = Fraction(1)
    sign = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign *= -1
        value = matrix[column][column]
        answer *= value
        for entry in range(column, size):
            matrix[column][entry] /= value
        for row in range(column + 1, size):
            multiple = matrix[row][column]
            for entry in range(column, size):
                matrix[row][entry] -= multiple * matrix[column][entry]
    return sign * answer


def factorial_term(n, factors=(), base=1):
    """base^n times product factorial(a*n+b)^exponent, exactly."""
    value = Fraction(base**n)
    for slope, shift, exponent in factors:
        argument = slope * n + shift
        require(argument >= 0, "factorial argument")
        value *= Fraction(factorial(argument)) ** exponent
    return value


def direct_ratio(n, factors=(), base=1):
    return factorial_term(n + 1, factors, base) / factorial_term(
        n, factors, base
    )


def product_ratio(n, factors=(), base=1):
    value = Fraction(base)
    for slope, shift, exponent in factors:
        block = 1
        for offset in range(1, slope + 1):
            block *= slope * n + shift + offset
        value *= Fraction(block) ** exponent
    return value


# Gamma/factorial monomials have the claimed positive rational shift ratio.
term_bank = (
    (((2, 0, 1), (1, 0, -2)), 1),
    (((3, 2, 1), (1, 0, -1), (1, 1, -2)), 1),
    (((4, 1, 2), (2, 0, -3), (1, 2, 1)), 3),
    (((5, 3, 1), (3, 1, -1), (1, 0, -2)), 2),
)
ratio_cells = 0
for factors, base in term_bank:
    for n in range(3, 18):
        ratio = direct_ratio(n, factors, base)
        require(ratio == product_ratio(n, factors, base), "shift quotient")
        require(ratio > 0, "positive shift quotient")
        ratio_cells += 1
require(ratio_cells == 60, "ratio census")


# Exact proportional terms represent the identity branch, not a flat ghost.
identity_cells = 0
for n in range(2, 27):
    central = factorial_term(n, ((2, 0, 1), (1, 0, -2)))
    half = factorial_term(n, ((2, -1, 1), (1, 0, -1), (1, -1, -1)))
    require(central == 2 * half, "exact Gamma proportionality")
    identity_cells += 1
require(identity_cells == 25, "identity census")


# A fixed-shift evaluation matrix is rational after column normalization.
# The three comparable terms are n!, (n+1)!, and (n+2)!.
shift_cells = 0
shifts = (0, 1, 2)
comparable = (
    ((1, 0, 1),),
    ((1, 1, 1),),
    ((1, 2, 1),),
)
for n in range(3, 23):
    matrix = []
    for shift in shifts:
        matrix.append(
            [
                factorial_term(n + shift, factors)
                / factorial_term(n, factors)
                for factors in comparable
            ]
        )
    require(determinant(matrix) != 0, "rational shift determinant")
    shift_cells += 1
require(shift_cells == 20, "shift determinant census")


# One good value does not preserve orientation: only eventual nonvanishing.
sign_values = [10 * 2**n - 3**n for n in range(1, 31)]
require(sign_values[0] > 0 and sign_values[4] > 0, "positive early hostile")
require(all(value < 0 for value in sign_values[5:]), "eventual negative hostile")


# If positive quotients are dropped, root-of-unity oscillation gives
# infinitely many zeros without an identity.
oscillation_cells = 0
for n in range(24):
    value = 1 + (-1) ** n
    require(value == (2 if n % 2 == 0 else 0), "oscillation boundary")
    oscillation_cells += 1
require(oscillation_cells == 24, "oscillation census")


# Every remote coefficient shape (k*C+b)!/(C!)^k has a positive rational
# quotient.  k=0 is the constant lower coefficient.
remote_coefficient_cells = 0
for degree in range(1, 7):
    for copies in range(degree + 1):
        for fixed_offset in (0, 1, 4, 11):
            factors = (
                ()
                if copies == 0
                else ((copies, fixed_offset, 1), (1, 0, -copies))
            )
            for terminal in (13, 19, 31):
                ratio = direct_ratio(terminal, factors)
                require(ratio == product_ratio(terminal, factors), "remote ratio")
                require(ratio > 0, "remote ratio positivity")
                remote_coefficient_cells += 1
require(remote_coefficient_cells == 324, "remote coefficient census")


# Exact two-slot physical remote resultant.  For support {a,C}, eliminating
# F1=x+y gives C(2a,a)-2C(a+C,a)+C(2C,C), which is a finite Gamma sum.
remote_resultant_cells = 0
for lower in range(1, 7):
    for terminal in range(lower + 1, 32):
        resultant = (
            comb(2 * lower, lower)
            - 2 * comb(lower + terminal, lower)
            + comb(2 * terminal, terminal)
        )
        require(resultant > 0, "two-slot physical resultant")
        remote_resultant_cells += 1
require(remote_resultant_cells == 165, "remote resultant census")


print("THM-3099 FINITE GAMMA-SUM QUASIANALYTICITY")
print(f"ratio_cells={ratio_cells} positive_rational_shift=PASS")
print(f"identity_cells={identity_cells} exact_dependency=PASS")
print(f"shift_cells={shift_cells} rational_evaluation_determinant=PASS")
print("sign_hostile=10*2^n-3^n;good_sample_does_not_preserve_sign")
print(f"oscillation_cells={oscillation_cells} positivity_boundary=PASS")
print(f"remote_coefficient_cells={remote_coefficient_cells} factorial_shapes=PASS")
print(f"remote_resultant_cells={remote_resultant_cells} two_slot_physical=PASS")
print("remote_dichotomy=exact_zero_sequence_or_eventual_fixed_nonzero_sign")
print("one_exact_good_extension=eventual_remote_nonvanishing")
print("boundary=one_parameter;no_effective_threshold;no_sign_preservation")
print("all_exact_checks=PASS")
