#!/usr/bin/env python3
"""Exact hostile probe for Khinchin-scalar transfers into LRC and JC(2).

The script separates three objects that are often conflated:

1. the finite geometric mean of continued-fraction digits;
2. the ordered continuant / constant SL2 word;
3. the target-specific sidecar (LRC midpoint phase and ties, or JC repair
   multiplier holonomy).

All arithmetic checks are exact.  The bounded universes are declared below;
the two main no-go families are symbolic and merely replayed on finite ranges.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, prod


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def continued_fraction_digits(numerator: int, denominator: int) -> tuple[int, ...]:
    """Canonical digits after the initial zero for 0 < numerator < denominator."""

    require(0 < numerator < denominator, "fraction must lie in (0,1)")
    require(gcd(numerator, denominator) == 1, "fraction must be reduced")
    digits: list[int] = []
    quotient_numerator = denominator
    quotient_denominator = numerator
    while quotient_denominator:
        digit, remainder = divmod(quotient_numerator, quotient_denominator)
        digits.append(digit)
        quotient_numerator, quotient_denominator = quotient_denominator, remainder
    require(digits[-1] > 1 or len(digits) == 1, "canonical final digit")
    return tuple(digits)


def fraction_from_digits(digits: tuple[int, ...]) -> Fraction:
    require(bool(digits), "continued fraction needs at least one digit")
    value = Fraction(0, 1)
    for digit in reversed(digits):
        value = Fraction(1, digit + value)
    return value


def finite_khinchin_signature(digits: tuple[int, ...]) -> tuple[int, int]:
    return len(digits), prod(digits)


def midpoint_phase_trace(digits: tuple[int, ...]) -> tuple[int, ...]:
    """THM-778 phase recursion s'=(a-s) mod 2, starting at s=1."""

    phase = 1
    trace = [phase]
    for digit in digits:
        phase = (digit - phase) % 2
        trace.append(phase)
    return tuple(trace)


def centered_endpoint_blocks(first_speed: int, second_speed: int) -> tuple[str, ...]:
    events: dict[Fraction, list[str]] = defaultdict(list)
    for owner, speed in (("u", first_speed), ("v", second_speed)):
        for index in range(speed):
            events[Fraction(2 * index + 1, 2 * speed)].append(owner)
    return tuple("".join(sorted(events[time])) for time in sorted(events))


def tie_count(first_speed: int, second_speed: int) -> int:
    return sum(len(block) == 2 for block in centered_endpoint_blocks(first_speed, second_speed))


def multiply_2x2(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def determinant_2x2(matrix: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def elementary_plus(parameter: int) -> tuple[tuple[int, int], tuple[int, int]]:
    return ((1, parameter), (0, 1))


def elementary_minus(parameter: int) -> tuple[tuple[int, int], tuple[int, int]]:
    return ((1, 0), (parameter, 1))


def alternating_constant_word(digits: tuple[int, ...]) -> tuple[tuple[int, int], tuple[int, int]]:
    matrix = ((1, 0), (0, 1))
    for index, digit in enumerate(digits):
        factor = elementary_plus(digit) if index % 2 == 0 else elementary_minus(digit)
        matrix = multiply_2x2(matrix, factor)
    return matrix


def determinant_fraction_matrix(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    determinant = Fraction(1, 1)
    size = len(work)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0, 1)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant *= -1
        pivot_value = work[column][column]
        determinant *= pivot_value
        for entry in range(column, size):
            work[column][entry] /= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column]
            if factor:
                for entry in range(column, size):
                    work[row][entry] -= factor * work[column][entry]
    return determinant


def cyclic_transport_determinant(multipliers: tuple[Fraction, ...]) -> Fraction:
    size = len(multipliers)
    matrix = [[Fraction(0, 1) for _ in range(size)] for _ in range(size)]
    for index, multiplier in enumerate(multipliers):
        matrix[index][index] = multiplier
        matrix[index][(index - 1) % size] += 1
    return determinant_fraction_matrix(matrix)


checks = 0

print("=== Rational finite-mean ambiguity ===")
for denominator in range(2, 65):
    canonical = (denominator,)
    alternate = (denominator - 1, 1)
    require(fraction_from_digits(canonical) == Fraction(1, denominator), "canonical rational CF")
    require(fraction_from_digits(alternate) == Fraction(1, denominator), "alternate rational CF")
    require(finite_khinchin_signature(canonical) != finite_khinchin_signature(alternate), "mean ambiguity")
    checks += 3
print("1/14=[0;14]=[0;13,1], with finite digit products 14 and 13.")
print("A rational finite Khinchin mean is convention-dependent before canonicalization.")


print("\n=== Infinite LRC endpoint hostile ===")
family_checks = 0
for odd_parameter in range(3, 102, 2):
    denominator = 2 * odd_parameter + 1
    even_numerator = 2 * odd_parameter
    odd_numerator = odd_parameter
    even_digits = (1, 2 * odd_parameter)
    odd_digits = (2, odd_parameter)
    require(fraction_from_digits(even_digits) == Fraction(even_numerator, denominator), "even family fraction")
    require(fraction_from_digits(odd_digits) == Fraction(odd_numerator, denominator), "odd family fraction")
    require(finite_khinchin_signature(even_digits) == finite_khinchin_signature(odd_digits), "same finite mean")
    require(midpoint_phase_trace(even_digits)[-1] == midpoint_phase_trace(odd_digits)[-1] == 0, "same final phase")
    require(tie_count(even_numerator, denominator) == 0, "even/odd no tie")
    require(tie_count(odd_numerator, denominator) == 1, "odd/odd one tie")
    family_checks += 6
checks += family_checks
print("For every odd m>=3:")
print("  [0;1,2m]=2m/(2m+1) and [0;2,m]=m/(2m+1)")
print("  share length 2, digit product 2m, denominator 2m+1, and final midpoint phase 0;")
print("  the first centered pair has no tie, while the second has exactly one tie.")
print(f"Bounded replay: {family_checks} exact gates for odd 3<=m<=101.")


print("\n=== Canonical bounded LRC census ===")
groups: dict[tuple[int, int, int, int], list[tuple[int, tuple[int, ...], int]]] = defaultdict(list)
for denominator in range(2, 201):
    for numerator in range(1, denominator):
        if gcd(numerator, denominator) != 1:
            continue
        digits = continued_fraction_digits(numerator, denominator)
        key = (
            len(digits),
            prod(digits),
            denominator,
            midpoint_phase_trace(digits)[-1],
        )
        groups[key].append((numerator, digits, tie_count(numerator, denominator)))
        checks += 2
mixed_tie_groups = [
    (key, values)
    for key, values in groups.items()
    if {value[2] for value in values} == {0, 1}
]
mixed_tie_groups.sort(key=lambda item: item[0][2:3] + item[0][:2])
require(bool(mixed_tie_groups), "expected mixed-tie signature collisions")
minimal_key, minimal_values = mixed_tie_groups[0]
require(minimal_key == (2, 6, 7, 0), "minimal hostile signature")
print(f"Reduced pairs with denominator<=200: {sum(len(values) for values in groups.values())}.")
print(f"Khinchin-length/product + denominator + final-phase groups mixing tie status: {len(mixed_tie_groups)}.")
print(f"Minimal key {minimal_key}: {minimal_values}.")


print("\n=== Current 182-skeleton ordered-word collision ===")
first_digits = continued_fraction_digits(43, 182)
second_digits = continued_fraction_digits(55, 182)
require(first_digits == (4, 4, 3, 3), "43/182 digits")
require(second_digits == (3, 3, 4, 4), "55/182 digits")
require(finite_khinchin_signature(first_digits) == finite_khinchin_signature(second_digits) == (4, 144), "same skeleton mean")
require(midpoint_phase_trace(first_digits)[-1] == midpoint_phase_trace(second_digits)[-1] == 1, "same skeleton phase")
checks += 4
print("43/182=[0;4,4,3,3] and 55/182=[0;3,3,4,4].")
print("They are distinct THM-3710 phases with the same length, product 144, and final phase 1.")


print("\n=== Constant continued-fraction tails in SL2 ===")
word_16 = alternating_constant_word((1, 6))
word_23 = alternating_constant_word((2, 3))
require(word_16 == ((7, 1), (6, 1)), "constant word (1,6)")
require(word_23 == ((7, 2), (3, 1)), "constant word (2,3)")
require(word_16 != word_23, "ordered words differ")
require(determinant_2x2(word_16) == determinant_2x2(word_23) == 1, "constant words lie in SL2")
checks += 4
print(f"E_+(1)E_-(6)={word_16}; E_+(2)E_-(3)={word_23}.")
print("The same digit product gives different constant SL2 words; THM-3736 closes the whole constant-right orbit.")


print("\n=== JC cyclic multiplier holonomy ===")
multiplier_values = tuple(Fraction(value, 1) for value in (-3, -2, -1, 1, 2, 3))
holonomy_checks = 0
for size in range(2, 6):
    for multipliers in product(multiplier_values, repeat=size):
        determinant = cyclic_transport_determinant(multipliers)
        expected = prod(multipliers) - Fraction((-1) ** size, 1)
        require(determinant == expected, "cyclic holonomy determinant")
        holonomy_checks += 1
checks += holonomy_checks

positive_control = (Fraction(2, 1), Fraction(1, 2))
sign_hostile = (Fraction(-2, 1), Fraction(1, 2))
require(cyclic_transport_determinant(positive_control) == 0, "reciprocal even control closes")
require(cyclic_transport_determinant(sign_hostile) != 0, "same absolute mean sign hostile")
require(prod(abs(value) for value in positive_control) == prod(abs(value) for value in sign_hostile) == 1, "same absolute product")
checks += 3

factorial_checks = 0
for size in range(1, 13):
    factorial_multipliers = tuple(Fraction(value, 1) for value in range(2, size + 2))
    determinant = cyclic_transport_determinant(factorial_multipliers)
    require(determinant == prod(factorial_multipliers) - Fraction((-1) ** size, 1), "factorial seam formula")
    require(determinant != 0, "factorial seam never closes")
    factorial_checks += 2
checks += factorial_checks
print(f"Exhaustive alpha_i in {{-3,-2,-1,1,2,3}}, lengths 2..5: {holonomy_checks} determinant gates.")
print("(2,1/2) closes, while (-2,1/2) has the same absolute geometric mean and does not.")
print(f"Factorial seams (2,3,...,n+1), 1<=n<=12: {factorial_checks} exact gates, all nonclosing.")


print("\n=== Verdict ===")
print("Khinchin's constant is a typical asymptotic scalar for irrational digit products.")
print("The exact finite objects needed here are ordered continuants plus target-specific sidecars:")
print("  LRC: common scale, midpoint phase trace, tie blocks, owners, and semantic arrival;")
print("  JC: the exact constant/polynomial SL2 word, multiplier holonomy, curl, and charge sector.")
print(f"CHECKS={checks}")
