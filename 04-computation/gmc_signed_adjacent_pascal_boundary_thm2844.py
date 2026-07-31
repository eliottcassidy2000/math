#!/usr/bin/env python3
"""Exact companion for THM-2844.

The computation checks the sharp signed boundary of the Pascal-ratio
orientation certificate.  It uses exact integers, rationals, and symbolic
polynomials only.  Every truth-bearing gate is an explicit ``require`` and
therefore remains active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial, gcd

import sympy as sp


s, n, m, t, u = sp.symbols("s n m t u")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def adjacent_difference(index: int) -> sp.Expr:
    return sp.expand(
        s ** (index + 1) / factorial(index + 1)
        - s**index / factorial(index)
    )


def factorial_functional(poly: sp.Expr) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * factorial(degree[0])
            for degree, coefficient in sp.Poly(sp.expand(poly), s).terms()
        )
    )


def pascal_symbolic(level: sp.Expr, index: int) -> sp.Expr:
    return sp.expand(sp.rf(level + 1, index) / factorial(index))


def triple_symbolic(level: sp.Expr, left: int, right: int) -> sp.Expr:
    total = left + right
    return sp.expand(
        sp.rf(level + 2, total) / (factorial(left) * factorial(right))
        + sp.rf(level + 1, total + 1)
        / (factorial(left + 1) * factorial(right))
        + sp.rf(level + 1, total + 1)
        / (factorial(left) * factorial(right + 1))
        - sp.rf(level + 1, total) / (factorial(left) * factorial(right))
    )


def symbolic_ab(
    level: sp.Expr,
    coefficients: dict[int, sp.Expr],
) -> tuple[sp.Expr, sp.Expr]:
    a_value = sum(
        coefficient * pascal_symbolic(level, index)
        for index, coefficient in coefficients.items()
    )
    b_value = sum(
        left_coefficient
        * right_coefficient
        * triple_symbolic(level, left, right)
        for left, left_coefficient in coefficients.items()
        for right, right_coefficient in coefficients.items()
    )
    return sp.expand(a_value), sp.expand(b_value)


def pascal_numeric(level: int, index: int) -> int:
    return comb(level + index, index)


def triple_numeric(level: int, left: int, right: int) -> int:
    total = level + left + right
    return (
        factorial(total + 1)
        // (factorial(level + 1) * factorial(left) * factorial(right))
        + factorial(total + 1)
        // (factorial(level) * factorial(left + 1) * factorial(right))
        + factorial(total + 1)
        // (factorial(level) * factorial(left) * factorial(right + 1))
        - factorial(total)
        // (factorial(level) * factorial(left) * factorial(right))
    )


def numeric_ab(
    level: int,
    coefficients: dict[int, int | Fraction],
) -> tuple[Fraction, Fraction]:
    a_value = sum(
        Fraction(coefficient) * pascal_numeric(level, index)
        for index, coefficient in coefficients.items()
    )
    b_value = sum(
        Fraction(left_coefficient)
        * Fraction(right_coefficient)
        * triple_numeric(level, left, right)
        for left, left_coefficient in coefficients.items()
        for right, right_coefficient in coefficients.items()
    )
    return a_value, b_value


def orientation_direct(upper: sp.Expr, lower: sp.Expr) -> sp.Expr:
    return sp.expand(
        2
        * factorial_functional(lower**3)
        * factorial_functional(upper * lower)
        - 3
        * factorial_functional(upper * lower**2)
        * factorial_functional(lower**2)
    )


def first_divisibility_invariant(
    upper: sp.Expr,
    lower: sp.Expr,
) -> sp.Expr:
    g11 = factorial_functional(upper**2)
    g12 = factorial_functional(upper * lower)
    g22 = factorial_functional(lower**2)
    t111 = factorial_functional(upper**3)
    t112 = factorial_functional(upper**2 * lower)
    t222 = factorial_functional(lower**3)
    return sp.expand(
        3 * t112 * g11 * g22
        - t222 * g11**2
        - 2 * t111 * g12 * g22
    )


def positive_coefficients(poly: sp.Expr, *variables: sp.Symbol) -> bool:
    return all(
        coefficient > 0
        for coefficient in sp.Poly(sp.expand(poly), *variables).coeffs()
    )


def boundary_polynomial(value: int | Fraction | sp.Expr) -> sp.Expr:
    return sp.expand(
        5 * value**3 + 30 * value**2 + 57 * value + 24
    )


# ---------------------------------------------------------------------------
# 1. Symbolic signed-adjacent boundary
# ---------------------------------------------------------------------------

adjacent_a, adjacent_b = symbolic_ab(n, {1: t, 2: sp.Integer(1)})
expected_a = sp.expand((n + 1) * (n + 2 * t + 2) / 2)
require(
    sp.expand(adjacent_a - expected_a) == 0,
    "signed-adjacent A_n identity",
)

base_b = (
    n**5
    + 10 * n**4
    + 47 * n**3
    + 122 * n**2
    + 168 * n
    + 96
)
linear_b = (
    5 * n**4
    + 38 * n**3
    + 121 * n**2
    + 184 * n
    + 108
)
quadratic_b = n**3 + 6 * n**2 + 13 * n + 10
expected_six_b = sp.expand(
    base_b + u * linear_b + 6 * u**2 * quadratic_b
)
require(
    sp.expand(6 * adjacent_b.subs(t, u - 1) - expected_six_b) == 0,
    "signed-adjacent B_n positive decomposition",
)
require(
    positive_coefficients(base_b, n)
    and positive_coefficients(linear_b, n)
    and positive_coefficients(quadratic_b, n),
    "signed-adjacent B_n coefficient positivity",
)

adjacent_delta = sp.expand(
    adjacent_a * adjacent_b.subs(n, n + 1)
    - adjacent_a.subs(n, n + 1) * adjacent_b
)
boundary_p = boundary_polynomial(t)
require(
    sp.expand(adjacent_delta.subs(n, 0) - 2 * boundary_p) == 0,
    "first signed-adjacent determinant",
)

shifted_delta_target = (
    m**6
    + 17 * m**5
    + 119 * m**4
    + 427 * m**3
    + 788 * m**2
    + 632 * m
    + 104
    + u
    * (
        6 * m**5
        + 86 * m**4
        + 494 * m**3
        + 1410 * m**2
        + 1980 * m
        + 1080
    )
    + u**2
    * (
        12 * m**4
        + 140 * m**3
        + 608 * m**2
        + 1152 * m
        + 792
    )
    + u**3 * (8 * m**3 + 72 * m**2 + 208 * m + 184)
)
require(
    sp.expand(
        4 * adjacent_delta.subs({n: m + 1, t: u - 1})
        - shifted_delta_target
    )
    == 0,
    "all later signed-adjacent determinants",
)
require(
    positive_coefficients(shifted_delta_target, m, u),
    "later signed-adjacent determinant coefficient positivity",
)

require(
    sp.expand(sp.diff(boundary_p, t) - (15 * (t + 2) ** 2 - 3)) == 0,
    "boundary derivative identity",
)
require(
    positive_coefficients(sp.diff(boundary_p, t).subs(t, u - 1), u),
    "boundary derivative positivity on t>-1",
)
require(
    boundary_polynomial(sp.Rational(-2, 3)) == sp.Rational(-58, 27)
    and boundary_polynomial(sp.Rational(-1, 2)) == sp.Rational(19, 8),
    "unique boundary root bracket",
)

direct_symbolic_checks = 0
symbolic_lower = sp.expand(t * adjacent_difference(1) + adjacent_difference(2))
for level in range(15):
    require(
        sp.expand(
            factorial_functional(adjacent_difference(level) * symbolic_lower)
            - adjacent_a.subs(n, level)
        )
        == 0,
        f"direct A_n check at n={level}",
    )
    require(
        sp.expand(
            factorial_functional(
                adjacent_difference(level) * symbolic_lower**2
            )
            - adjacent_b.subs(n, level)
        )
        == 0,
        f"direct B_n check at n={level}",
    )
    direct_symbolic_checks += 2

direct_orientation = orientation_direct(adjacent_difference(0), symbolic_lower)
require(
    sp.expand(direct_orientation - 12 * boundary_p) == 0,
    "orientation/boundary identity",
)


# ---------------------------------------------------------------------------
# 2. Singleton support-cut law and Cartesian maximality
# ---------------------------------------------------------------------------

singleton_positive = 0
singleton_zero = 0
singleton_negative = 0

for upper_index in range(13):
    for lower_index in range(1, 13):
        direct_value = orientation_direct(
            adjacent_difference(upper_index),
            adjacent_difference(lower_index),
        )
        a_upper, b_upper = numeric_ab(
            upper_index,
            {lower_index: 1},
        )
        a_predecessor, b_predecessor = numeric_ab(
            lower_index - 1,
            {lower_index: 1},
        )
        ratio_value = (
            6
            * a_upper
            * a_predecessor
            * (
                b_predecessor / a_predecessor
                - b_upper / a_upper
            )
        )
        require(
            direct_value == ratio_value,
            f"singleton ratio identity at {(upper_index, lower_index)}",
        )
        expected_sign = lower_index - upper_index - 1
        require(
            (direct_value > 0) == (expected_sign > 0)
            and (direct_value == 0) == (expected_sign == 0)
            and (direct_value < 0) == (expected_sign < 0),
            f"singleton support-cut sign at {(upper_index, lower_index)}",
        )
        if direct_value > 0:
            singleton_positive += 1
        elif direct_value == 0:
            singleton_zero += 1
        else:
            singleton_negative += 1

require(
    (singleton_positive, singleton_zero, singleton_negative) == (66, 12, 78),
    "singleton sign census",
)
require(
    orientation_direct(adjacent_difference(2), adjacent_difference(1))
    == -228,
    "reversed singleton hostile",
)
require(
    orientation_direct(adjacent_difference(1), adjacent_difference(2)) == 0,
    "ordered adjacent equality",
)

cartesian_support_pairs = 0
cartesian_legal_pairs = 0
for upper_mask in range(1, 1 << 6):
    upper_support = {
        index for index in range(6) if (upper_mask >> index) & 1
    }
    for lower_mask in range(1, 1 << 6):
        lower_support = {
            index + 1 for index in range(6) if (lower_mask >> index) & 1
        }
        all_singletons_nonnegative = all(
            lower_index > upper_index
            for upper_index in upper_support
            for lower_index in lower_support
        )
        support_cut = max(upper_support) < min(lower_support)
        require(
            all_singletons_nonnegative == support_cut,
            "Cartesian support-cut maximality census",
        )
        cartesian_support_pairs += 1
        cartesian_legal_pairs += int(support_cut)


# ---------------------------------------------------------------------------
# 3. Exact rational bank and primitive adjacent integer scan
# ---------------------------------------------------------------------------

rational_bank = sorted(
    {
        Fraction(numerator, denominator)
        for denominator in range(1, 31)
        for numerator in range(-denominator + 1, 3 * denominator + 1)
    }
)
rational_below = []
rational_above = []
rational_bank_checks = 0

for rational_t in rational_bank:
    coefficients = {1: rational_t, 2: Fraction(1)}
    values = [numeric_ab(level, coefficients) for level in range(41)]
    require(
        all(a_value > 0 and b_value > 0 for a_value, b_value in values),
        f"rational-bank positive transforms at t={rational_t}",
    )
    determinants = [
        values[level][0] * values[level + 1][1]
        - values[level + 1][0] * values[level][1]
        for level in range(40)
    ]
    p_value = Fraction(
        boundary_polynomial(
            sp.Rational(rational_t.numerator, rational_t.denominator)
        )
    )
    require(p_value != 0, f"unexpected rational boundary root {rational_t}")
    require(
        determinants[0] == 2 * p_value,
        f"rational-bank first determinant at t={rational_t}",
    )
    require(
        all(value > 0 for value in determinants[1:]),
        f"rational-bank later determinants at t={rational_t}",
    )
    require(
        all(value >= 0 for value in determinants) == (p_value > 0),
        f"rational-bank sharp classification at t={rational_t}",
    )
    if p_value < 0:
        rational_below.append(rational_t)
    else:
        rational_above.append(rational_t)
    rational_bank_checks += 1

require(
    len(rational_bank) == 1112
    and len(rational_below) == 116
    and len(rational_above) == 996,
    "rational-bank census",
)
require(
    rational_below[-1] == Fraction(-7, 12)
    and rational_above[0] == Fraction(-11, 19),
    "rational-bank nearest boundary points",
)

primitive_adjacent_pairs = []
primitive_adjacent_hostiles = []

for coefficient_height in range(2, 51):
    for lower_coefficient in range(1, coefficient_height):
        upper_coefficient = -(coefficient_height - lower_coefficient)
        if (
            upper_coefficient + lower_coefficient <= 0
            or gcd(-upper_coefficient, lower_coefficient) != 1
        ):
            continue
        pair = (upper_coefficient, lower_coefficient)
        primitive_adjacent_pairs.append(pair)
        coefficients = {
            1: upper_coefficient,
            2: lower_coefficient,
        }
        values = [numeric_ab(level, coefficients) for level in range(31)]
        require(
            all(a_value > 0 and b_value > 0 for a_value, b_value in values),
            f"primitive adjacent positivity at pair={pair}",
        )
        determinants = [
            values[level][0] * values[level + 1][1]
            - values[level + 1][0] * values[level][1]
            for level in range(30)
        ]
        projective_t = Fraction(upper_coefficient, lower_coefficient)
        p_value = Fraction(
            boundary_polynomial(
                sp.Rational(projective_t.numerator, projective_t.denominator)
            )
        )
        require(
            determinants[0] == 2 * lower_coefficient**3 * p_value,
            f"primitive adjacent first determinant at pair={pair}",
        )
        require(
            all(value > 0 for value in determinants[1:]),
            f"primitive adjacent later determinants at pair={pair}",
        )
        if p_value < 0:
            primitive_adjacent_hostiles.append(pair)

require(
    len(primitive_adjacent_pairs) == 386
    and len(primitive_adjacent_hostiles) == 102,
    "primitive adjacent scan census",
)
minimum_hostile_height = min(
    -upper_coefficient + lower_coefficient
    for upper_coefficient, lower_coefficient in primitive_adjacent_hostiles
)
minimum_hostiles = [
    pair
    for pair in primitive_adjacent_hostiles
    if -pair[0] + pair[1] == minimum_hostile_height
]
require(
    minimum_hostile_height == 5 and minimum_hostiles == [(-2, 3)],
    "primitive adjacent minimal hostile",
)


# ---------------------------------------------------------------------------
# 4. Two sharp primitive hostiles
# ---------------------------------------------------------------------------

star_a, star_b = symbolic_ab(n, {1: -2, 2: 3})
star_delta = sp.expand(
    star_a * star_b.subs(n, n + 1)
    - star_a.subs(n, n + 1) * star_b
)
star_a_target = sp.expand((n + 1) * (3 * n + 2) / 2)
star_b_target = sp.expand(
    (n + 2)
    * (
        3 * n**4
        + 29 * n**3
        + 123 * n**2
        + 253 * n
        + 208
    )
    / 2
)
star_four_delta_target = (
    27 * n**6
    + 351 * n**5
    + 1863 * n**4
    + 4901 * n**3
    + 6066 * n**2
    + 2344 * n
    - 464
)
star_shifted_target = (
    27 * m**6
    + 513 * m**5
    + 4023 * m**4
    + 16403 * m**3
    + 35862 * m**2
    + 38548 * m
    + 15088
)
require(
    sp.expand(star_a - star_a_target) == 0
    and sp.expand(star_b - star_b_target) == 0,
    "V-star transform formulas",
)
require(
    sp.expand(4 * star_delta - star_four_delta_target) == 0
    and sp.expand(
        star_four_delta_target.subs(n, m + 1) - star_shifted_target
    )
    == 0
    and positive_coefficients(star_shifted_target, m),
    "V-star determinant formulas",
)
star_lower = sp.expand(
    3 * adjacent_difference(2) - 2 * adjacent_difference(1)
)
require(
    star_delta.subs(n, 0) == -116
    and orientation_direct(adjacent_difference(0), star_lower) == -2088,
    "V-star exact hostile",
)

dagger_a, dagger_b = symbolic_ab(n, {1: -1, 3: 2})
dagger_delta = sp.expand(
    dagger_a * dagger_b.subs(n, n + 1)
    - dagger_a.subs(n, n + 1) * dagger_b
)
dagger_a_target = sp.expand((n + 1) * (n**2 + 5 * n + 3) / 3)
dagger_b_target = sp.expand(
    (n + 2)
    * (
        n**6
        + 26 * n**5
        + 273 * n**4
        + 1519 * n**3
        + 4796 * n**2
        + 8151 * n
        + 5814
    )
    / 18
)
dagger_shift_polynomial = (
    m**9
    + 42 * m**8
    + 765 * m**7
    + 7920 * m**6
    + 51198 * m**5
    + 213129 * m**4
    + 566273 * m**3
    + 912909 * m**2
    + 791487 * m
    + 269703
)
require(
    sp.expand(dagger_a - dagger_a_target) == 0
    and sp.expand(dagger_b - dagger_b_target) == 0,
    "V-dagger transform formulas",
)
require(
    sp.expand(
        dagger_delta.subs(n, m + 1)
        - sp.Rational(2, 27) * dagger_shift_polynomial
    )
    == 0
    and positive_coefficients(dagger_shift_polynomial, m),
    "V-dagger later determinant formula",
)
dagger_lower = sp.expand(
    2 * adjacent_difference(3) - adjacent_difference(1)
)
require(
    dagger_delta.subs(n, 0) == -446
    and dagger_b.subs(n, 0) / dagger_a.subs(n, 0) == 646
    and dagger_b.subs(n, 1) / dagger_a.subs(n, 1)
    == sp.Rational(1715, 3)
    and orientation_direct(adjacent_difference(0), dagger_lower) == 24792,
    "V-dagger exact hostile and long-secant survivor",
)

height_two_controls = 0
for left_index in range(9):
    for right_index in range(left_index + 1, 9):
        height_two = sp.expand(
            adjacent_difference(right_index)
            - adjacent_difference(left_index)
        )
        require(
            factorial_functional(adjacent_difference(0) * height_two) == 0,
            "height-two signed boundary has A_0=0",
        )
        height_two_controls += 1


# ---------------------------------------------------------------------------
# 5. Failure of orientation is not failure of moment-three detection
# ---------------------------------------------------------------------------

i1_direct = first_divisibility_invariant(
    adjacent_difference(0),
    symbolic_lower,
)
i1_target = -2 * (7 * u**3 + 43 * u**2 + 97 * u + 81)
require(
    sp.expand(i1_direct.subs(t, u - 1) - i1_target) == 0,
    "first divisibility invariant identity",
)
require(
    positive_coefficients(-i1_target / 2, u),
    "first divisibility invariant strict negativity",
)
require(
    i1_direct.subs(t, sp.Rational(-2, 3)) == sp.Rational(-6392, 27),
    "V-star moment-three survivor",
)


print("THM-2844 SIGNED ADJACENT PASCAL BOUNDARY - exact referee")
print("symbolic direct A/B checks:", direct_symbolic_checks)
print(
    "singleton signs positive/zero/negative:",
    singleton_positive,
    "/",
    singleton_zero,
    "/",
    singleton_negative,
)
print(
    "Cartesian support pairs / support-cut legal:",
    cartesian_support_pairs,
    "/",
    cartesian_legal_pairs,
)
print(
    "rational t bank / below / above:",
    rational_bank_checks,
    "/",
    len(rational_below),
    "/",
    len(rational_above),
)
print(
    "nearest rational bank bracket:",
    rational_below[-1],
    "< alpha <",
    rational_above[0],
)
print(
    "primitive adjacent pairs / hostiles / first:",
    len(primitive_adjacent_pairs),
    "/",
    len(primitive_adjacent_hostiles),
    "/",
    minimum_hostiles[0],
)
print("V-star Delta0 / orientation:", star_delta.subs(n, 0), "/", -2088)
print(
    "V-dagger Delta0 / R0 / R1 / long-secant orientation:",
    dagger_delta.subs(n, 0),
    "/",
    646,
    "/",
    sp.Rational(1715, 3),
    "/",
    24792,
)
print("height-two A0-zero controls:", height_two_controls)
print("V-star I1 survivor:", sp.Rational(-6392, 27))
print("all arithmetic exact: integers/rationals/symbolic polynomials only")
print("THM-2844 exact companion: PASS")
