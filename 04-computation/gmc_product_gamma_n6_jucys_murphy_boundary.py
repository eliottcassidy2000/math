#!/usr/bin/env python3
"""Exact N=6 Farahat--Higman/Jucys--Murphy boundary for THM-3110.

The script derives the degree-six power-sum class coefficients from the
actual 24/25 response banks, reconstructs the corresponding central
S_6 element by Murnaghan--Nakayama character arithmetic, and expresses
its spectrum in the first three power sums of Young-diagram contents.

All calculations are symbolic over Q[a,b].
"""

import importlib.util
from collections import Counter
from functools import lru_cache
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
SPEC = importlib.util.spec_from_file_location("thm3110_schur", UPSTREAM)
THM3110 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(THM3110)
BANKS = THM3110.BANKS

a, b, summation_index = sp.symbols("a b summation_index")


def linear(form):
    return form[0] * a + form[1] * b


@lru_cache(maxsize=None)
def interval_power(power, form):
    endpoint = linear(form)
    return sp.expand(
        sp.summation(
            summation_index ** power,
            (summation_index, 1, endpoint - 1),
        )
    )


@lru_cache(maxsize=None)
def integer_partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    answer = []
    for first in range(largest, 0, -1):
        for rest in integer_partitions(total - first, first):
            answer.append((first,) + rest)
    return tuple(answer)


def class_quotients(total, invariant, chamber):
    middle = (0, 1) if chamber == "wide" else (2, 0)
    common = ((1, 0), (1, 0), (1, 0), middle)
    if invariant == 1:
        common += (middle,)
    divisor = (
        a ** (4 if invariant == 0 else 3)
        * b**2
        * (b - a) ** (3 if invariant == 0 else 4)
    )

    atom_powers = []
    for coefficient, row in BANKS[invariant]:
        powers = {
            degree: sp.expand(
                sum(interval_power(degree, form) for form in row)
                - sum(interval_power(degree, form) for form in common)
            )
            for degree in range(1, total + 1)
        }
        atom_powers.append((coefficient, powers))

    answer = {}
    for cycle_type in integer_partitions(total):
        numerator = sp.expand(
            sum(
                coefficient
                * sp.prod(powers[part] for part in cycle_type)
                for coefficient, powers in atom_powers
            )
        )
        quotient_numerator, quotient_denominator = sp.fraction(
            sp.cancel(numerator / divisor)
        )
        require(
            quotient_denominator == 1,
            "collision divisor failed in class quotient",
        )
        answer[cycle_type] = sp.factor(quotient_numerator)
    return answer


def z_value(cycle_type):
    answer = 1
    for part, multiplicity in Counter(cycle_type).items():
        answer *= part**multiplicity * factorial(multiplicity)
    return answer


def subpartitions(shape, target):
    def recurse(row, previous, remaining, prefix):
        if row == len(shape):
            if remaining == 0:
                yield tuple(part for part in prefix if part)
            return
        for part in range(min(shape[row], previous, remaining), -1, -1):
            yield from recurse(
                row + 1, part, remaining - part, prefix + (part,)
            )

    yield from recurse(0, sum(shape), target, ())


def border_strip_height(shape, subshape, size):
    padded = subshape + (0,) * (len(shape) - len(subshape))
    cells = {
        (row, column)
        for row in range(len(shape))
        for column in range(padded[row] + 1, shape[row] + 1)
    }
    if len(cells) != size or not cells:
        return None

    seen = {next(iter(cells))}
    stack = list(seen)
    while stack:
        row, column = stack.pop()
        for neighbour in (
            (row - 1, column),
            (row + 1, column),
            (row, column - 1),
            (row, column + 1),
        ):
            if neighbour in cells and neighbour not in seen:
                seen.add(neighbour)
                stack.append(neighbour)
    if seen != cells:
        return None
    if any(
        (row, column + 1) in cells
        and (row + 1, column) in cells
        and (row + 1, column + 1) in cells
        for row, column in cells
    ):
        return None
    return len({row for row, _ in cells}) - 1


@lru_cache(maxsize=None)
def character(shape, cycle_type):
    if not cycle_type:
        return int(not shape)
    strip_size = cycle_type[0]
    answer = 0
    for subshape in subpartitions(shape, sum(shape) - strip_size):
        height = border_strip_height(shape, subshape, strip_size)
        if height is not None:
            answer += (-1) ** height * character(subshape, cycle_type[1:])
    return answer


def contents(shape):
    return tuple(
        column - row
        for row, length in enumerate(shape, 1)
        for column in range(1, length + 1)
    )


PARTITIONS_SIX = integer_partitions(6)
CONTENT_BASIS = []
for shape in PARTITIONS_SIX:
    content = contents(shape)
    p1 = sum(content)
    p2 = sum(entry**2 for entry in content)
    p3 = sum(entry**3 for entry in content)
    CONTENT_BASIS.append((1, p1, p2, p1**2, p3, p1 * p2, p1**3))
CONTENT_MATRIX = sp.Matrix(CONTENT_BASIS)
require(CONTENT_MATRIX.rank() == 7, "content basis rank drift")


EXPECTED_COMMON_SIX = {
    (6,): 0,
    (5, 1): 0,
    (4, 2): 0,
    (3, 3): 0,
    (4, 1, 1): 36,
    (3, 2, 1): -21,
    (3, 1, 1, 1): 9 * (18 * a + 18 * b - 5),
    (2, 2, 2): 48,
}


def central_spectrum(class_values):
    eigenvalues = []
    for shape in PARTITIONS_SIX:
        dimension = character(shape, (1,) * 6)
        require(dimension > 0, "character dimension drift")
        eigenvalues.append(
            sp.factor(
                sum(
                    class_values[cycle_type]
                    * sp.Rational(character(shape, cycle_type), z_value(cycle_type))
                    for cycle_type in PARTITIONS_SIX
                )
                / dimension
            )
        )
    solution = tuple(next(iter(sp.linsolve((CONTENT_MATRIX, sp.Matrix(eigenvalues))))))
    reconstructed = CONTENT_MATRIX * sp.Matrix(solution)
    require(
        all(
            sp.expand(reconstructed[index] - eigenvalues[index]) == 0
            for index in range(len(eigenvalues))
        ),
        "content reconstruction drift",
    )
    return eigenvalues, solution


records = {}
for invariant in (0, 1):
    for chamber in ("tight", "wide"):
        q5 = class_quotients(5, invariant, chamber)
        linear_parameter = 3 * a + (5 if invariant == 0 else 4) * b
        expected_five = {
                (5,): 0,
                (4, 1): 0,
                (3, 2): 0,
                (3, 1, 1): 0,
                (2, 2, 1): 0,
                (2, 1, 1, 1): 30,
                (1, 1, 1, 1, 1): 60 * linear_parameter,
            }
        require(
            all(
                sp.expand(q5[cycle_type] - expected) == 0
                for cycle_type, expected in expected_five.items()
            ),
            "degree-five class table drift",
        )

        q6 = class_quotients(6, invariant, chamber)
        for cycle_type, expected in EXPECTED_COMMON_SIX.items():
            require(
                sp.expand(q6[cycle_type] - expected) == 0,
                "universal degree-six class drift",
            )
        require(
            all(
                q6[cycle_type] == 0
                for cycle_type in ((6,), (5, 1), (4, 2), (3, 3))
            ),
            "high-defect support drift",
        )

        eigenvalues, coefficients = central_spectrum(q6)
        require(
            coefficients[-3:]
            == (sp.Rational(7, 18), sp.Rational(-31, 240), sp.Rational(1, 90)),
            "universal cubic content drift",
        )

        first = PARTITIONS_SIX.index((4, 1, 1))
        second = PARTITIONS_SIX.index((3, 3))
        expected_difference = (
            (4 * a + 18 * b + 31) / 10
            if invariant == 0
            else (4 * a + 11 * b + 31) / 10
        )
        require(
            sp.expand(
                eigenvalues[first] - eigenvalues[second] - expected_difference
            )
            == 0,
            "equal-p1 collision difference drift",
        )
        records[invariant, chamber] = (q5, q6, eigenvalues, coefficients)


def top_cubic(shape):
    content = contents(shape)
    p1 = sum(content)
    p2 = sum(entry**2 for entry in content)
    p3 = sum(entry**3 for entry in content)
    return sp.Rational(8 * p1**3 - 93 * p1 * p2 + 280 * p3, 720)


require(top_cubic((6,)) == sp.Rational(295, 16), "top cubic positive pole")
require(top_cubic((3, 2, 1)) == 0, "top cubic self-conjugate wall")
require(top_cubic((1, 1, 1, 1, 1, 1)) == sp.Rational(-295, 16),
        "top cubic negative pole")

# The N=5 central factor on a shape mu is L/2+p1(C(mu))/4.  In the
# universal filtered polynomial ring it cannot divide N=6: after setting
# p1=-2L, the coefficient of p3 remains 7/18.
p1_symbol, p2_symbol, p3_symbol, L_symbol = sp.symbols("p1 p2 p3 L")
universal_top = (
    sp.Rational(1, 90) * p1_symbol**3
    - sp.Rational(31, 240) * p1_symbol * p2_symbol
    + sp.Rational(7, 18) * p3_symbol
)
require(
    sp.Poly(universal_top.subs(p1_symbol, -2 * L_symbol), p3_symbol).coeff_monomial(p3_symbol)
    == sp.Rational(7, 18),
    "degree-five divisibility hostile drift",
)

print("normalization=q_rho=Phi(p_rho)/base;Z6=(1/720)*sum_sigma q_type(sigma)*sigma")
print("class_support_N6=111111,21111,2211,3111,222,321,411")
print("universal_classes=q411:36;q321:-21;q222:48;q311:9*(18a+18b-5)")
print("content_basis=1,p1,p2,p1^2,p3,p1*p2,p1^3;rank=7")
print("universal_cubic=(8*p1^3-93*p1*p2+280*p3)/720")
print("same_p1_pair=(411),(33);delta_I1=(4a+18b+31)/10;delta_I2=(4a+11b+31)/10")
print("top_cubic_values=295/16,0,-295/16:on:(6),(321),(111111)")
print("no_go=not_p1_only;not_N5_linear_multiple;universal_cubic_not_PSD")
print("needed_sidecar=p2,p3_or_rooted_NBC_flag")
print("all_exact_checks=PASS")
