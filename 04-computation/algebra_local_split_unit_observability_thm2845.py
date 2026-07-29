#!/usr/bin/env python3
"""Exact finite-universe companion for THM-2845.

The script checks the theorem's smallest local, nonsplit, matrix, split,
finite-field, and characteristic-zero boundaries:

1. in truncated local algebras F_q[e]/(e^n), the unit-safe and exact linear
   functionals are precisely the nonzero residue multiples;
2. every radical-sensitive functional kills an explicitly enumerated unit;
3. F_4 over F_2 and M_2(F_q) admit no unit-safe scalar in the tested universes;
4. on F_q^D, direct enumeration agrees with the support classification,
   including the odd-support F_2 exception;
5. the zero-sum unit census agrees with
      N_D(q)=((q-1)^D+(q-1)(-1)^D)/q;
6. the full group sum in Q[C_p] has nonzero augmentation but singular regular
   representation, while the delta mask is the unit positive control.

All truth-bearing gates use ``require``.  Arithmetic is finite-field integer
arithmetic or exact ``Fraction`` arithmetic; there are no floats and no Python
``assert`` statements, so ``python -O`` performs the same audit.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot_mod(coefficients: tuple[int, ...], value: tuple[int, ...], q: int) -> int:
    return sum(c * x for c, x in zip(coefficients, value)) % q


def all_vectors(q: int, dimension: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product(range(q), repeat=dimension))


def classify_functionals(
    elements: tuple[tuple[int, ...], ...],
    units: frozenset[tuple[int, ...]],
    q: int,
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    """Return all unit-safe and exact coefficient vectors."""
    dimension = len(elements[0])
    safe: list[tuple[int, ...]] = []
    exact: list[tuple[int, ...]] = []
    for coefficients in all_vectors(q, dimension):
        if all(dot_mod(coefficients, unit, q) != 0 for unit in units):
            safe.append(coefficients)
        if all(
            (dot_mod(coefficients, value, q) != 0) == (value in units)
            for value in elements
        ):
            exact.append(coefficients)
    return tuple(safe), tuple(exact)


# ---------------------------------------------------------------------------
# 1. Truncated local algebras: residue is the unique unit-safe direction.
# ---------------------------------------------------------------------------

local_configurations = ((2, 2), (2, 3), (3, 2), (3, 3), (5, 2))
local_elements = 0
local_functionals = 0
local_safe = 0
local_exact = 0
radical_sensitive = 0
radical_killed_unit = 0

for q, nilpotence_degree in local_configurations:
    elements = all_vectors(q, nilpotence_degree)
    units = frozenset(value for value in elements if value[0] != 0)
    safe, exact = classify_functionals(elements, units, q)
    predicted = tuple(
        (residue_scalar,) + (0,) * (nilpotence_degree - 1)
        for residue_scalar in range(1, q)
    )
    require(
        safe == predicted,
        f"local unit-safe classification failed for q={q}, n={nilpotence_degree}",
    )
    require(
        exact == predicted,
        f"local exact classification failed for q={q}, n={nilpotence_degree}",
    )

    for coefficients in all_vectors(q, nilpotence_degree):
        if any(coefficients[index] != 0 for index in range(1, nilpotence_degree)):
            radical_sensitive += 1
            witness = next(
                (
                    unit
                    for unit in units
                    if dot_mod(coefficients, unit, q) == 0
                ),
                None,
            )
            require(
                witness is not None,
                "radical-sensitive functional unexpectedly survived all units "
                f"for q={q}, n={nilpotence_degree}, ell={coefficients}",
            )
            radical_killed_unit += 1

    local_elements += len(elements)
    local_functionals += q**nilpotence_degree
    local_safe += len(safe)
    local_exact += len(exact)

require(local_elements == 73, "local element census drift")
require(local_functionals == 73, "local functional census drift")
require(local_safe == 10, "local unit-safe census drift")
require(local_exact == 10, "local exact census drift")
require(radical_sensitive == 58, "radical-sensitive census drift")
require(
    radical_killed_unit == radical_sensitive,
    "not every radical-sensitive functional acquired a killed unit",
)


# ---------------------------------------------------------------------------
# 2. The residue-degree hostile F_4/F_2.
# ---------------------------------------------------------------------------

def f4_add(left: int, right: int) -> int:
    return left ^ right


def f4_multiply(left: int, right: int) -> int:
    """Multiply a+b*alpha modulo alpha^2+alpha+1 over F_2."""
    a0, a1 = left & 1, (left >> 1) & 1
    b0, b1 = right & 1, (right >> 1) & 1
    constant = (a0 * b0 + a1 * b1) % 2
    alpha = (a0 * b1 + a1 * b0 + a1 * b1) % 2
    return constant | (alpha << 1)


for left, middle, right in product(range(4), repeat=3):
    require(
        f4_multiply(f4_multiply(left, middle), right)
        == f4_multiply(left, f4_multiply(middle, right)),
        "F4 associativity failed",
    )
    require(
        f4_multiply(left, f4_add(middle, right))
        == f4_add(f4_multiply(left, middle), f4_multiply(left, right)),
        "F4 distributivity failed",
    )

for value in range(1, 4):
    require(
        any(f4_multiply(value, inverse) == 1 for inverse in range(1, 4)),
        f"F4 nonzero element {value} has no inverse",
    )

f4_safe = 0
f4_exact = 0
f4_kernel_unit_hostiles = 0
for coefficients in all_vectors(2, 2):
    values = tuple(
        (coefficients[0] * (element & 1)
         + coefficients[1] * ((element >> 1) & 1))
        % 2
        for element in range(4)
    )
    safe = all(values[element] != 0 for element in range(1, 4))
    exact = all((values[element] != 0) == (element != 0) for element in range(4))
    f4_safe += int(safe)
    f4_exact += int(exact)
    if coefficients != (0, 0):
        require(
            any(values[element] == 0 for element in range(1, 4)),
            f"nonzero F2-linear F4 functional {coefficients} has trivial kernel",
        )
        f4_kernel_unit_hostiles += 1

require(f4_safe == 0, "F4/F2 unexpectedly has a unit-safe functional")
require(f4_exact == 0, "F4/F2 unexpectedly has an exact functional")
require(f4_kernel_unit_hostiles == 3, "F4/F2 kernel-hostile census drift")


# ---------------------------------------------------------------------------
# 3. Matrix hostiles: M_2(F_q) has no unit-safe scalar.
# ---------------------------------------------------------------------------

matrix_fields = (2, 3, 5)
matrix_elements = 0
matrix_units = 0
matrix_functionals = 0
matrix_safe = 0
matrix_exact = 0

for q in matrix_fields:
    elements = all_vectors(q, 4)
    units = frozenset(
        matrix
        for matrix in elements
        if (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % q != 0
    )
    safe, exact = classify_functionals(elements, units, q)
    require(safe == (), f"M2(F_{q}) unexpectedly has a unit-safe scalar")
    require(exact == (), f"M2(F_{q}) unexpectedly has an exact scalar")
    require(
        len(units) == (q**2 - 1) * (q**2 - q),
        f"GL2(F_{q}) census drift",
    )
    matrix_elements += len(elements)
    matrix_units += len(units)
    matrix_functionals += q**4
    matrix_safe += len(safe)
    matrix_exact += len(exact)

require(matrix_elements == 722, "matrix element census drift")
require(matrix_units == 534, "matrix unit census drift")
require(matrix_functionals == 722, "matrix functional census drift")
require(matrix_safe == 0 and matrix_exact == 0, "matrix classification drift")


# ---------------------------------------------------------------------------
# 4. Split products: support-one, odd-support F_2, exact only in D=1.
# ---------------------------------------------------------------------------

split_dimensions = {
    2: range(1, 7),
    3: range(1, 6),
    5: range(1, 5),
}
split_configurations = 0
split_weight_vectors = 0
split_safe = 0
split_exact = 0

for q, dimensions in split_dimensions.items():
    for dimension in dimensions:
        elements = all_vectors(q, dimension)
        units = frozenset(
            value for value in elements if all(coordinate != 0 for coordinate in value)
        )
        safe, exact = classify_functionals(elements, units, q)
        predicted_safe = tuple(
            weights
            for weights in all_vectors(q, dimension)
            if (
                sum(weight != 0 for weight in weights) == 1
                if q != 2
                else sum(weight != 0 for weight in weights) % 2 == 1
            )
        )
        predicted_exact = (
            tuple((weight,) for weight in range(1, q))
            if dimension == 1
            else ()
        )
        require(
            safe == predicted_safe,
            f"split unit-safe support law failed for q={q}, D={dimension}",
        )
        require(
            exact == predicted_exact,
            f"split exact law failed for q={q}, D={dimension}",
        )

        split_configurations += 1
        split_weight_vectors += q**dimension
        split_safe += len(safe)
        split_exact += len(exact)

require(split_configurations == 15, "split configuration census drift")
require(split_weight_vectors == 1269, "split weight-vector census drift")
require(split_safe == 133, "split unit-safe census drift")
require(split_exact == 7, "split exact census drift")

# Name the smallest parity boundary explicitly.
f2_square_trace = lambda value: sum(value) % 2
require(
    f2_square_trace((1, 1)) == 0,
    "F2^2 full trace failed to kill its sole unit",
)
require(
    f2_square_trace((1, 1, 1)) == 1,
    "F2^3 full trace failed to preserve its sole unit",
)
require(
    f2_square_trace((1, 0, 0)) == 1,
    "F2^3 full trace failed to expose its nonexact boundary",
)


# ---------------------------------------------------------------------------
# 5. Direct zero-sum census for the finite-field formula N_D(q).
# ---------------------------------------------------------------------------

zero_sum_fields = (2, 3, 5, 7)
zero_sum_degrees = range(1, 7)
zero_sum_configurations = 0
zero_sum_units = 0
zero_sum_solutions = 0

for q in zero_sum_fields:
    for dimension in zero_sum_degrees:
        nonzero_tuples = tuple(product(range(1, q), repeat=dimension))
        direct = sum(sum(value) % q == 0 for value in nonzero_tuples)
        numerator = (q - 1) ** dimension + (q - 1) * (-1) ** dimension
        require(numerator % q == 0, "N_D(q) formula lost integrality")
        formula = numerator // q
        require(
            direct == formula,
            f"N_D(q) formula failed for q={q}, D={dimension}",
        )
        predicted_zero = dimension == 1 or (q == 2 and dimension % 2 == 1)
        require(
            (direct == 0) == predicted_zero,
            f"zero-sum exception classification failed for q={q}, D={dimension}",
        )
        zero_sum_configurations += 1
        zero_sum_units += len(nonzero_tuples)
        zero_sum_solutions += direct

require(zero_sum_configurations == 24, "zero-sum configuration census drift")
require(zero_sum_units == 61578, "zero-sum unit census drift")


# ---------------------------------------------------------------------------
# 6. Q[C_p]: nonzero augmentation alone is not a unit detector.
# ---------------------------------------------------------------------------

def cyclic_convolution(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = len(left)
    require(len(right) == size, "cyclic convolution size mismatch")
    return tuple(
        sum(left[index] * right[(output - index) % size] for index in range(size))
        for output in range(size)
    )


def circulant_columns(mask: tuple[int, ...]) -> list[list[Fraction]]:
    size = len(mask)
    return [
        [
            Fraction(mask[(row - column) % size])
            for column in range(size)
        ]
        for row in range(size)
    ]


def rational_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            multiplier = work[row][column]
            if multiplier != 0:
                work[row] = [
                    entry - multiplier * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


cyclic_primes = (2, 3, 5, 7)
cyclic_hostiles = 0
cyclic_controls = 0

for prime in cyclic_primes:
    group_sum = (1,) * prime
    difference = (-1, 1) + (0,) * (prime - 2)
    require(
        sum(group_sum) == prime,
        f"Q[C_{prime}] group-sum augmentation drift",
    )
    require(
        cyclic_convolution(difference, group_sum) == (0,) * prime,
        f"Q[C_{prime}] group sum lost its zero-divisor witness",
    )
    require(
        rational_rank(circulant_columns(group_sum)) == 1,
        f"Q[C_{prime}] group-sum regular rank drift",
    )
    require(
        sum(group_sum) % prime == 0,
        f"F_{prime}[C_{prime}] group-sum augmentation should vanish",
    )
    cyclic_hostiles += 1

    delta = (1,) + (0,) * (prime - 1)
    require(sum(delta) == 1, f"Q[C_{prime}] delta augmentation drift")
    require(
        rational_rank(circulant_columns(delta)) == prime,
        f"Q[C_{prime}] delta regular representation is not invertible",
    )
    cyclic_controls += 1

require(cyclic_hostiles == 4, "cyclic hostile census drift")
require(cyclic_controls == 4, "cyclic positive-control census drift")


print("THM-2845 LOCAL-VS-SPLIT UNIT OBSERVABILITY - exact audit")
print(
    "truncated locals (configs/elements/functionals):",
    f"{len(local_configurations)}/{local_elements}/{local_functionals}",
)
print(
    "local residue functionals (unit-safe/exact):",
    f"{local_safe}/{local_exact}",
)
print(
    "radical-sensitive killed-unit hostiles:",
    f"{radical_killed_unit}/{radical_sensitive}",
)
print(
    "F4/F2 (linear functionals/kernel-unit hostiles/safe/exact):",
    f"4/{f4_kernel_unit_hostiles}/{f4_safe}/{f4_exact}",
)
print(
    "M2 finite fields (configs/elements/units/functionals/safe/exact):",
    f"{len(matrix_fields)}/{matrix_elements}/{matrix_units}/"
    f"{matrix_functionals}/{matrix_safe}/{matrix_exact}",
)
print(
    "split products (configs/weight-vectors/unit-safe/exact):",
    f"{split_configurations}/{split_weight_vectors}/{split_safe}/{split_exact}",
)
print("F2 parity boundary (D2-unit/D3-unit/D3-nonunit): 0/1/1")
print(
    "zero-sum census (configs/units/solutions):",
    f"{zero_sum_configurations}/{zero_sum_units}/{zero_sum_solutions}",
)
print("N_D(q) direct/formula checks:", zero_sum_configurations)
print(
    "Q[C_p] augmentation hostiles/delta controls:",
    f"{cyclic_hostiles}/{cyclic_controls}",
)
print("ALL EXACT CHECKS PASSED")
