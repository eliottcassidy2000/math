#!/usr/bin/env python3
"""Exact companion for THM-2909.

The proved theorem is the weighted-chord upgrade of THM-2873 from
consecutive response rows to arbitrary triples of rows at least one.
This companion checks the exact chord identity, the two sharp hostiles,
and a separately scoped finite scan beyond the proved support gaps.

All truth-bearing gates use ``require`` rather than ``assert``, so normal
and optimized Python execute the same checks.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def tensor(indices):
    """Return L(prod_i d_(indices[i])) by exact inclusion-exclusion."""

    block_sizes = [index + 1 for index in indices]
    elementary = [1]
    for block_size in block_sizes:
        elementary.append(0)
        for degree in range(len(elementary) - 1, 0, -1):
            elementary[degree] += (
                block_size * elementary[degree - 1]
            )

    total = sum(block_sizes)
    numerator = sum(
        (-1) ** degree
        * elementary[degree]
        * factorial(total - degree)
        for degree in range(len(elementary))
    )
    denominator = 1
    for block_size in block_sizes:
        denominator *= factorial(block_size)
    return Fraction(numerator, denominator)


def determinant(matrix):
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    )


def two_ray_moment(row, power, pivot, gap, low_weight, high_weight):
    return sum(
        Fraction(comb(power, high_count))
        * low_weight ** (power - high_count)
        * high_weight ** high_count
        * tensor(
            tuple(
                sorted(
                    (row,)
                    + (pivot,) * (power - high_count)
                    + (pivot + gap,) * high_count
                )
            )
        )
        for high_count in range(power + 1)
    )


def signed_moment(row, power, terms):
    """Expand L(d_row W^power) for a short signed list (coefficient,index)."""

    total = Fraction(0)

    def expand(slot, coefficient, indices):
        nonlocal total
        if slot == power:
            total += coefficient * tensor(
                tuple(sorted((row,) + tuple(indices)))
            )
            return
        for weight, index in terms:
            expand(
                slot + 1,
                coefficient * weight,
                indices + [index],
            )

    expand(0, Fraction(1), [])
    return total


def row_matrix(rows, pivot, gap, low_weight, high_weight):
    return [
        [
            two_ray_moment(
                row,
                power,
                pivot,
                gap,
                low_weight,
                high_weight,
            )
            for power in (1, 2, 3)
        ]
        for row in rows
    ]


def weighted_chord_side(matrix_by_row, rows):
    """Return the normalized determinant and its exact double-sum form."""

    r0, r1, r2 = rows
    first = {
        row: matrix_by_row[row][0]
        for row in range(r0, r2 + 1)
    }
    response_r = {
        row: matrix_by_row[row][1] / first[row]
        for row in range(r0, r2 + 1)
    }
    response_s = {
        row: matrix_by_row[row][2] / first[row]
        for row in range(r0, r2 + 1)
    }
    delta_r = {
        row: response_r[row + 1] - response_r[row]
        for row in range(r0, r2)
    }
    require(
        all(value > 0 for value in delta_r.values()),
        "response R must increase in every checked cell",
    )
    slope = {
        row: (
            response_s[row + 1] - response_s[row]
        ) / delta_r[row]
        for row in range(r0, r2)
    }

    normalized = determinant(
        [matrix_by_row[row] for row in rows]
    ) / (first[r0] * first[r1] * first[r2])
    double_sum = sum(
        delta_r[left]
        * delta_r[right]
        * (slope[right] - slope[left])
        for left in range(r0, r1)
        for right in range(r1, r2)
    )
    require(normalized == double_sum, "weighted-chord identity")
    return normalized


def exact_scan():
    """Finite evidence, with gaps 3..8 kept outside the proved theorem."""

    theorem_cells = 0
    evidence_cells = 0
    minimum = None
    minimum_cell = None
    weights = [
        (Fraction(low), Fraction(high))
        for low in range(4)
        for high in range(4)
        if low or high
    ]

    for pivot in range(7):
        for gap in range(1, 9):
            for low_weight, high_weight in weights:
                matrix_by_row = {
                    row: row_matrix(
                        (row,),
                        pivot,
                        gap,
                        low_weight,
                        high_weight,
                    )[0]
                    for row in range(1, 10)
                }
                for rows in combinations(range(1, 10), 3):
                    normalized = weighted_chord_side(
                        matrix_by_row,
                        rows,
                    )
                    require(
                        normalized > 0,
                        "nonpositive scanned response determinant",
                    )
                    raw = determinant(
                        [matrix_by_row[row] for row in rows]
                    )
                    require(raw > 0, "nonpositive raw scanned determinant")
                    if minimum is None or raw < minimum:
                        minimum = raw
                        minimum_cell = (
                            pivot,
                            gap,
                            int(low_weight),
                            int(high_weight),
                            rows,
                        )
                    if gap in (1, 2):
                        theorem_cells += 1
                    else:
                        evidence_cells += 1

    return theorem_cells, evidence_cells, minimum, minimum_cell


def main():
    print("THM-2909 exact nonconsecutive response-row TP3 arc")

    d4_rows = (0, 1, 2)
    d4_matrix = [
        [
            signed_moment(
                row,
                power,
                [(Fraction(1), 4)],
            )
            for power in (1, 2, 3)
        ]
        for row in d4_rows
    ]
    require(
        d4_matrix
        == [
            [Fraction(1), Fraction(812), Fraction(3854466)],
            [Fraction(5), Fraction(5040), Fraction(33243210)],
            [Fraction(15), Fraction(22260), Fraction(201170970)],
        ],
        "d4 row-zero matrix",
    )
    d4_determinant = determinant(d4_matrix)
    d4_curvature = d4_determinant / Fraction(1 * 5 * 15)
    require(d4_determinant == -339570000, "d4 row-zero determinant")
    require(d4_curvature == -4527600, "d4 row-zero curvature")

    signed_terms = [
        (Fraction(1), 0),
        (Fraction(-1, 5), 1),
    ]
    signed_rows = (1, 2, 3)
    signed_matrix = [
        [
            signed_moment(row, power, signed_terms)
            for power in (1, 2, 3)
        ]
        for row in signed_rows
    ]
    require(
        signed_matrix
        == [
            [Fraction(3, 5), Fraction(6, 5), Fraction(351, 125)],
            [Fraction(2, 5), Fraction(28, 25), Fraction(66, 25)],
            [Fraction(1, 5), Fraction(4, 5), Fraction(33, 25)],
        ],
        "signed hostile matrix",
    )
    signed_determinant = determinant(signed_matrix)
    require(
        signed_determinant == Fraction(-1728, 15625),
        "signed hostile determinant",
    )

    theorem_cells, evidence_cells, minimum, minimum_cell = exact_scan()
    require(theorem_cells == 17640, "proved-gap scan cell count")
    require(evidence_cells == 52920, "unproved-gap scan cell count")
    require(minimum == 12, "finite scan minimum determinant")
    require(
        minimum_cell == (0, 1, 1, 0, (1, 2, 3)),
        "finite scan minimum cell",
    )

    print("row_zero_d4_matrix=1,812,3854466;5,5040,33243210;15,22260,201170970")
    print(f"row_zero_d4_determinant={d4_determinant}")
    print(f"row_zero_d4_normalized_curvature={d4_curvature}")
    print("signed_hostile_A=3/5,2/5,1/5")
    print(f"signed_hostile_determinant={signed_determinant}")
    print(f"proved_gap_weighted_chord_cells={theorem_cells}")
    print(f"unproved_gap_finite_evidence_cells={evidence_cells}")
    print("finite_scan_minimum=12;cell=p0,h1,x1,y0,rows1-2-3")
    print("scope=h in {1,2} proved; h=3..8 is finite-exact evidence only")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
