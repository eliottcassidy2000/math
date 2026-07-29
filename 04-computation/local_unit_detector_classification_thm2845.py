#!/usr/bin/env python3
"""Finite exact controls for THM-2845.

The theorem itself is proved abstractly.  This companion exhausts the small
hostile algebras where a false strengthening is most tempting:

* split products over F_2, F_3, and F_5;
* dual numbers, checking that a detector kills the radical;
* the noncommutative upper-triangular algebras UT_2(F_3) and UT_3(F_2);
* M_2(F_2), whose every nonzero functional misses a unit; and
* F_4 as a two-dimensional F_2-algebra, where every nonzero functional
  has a nonzero field element in its kernel.

All gates use explicit exceptions and therefore survive ``python -O``.
"""

from __future__ import annotations

from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(coefficients: tuple[int, ...], vector: tuple[int, ...], q: int) -> int:
    return sum(c * x for c, x in zip(coefficients, vector)) % q


def detector_coefficients(
    q: int,
    dimension: int,
    units: list[tuple[int, ...]],
) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []
    for coefficients in product(range(q), repeat=dimension):
        if not any(coefficients):
            continue
        if all(dot(coefficients, unit, q) != 0 for unit in units):
            result.append(coefficients)
    return result


def split_controls() -> list[str]:
    lines: list[str] = []
    for q in (2, 3, 5):
        for rank in range(1, 6):
            units = list(product(range(1, q), repeat=rank))
            actual = detector_coefficients(q, rank, units)
            if q == 2:
                expected = [
                    coefficients
                    for coefficients in product(range(q), repeat=rank)
                    if sum(coefficients) % 2 == 1
                ]
            else:
                expected = [
                    coefficients
                    for coefficients in product(range(q), repeat=rank)
                    if sum(c != 0 for c in coefficients) == 1
                ]
            require(
                actual == expected,
                f"split classification failed at q={q}, rank={rank}",
            )
        last_count = len(
            detector_coefficients(
                q,
                5,
                list(product(range(1, q), repeat=5)),
            )
        )
        lines.append(f"split F_{q} ranks 1..5: exact; rank5 detectors={last_count}")
    return lines


def dual_number_controls() -> list[str]:
    lines: list[str] = []
    for q in (2, 3, 5):
        # Coordinates are (constant, epsilon); a unit has nonzero constant.
        units = [
            (constant, epsilon)
            for constant in range(1, q)
            for epsilon in range(q)
        ]
        actual = detector_coefficients(q, 2, units)
        expected = [(constant_weight, 0) for constant_weight in range(1, q)]
        require(actual == expected, f"dual-number radical gate failed at q={q}")
        lines.append(
            f"dual F_{q}[e]/(e^2): detectors={len(actual)}; "
            "every detector kills (e)"
        )
    return lines


def upper_triangular_controls() -> list[str]:
    lines: list[str] = []

    # UT_2(F_3), coordinates (a_11,a_12,a_22).
    q = 3
    units_2 = [
        (a, b, d)
        for a in range(1, q)
        for b in range(q)
        for d in range(1, q)
    ]
    actual_2 = detector_coefficients(q, 3, units_2)
    expected_2 = [
        (c, 0, 0) for c in range(1, q)
    ] + [
        (0, 0, c) for c in range(1, q)
    ]
    require(
        set(actual_2) == set(expected_2),
        "UT_2(F_3) classification failed",
    )
    lines.append(
        f"UT_2(F_3): units={len(units_2)} detectors={len(actual_2)}; "
        "off-diagonal radical killed"
    )

    # UT_3(F_2), coordinates (a_11,a_12,a_13,a_22,a_23,a_33).
    units_3 = [
        (1, a12, a13, 1, a23, 1)
        for a12, a13, a23 in product(range(2), repeat=3)
    ]
    actual_3 = detector_coefficients(2, 6, units_3)
    expected_3 = [
        (c1, 0, 0, c2, 0, c3)
        for c1, c2, c3 in product(range(2), repeat=3)
        if (c1 + c2 + c3) % 2 == 1
    ]
    require(actual_3 == expected_3, "UT_3(F_2) parity classification failed")
    require(
        (1, 0, 0, 1, 0, 1) in actual_3,
        "UT_3(F_2) full odd trace was not detected",
    )
    lines.append(
        f"UT_3(F_2): units={len(units_3)} detectors={len(actual_3)}; "
        "full three-diagonal trace works"
    )
    return lines


def determinant_2_by_2(
    matrix: tuple[int, int, int, int],
    q: int,
) -> int:
    a, b, c, d = matrix
    return (a * d - b * c) % q


def nonscalar_simple_controls() -> list[str]:
    lines: list[str] = []

    q = 2
    matrices = list(product(range(q), repeat=4))
    matrix_units = [
        matrix for matrix in matrices if determinant_2_by_2(matrix, q) != 0
    ]
    matrix_detectors = detector_coefficients(q, 4, matrix_units)
    require(matrix_detectors == [], "M_2(F_2) acquired a detector")
    require(len(matrix_units) == 6, "M_2(F_2) unit count changed")
    lines.append(
        f"M_2(F_2): units={len(matrix_units)} nonzero forms=15 detectors=0"
    )

    # F_4=F_2[a]/(a^2+a+1).  Every nonzero coordinate pair is a unit.
    field_units = [
        vector for vector in product(range(2), repeat=2) if any(vector)
    ]
    field_detectors = detector_coefficients(2, 2, field_units)
    require(field_detectors == [], "F_4/F_2 acquired a detector")
    lines.append("F_4 over F_2: units=3 nonzero forms=3 detectors=0")
    return lines


def main() -> None:
    lines = [
        "THM-2845 finite exact controls",
        *split_controls(),
        *dual_number_controls(),
        *upper_triangular_controls(),
        *nonscalar_simple_controls(),
        "all exact controls passed",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    main()
