#!/usr/bin/env python3
"""Exact all-width certificate for the first-gap transfer-one signs.

For the cardinal response space

    span{1-x^M, x-x^M, x^2-x^M},

the six off-diagonal values

    (E^2+E)(lambda_k/lambda_j)|_{x=j},  j,k in {2,3,4},

have a common positive squared generalized-Vandermonde denominator. This
companion checks the explicit rational-exponential numerators against direct
rational response inversion on M=3..50, then verifies elementary dominant-base
tail certificates covering every M beyond the finite prefix.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from flint import fmpq_mat


# A polynomial tuple means a*M^2+b*M+c. The sign field multiplies the raw
# numerator of S_jk so that the displayed exponential sum must be positive.
CASES = {
    (2, 3): {
        "sign": 1,
        "den_factor": 1,
        "terms": {
            16: (0, 0, 6),
            12: (0, 0, -12),
            8: (-3, -9, 12),
            6: (9, -9, 36),
            4: (9, -27, -30),
            3: (0, 0, -24),
            2: (3, 45, -12),
            1: (0, 0, 24),
        },
        "tail": 25,
        "margin": Fraction(
            14081865922424343364772654133,
            39614081257132168796771975168,
        ),
    },
    (2, 4): {
        "sign": -1,
        "den_factor": 1,
        "terms": {
            12: (0, 0, 6),
            9: (0, 0, -12),
            8: (-1, -1, -8),
            6: (3, -9, 42),
            4: (3, -3, -22),
            3: (0, 0, -24),
            2: (1, 13, 14),
            1: (0, 0, 4),
        },
        "tail": 25,
        "margin": Fraction(
            1423058678816606772702161,
            13249474533898474022240256,
        ),
    },
    (3, 2): {
        "sign": 1,
        "den_factor": 2,
        "terms": {
            16: (0, 0, 48),
            12: (6, -120, 171),
            9: (18, 36, -216),
            8: (0, 0, -288),
            6: (18, 72, 297),
            4: (0, 0, 21),
            3: (-6, 12, -36),
            2: (0, 0, -9),
            1: (0, 0, 12),
        },
        "tail": 30,
        "margin": Fraction(
            1334504361243610627652012697506433,
            5192296858534827628530496329220096,
        ),
    },
    (3, 4): {
        "sign": 1,
        "den_factor": 2,
        "terms": {
            12: (2, -16, 45),
            9: (6, -60, 108),
            8: (0, 0, -48),
            6: (6, 96, -369),
            4: (0, 0, 291),
            3: (-2, -20, 108),
            2: (0, 0, -159),
            1: (0, 0, 24),
        },
        "tail": 18,
        "margin": Fraction(
            86807696439009361,
            369768517790072832,
        ),
    },
    (4, 2): {
        "sign": -1,
        "den_factor": 1,
        "terms": {
            16: (3, -35, 64),
            12: (9, -15, 44),
            9: (0, 0, -108),
            8: (-9, 87, -320),
            6: (0, 0, 324),
            4: (3, -37, 148),
            3: (0, 0, -152),
            2: (0, 0, -4),
            1: (0, 0, 4),
        },
        "tail": 15,
        "margin": Fraction(
            333379720232789933,
            32425917317067571200,
        ),
    },
    (4, 3): {
        "sign": 1,
        "den_factor": 1,
        "terms": {
            16: (3, -23, 36),
            12: (9, -51, 92),
            8: (-9, 123, -296),
            6: (0, 0, -108),
            4: (3, -49, 456),
            3: (0, 0, 16),
            2: (0, 0, -244),
            1: (0, 0, 48),
        },
        "tail": 16,
        "margin": Fraction(
            49905085875636783,
            288230376151711744,
        ),
    },
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def poly_value(poly: tuple[int, int, int], width: int) -> int:
    a, b, c = poly
    return a * width * width + b * width + c


def exponential_sum(case: dict, width: int) -> int:
    return sum(
        poly_value(poly, width) * base**width
        for base, poly in case["terms"].items()
    )


def as_fraction(value) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def direct_cardinal_value(width: int, j: int, k: int) -> Fraction:
    nodes = (2, 3, 4)
    response = fmpq_mat(
        [
            [1 - node**width, node - node**width, node * node - node**width]
            for node in nodes
        ]
    )
    inverse = response.inv()
    j_power = j**width
    first_row = [
        -width * j_power // j,
        1 - width * j_power // j,
        2 * j - width * j_power // j,
    ]
    second_row = [
        -width * (width - 1) * j_power // (j * j),
        -width * (width - 1) * j_power // (j * j),
        2 - width * (width - 1) * j_power // (j * j),
    ]
    j_index = nodes.index(j)
    k_index = nodes.index(k)
    lambda_k_1 = sum(
        (as_fraction(inverse[row, k_index]) * first_row[row] for row in range(3)),
        start=Fraction(0),
    )
    lambda_k_2 = sum(
        (as_fraction(inverse[row, k_index]) * second_row[row] for row in range(3)),
        start=Fraction(0),
    )
    lambda_j_1 = sum(
        (as_fraction(inverse[row, j_index]) * first_row[row] for row in range(3)),
        start=Fraction(0),
    )
    return (
        j * j * lambda_k_2
        + 2 * j * lambda_k_1
        - 2 * j * j * lambda_k_1 * lambda_j_1
    )


def denominator_core(width: int) -> int:
    return 3 * 2**width - 3 * 3**width + 4**width - 1


def coefficient_bound(poly: tuple[int, int, int]) -> int:
    return sum(abs(value) for value in poly)


def certify_tail(case: dict) -> Fraction:
    terms = case["terms"]
    dominant = max(terms)
    leading_poly = terms[dominant]
    threshold = case["tail"]
    lower = {base: poly for base, poly in terms.items() if base != dominant}
    if leading_poly[:2] == (0, 0):
        leading_floor = Fraction(leading_poly[2])
        tail_bound = sum(
            (
                Fraction(coefficient_bound(poly) * threshold * threshold)
                * Fraction(base, dominant) ** threshold
                for base, poly in lower.items()
            ),
            start=Fraction(0),
        )
        # Every M^2(base/dominant)^M upper bound decreases from threshold on.
        for base in lower:
            require(
                (threshold + 1) ** 2 * base < threshold**2 * dominant,
                ("tail monotonicity", threshold, base, dominant),
            )
    else:
        _a, b, c = leading_poly
        require(b < 0 and c > 0, ("leading polynomial", leading_poly))
        leading_floor = Fraction(
            poly_value(leading_poly, threshold), threshold * threshold
        )
        # p(M)/M^2 increases: -b M(M+1)-c(2M+1)>0, then grows.
        require(
            -b * threshold * (threshold + 1) > c * (2 * threshold + 1),
            ("leading monotonicity", threshold),
        )
        require(-b * (threshold + 1) > c, ("leading growth", threshold))
        tail_bound = sum(
            (
                Fraction(coefficient_bound(poly))
                * Fraction(base, dominant) ** threshold
                for base, poly in lower.items()
            ),
            start=Fraction(0),
        )
    margin = leading_floor - tail_bound
    require(margin == case["margin"] and margin > 0, ("tail margin", margin))
    return margin


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for (j, k), case in CASES.items():
        for width in range(3, 51):
            direct = direct_cardinal_value(width, j, k)
            predicted = Fraction(
                case["sign"] * exponential_sum(case, width),
                case["den_factor"] * denominator_core(width) ** 2,
            )
            require(direct == predicted, ("formula", width, j, k))
        require(
            all(
                exponential_sum(case, width) > 0
                for width in range(3, case["tail"])
            ),
            ("finite prefix", j, k),
        )
        margin = certify_tail(case)
        records.append(
            (
                j,
                k,
                case["sign"],
                case["tail"],
                margin.numerator,
                margin.denominator,
            )
        )
    lines = [
        "FIRST-GAP TRANSFER-ONE CHECKERBOARD ALL-WIDTH CERTIFICATE",
        "response_space=span{1-x^M,x-x^M,x^2-x^M}",
        "direct_formula_match=M3..50;cases=6;equalities=288",
    ]
    for j, k, sign, threshold, numerator, denominator in records:
        lines.append(
            f"pair={j}{k};signed_numerator_sign={sign};tail={threshold};"
            f"margin={numerator}/{denominator}"
        )
    digest = sha256("\n".join(lines).encode()).hexdigest()
    lines.extend(
        [
            f"record_digest={digest}",
            "singular_boundary=M1,M2;strict_cardinal_sign=ALL_M_GE_3",
            "first_gap_checkerboard_sign=PROVED_FOR_ALL_M_GE_6",
            "all_exact_checks=PASS",
        ]
    )
    transcript = "\n".join(lines) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
