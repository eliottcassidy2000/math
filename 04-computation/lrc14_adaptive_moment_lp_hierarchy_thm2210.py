#!/usr/bin/env python3
"""Exact finite referee for THM-2210's Boolean moment LP hierarchy.

For every nonzero packet p in {0,1,2}^{j+1}, 0 <= j <= 5, this script
independently enumerates primal basic feasible solutions and dual binomial
minorants at every degree.  It also checks the sharp fixed-prefix coefficient
through j=20 and freezes the two minimal structural examples used in the
theorem.  All arithmetic is rational.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


def solve_square(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> list[Fraction]:
    n = len(rhs)
    aug = [row[:] + [value] for row, value in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((row for row in range(col, n) if aug[row][col]), None)
        require(pivot is not None, "singular binomial evaluation matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [value / scale for value in aug[col]]
        for row in range(n):
            if row == col:
                continue
            scale = aug[row][col]
            if scale:
                aug[row] = [
                    left - scale * right
                    for left, right in zip(aug[row], aug[col])
                ]
    return [aug[row][-1] for row in range(n)]


def moments(packet: tuple[int | Fraction, ...]) -> tuple[Fraction, ...]:
    j = len(packet) - 1
    return tuple(
        sum(Fraction(choose(k, r)) * packet[k] for k in range(j + 1))
        for r in range(j + 1)
    )


def primal_value(
    moment: tuple[Fraction, ...], j: int, d: int
) -> Fraction:
    best = None
    for support in combinations(range(j + 1), d + 1):
        matrix = [
            [Fraction(choose(k, r)) for k in support]
            for r in range(d + 1)
        ]
        solution = solve_square(matrix, list(moment[: d + 1]))
        if any(value < 0 for value in solution):
            continue
        value = solution[support.index(0)] if 0 in support else Fraction(0)
        best = value if best is None else min(best, value)
    require(best is not None, "no primal basic feasible solution")
    return best


def dual_value(
    moment: tuple[Fraction, ...], j: int, d: int
) -> Fraction:
    best = None
    for active_nodes in combinations(range(j + 1), d + 1):
        matrix = [
            [Fraction(choose(k, r)) for r in range(d + 1)]
            for k in active_nodes
        ]
        rhs = [Fraction(1 if k == 0 else 0) for k in active_nodes]
        coefficient = solve_square(matrix, rhs)

        def q(k: int) -> Fraction:
            return sum(
                coefficient[r] * choose(k, r) for r in range(d + 1)
            )

        if any(q(k) > (1 if k == 0 else 0) for k in range(j + 1)):
            continue
        value = sum(coefficient[r] * moment[r] for r in range(d + 1))
        best = value if best is None else max(best, value)
    require(best is not None, "no dual vertex")
    return best


def sharp_prefix_coefficient(j: int, d: int) -> Fraction:
    require(0 <= d <= j, "invalid prefix degree")
    if d == 0:
        return Fraction(1 if j == 0 else 0)
    return Fraction(d, j) if d % 2 == 0 else Fraction(-1)


def prefix_value(
    moment: tuple[Fraction, ...], j: int, d: int
) -> Fraction:
    coefficient = sharp_prefix_coefficient(j, d)
    return sum(
        Fraction((-1) ** r) * moment[r] for r in range(d)
    ) + coefficient * moment[d]


def prefix_polynomial(j: int, d: int, c: Fraction, k: int) -> Fraction:
    if d == 0:
        return c
    return sum(Fraction((-1) ** r * choose(k, r)) for r in range(d)) + (
        c * choose(k, d)
    )


def check_sharp_prefixes() -> int:
    cases = 0
    epsilon = Fraction(1, 100)
    for j in range(21):
        for d in range(j + 1):
            cases += 1
            c = sharp_prefix_coefficient(j, d)
            require(
                all(
                    prefix_polynomial(j, d, c, k)
                    <= (1 if k == 0 else 0)
                    for k in range(j + 1)
                ),
                "sharp prefix is infeasible",
            )
            larger = c + epsilon
            require(
                any(
                    prefix_polynomial(j, d, larger, k)
                    > (1 if k == 0 else 0)
                    for k in range(j + 1)
                ),
                "prefix coefficient was not sharp",
            )
    require(cases == 231, "prefix census drift")
    return cases


def check_packet_census() -> int:
    cases = 0
    for j in range(6):
        for raw_packet in product(range(3), repeat=j + 1):
            if not any(raw_packet):
                continue
            cases += 1
            packet = tuple(Fraction(value) for value in raw_packet)
            moment = moments(packet)
            values = []
            for d in range(j + 1):
                primal = primal_value(moment, j, d)
                dual = dual_value(moment, j, d)
                require(primal == dual, "primal/dual mismatch")
                require(Fraction(0) <= primal <= packet[0], "LP range failure")
                require(
                    prefix_value(moment, j, d) <= primal,
                    "fixed prefix exceeds adaptive LP",
                )
                values.append(primal)
            require(
                all(values[d] <= values[d + 1] for d in range(j)),
                "adaptive hierarchy is not nested",
            )
            require(values[-1] == packet[0], "full-depth inversion failure")
            if j >= 1:
                require(
                    values[1] == max(Fraction(0), moment[0] - moment[1]),
                    "degree-one closed form failure",
                )
            if j >= 2:
                pair_bound = (
                    moment[0]
                    - moment[1]
                    + Fraction(2, j) * moment[2]
                )
                require(
                    values[2] == max(Fraction(0), pair_bound),
                    "degree-two closed form failure",
                )
    require(cases == 1086, "packet census drift")
    return cases


def check_structural_examples() -> None:
    nonnested = (
        Fraction(1, 2),
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(1, 2),
    )
    moment = moments(nonnested)
    require(moment == (1, 2, 3, 2, Fraction(1, 2)), "moment drift")
    require(prefix_value(moment, 4, 2) == Fraction(1, 2), "B2 drift")
    require(prefix_value(moment, 4, 3) == 0, "B3 drift")
    require(primal_value(moment, 4, 2) == Fraction(1, 2), "L2 drift")
    require(primal_value(moment, 4, 3) == Fraction(1, 2), "L3 drift")

    cubic_gain = (
        Fraction(1, 2),
        Fraction(0),
        Fraction(0),
        Fraction(1, 2),
        Fraction(0),
    )
    moment = moments(cubic_gain)
    require(
        moment == (1, Fraction(3, 2), Fraction(3, 2), Fraction(1, 2), 0),
        "cubic example moment drift",
    )
    require(primal_value(moment, 4, 2) == Fraction(1, 4), "cubic L2 drift")
    require(primal_value(moment, 4, 3) == Fraction(1, 2), "cubic L3 drift")


def main() -> None:
    packet_cases = check_packet_census()
    prefix_cases = check_sharp_prefixes()
    check_structural_examples()
    print(f"nonzero_packet_cases_j0_to_j5={packet_cases}")
    print(f"sharp_prefix_cases_j_to_20={prefix_cases}")
    print("exact_primal_dual_values=PASS")
    print("adaptive_lp_nesting_and_full_inversion=PASS")
    print("sharp_fixed_prefix_coefficients=PASS")
    print("degree1_degree2_closed_forms=PASS")
    print("fixed_prefix_nonnesting_witness=PASS j=4 support={0,4}")
    print("strict_cubic_gain_witness=PASS j=4 support={0,3}")
    print("status=THM-2210_EXACT_FINITE_REFEREE_ONLY")


if __name__ == "__main__":
    main()
