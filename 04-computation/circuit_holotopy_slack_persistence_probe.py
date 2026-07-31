#!/usr/bin/env python3
"""Exact controls for metric slack holotopy and weak ray compactification."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def two_grade() -> None:
    grid = [F(k, 2) for k in range(3)]
    values = [(max(1 - x, x), x) for x in grid]
    sigma, point = min(values)
    require((sigma, point) == (F(1, 2), F(1, 2)), values)
    for test in (F(0), F(1, 3), F(1, 2)):
        full = 1 - test <= test
        require(full == (test >= F(1, 2)), (test, full))
    print(
        "two_grade;nerve_before=S0;death=1/2;"
        "dual=(1/2,1/2);nerve_after=filled_edge"
    )

    scaled_point = F(100, 101)
    scaled_death = max(100 * (1 - scaled_point), scaled_point)
    require(scaled_death == F(100, 101), scaled_death)
    print("row_scaling_hostile;same_sign_type=1;scaled_death=100/101")


def three_grade() -> None:
    pair_witnesses = {
        (0, 1): (F(1), F(1)),
        (0, 2): (F(1), F(0)),
        (1, 2): (F(0), F(1)),
    }
    for pair in combinations(range(3), 2):
        x, y = pair_witnesses[pair]
        tests = (x >= 1, y >= 1, x + y <= 1)
        require(all(tests[j] for j in pair), (pair, tests))

    candidates = []
    for i in range(13):
        for j in range(13):
            x, y = F(i, 6), F(j, 6)
            candidates.append((max(1 - x, 1 - y, x + y - 1), x, y))
    sigma, x, y = min(candidates)
    require((sigma, x, y) == (F(1, 3), F(2, 3), F(2, 3)), (sigma, x, y))
    for x, y in ((F(0), F(0)), (F(2), F(0)), (F(3, 5), F(7, 6))):
        dual = ((1 - x) + (1 - y) + (x + y - 1)) / 3
        require(dual == F(1, 3), (x, y, dual))
    print(
        "three_grade;nerve_before=S1;death=1/3;"
        "dual=(1/3,1/3,1/3);nerve_after=filled_triangle"
    )


def nonminimal_hostile() -> None:
    vertices = {0, 1, 2}
    edges = {(0, 2), (1, 2)}
    require((0, 1) not in edges, edges)
    require(len(vertices) - len(edges) == 1, "nerve is a path, not a circle")
    print("nonminimal_hostile;missing_pair=(0,1);nerve=contractible_path;S1=0")


def projective_ray_compactification() -> None:
    period, residue, amplitude, target = 10, 3, F(12), F(1)
    finite_residual = target - F(3, 4)
    finite = []
    for m in range(20):
        z = residue + m * period
        if amplitude / z >= finite_residual:
            finite.append(m)
    require(finite == list(range(5)), finite)
    require(amplitude / (residue + 4 * period) >= finite_residual, "last point")
    require(amplitude / (residue + 5 * period) < finite_residual, "first omission")

    infinite_residual = target - F(1)
    require(infinite_residual == 0, infinite_residual)
    for m in (0, 1, 7, 10**4):
        score = amplitude / (residue + m * period)
        require(score >= infinite_residual and score > infinite_residual, m)
    require(F(0) >= infinite_residual, "weak infinity endpoint")
    require(not (F(0) > infinite_residual), "strict infinity hostile")

    divisor = gcd(period, residue)
    denominator, unit = period // divisor, residue // divisor
    for m in (0, 1, 17, 10**4):
        z = residue + m * period
        d = gcd(period, z)
        require((period // d, (z % period) // d) == (denominator, unit), m)
    print(
        "ray_compactification;finite_bar=0..4;infinite_core=1;"
        f"phase=(d:{denominator},u:{unit});horizon_free=1;"
        "strict_R0_infinity=EXCLUDED"
    )


def main() -> None:
    two_grade()
    three_grade()
    nonminimal_hostile()
    projective_ray_compactification()
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
