#!/usr/bin/env python3
"""Exact referee for THM-1134's five-residue and step-two certificates."""

from fractions import Fraction as F
from itertools import combinations


def gaps(values: tuple[int, ...]) -> tuple[int, ...]:
    values = tuple(sorted(values))
    return tuple(values[i + 1] - values[i] for i in range(len(values) - 1)) + (
        13 - values[-1] + values[0],
    )


def orbit(seed: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {
        tuple(sorted((unit * x + shift) % 13 for x in seed))
        for unit in range(1, 13)
        for shift in range(13)
    }


def distance(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def main() -> None:
    remaining = set(combinations(range(13), 5))
    rows = []
    while remaining:
        representative = min(remaining)
        affine_orbit = orbit(representative)
        assert affine_orbit <= set(combinations(range(13), 5))
        assert affine_orbit <= remaining
        remaining -= affine_orbit
        candidates = []
        for unit in range(1, 13):
            word = gaps(tuple((unit * x) % 13 for x in representative))
            candidates.append((max(word), -unit, word))
        maximum, negative_unit, word = max(candidates)
        rows.append(
            (representative, len(affine_orbit), -negative_unit, word, maximum)
        )

    expected = [
        ((0, 1, 2, 3, 4), 78, 1, (1, 1, 1, 1, 9), 9),
        ((0, 1, 2, 3, 5), 156, 1, (1, 1, 1, 2, 8), 8),
        ((0, 1, 2, 3, 6), 156, 1, (1, 1, 1, 3, 7), 7),
        ((0, 1, 2, 3, 7), 156, 2, (1, 1, 2, 2, 7), 7),
        ((0, 1, 2, 3, 8), 78, 2, (2, 1, 1, 2, 7), 7),
        ((0, 1, 2, 4, 5), 156, 1, (1, 1, 2, 1, 8), 8),
        ((0, 1, 2, 4, 7), 156, 1, (1, 1, 2, 3, 6), 6),
        ((0, 1, 2, 4, 10), 156, 1, (1, 1, 2, 6, 3), 6),
        ((0, 1, 2, 5, 6), 156, 1, (1, 1, 3, 1, 7), 7),
        ((0, 1, 2, 6, 9), 39, 2, (2, 2, 1, 7, 1), 7),
    ]
    assert rows == expected
    assert sum(row[1] for row in rows) == 1287

    all_subsets = list(combinations(range(13), 5))
    maxima = []
    for subset in all_subsets:
        chart_maximum = max(
            max(gaps(tuple((unit * x) % 13 for x in subset)))
            for unit in range(1, 13)
        )
        maxima.append(chart_maximum)
    assert min(maxima) == 6
    assert sum(value == 6 for value in maxima) == 312

    minimum_by_size = {}
    for size in range(1, 6):
        minimum_by_size[size] = min(
            max(
                max(gaps(tuple((unit * x) % 13 for x in subset)))
                for unit in range(1, 13)
            )
            for subset in combinations(range(13), size)
        )
    assert minimum_by_size == {1: 13, 2: 12, 3: 9, 4: 7, 5: 6}

    core = (1, 2, 4, 7, 9, 11, 12)
    left, right = F(15, 98), F(9, 56)
    assert right - left == F(3, 392)
    breakpoints = {left, right}
    for p in core:
        for tooth in range(p + 1):
            for sign in (-1, 1):
                point = F(14 * tooth + sign, 14 * p)
                if left <= point <= right:
                    breakpoints.add(point)
    ordered = sorted(breakpoints)
    for a, b in zip(ordered, ordered[1:]):
        midpoint = (a + b) / 2
        assert all(distance(p * midpoint) >= F(1, 14) for p in core)
    assert 2 * left == F(15, 49)
    assert 2 * right == F(9, 28)
    vertical = 2 * left - F(1, 7)
    assert vertical == F(8, 49)
    assert (148 + 4) * (right - left) == 1 + vertical
    assert F(8, 49 * (148 + 4)) > F(1, 7 * (148 + 8))

    print("THM-1134 five-residue multiplier Kakeya certificate")
    print(f"five_subsets={len(all_subsets)}")
    print(f"affine_orbits={len(rows)}")
    for representative, size, unit, word, maximum in rows:
        print(
            f"orbit rep={representative} size={size} "
            f"unit={unit} gaps={word} maximum={maximum}"
        )
    print(f"min_over_subsets_max_over_charts={min(maxima)}")
    print(f"sharp_subsets={sum(value == 6 for value in maxima)}")
    print(f"min_by_distinct_cardinality={minimum_by_size}")
    print(f"step2_core={core}")
    print(f"step2_J=[{left},{right}]")
    print(f"step2_J_length={right-left}")
    print(f"step2_fixed_vertical_width={vertical}")
    print("step2_first_geometric_base=148")
    print("step2_first_legal_base=157")
    print("tournament=affine residue orbit plus labelled cyclic-gap sidecar")
    print("destroyed_by_single_chart=maximum-over-multipliers predicate")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
