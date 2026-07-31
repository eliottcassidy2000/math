#!/usr/bin/env python3
"""Exact referee for THM-2789.

The script checks:

1. every family consisting of one interval of each length 1,...,m on an
   m-point line, through m=8;
2. the quantitative one-dimensional Helly identity;
3. Boolean Mobius reconstruction of every incidence-column multiplicity
   from the interval Gram matrix;
4. the exact suffix/tail-column census;
5. sharp arbitrary-set, ordered-column, and circular-arc hostiles.

All truth-bearing checks use explicit exceptions, so ``python -O`` performs
the same verification.
"""

from __future__ import annotations

from itertools import product
from math import factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def interval_mask(start: int, length: int) -> int:
    return ((1 << length) - 1) << start


def gram_from_rows(rows: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple((row_i & row_j).bit_count() for row_j in rows)
        for row_i in rows
    )


def column_patterns(rows: tuple[int, ...], ambient_size: int) -> tuple[int, ...]:
    patterns = []
    for point in range(ambient_size):
        pattern = 0
        for index, row in enumerate(rows):
            if (row >> point) & 1:
                pattern |= 1 << index
        patterns.append(pattern)
    return tuple(patterns)


def pattern_counts(
    rows: tuple[int, ...], ambient_size: int
) -> tuple[int, ...]:
    counts = [0] * (1 << len(rows))
    for pattern in column_patterns(rows, ambient_size):
        counts[pattern] += 1
    return tuple(counts)


def intersection_counts(
    rows: tuple[int, ...], ambient_size: int
) -> tuple[int, ...]:
    rank = len(rows)
    full = (1 << ambient_size) - 1
    intersections = [0] * (1 << rank)
    intersections[0] = ambient_size
    masks = [full] * (1 << rank)
    for subset in range(1, 1 << rank):
        low_bit = subset & -subset
        index = low_bit.bit_length() - 1
        masks[subset] = masks[subset ^ low_bit] & rows[index]
        intersections[subset] = masks[subset].bit_count()
    return tuple(intersections)


def helly_counts_from_gram(
    gram: tuple[tuple[int, ...], ...], ambient_size: int
) -> tuple[int, ...]:
    rank = len(gram)
    row_minima = []
    for row in range(rank):
        minima = [ambient_size + 1] * (1 << rank)
        for subset in range(1, 1 << rank):
            low_bit = subset & -subset
            column = low_bit.bit_length() - 1
            minima[subset] = min(
                minima[subset ^ low_bit], gram[row][column]
            )
        row_minima.append(minima)

    counts = [ambient_size] + [0] * ((1 << rank) - 1)
    for subset in range(1, 1 << rank):
        low_bit = subset & -subset
        row = low_bit.bit_length() - 1
        rest = subset ^ low_bit
        counts[subset] = (
            gram[row][row]
            if rest == 0
            else min(counts[rest], row_minima[row][subset])
        )
    return tuple(counts)


def mobius_exact_patterns(
    upper_counts: tuple[int, ...], rank: int
) -> tuple[int, ...]:
    exact = list(upper_counts)
    for index in range(rank):
        bit = 1 << index
        for subset in range(1 << rank):
            if not subset & bit:
                exact[subset] -= exact[subset | bit]
    return tuple(exact)


def tail_mask(rank: int, first_difference: int) -> int:
    """The edge-difference suffix {first_difference,...,rank}."""
    require(
        1 <= first_difference <= rank,
        "tail start is outside the edge-difference range",
    )
    return ((1 << rank) - 1) ^ ((1 << (first_difference - 1)) - 1)


def interval_family_audit(max_rank: int = 8):
    summary = {}
    total_families = 0
    total_subsets = 0
    for rank in range(1, max_rank + 1):
        starts_by_length = [
            range(rank - length + 1) for length in range(1, rank + 1)
        ]
        family_count = 0
        tail_count = 0
        gram_to_patterns = {}
        mixed_gram_fibres = 0
        for starts in product(*starts_by_length):
            rows = tuple(
                interval_mask(start, length)
                for length, start in enumerate(starts, start=1)
            )
            gram = gram_from_rows(rows)
            literal_intersections = intersection_counts(rows, rank)
            helly_intersections = helly_counts_from_gram(gram, rank)
            require(
                literal_intersections == helly_intersections,
                "quantitative interval Helly identity failed",
            )
            reconstructed = mobius_exact_patterns(
                helly_intersections, rank
            )
            actual = pattern_counts(rows, rank)
            require(
                reconstructed == actual,
                "Mobius reconstruction lost an interval column",
            )
            require(
                all(value >= 0 for value in reconstructed),
                "reconstructed interval multiplicity became negative",
            )

            previous = gram_to_patterns.setdefault(gram, actual)
            if previous != actual:
                mixed_gram_fibres += 1

            if any(
                actual[tail_mask(rank, first)] > 0
                for first in range(1, rank + 1)
            ):
                tail_count += 1
            family_count += 1
            total_subsets += (1 << rank) - 1

        require(
            family_count == factorial(rank),
            "one-of-each-length interval-family count changed",
        )
        require(
            mixed_gram_fibres == 0,
            "one Gram matrix reconstructed two column multisets",
        )
        summary[rank] = (
            family_count,
            len(gram_to_patterns),
            tail_count,
        )
        total_families += family_count

    expected_tail_counts = {
        1: 1,
        2: 2,
        3: 6,
        4: 24,
        5: 120,
        6: 718,
        7: 4990,
        8: 39744,
    }
    require(
        {rank: values[2] for rank, values in summary.items()}
        == expected_tail_counts,
        "tail-column census changed",
    )
    return summary, total_families, total_subsets


def arbitrary_set_hostile():
    # Columns {12},{13},{23},empty versus {123},{1},{2},{3}.
    first = (0b011, 0b101, 0b110, 0b000)
    second = (0b111, 0b001, 0b010, 0b100)

    def rows_from_columns(columns):
        return tuple(
            sum(
                ((column >> row) & 1) << point
                for point, column in enumerate(columns)
            )
            for row in range(3)
        )

    rows_first = rows_from_columns(first)
    rows_second = rows_from_columns(second)
    gram_first = gram_from_rows(rows_first)
    gram_second = gram_from_rows(rows_second)
    require(gram_first == gram_second, "set hostile lost its common Gram")
    require(
        gram_first == ((2, 1, 1), (1, 2, 1), (1, 1, 2)),
        "set hostile Gram changed",
    )
    require(
        (rows_first[0] & rows_first[1] & rows_first[2]).bit_count() == 0,
        "first set hostile gained a triple point",
    )
    require(
        (rows_second[0] & rows_second[1] & rows_second[2]).bit_count() == 1,
        "second set hostile lost its triple point",
    )
    return gram_first


def ordered_column_hostile():
    # One interval of each length 1,2,3 on a three-point line.
    first_rows = (0b001, 0b011, 0b111)
    second_rows = (0b010, 0b011, 0b111)
    first_columns = column_patterns(first_rows, 3)
    second_columns = column_patterns(second_rows, 3)
    require(
        gram_from_rows(first_rows) == gram_from_rows(second_rows),
        "ordered hostile lost its common Gram",
    )
    require(
        sorted(first_columns) == sorted(second_columns),
        "ordered hostile lost its common column multiset",
    )
    require(
        first_columns != second_columns,
        "ordered hostile unexpectedly retained column order",
    )
    require(
        tuple(reversed(first_columns)) != second_columns,
        "ordered hostile collapsed to reflection",
    )
    return first_columns, second_columns


def circular_arc_hostile():
    # Three length-three arcs covering C6 cyclically.
    arcs = (
        (1 << 0) | (1 << 1) | (1 << 2),
        (1 << 2) | (1 << 3) | (1 << 4),
        (1 << 4) | (1 << 5) | (1 << 0),
    )
    gram = gram_from_rows(arcs)
    triple = (arcs[0] & arcs[1] & arcs[2]).bit_count()
    pair_min = min(gram[i][j] for i in range(3) for j in range(3))
    require(
        gram == ((3, 1, 1), (1, 3, 1), (1, 1, 3)),
        "circular-arc hostile Gram changed",
    )
    require(
        (triple, pair_min) == (0, 1),
        "circular-arc hostile no longer breaks interval Helly tomography",
    )
    require(
        (arcs[0] | arcs[1] | arcs[2]) == (1 << 6) - 1,
        "circular hostile gained a common omitted cut point",
    )
    return gram, triple, pair_min


def main() -> None:
    summary, total_families, total_subsets = interval_family_audit()
    set_gram = arbitrary_set_hostile()
    ordered_first, ordered_second = ordered_column_hostile()
    circle_gram, circle_triple, circle_pair_min = circular_arc_hostile()

    print("THM2789 INTERVAL GRAM / GAP-TAIL EXACT REFEREE")
    for rank in sorted(summary):
        families, gram_classes, tails = summary[rank]
        print(
            f"m={rank}: families={families} "
            f"gram_classes={gram_classes} tail_present={tails}"
        )
    print(
        "totals:",
        total_families,
        "families;",
        total_subsets,
        "nonempty subset intersections",
    )
    print("arbitrary-set common Gram:", set_gram, "triple=0/1")
    print("ordered interval columns:", ordered_first, ordered_second)
    print(
        "circular C6 Gram:",
        circle_gram,
        f"triple={circle_triple} pair_min={circle_pair_min}",
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
