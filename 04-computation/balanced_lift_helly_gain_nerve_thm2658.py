#!/usr/bin/env python3
"""Exact referee for THM-2658's balanced lift-Helly gain nerve."""

from fractions import Fraction
from itertools import combinations_with_replacement, product


ZERO = Fraction(0)
ONE = Fraction(1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def overlap_length(first, second):
    return max(ZERO, min(first[1], second[1]) - max(first[0], second[0]))


def shifted(interval, amount):
    return interval[0] + amount, interval[1] + amount


def gains(first, second):
    """All integer n for which first meets second+n in positive length."""
    low = first[0] - second[1]
    high = first[1] - second[0]
    return tuple(
        integer
        for integer in range(floor_fraction(low) + 1, ceil_fraction(high))
        if overlap_length(first, shifted(second, integer)) > 0
    )


def gain_matrix(intervals):
    size = len(intervals)
    matrix = [[tuple() for _ in range(size)] for _ in range(size)]
    for first in range(size):
        for second in range(first + 1, size):
            values = gains(intervals[first], intervals[second])
            matrix[first][second] = values
            matrix[second][first] = tuple(sorted(-value for value in values))
            require(matrix[second][first] == gains(intervals[second], intervals[first]),
                    "gain antisymmetry")
    return matrix


def balanced_sections(intervals):
    """Return (potentials, common length) for every balanced gain section."""
    require(len(intervals) >= 2, "at least two arcs")
    require(all(0 < right - left < 1 for left, right in intervals),
            "proper open arc lifts")
    matrix = gain_matrix(intervals)
    result = []
    for tail in product(*(matrix[0][index] for index in range(1, len(intervals)))):
        potentials = (0,) + tail
        if any(
            potentials[second] - potentials[first] not in matrix[first][second]
            for first in range(len(intervals))
            for second in range(first + 1, len(intervals))
        ):
            continue
        lifts = [shifted(interval, potential)
                 for interval, potential in zip(intervals, potentials)]
        left = max(interval[0] for interval in lifts)
        right = min(interval[1] for interval in lifts)
        require(left < right, "pairwise interval Helly failure")
        common = right - left
        pair_lengths = [
            overlap_length(
                intervals[first],
                shifted(intervals[second],
                        potentials[second] - potentials[first]),
            )
            for first in range(len(intervals))
            for second in range(first + 1, len(intervals))
        ]
        require(common == min(pair_lengths), "minimum pair length identity")
        result.append((potentials, common))
    return tuple(result)


def circle_pieces(interval):
    """Cut the projection of one proper lift into disjoint [0,1] pieces."""
    left, right = interval
    pieces = []
    for integer in range(floor_fraction(left) - 1, ceil_fraction(right) + 1):
        piece = (max(ZERO, left - integer), min(ONE, right - integer))
        if piece[0] < piece[1]:
            pieces.append(piece)
    pieces.sort()
    return tuple(pieces)


def circle_intersection_pieces(intervals):
    pieces = ((ZERO, ONE),)
    for interval in intervals:
        next_pieces = []
        for first in pieces:
            for second in circle_pieces(interval):
                left = max(first[0], second[0])
                right = min(first[1], second[1])
                if left < right:
                    next_pieces.append((left, right))
        pieces = tuple(sorted(next_pieces))
    return pieces


def circle_intersection_measure(intervals):
    return sum((right - left for left, right in circle_intersection_pieces(intervals)),
               ZERO)


def check_family(intervals):
    sections = balanced_sections(intervals)
    predicted = sum((length for _, length in sections), ZERO)
    actual = circle_intersection_measure(intervals)
    require(predicted == actual, "balanced-section measure formula")
    return len(sections), actual


def union_intersection_measure(charts):
    """Circle measure of a common intersection of finite disjoint arc unions."""
    choices = product(*charts)
    return sum((circle_intersection_measure(choice) for choice in choices), ZERO)


def main():
    # Complete fifth-grid bank up to relabelling: every proper arc and every
    # multiset family of sizes 2--4.
    arcs = tuple(
        (Fraction(start, 5), Fraction(start + length, 5))
        for start in range(5)
        for length in range(1, 5)
    )
    by_size = []
    family_checks = 0
    positive_families = 0
    section_count = 0
    measure_checksum = ZERO
    switch_checks = 0
    for size in range(2, 5):
        count = 0
        for indices in combinations_with_replacement(range(len(arcs)), size):
            intervals = tuple(arcs[index] for index in indices)
            sections, measure = check_family(intervals)
            family_checks += 1
            count += 1
            section_count += sections
            positive_families += measure > 0
            measure_checksum += measure
            if size == 3:
                shifts = (0, 1, -1)
                switched = tuple(shifted(interval, amount)
                                 for interval, amount in zip(intervals, shifts))
                switched_lengths = sorted(length for _, length
                                          in balanced_sections(switched))
                original_lengths = sorted(length for _, length
                                          in balanced_sections(intervals))
                require(switched_lengths == original_lengths,
                        "reference-lift switching changed decorated nerve")
                switch_checks += 1
        by_size.append((size, count))

    # Integer winding is essential: all three pair edges exist, but their
    # unique gains have holonomy one and the triple intersection is empty.
    wrap = (
        (Fraction(-1, 10), Fraction(2, 5)),
        (Fraction(1, 4), Fraction(3, 4)),
        (Fraction(3, 5), Fraction(11, 10)),
    )
    wrap_gains = (gains(wrap[0], wrap[1]),
                  gains(wrap[1], wrap[2]),
                  gains(wrap[2], wrap[0]))
    wrap_pair_lengths = tuple(
        sum((overlap_length(first, shifted(second, amount))
             for amount in gains(first, second)), ZERO)
        for first, second in ((wrap[0], wrap[1]),
                              (wrap[1], wrap[2]),
                              (wrap[2], wrap[0]))
    )
    require(wrap_gains == ((0,), (0,), (1,))
            and wrap_pair_lengths == (Fraction(3, 20),
                                      Fraction(3, 20), Fraction(1, 5))
            and not balanced_sections(wrap)
            and circle_intersection_measure(wrap) == 0,
            "integer-winding hostile")

    # Disconnected chart labels can use incompatible components on their
    # three edges.  Componentwise refinement is therefore load-bearing.
    x = ((Fraction(0), Fraction(1, 10)),)
    y = ((Fraction(3, 10), Fraction(2, 5)),)
    z = ((Fraction(3, 5), Fraction(7, 10)),)
    disconnected = (x + z, x + y, y + z)
    disconnected_pairs = tuple(
        union_intersection_measure((disconnected[first], disconnected[second]))
        for first, second in ((0, 1), (1, 2), (2, 0))
    )
    disconnected_triple = union_intersection_measure(disconnected)
    require(disconnected_pairs == (Fraction(1, 10),) * 3
            and disconnected_triple == 0,
            "disconnected-component hostile")

    # Endpoint-only and incomplete-graph boundaries.
    endpoint_gain = gains((ZERO, Fraction(1, 3)),
                          (Fraction(1, 3), Fraction(2, 3)))
    incomplete = ((ZERO, Fraction(2, 5)),
                  (Fraction(3, 10), Fraction(7, 10)),
                  (Fraction(3, 5), ONE))
    incomplete_path = (gains(incomplete[0], incomplete[1]),
                       gains(incomplete[1], incomplete[2]))
    require(endpoint_gain == tuple()
            and incomplete_path == ((0,), (0,))
            and not balanced_sections(incomplete)
            and circle_intersection_measure(incomplete) == 0,
            "endpoint or incomplete-graph hostile")

    print("THM2658 BALANCED LIFT HELLY GAIN NERVE EXACT REFEREE")
    print(f"grid_arcs={len(arcs)} unordered_family_checks={family_checks} "
          f"by_size={tuple(by_size)}")
    print(f"positive_families={positive_families} "
          f"balanced_sections={section_count} "
          f"measure_checksum={measure_checksum}")
    print(f"switch_checks={switch_checks}")
    print(f"wrap_gains={wrap_gains} pair_lengths={wrap_pair_lengths} "
          "balanced_sections=0")
    print(f"disconnected_pair_measures={disconnected_pairs} "
          f"triple={disconnected_triple}")
    print(f"endpoint_open_gain={endpoint_gain} "
          f"incomplete_path_gains={incomplete_path} full_balanced=0")
    print("PASS")


if __name__ == "__main__":
    main()
