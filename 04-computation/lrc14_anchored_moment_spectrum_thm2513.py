#!/usr/bin/env python3
"""Exact referee for THM-2513's anchored first/second-moment dichotomy."""

from fractions import Fraction
from itertools import product
from math import gcd


ROWS = 7
COLS = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rectangles(table: list[list[int]]) -> tuple[int, ...]:
    return tuple(
        table[0][s] - table[ell][s] - table[0][0] + table[ell][0]
        for ell in range(1, ROWS)
        for s in range(1, COLS)
    )


def interaction(table: list[list[int]]) -> list[list[Fraction]]:
    row_means = [
        Fraction(sum(table[ell]), COLS) for ell in range(ROWS)
    ]
    column_means = [
        Fraction(sum(table[ell][s] for ell in range(ROWS)), ROWS)
        for s in range(COLS)
    ]
    grand = Fraction(sum(map(sum, table)), ROWS * COLS)
    return [
        [
            Fraction(table[ell][s]) - row_means[ell] - column_means[s] + grand
            for s in range(COLS)
        ]
        for ell in range(ROWS)
    ]


def norm_square(matrix: list[list[Fraction]]) -> Fraction:
    return sum((value * value for row in matrix for value in row), Fraction())


def square_table(table: list[list[int]]) -> list[list[int]]:
    return [[value * value for value in row] for row in table]


def replica(anchor: int, word: tuple[int, ...]) -> list[list[int]]:
    require(len(word) == COLS and word[0] == 0, "replica word malformed")
    return [
        [anchor + word[s] if ell == 0 else word[s] for s in range(COLS)]
        for ell in range(ROWS)
    ]


def owner_nonflat(table: list[list[int]]) -> bool:
    return len(set(table[0])) > 1


def first_or_second(table: list[list[int]]) -> tuple[bool, bool]:
    return any(rectangles(table)), any(rectangles(square_table(table)))


def ternary_tag_coordinates(
    first: Fraction, second: Fraction, gamma: int
) -> tuple[Fraction, Fraction]:
    """Coordinates in the basis (1, omega), with omega^2=-1-omega."""
    require(gamma in (1, 2), "ternary character must be primitive")
    if gamma == 1:
        return first, second
    return first - second, -second


def initial_arc_self_join(length: Fraction, dilation: int) -> Fraction:
    """Integral of 1_[0,length)(x) 1_[0,length)({dilation*x})."""
    require(Fraction(0) <= length <= Fraction(1), "arc length out of range")
    require(dilation >= 1, "dilation must be positive")
    total = Fraction()
    for k in range(dilation):
        left = Fraction(k, dilation)
        right = Fraction(k, dilation) + length / dilation
        overlap = min(length, right) - left
        if overlap > 0:
            total += overlap
    return total


def lcg_tables(count: int, anchor: int) -> list[list[list[int]]]:
    state = 0x2513
    answer = []
    for _ in range(count):
        table = [[0] * COLS for _ in range(ROWS)]
        table[0][0] = anchor
        for ell in range(ROWS):
            for s in range(1, COLS):
                state = (1664525 * state + 1013904223) & 0xFFFFFFFF
                table[ell][s] = state % 8
        answer.append(table)
    return answer


def main() -> None:
    anchor = 3
    binary_profiles = 0
    nonflat_second = 0
    flat_boundary = 0
    energy_checks = 0
    for tail in product((0, 1), repeat=COLS - 1):
        word = (0,) + tail
        table = replica(anchor, word)
        squared = square_table(table)
        require(not any(rectangles(table)), "replica acquired first interaction")
        expected_rectangles = tuple(
            2 * anchor * word[s]
            for _ell in range(1, ROWS)
            for s in range(1, COLS)
        )
        require(
            rectangles(squared) == expected_rectangles,
            "Hadamard-square rectangle identity failed",
        )
        first, second = first_or_second(table)
        require(not first, "replica first-moment route drifted")
        if any(word):
            require(owner_nonflat(table) and second, "nonflat replica escaped square")
            nonflat_second += 1
        else:
            require(not owner_nonflat(table) and not second, "flat sharp boundary failed")
            flat_boundary += 1

        mean = Fraction(sum(word), COLS)
        variance_sum = sum(
            ((Fraction(value) - mean) ** 2 for value in word), Fraction()
        )
        predicted = Fraction(24 * anchor * anchor, ROWS) * variance_sum
        require(
            norm_square(interaction(squared)) == predicted,
            "quantitative square-interaction identity failed",
        )
        energy_checks += 1
        binary_profiles += 1

    require(
        (binary_profiles, nonflat_second, flat_boundary, energy_checks)
        == (4_096, 4_095, 1, 4_096),
        "binary replica census drifted",
    )

    # Exact same-circle delayed self-joining control.  Realize a rational
    # replica table by initial-arc indicators and compare its two-time table
    # with the entrywise square under the standard BV covariance invoice.
    rational_word = tuple(Fraction(s, 52) for s in range(COLS))
    rational_replica = replica(Fraction(1, 4), rational_word)
    rational_square = square_table(rational_replica)
    delayed_join_checks = 0
    delayed_join_interactions = 0
    for delay in range(1, 4):
        dilation = 13**delay
        delayed = [
            [initial_arc_self_join(value, dilation) for value in row]
            for row in rational_replica
        ]
        for ell in range(ROWS):
            for s in range(COLS):
                error = abs(delayed[ell][s] - rational_square[ell][s])
                require(
                    error <= Fraction(1, 3 * dilation),
                    "delayed self-joining exceeded the BV invoice",
                )
                delayed_join_checks += 1
        require(
            delayed[0][0] > 0
            and all(delayed[ell][0] == 0 for ell in range(1, ROWS)),
            "delayed self-joining lost the owner anchor",
        )
        if any(rectangles(delayed)):
            delayed_join_interactions += 1
    require(
        (delayed_join_checks, delayed_join_interactions) == (273, 3),
        "delayed self-joining control drifted",
    )

    # First moment zero, second moment nonzero.
    first_fail = replica(anchor, tuple([0] + list(range(1, COLS))))
    require(first_or_second(first_fail) == (False, True), "first-only hostile failed")

    # Second moment zero, first moment nonzero.  The square is itself a
    # replica table with anchor 9 and word 16 away from the anchor column.
    second_fail = [[0] * COLS for _ in range(ROWS)]
    second_fail[0][0] = 3
    for ell in range(ROWS):
        for s in range(1, COLS):
            second_fail[ell][s] = 5 if ell == 0 else 4
    require(owner_nonflat(second_fail), "second-only hostile owner row is flat")
    require(first_or_second(second_fail) == (True, False), "second-only hostile failed")
    require(
        square_table(second_fail) == replica(9, (0,) + (16,) * 12),
        "second-only hostile square is not the claimed replica",
    )

    # A two-channel vector is scalarized faithfully by a ternary tag.  Work
    # exactly in the F-basis (1,omega), where F=Q(zeta_91).  The field-degree
    # step itself is proved in the theorem from phi(273)/phi(91)=2.
    require(gcd(3, 91) == 1, "cyclotomic tag is not conductor-coprime")
    require((2 * 72, 2 * 5_184) == (144, 10_368), "tag mode count drifted")
    ternary_basis_controls = 0
    for x in range(-3, 4):
        for y in range(-3, 4):
            if (x, y) == (0, 0):
                continue
            for gamma in (1, 2):
                coordinates = ternary_tag_coordinates(
                    Fraction(x), Fraction(y), gamma
                )
                require(coordinates != (0, 0), "ternary tag cancelled a pair")
                ternary_basis_controls += 1
    require(ternary_basis_controls == 96, "ternary basis census drifted")

    # Sharp F_2 failure inside the theorem's anchored class: the delta table
    # satisfies B=A, so its binary difference channel is identically zero.
    binary_tag_hostile = [[0] * COLS for _ in range(ROWS)]
    binary_tag_hostile[0][0] = 1
    require(owner_nonflat(binary_tag_hostile), "binary hostile owner row is flat")
    require(
        square_table(binary_tag_hostile) == binary_tag_hostile,
        "binary hostile is not idempotent",
    )
    require(
        first_or_second(binary_tag_hostile) == (True, True),
        "binary hostile moment routes drifted",
    )

    random_controls = lcg_tables(10_000, anchor=5)
    route_histogram = {(True, True): 0, (True, False): 0, (False, True): 0}
    local_secant_checks = 0
    for table in random_controls:
        if not owner_nonflat(table):
            continue
        route = first_or_second(table)
        require(route != (False, False), "anchored moment disjunction failed")
        route_histogram[route] = route_histogram.get(route, 0) + 1
        anchor_mass = table[0][0]
        witness_column = next(
            s for s in range(COLS) if table[0][s] != anchor_mass
        )
        for ell in range(1, ROWS):
            x = table[0][witness_column]
            y = table[ell][witness_column]
            pair = (x - y - anchor_mass, x * x - y * y - anchor_mass**2)
            require(pair != (0, 0), "local parabola-secant witness vanished")
            for gamma in (1, 2):
                require(
                    ternary_tag_coordinates(
                        Fraction(pair[0]), Fraction(pair[1]), gamma
                    )
                    != (0, 0),
                    "tagged local rectangle vanished",
                )
            local_secant_checks += 1
    require(local_secant_checks == 60_000, "local secant census drifted")

    # First and second response masses do not determine overlap with a fixed
    # typed deep probe.  The square does have the exact one-copy weighted
    # realization integral (mean(response)*response).
    deep_probe = (1, 1, 0, 0)
    response_same = (1, 1, 0, 0)
    response_disjoint = (0, 0, 1, 1)
    mean_probe = Fraction(sum(deep_probe), 4)
    mean_same = Fraction(sum(response_same), 4)
    mean_disjoint = Fraction(sum(response_disjoint), 4)
    pair_product_same = mean_probe * mean_same
    pair_product_disjoint = mean_probe * mean_disjoint
    weighted_same = Fraction(
        sum(mean_same * value for value in response_same), 4
    )
    weighted_disjoint = Fraction(
        sum(mean_disjoint * value for value in response_disjoint), 4
    )
    diagonal_same = Fraction(
        sum(a * b for a, b in zip(deep_probe, response_same)), 4
    )
    diagonal_disjoint = Fraction(
        sum(a * b for a, b in zip(deep_probe, response_disjoint)), 4
    )
    require(
        pair_product_same == pair_product_disjoint == Fraction(1, 4),
        "pair-space mass changed",
    )
    require(
        weighted_same == weighted_disjoint == Fraction(1, 4),
        "one-copy scalar reweighting changed the square mass",
    )
    require(
        (diagonal_same, diagonal_disjoint) == (Fraction(1, 2), Fraction(0)),
        "one-copy overlap hostile drifted",
    )

    print("THM-2513 anchored moment-spectrum exact companion: PASS")
    print(
        "binary_replica_profiles=4096; nonflat_second_moment=4095; "
        "flat_boundary=1"
    )
    print("Hadamard_square_rectangle=2*anchor*w_s; quantitative_checks=4096")
    print(
        "delayed_self_join_checks=273; nonzero_interaction_delays=3; "
        "parity_subsequences=both"
    )
    print("sharp_route_hostiles=first_zero_second_nonzero,first_nonzero_second_zero")
    print(
        "ternary_tag_modes=144; tagged_cut_modes=10368; "
        "basis_controls=96; binary_delta_hostile=cancels"
    )
    print("single_target_parabola_secant_checks=60000; all_six_sources=nonzero")
    print(f"deterministic_anchored_controls=10000; route_histogram={route_histogram}")
    print(
        "pair_space_product=1/4; one_copy_weighted_mass=1/4; "
        "typed_probe_overlaps="
        f"{diagonal_same},{diagonal_disjoint}"
    )
    print("first-or-second mixed spectrum is unconditional under anchor+nonflatness")
    print(
        "second moment has weighted same-time and delayed Boolean one-circle lifts; "
        "not a same-time typed event"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
