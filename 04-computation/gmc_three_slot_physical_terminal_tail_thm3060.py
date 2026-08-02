#!/usr/bin/env python3
"""Exact referee for THM-3060.

The proof is asymptotic and works for every fixed depth and terminal gap.  This
companion reconstructs the *physical* signed inclusion forms over Q, checks the
two-line initial resultant, exercises the positive physical resultant bank,
and tests finite generalized-Hankel recovery along rational points of the
carrier-to-physical homotopy.  It also checks the rank-two collapse of the
same initial face from four slots onward.
"""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(value, length):
    answer = Fraction(1)
    for offset in range(length):
        answer *= value + offset
    return answer


def determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    size = len(work)
    sign = 1
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        diagonal = work[column][column]
        answer *= diagonal
        for row in range(column + 1, size):
            multiplier = work[row][column] / diagonal
            for index in range(column + 1, size):
                work[row][index] -= multiplier * work[column][index]
    return sign * answer


def sylvester_resultant(first, second):
    """Resultant of ascending univariate coefficient lists."""
    degree_first = len(first) - 1
    degree_second = len(second) - 1
    first_descending = list(reversed(first))
    second_descending = list(reversed(second))
    rows = []
    for shift in range(degree_second):
        rows.append(
            [Fraction(0)] * shift
            + first_descending
            + [Fraction(0)] * (degree_second - 1 - shift)
        )
    for shift in range(degree_first):
        rows.append(
            [Fraction(0)] * shift
            + second_descending
            + [Fraction(0)] * (degree_first - 1 - shift)
        )
    return determinant(rows)


def normalized_tensor(depth, offsets):
    order = len(offsets)
    return rising(order * depth + 1, sum(offsets)) / reduce(
        lambda value, offset: value * rising(depth + 1, offset),
        offsets,
        Fraction(1),
    )


def signed_inclusion(depth, width, directions):
    answer = Fraction(0)
    for mask in range(1 << len(directions)):
        offsets = tuple(
            width if mask & (1 << index) else direction
            for index, direction in enumerate(directions)
        )
        answer += (-1) ** mask.bit_count() * normalized_tensor(depth, offsets)
    return answer


def physical_form(order, depth, moving, gap):
    width = moving + gap
    return [
        comb(order, moving_count)
        * signed_inclusion(
            depth,
            width,
            (moving,) * moving_count + (0,) * (order - moving_count),
        )
        for moving_count in range(order + 1)
    ]


def physical_resultant(depth, moving, gap):
    return sylvester_resultant(
        physical_form(2, depth, moving, gap),
        physical_form(3, depth, moving, gap),
    )


def width_flag(depth, width):
    # For three slots THM-3047 has (A,B,I)=(5,2,7).
    return (
        rising(depth + 1, width) ** 5
        * rising(depth + 2, width) ** 2
        / depth ** (7 * width)
    )


def beta_carrier(moving):
    # Moments of 46656*Beta(1,7/2).
    return (
        46656**moving
        * rising(Fraction(1), moving)
        / rising(Fraction(9, 2), moving)
    )


def generalized_minor(sequence, base, rows, columns):
    return determinant(
        [[sequence[base + row + column] for column in columns] for row in rows]
    )


# The all-large inclusion face is an exact binomial power.
face_cells = 0
for order in (2, 3):
    for gap in range(1, 9):
        for moving_count in range(order + 1):
            direct = sum(
                (-1) ** (order - moving_count + replaced)
                * comb(moving_count, replaced)
                * order ** ((order - moving_count + replaced) * gap)
                for replaced in range(moving_count + 1)
            )
            expected = (
                (-order**gap) ** (order - moving_count)
                * (1 - order**gap) ** moving_count
            )
            require(direct == expected, "all-large finite difference changed")
            face_cells += 1


# The K3 leading resultant is the sixth power of the line determinant.
line_cells = 0
for gap in range(1, 9):
    line_two = [-2**gap, 1 - 2**gap]
    line_three = [-3**gap, 1 - 3**gap]
    square = [line_two[0] ** 2, 2 * line_two[0] * line_two[1], line_two[1] ** 2]
    cube = [
        line_three[0] ** 3,
        3 * line_three[0] ** 2 * line_three[1],
        3 * line_three[0] * line_three[1] ** 2,
        line_three[1] ** 3,
    ]
    direct = sylvester_resultant(square, cube)
    require(direct == (3**gap - 2**gap) ** 6, "line resultant changed")
    line_cells += 1


# Exact physical Q-resultants.  This is evidence, not the all-parameter proof.
physical_cells = 0
physical_digest_rows = []
for depth in range(1, 6):
    for gap in range(1, 6):
        row = []
        for moving in range(1, 16):
            value = physical_resultant(depth, moving, gap)
            require(value > 0, "positive physical-resultant control failed")
            row.append(f"{value.numerator}/{value.denominator}")
            physical_cells += 1
        physical_digest_rows.append(f"{depth}:{gap}:" + ",".join(row))
physical_digest = sha256("\n".join(physical_digest_rows).encode("ascii")).hexdigest()


# The normalized physical ratio has the predicted C^{-7/2} logarithmic slope.
ratio_cells = 0
for depth in range(1, 5):
    for gap in range(1, 5):
        moving = 80
        ratio = (
            physical_resultant(depth, moving + 1, gap)
            / physical_resultant(depth, moving, gap)
            / 46656
        )
        slope = moving * (1 - ratio)
        require(Fraction(3) < slope < Fraction(4), "physical slope left -7/2 window")
        ratio_cells += 1


# Exact carrier and finite physical/homotopy generalized-Hankel controls.
carrier_cells = 0
homotopy_cells = 0
threshold_records = []
theta_grid = tuple(Fraction(index, 8) for index in range(9))
for depth in range(1, 5):
    for gap in range(1, 4):
        maximum = 28
        actual = {
            moving: width_flag(depth, moving + gap)
            * physical_resultant(depth, moving, gap)
            for moving in range(1, maximum + 7)
        }
        carrier = {
            moving: width_flag(depth, moving + gap) * beta_carrier(moving)
            for moving in range(1, maximum + 7)
        }
        last_bad = 0
        for order in (1, 2, 3):
            offsets = tuple(range(order))
            for base in range(1, maximum + 1):
                carrier_minor = generalized_minor(carrier, base, offsets, offsets)
                require(carrier_minor > 0, "strict Beta-Gamma carrier minor failed")
                carrier_cells += 1
                for theta in theta_grid:
                    mixed = {
                        index: (1 - theta) * carrier[index] + theta * actual[index]
                        for index in actual
                    }
                    value = generalized_minor(mixed, base, offsets, offsets)
                    if value <= 0:
                        last_bad = max(last_bad, base)
                    homotopy_cells += 1
        require(last_bad <= 2, "finite physical homotopy recovered too late")
        threshold_records.append((depth, gap, last_bad + 1))


# From K4 onward the first all-large face has rank at most two: all frozen
# columns coincide.  This is the exact boundary handed to THM-3058.
rank_collapse_cells = 0
for slot_count in range(4, 10):
    gap = slot_count + 1
    rows = [
        [-order**gap] * (slot_count - 2) + [1 - order**gap]
        for order in range(2, slot_count + 1)
    ]
    require(determinant(rows) == 0, "higher-slot initial-face rank did not collapse")
    rank_collapse_cells += 1


print("THM-3060 THREE-SLOT PHYSICAL TERMINAL FACE AND TAIL HOLOTOPY")
print(f"all_large_face_cells={face_cells} line_resultant_cells={line_cells}")
print(f"physical_resultant_cells={physical_cells} all_positive=1")
print(f"physical_resultant_digest={physical_digest}")
print(f"ratio_slope_cells={ratio_cells} target_exponent=-7/2")
print(f"carrier_minor_cells={carrier_cells}")
print(f"homotopy_minor_cells={homotopy_cells}")
print(f"finite_control_thresholds={threshold_records}")
print(f"higher_slot_rank_collapse_cells={rank_collapse_cells}")
print("scope=normalized intrinsic physical inclusion resultant;not raw chart or wall-stripped core")
print("all_exact_checks=PASS")
