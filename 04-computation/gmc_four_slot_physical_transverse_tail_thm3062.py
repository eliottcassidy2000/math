#!/usr/bin/env python3
"""Exact referee for THM-3062's four-slot physical augmentation."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from math import comb, factorial

from flint import fmpq_mat


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(value, length):
    answer = Fraction(1)
    for offset in range(length):
        answer *= value + offset
    return answer


def compositions(total, variables=3):
    if variables == 1:
        return [(total,)]
    return [
        (first,) + suffix
        for first in range(total + 1)
        for suffix in compositions(total - first, variables - 1)
    ]


MONOMIALS = {degree: compositions(degree) for degree in range(8)}
SELECTED_ROWS = (
    tuple(range(20))
    + tuple(range(21, 30))
    + (35,)
    + tuple(range(36, 42))
)


def multinomial(alpha):
    answer = factorial(sum(alpha))
    for entry in alpha:
        answer //= factorial(entry)
    return answer


def determinant_fraction(matrix):
    encoded = [
        [f"{value.numerator}/{value.denominator}" for value in row]
        for row in matrix
    ]
    value = fmpq_mat(encoded).det()
    return Fraction(int(value.numer()), int(value.denom()))


def normalized_tensor(depth, offsets):
    return rising(len(offsets) * depth + 1, sum(offsets)) / reduce(
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


def physical_forms(depth, first, second, moving, gap):
    width = moving + gap
    directions = (first, second, moving)
    forms = {}
    for order in (2, 3, 4):
        form = {}
        for alpha in MONOMIALS[order]:
            entries = tuple(
                directions[index]
                for index, count in enumerate(alpha)
                for _ in range(count)
            )
            form[alpha] = multinomial(alpha) * signed_inclusion(
                depth, width, entries
            )
        forms[order] = form
    return forms


def physical_layers(order, depth, first, second, moving, gap):
    """Split the exact form by the number of C/M occurrences."""
    width = moving + gap
    directions = (first, second, moving)
    layers = {count: {alpha: Fraction(0) for alpha in MONOMIALS[order]}
              for count in range(order + 1)}
    for alpha in MONOMIALS[order]:
        entries = tuple(
            (directions[index], index)
            for index, count in enumerate(alpha)
            for _ in range(count)
        )
        for mask in range(1 << order):
            offsets = tuple(
                width if mask & (1 << index) else entry[0]
                for index, entry in enumerate(entries)
            )
            high_count = sum(
                1
                for index, entry in enumerate(entries)
                if entry[1] == 2 or mask & (1 << index)
            )
            layers[high_count][alpha] += (
                multinomial(alpha)
                * (-1) ** mask.bit_count()
                * normalized_tensor(depth, offsets)
            )
    return layers


def normal_order(alpha, chosen_w_power):
    # x0=y, x1=w-y, x2=z.  A term chosen_w_power from x1 has normal
    # degree chosen_w_power+alpha[2].
    return chosen_w_power + alpha[2]


def check_layer_ideal(layer, required_order):
    # Expand x1^j=(w-y)^j and collect y,w,z.  Cancellation is checked exactly.
    expanded = {}
    for alpha, coefficient in layer.items():
        for w_power in range(alpha[1] + 1):
            key = (
                alpha[0] + alpha[1] - w_power,
                w_power,
                alpha[2],
            )
            expanded[key] = expanded.get(key, Fraction(0)) + (
                coefficient
                * comb(alpha[1], w_power)
                * (-1) ** (alpha[1] - w_power)
            )
    return all(
        coefficient == 0 or w_power + z_power >= required_order
        for (_y_power, w_power, z_power), coefficient in expanded.items()
    )


def intrinsic_resultant(forms):
    columns = {
        monomial: index for index, monomial in enumerate(MONOMIALS[7])
    }
    rows = []
    for order in (2, 3, 4):
        for beta in MONOMIALS[7 - order]:
            row = [Fraction(0)] * len(MONOMIALS[7])
            for alpha, coefficient in forms[order].items():
                target = tuple(alpha[index] + beta[index] for index in range(3))
                row[columns[target]] += coefficient
            rows.append(row)
    raw = determinant_fraction([rows[index] for index in SELECTED_ROWS])
    quadratic, cubic = forms[2], forms[3]
    q200 = quadratic[(2, 0, 0)]
    q110 = quadratic[(1, 1, 0)]
    q020 = quadratic[(0, 2, 0)]
    c300 = cubic[(3, 0, 0)]
    c210 = cubic[(2, 1, 0)]
    c120 = cubic[(1, 2, 0)]
    flag = (
        c120 * q200**2
        - c210 * q110 * q200
        - c300 * q020 * q200
        + c300 * q110**2
    )
    require(q200 and c300 and flag, "selected Macaulay chart hit its flag wall")
    return raw / (q200**6 * c300 * flag)


def evaluate(form, point):
    return sum(
        coefficient
        * reduce(
            lambda value, item: value * item[0] ** item[1],
            zip(point, alpha),
            Fraction(1),
        )
        for alpha, coefficient in form.items()
    )


def fixed_difference_moment(depth, first, second, order):
    return sum(
        Fraction((-1) ** count * comb(order, count))
        * normalized_tensor(
            depth,
            (first,) * (order - count) + (second,) * count,
        )
        for count in range(order + 1)
    )


def top_tensor(depth, order, moving):
    return normalized_tensor(depth, (moving,) * order)


def line_resultant(first, second):
    # Res((a+bt)^3,(c+dt)^4)=det^12.
    determinant = first[0] * second[1] - second[0] * first[1]
    return determinant**12


def carrier(depth, first, second, moving, gap):
    g2 = fixed_difference_moment(depth, first, second, 2)
    line_det = 4**gap - 3**gap
    return (
        g2**12
        * line_det**24
        * top_tensor(depth, 3, moving) ** 8
        * top_tensor(depth, 4, moving) ** 6
    )


def width_flag(depth, width):
    return (
        rising(depth, width) ** 26
        * rising(depth + 1, width) ** 20
        / depth ** (46 * width)
    )


def determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    sign = 1
    answer = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        diagonal = work[column][column]
        answer *= diagonal
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] / diagonal
            for index in range(column + 1, len(work)):
                work[row][index] -= multiplier * work[column][index]
    return sign * answer


# Exact normal-ideal layer typing and fixed transverse evaluation.
layer_cells = 0
transverse_cells = 0
for depth, first, second, moving, gap in (
    (1, 0, 1, 3, 1),
    (1, 0, 2, 5, 2),
    (2, 1, 3, 6, 1),
    (3, 0, 2, 7, 3),
):
    forms = physical_forms(depth, first, second, moving, gap)
    for order in (2, 3, 4):
        layers = physical_layers(order, depth, first, second, moving, gap)
        for high_count, layer in layers.items():
            require(
                check_layer_ideal(layer, high_count),
                "high layer escaped its normal-ideal power",
            )
            layer_cells += 1
        direct = evaluate(forms[order], (1, -1, 0))
        expected = fixed_difference_moment(depth, first, second, order)
        require(direct == expected, "transverse fixed form changed")
        transverse_cells += 1


# Exact powered-line and carrier factorizations.
line_cells = 0
beta_cells = 0
for gap in range(1, 9):
    line_three = (-3**gap, 1 - 3**gap)
    line_four = (-4**gap, 1 - 4**gap)
    require(
        line_resultant(line_three, line_four) == (4**gap - 3**gap) ** 12,
        "top-line resultant changed",
    )
    line_cells += 1
for depth in range(1, 5):
    for moving in range(1, 9):
        for order in (3, 4):
            product = Fraction(order ** (order * moving))
            for residue in range(1, order):
                product *= (
                    rising(Fraction(depth * order + residue, order), moving)
                    / rising(depth + 1, moving)
                )
            require(product == top_tensor(depth, order, moving), "Gauss/Beta split failed")
            beta_cells += 1


# Physical resultants and their positive transverse carriers.
physical_records = []
physical_cells = 0
families = (
    (1, 0, 1, 1, range(2, 17)),
    (1, 0, 2, 2, range(3, 13)),
    (2, 1, 3, 1, range(4, 12)),
)
cached = {}
for depth, first, second, gap, moving_values in families:
    for moving in moving_values:
        forms = physical_forms(depth, first, second, moving, gap)
        value = intrinsic_resultant(forms)
        comparison = carrier(depth, first, second, moving, gap)
        require(value > 0 and comparison > 0, "physical/carrier positivity control failed")
        ratio = value / comparison
        cached[(depth, first, second, moving, gap)] = value
        physical_records.append(
            f"{depth}:{first}:{second}:{moving}:{gap}:"
            f"{ratio.numerator}/{ratio.denominator}"
        )
        physical_cells += 1
physical_digest = sha256("\n".join(physical_records).encode("ascii")).hexdigest()


# Far-tail exact controls for the predicted limit one.
tail_ratios = []
for moving in (40, 60, 80):
    value = intrinsic_resultant(physical_forms(1, 0, 1, moving, 1))
    ratio = value / carrier(1, 0, 1, moving, 1)
    require(ratio > 1, "tail ratio crossed its positive limit")
    tail_ratios.append(ratio)
require(tail_ratios[0] > tail_ratios[1] > tail_ratios[2], "tail ratios lost descent")
tail_ratio_scaled = tuple(int(value * 10**9) for value in tail_ratios)
require(
    tail_ratio_scaled == (25913835936, 11091131268, 6782003397),
    "exact tail-ratio interval controls changed",
)


# Exact finite carrier and straight-homotopy Hankel controls.
homotopy_cells = 0
carrier_minor_cells = 0
depth, first, second, gap = 1, 0, 1, 1
values = {}
carriers = {}
for moving in range(5, 24):
    result = cached.get((depth, first, second, moving, gap))
    if result is None:
        result = intrinsic_resultant(
            physical_forms(depth, first, second, moving, gap)
        )
    values[moving] = width_flag(depth, moving + gap) * result
    carriers[moving] = width_flag(depth, moving + gap) * carrier(
        depth, first, second, moving, gap
    )
for order in (1, 2, 3):
    offsets = tuple(range(order))
    for base in range(5, 18):
        carrier_minor = determinant(
            [[carriers[base + row + column] for column in offsets] for row in offsets]
        )
        require(carrier_minor > 0, "K4 carrier minor failed")
        carrier_minor_cells += 1
        for theta in (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
            mixed = {
                index: (1 - theta) * carriers[index] + theta * values[index]
                for index in values
            }
            minor = determinant(
                [[mixed[base + row + column] for column in offsets] for row in offsets]
            )
            require(minor > 0, "finite K4 homotopy minor failed")
            homotopy_cells += 1


print("THM-3062 FOUR-SLOT PHYSICAL TRANSVERSE AUGMENTATION")
print(f"normal_ideal_layer_cells={layer_cells} transverse_cells={transverse_cells}")
print(f"powered_line_cells={line_cells} beta_factor_cells={beta_cells}")
print(f"physical_resultant_cells={physical_cells} all_positive=1")
print(f"physical_ratio_digest={physical_digest}")
print("tail_ratio_floor_times_1e9=" + ",".join(map(str, tail_ratio_scaled)))
print(f"carrier_minor_cells={carrier_minor_cells} homotopy_minor_cells={homotopy_cells}")
print("asymptotic=G2^12*(4^h-3^h)^24*U3^8*U4^6;relative_limit=1")
print("scope=normalized_intrinsic_physical_resultant;fixed_two_lows_and_gap;finite_order_tail")
print("all_exact_checks=PASS")
