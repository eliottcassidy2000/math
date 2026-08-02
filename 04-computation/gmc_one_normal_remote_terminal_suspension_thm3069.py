#!/usr/bin/env python3
"""Exact companion for THM-3069's one-normal physical suspension."""

from decimal import Decimal, getcontext
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import comb, factorial
import sys

from flint import fmpq_mat


sys.set_int_max_str_digits(0)
getcontext().prec = 80


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
    """L(prod f_(n+a))/L(f_n^r)."""
    return rising(len(offsets) * depth + 1, sum(offsets)) / reduce(
        lambda value, offset: value * rising(depth + 1, offset),
        offsets,
        Fraction(1),
    )


def signed_terminal_inclusion(depth, terminal, directions):
    """Moment of prod_i(f_(n+directions_i)-f_(n+terminal))."""
    answer = Fraction(0)
    for mask in range(1 << len(directions)):
        offsets = tuple(
            terminal if mask & (1 << index) else direction
            for index, direction in enumerate(directions)
        )
        answer += (-1) ** mask.bit_count() * normalized_tensor(depth, offsets)
    return answer


def physical_forms(depth, lows, terminal):
    forms = {}
    for order in (2, 3, 4):
        form = {}
        for alpha in MONOMIALS[order]:
            entries = tuple(
                lows[index]
                for index, count in enumerate(alpha)
                for _ in range(count)
            )
            form[alpha] = multinomial(alpha) * signed_terminal_inclusion(
                depth, terminal, entries
            )
        forms[order] = form
    return forms


def transform_to_yw(form):
    """Substitute x_2=w-y_0-y_1 in a ternary homogeneous form."""
    transformed = {}
    for alpha, coefficient in form.items():
        for split in compositions(alpha[2], 3):
            w_power, y0_extra, y1_extra = split
            target = (
                alpha[0] + y0_extra,
                alpha[1] + y1_extra,
                w_power,
            )
            transformed[target] = transformed.get(target, Fraction(0)) + (
                coefficient
                * multinomial(split)
                * (-1) ** (y0_extra + y1_extra)
            )
    return transformed


def remote_layer_form(depth, lows, terminal, order, remote_count):
    """The exact ell-remote part after Y=E(y)+w(f_c-f_C)."""
    first, second, last = lows
    layer = {alpha: Fraction(0) for alpha in MONOMIALS[order]}
    for y0_count, y1_count, w_count in MONOMIALS[order]:
        if remote_count > w_count:
            continue
        fixed_bases = (first,) * y0_count + (second,) * y1_count
        subtotal = Fraction(0)
        for mask in range(1 << (y0_count + y1_count)):
            offsets = tuple(
                last if mask & (1 << index) else base
                for index, base in enumerate(fixed_bases)
            )
            offsets += (last,) * (w_count - remote_count)
            offsets += (terminal,) * remote_count
            subtotal += (-1) ** mask.bit_count() * normalized_tensor(
                depth, offsets
            )
        subtotal *= (-1) ** remote_count * comb(w_count, remote_count)
        layer[(y0_count, y1_count, w_count)] = (
            multinomial((y0_count, y1_count, w_count)) * subtotal
        )
    return layer


def intrinsic_resultant(forms):
    """Standard ternary resultant through the proved THM-2942 chart."""
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


def binary_resultant(form_f, degree_f, form_g, degree_g):
    coeff_f = [
        form_f.get((degree_f - index, index), Fraction(0))
        for index in range(degree_f + 1)
    ]
    coeff_g = [
        form_g.get((degree_g - index, index), Fraction(0))
        for index in range(degree_g + 1)
    ]
    descending_f = list(reversed(coeff_f))
    descending_g = list(reversed(coeff_g))
    size = degree_f + degree_g
    rows = []
    for shift in range(degree_g):
        rows.append(
            [Fraction(0)] * shift
            + descending_f
            + [Fraction(0)] * (degree_g - 1 - shift)
        )
    for shift in range(degree_f):
        rows.append(
            [Fraction(0)] * shift
            + descending_g
            + [Fraction(0)] * (degree_f - 1 - shift)
        )
    require(all(len(row) == size for row in rows), "Sylvester shape changed")
    return determinant_fraction(rows)


def low_binary_form(depth, lows, order):
    first, second, last = lows
    form = {}
    for y0_count in range(order + 1):
        y1_count = order - y0_count
        bases = (first,) * y0_count + (second,) * y1_count
        total = Fraction(0)
        for mask in range(1 << order):
            offsets = tuple(
                last if mask & (1 << index) else base
                for index, base in enumerate(bases)
            )
            total += (-1) ** mask.bit_count() * normalized_tensor(depth, offsets)
        form[(y0_count, y1_count)] = comb(order, y0_count) * total
    return form


def transverse_resultant(depth, lows):
    return binary_resultant(
        low_binary_form(depth, lows, 2),
        2,
        low_binary_form(depth, lows, 3),
        3,
    )


def separated_forms(depth, lows, normal_sign):
    forms = {}
    for order in (2, 3, 4):
        binary = low_binary_form(depth, lows, order)
        form = {
            (y0_power, y1_power, 0): coefficient
            for (y0_power, y1_power), coefficient in binary.items()
        }
        if order == 4:
            form[(0, 0, 4)] = Fraction(normal_sign)
        forms[order] = form
    return forms


def top_tensor(depth, order, terminal):
    return normalized_tensor(depth, (terminal,) * order)


def gauss_top_tensor(depth, order, terminal):
    answer = Fraction(order ** (order * terminal))
    for residue in range(1, order):
        answer *= rising(Fraction(order * depth + residue, order), terminal)
        answer /= rising(depth + 1, terminal)
    return answer


def fraction_hash(value):
    payload = f"{value.numerator}/{value.denominator}".encode("ascii")
    return sha256(payload).hexdigest()


def decimal_text(value, places=15):
    decimal = Decimal(value.numerator) / Decimal(value.denominator)
    return f"{decimal:.{places}f}"


def carrier_value(depth, lows, terminal):
    sidecar = transverse_resultant(depth, lows)
    return sidecar**4 * top_tensor(depth, 4, terminal) ** 6


def maximal_partition_product(total, maximum_part):
    """Maximize prod part^part over positive parts summing to total."""
    best = [0] * (total + 1)
    best[0] = 1
    for weight in range(1, total + 1):
        best[weight] = max(
            part**part * best[weight - part]
            for part in range(1, min(maximum_part, weight) + 1)
        )
    return best[total]


# General covariance, sign, isotropic gap, and resultant-character gap ledger.
general_cells = 0
layer_base_cells = 0
partition_cells = 0
for slots in range(3, 13):
    mu = factorial(slots - 1)
    require(slots * mu == factorial(slots), "one-normal covariance changed")
    require((-1) ** (slots * mu) == 1, "surviving normal sign changed")
    isotropic_gap = Fraction(slots - 1, slots) ** (slots - 1)
    resultant_gap = isotropic_gap / slots
    require(resultant_gap < isotropic_gap < 1, "gap ordering changed")
    equality_cells = []
    for w_count in range(1, slots + 1):
        for remote_count in range(w_count + 1):
            if w_count == remote_count == slots:
                continue
            remote_power = 1 if remote_count == 0 else remote_count**remote_count
            base = Fraction(remote_power, slots**w_count)
            require(base <= isotropic_gap, "physical layer gap failed")
            if base == isotropic_gap:
                equality_cells.append((w_count, remote_count))
            layer_base_cells += 1
    require(
        equality_cells == [(slots - 1, slots - 1)],
        "sharp physical layer boundary changed",
    )
    for removed_top_cells in range(1, 2 * slots + 1):
        redistributed = maximal_partition_product(
            slots * removed_top_cells, slots - 1
        )
        ratio = Fraction(
            redistributed, slots ** (slots * removed_top_cells)
        )
        require(ratio <= resultant_gap, "resultant partition gap failed")
        require(
            (ratio == resultant_gap) == (removed_top_cells == 1),
            "resultant partition equality boundary changed",
        )
        partition_cells += 1
    general_cells += 1


# Exact coordinate identity and remote-layer partition on physical K4 cells.
LOW_TRIPLES = ((0, 1, 2), (0, 1, 3), (0, 2, 5), (1, 4, 7))
coordinate_cells = 0
for depth, lows, terminal in (
    (0, (0, 1, 2), 7),
    (1, (0, 1, 3), 8),
    (2, (0, 2, 5), 10),
    (3, (1, 4, 7), 13),
):
    source = physical_forms(depth, lows, terminal)
    for order in (2, 3, 4):
        transformed = transform_to_yw(source[order])
        rebuilt = {alpha: Fraction(0) for alpha in MONOMIALS[order]}
        for remote_count in range(order + 1):
            layer = remote_layer_form(
                depth, lows, terminal, order, remote_count
            )
            for alpha, coefficient in layer.items():
                rebuilt[alpha] += coefficient
        require(transformed == rebuilt, "one-normal layer partition changed")
        coordinate_cells += 1


# Low sidecars and exact separated-resultant sign controls.
sidecar_records = []
separated_cells = 0
for lows in LOW_TRIPLES:
    sidecar = transverse_resultant(0, lows)
    require(sidecar > 0, "three-slot sidecar lost positivity")
    sidecar_records.append([list(lows), str(sidecar)])
    for normal_sign in (-1, 1):
        separated = intrinsic_resultant(separated_forms(0, lows, normal_sign))
        require(separated == sidecar**4, "separated resultant scalar changed")
        separated_cells += 1
sidecar_digest = sha256(
    json.dumps(sidecar_records, separators=(",", ":")).encode("ascii")
).hexdigest()


# Independent standard-orientation control for the intrinsic resultant chart.
orientation_forms = {
    2: {(2, 0, 0): Fraction(1), (0, 2, 0): Fraction(1)},
    3: {(3, 0, 0): Fraction(1), (0, 0, 3): Fraction(1)},
    4: {(0, 0, 4): Fraction(1)},
}
for degree in (2, 3, 4):
    orientation_forms[degree] = {
        alpha: orientation_forms[degree].get(alpha, Fraction(0))
        for alpha in MONOMIALS[degree]
    }
require(intrinsic_resultant(orientation_forms) == 1, "resultant orientation changed")


# Exact intrinsic physical bank and normalized convergence controls.
physical_records = []
ratio_rows = []
for lows in LOW_TRIPLES:
    sidecar = transverse_resultant(0, lows)
    row = []
    for gap in (3, 6, 9):
        terminal = lows[-1] + gap
        resultant = intrinsic_resultant(physical_forms(0, lows, terminal))
        carrier = sidecar**4 * top_tensor(0, 4, terminal) ** 6
        ratio = resultant / carrier
        require(resultant > 0 and ratio > 0, "physical resultant lost positivity")
        row.append(ratio)
        physical_records.append(
            [
                list(lows),
                terminal,
                str(sidecar),
                fraction_hash(resultant),
                fraction_hash(ratio),
                decimal_text(ratio),
            ]
        )
    require(row[0] < row[1] < row[2] < 1, "physical ratio control changed")
    ratio_rows.append(row)
physical_digest = sha256(
    json.dumps(physical_records, separators=(",", ":")).encode("ascii")
).hexdigest()


# Gauss/Beta factorization and finite strict carrier-Hankel controls.
beta_cells = 0
for depth in range(4):
    for terminal in range(1, 10):
        require(
            top_tensor(depth, 4, terminal)
            == gauss_top_tensor(depth, 4, terminal),
            "Gauss/Beta factorization changed",
        )
        beta_cells += 1

carrier_cells = 0
depth = 0
lows = (0, 1, 2)
for size in range(1, 5):
    offset_banks = (
        (tuple(range(size)), tuple(range(size))),
        (tuple(2 * index for index in range(size)), tuple(range(size))),
        (tuple(range(size)), tuple(3 * index for index in range(size))),
    )
    for row_offsets, column_offsets in offset_banks:
        base = 5
        matrix = [
            [
                carrier_value(depth, lows, base + row_offset + column_offset)
                for column_offset in column_offsets
            ]
            for row_offset in row_offsets
        ]
        require(determinant_fraction(matrix) > 0, "carrier Hankel minor lost positivity")
        carrier_cells += 1


print("THM-3069 ONE-NORMAL REMOTE-TERMINAL SUSPENSION")
print(
    f"general_k_cells={general_cells} layer_base_cells={layer_base_cells} "
    f"partition_cells={partition_cells} sharp_boundary=(k-1,k-1)"
)
print("isotropic_gap=((k-1)/k)^(k-1);resultant_gap=isotropic_gap/k;distinct=1")
print(f"coordinate_layer_cells={coordinate_cells} exact_partition=1")
print(f"sidecar_cells={len(sidecar_records)} all_positive=1")
print(f"sidecar_digest={sidecar_digest}")
print(f"separated_sign_cells={separated_cells} exact_S4=1 orientation=1")
print(f"physical_resultant_cells={len(physical_records)} all_positive=1")
print(f"physical_digest={physical_digest}")
print(
    "largest_gap_ratios="
    + ",".join(decimal_text(row[-1]) for row in ratio_rows)
)
print(f"beta_factor_cells={beta_cells} carrier_hankel_cells={carrier_cells}")
print("k4_asymptotic=S3^4*U4^6;base=4^(24C);power=-9;gaps=27/64,27/256")
print("scope=fixed_lower_support;remote_terminal;finite_order_holotopy")
print("all_exact_checks=PASS")
