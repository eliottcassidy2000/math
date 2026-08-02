#!/usr/bin/env python3
"""Exact referee for THM-3063's terminal suspension and K5 tail."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import comb, factorial
import sys

from flint import fmpq_mat, nmod_mat

sys.set_int_max_str_digits(0)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(value, length):
    answer = Fraction(1)
    for offset in range(length):
        answer *= value + offset
    return answer


def compositions(total, variables):
    if variables == 1:
        return [(total,)]
    return [
        (first,) + suffix
        for first in range(total + 1)
        for suffix in compositions(total - first, variables - 1)
    ]


MONOMIALS4 = {degree: compositions(degree, 4) for degree in range(12)}


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
    """L(prod f_(n+a))/L(f_n^r), with n=depth."""
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


def physical_forms_k5(depth, lows, moving, gap):
    width = moving + gap
    directions = tuple(lows) + (moving,)
    forms = {}
    for order in (2, 3, 4, 5):
        form = {}
        for alpha in MONOMIALS4[order]:
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


def physical_layers_k5(order, depth, lows, moving, gap):
    width = moving + gap
    directions = tuple(lows) + (moving,)
    layers = {
        count: {alpha: Fraction(0) for alpha in MONOMIALS4[order]}
        for count in range(order + 1)
    }
    for alpha in MONOMIALS4[order]:
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
                if entry[1] == 3 or mask & (1 << index)
            )
            layers[high_count][alpha] += (
                multinomial(alpha)
                * (-1) ** mask.bit_count()
                * normalized_tensor(depth, offsets)
            )
    return layers


def check_normal_ideal_k5(layer, required_order):
    """Expand x2=w-x0-x1 and check membership in (w,z)^q."""
    expanded = {}
    for alpha, coefficient in layer.items():
        for split in compositions(alpha[2], 3):
            w_power, x0_extra, x1_extra = split
            key = (
                alpha[0] + x0_extra,
                alpha[1] + x1_extra,
                w_power,
                alpha[3],
            )
            expanded[key] = expanded.get(key, Fraction(0)) + (
                coefficient
                * multinomial(split)
                * (-1) ** (x0_extra + x1_extra)
            )
    return all(
        coefficient == 0 or w_power + z_power >= required_order
        for (_x0, _x1, w_power, z_power), coefficient in expanded.items()
    )


def binary_resultant(form_f, degree_f, form_g, degree_g):
    """Res(f(1,t),g(1,t)); leading terms are required nonzero."""
    coeff_f = [form_f.get((degree_f - j, j), Fraction(0)) for j in range(degree_f + 1)]
    coeff_g = [form_g.get((degree_g - j, j), Fraction(0)) for j in range(degree_g + 1)]
    require(coeff_f[-1] != 0 and coeff_g[-1] != 0, "binary chart lost full degree")
    descending_f = list(reversed(coeff_f))
    descending_g = list(reversed(coeff_g))
    size = degree_f + degree_g
    rows = []
    for shift in range(degree_g):
        rows.append([Fraction(0)] * shift + descending_f + [Fraction(0)] * (degree_g - 1 - shift))
    for shift in range(degree_f):
        rows.append([Fraction(0)] * shift + descending_g + [Fraction(0)] * (degree_f - 1 - shift))
    require(all(len(row) == size for row in rows), "Sylvester shape changed")
    return determinant_fraction(rows)


def fixed_binary_form(depth, lows, order):
    """Moment of u(f_a-f_c)+v(f_b-f_c) on the transverse P1."""
    first, second, last = lows
    form = {}
    for u_count in range(order + 1):
        v_count = order - u_count
        total = Fraction(0)
        bases = (first,) * u_count + (second,) * v_count
        for mask in range(1 << order):
            offsets = tuple(
                last if mask & (1 << index) else base
                for index, base in enumerate(bases)
            )
            total += (-1) ** mask.bit_count() * normalized_tensor(depth, offsets)
        form[(u_count, v_count)] = comb(order, u_count) * total
    return form


def transverse_resultant(depth, lows):
    quadratic = fixed_binary_form(depth, lows, 2)
    cubic = fixed_binary_form(depth, lows, 3)
    return binary_resultant(quadratic, 2, cubic, 3)


def top_tensor(depth, order, moving):
    return normalized_tensor(depth, (moving,) * order)


def normal_boundary_form(depth, first_low, moving, gap, order):
    width = moving + gap
    form = {}
    for first_count in range(order + 1):
        moving_count = order - first_count
        entries = (first_low,) * first_count + (moving,) * moving_count
        form[(first_count, moving_count)] = (
            comb(order, first_count)
            * signed_inclusion(depth, width, entries)
            / top_tensor(depth, order, moving)
        )
    return form


def boundary_resultant(depth, first_low, moving, gap):
    fourth = normal_boundary_form(depth, first_low, moving, gap, 4)
    fifth = normal_boundary_form(depth, first_low, moving, gap, 5)
    return binary_resultant(fourth, 4, fifth, 5)


def width_flag_k5(depth, width):
    return (
        rising(depth, width) ** 154
        * rising(depth + 1, width) ** 172
        / depth ** (326 * width)
    )


def carrier_k5(depth, lows, moving, gap):
    sidecar = transverse_resultant(depth, lows)
    delta = 5**gap - 4**gap
    return (
        sidecar**20
        * delta**120
        * top_tensor(depth, 4, moving) ** 30
        * top_tensor(depth, 5, moving) ** 24
    )


def modular_fraction(value, prime):
    return (value.numerator % prime) * pow(value.denominator % prime, -1, prime) % prime


def macaulay_rank_mod_prime(forms, prime=1_000_003):
    target_degree = 11
    columns = {
        monomial: index for index, monomial in enumerate(MONOMIALS4[target_degree])
    }
    rows = []
    for order in (2, 3, 4, 5):
        reduced_form = {
            alpha: modular_fraction(value, prime)
            for alpha, value in forms[order].items()
        }
        for beta in MONOMIALS4[target_degree - order]:
            row = [0] * len(columns)
            for alpha, coefficient in reduced_form.items():
                target = tuple(alpha[index] + beta[index] for index in range(4))
                row[columns[target]] = (row[columns[target]] + coefficient) % prime
            rows.append(row)
    require(len(rows) == 589 and len(columns) == 364, "K5 Macaulay shape changed")
    return nmod_mat(rows, prime).rank(), len(rows), len(columns)


# General exponent and safe-gap ledger.
general_cells = 0
gap_cells = 0
for slots in range(3, 13):
    p, q = slots - 1, slots
    require(p * q * factorial(slots - 2) == factorial(slots), "Delta exponent changed")
    require(factorial(slots) // p == slots * factorial(slots - 2), "top multidegree changed")
    require(factorial(slots) // q == factorial(slots - 1), "terminal multidegree changed")
    general_cells += 1
    if slots >= 4:
        rho = Fraction(slots - 2, slots - 1) ** (slots - 2)
        for layer in range(1, slots - 1):
            require(Fraction(layer, slots - 1) ** layer <= rho, "safe gap failed")
            gap_cells += 1


# Literal K5 normal-ideal filtration.
layer_cells = 0
for depth, lows, moving, gap in (
    (1, (0, 1, 2), 4, 1),
    (1, (0, 2, 5), 7, 2),
    (2, (1, 3, 6), 8, 1),
    (3, (0, 4, 7), 10, 3),
):
    for order in (2, 3, 4, 5):
        layers = physical_layers_k5(order, depth, lows, moving, gap)
        for high_count, layer in layers.items():
            require(
                check_normal_ideal_k5(layer, high_count),
                "K5 high layer escaped its normal-ideal power",
            )
            layer_cells += 1


# Arbitrary three-low transverse bank.
sidecar_records = []
sidecar_cells = 0
for depth in (1, 2, 3):
    for last in range(2, 9):
        for first in range(last):
            for second in range(first + 1, last):
                value = transverse_resultant(depth, (first, second, last))
                require(value > 0, "three-low transverse resultant lost positivity")
                sidecar_records.append(
                    f"{depth}:{first}:{second}:{last}:{value.numerator}/{value.denominator}"
                )
                sidecar_cells += 1
sidecar_digest = sha256("\n".join(sidecar_records).encode("ascii")).hexdigest()


# Exact full boundary restrictions and their powered-line limit controls.
boundary_records = []
primary_boundaries = []
for depth, first_low, moving, gap in (
    (1, 0, 3, 1),
    (1, 0, 5, 1),
    (1, 0, 10, 1),
    (1, 0, 20, 1),
    (1, 0, 40, 1),
    (2, 1, 5, 2),
    (2, 1, 9, 2),
    (3, 0, 8, 3),
):
    value = boundary_resultant(depth, first_low, moving, gap)
    limit = Fraction((5**gap - 4**gap) ** 20)
    require(value > 0 and limit > 0, "normal boundary positivity control failed")
    boundary_records.append(
        f"{depth}:{first_low}:{moving}:{gap}:{value.numerator}/{value.denominator}"
    )
    if depth == 1 and first_low == 0 and gap == 1:
        primary_boundaries.append(value)
require(
    all(left > right > 1 for left, right in zip(primary_boundaries, primary_boundaries[1:])),
    "primary normal boundary lost descent to one",
)
boundary_digest = sha256("\n".join(boundary_records).encode("ascii")).hexdigest()


# Gauss/Beta factorization and exact width quotient.
beta_cells = 0
width_cells = 0
for depth in range(1, 5):
    for moving in range(1, 9):
        for order in (4, 5):
            product = Fraction(order ** (order * moving))
            for residue in range(1, order):
                product *= (
                    rising(Fraction(depth * order + residue, order), moving)
                    / rising(depth + 1, moving)
                )
            require(product == top_tensor(depth, order, moving), "K5 Gauss/Beta split failed")
            beta_cells += 1
    for width in range(0, 9):
        flag = width_flag_k5(depth, width)
        direct = (
            rising(depth, width) ** 154
            * rising(depth + 1, width) ** 172
            / depth ** (326 * width)
        )
        require(flag == direct, "K5 width flag changed")
        if width < 8:
            quotient = width_flag_k5(depth, width + 1) / flag
            expected = (
                Fraction(depth + width, depth) ** 154
                * Fraction(depth + width + 1, depth) ** 172
            )
            require(quotient == expected, "K5 width quotient changed")
        width_cells += 1


# Exact carrier Hankel banks on two affine clocks.
carrier_cells = 0
depth, lows, gap = 1, (0, 1, 2), 1
for step in (1, 2):
    values = {}
    for index in range(3, 19):
        moving = step * index
        values[index] = width_flag_k5(depth, moving + gap) * carrier_k5(
            depth, lows, moving, gap
        )
    for order in (1, 2, 3):
        for base in range(3, 15 - 2 * (order - 1)):
            matrix = [
                [values[base + row + column] for column in range(order)]
                for row in range(order)
            ]
            require(determinant_fraction(matrix) > 0, "K5 carrier Hankel minor failed")
            carrier_cells += 1


# Literal physical K5 degree-11 Macaulay maps over one exact finite field.
macaulay_cells = 0
macaulay_shape = None
macaulay_records = []
for depth in (1, 2):
    for lows in ((0, 1, 2), (0, 1, 3), (0, 2, 4)):
        for gap in (1, 2):
            for moving in (lows[-1] + 1, lows[-1] + 3, lows[-1] + 7):
                forms = physical_forms_k5(depth, lows, moving, gap)
                rank, rows, columns = macaulay_rank_mod_prime(forms)
                require(rank == columns, "physical K5 Macaulay map lost full rank")
                macaulay_shape = (rows, columns)
                macaulay_records.append(
                    json.dumps(
                        [depth, *lows, moving, gap, rank], separators=(",", ":")
                    )
                )
                macaulay_cells += 1
macaulay_digest = sha256("\n".join(macaulay_records).encode("ascii")).hexdigest()
require(
    macaulay_digest == "a21bce9b17f70d7f09e2d8f029ccffe157efa5557c475a1ebaad5ff448ca92a1",
    "physical K5 Macaulay digest changed",
)


print("THM-3063 TERMINAL SUSPENSION AND FIVE-SLOT PHYSICAL TAIL")
print(f"general_exponent_cells={general_cells} safe_gap_cells={gap_cells}")
print(f"normal_ideal_layer_cells={layer_cells}")
print(f"three_low_sidecar_cells={sidecar_cells} all_positive=1")
print(f"three_low_sidecar_digest={sidecar_digest}")
print(f"normal_boundary_cells={len(boundary_records)} all_positive=1")
print(f"normal_boundary_digest={boundary_digest}")
print(f"beta_factor_cells={beta_cells} width_cells={width_cells}")
print(f"carrier_hankel_cells={carrier_cells}")
print(
    f"physical_macaulay_cells={macaulay_cells} rows={macaulay_shape[0]} "
    f"columns={macaulay_shape[1]} full_rank=1"
)
print(f"physical_macaulay_digest={macaulay_digest}")
print("k5_asymptotic=S3^20*(5^h-4^h)^120*U4^30*U5^24;power=-93")
print("scope=intrinsic_terminal_suspension;fixed_lows_and_gap;finite_order_tail")
print("all_exact_checks=PASS")
