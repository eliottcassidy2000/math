#!/usr/bin/env python3
"""Exact finite companion for provisional THM-3488.

The universal arguments are in the theorem.  This script independently
checks the physical Rule 30 source, two Green engines, ballistic faces,
targetwise monicity, the parity-Cartier atlas, and the 64-class shallow
source-Taylor certificate.  Every gate is explicit and remains active under
``python -O``.
"""

from __future__ import annotations

from functools import lru_cache


ROW_HORIZON = 320
TARGET_CAP = 192
GREEN_CROSSCHECK_CAP = 160
GREEN_CENTRAL_CAP = 512
GLOBAL_JET_CAP = 96
CARTIER_TARGET_CAP = 128
TAYLOR_K_CAP = 255
Z_TRUNCATION = 24


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def centered_rows(horizon: int) -> list[frozenset[int]]:
    rows: list[frozenset[int]] = []
    row = frozenset({0})
    for time in range(horizon + 1):
        rows.append(row)
        row = frozenset(
            site
            for site in range(-time - 1, time + 2)
            if rule30(
                int(site - 1 in row),
                int(site in row),
                int(site + 1 in row),
            )
        )
    return rows


def edge_step(
    right: tuple[int, ...],
    left: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    next_right = []
    next_left = []
    for index in range(len(right)):
        right_one = right[index - 1] if index >= 1 else 0
        right_two = right[index - 2] if index >= 2 else 0
        next_right.append(right[index] ^ (right_one | right_two))
    for index in range(len(left)):
        left_one = left[index - 1] if index >= 1 else 0
        left_two = left[index - 2] if index >= 2 else 0
        next_left.append(
            left_two ^ left_one ^ ((1 ^ left_one) & left[index])
        )
    return tuple(next_right), tuple(next_left)


def edge_states(
    source_depth: int,
    count: int,
) -> tuple[list[int], list[tuple[tuple[int, ...], tuple[int, ...]]]]:
    """Return alpha_u(d) and the closed prefix state at each distance."""

    right = (1,) + (0,) * source_depth
    left = (1,) + (0,) * (source_depth + 1)
    for _ in range(source_depth):
        right, left = edge_step(right, left)

    word: list[int] = []
    states = []
    for distance in range(count):
        states.append((right, left))
        source = right[source_depth] & right[source_depth - 1]
        if distance:
            source ^= left[source_depth] & left[source_depth + 1]
        word.append(source)
        right, left = edge_step(right, left)
    return word, states


def folded_source(
    rows: list[frozenset[int]],
    source_depth: int,
    distance: int,
) -> int:
    source_time = source_depth + distance
    if source_time < 2:
        return 0
    row = rows[source_time]

    def edge(site: int) -> int:
        return int(site in row and site + 1 in row)

    if distance == 0:
        return edge(0)
    return edge(distance) ^ edge(-distance)


SOURCE_CYCLES = {
    3: ("1", "0100"),
    4: ("", "01000111"),
    5: ("1", "00010100"),
    6: ("", "0000000001000101"),
    7: ("", "01010101000010000010000010000010"),
    8: ("1", "11011100110110110011001110001000"),
    9: ("", "0111101001100011101110110010011000110111001110010001000100110001"),
    10: ("", "0001010010000001001000010100101001000000000001001000010000010101"),
    11: ("", "0101101100010001010000010000000110111001000110010011100100011001"),
    12: ("", "0111011011001100111001100110010000100110011001100111011001000110"),
}


def source_from_cycle(source_depth: int, distance: int) -> int:
    head, cycle = SOURCE_CYCLES[source_depth]
    if distance < len(head):
        return int(head[distance])
    return int(cycle[(distance - len(head)) % len(cycle)])


def ternary_powers(cap: int) -> list[int]:
    powers = [1]
    for _ in range(cap):
        current = powers[-1]
        powers.append(current ^ (current << 1) ^ (current << 2))
    return powers


def ternary_coefficient_polynomial(
    powers: list[int], exponent: int, degree: int
) -> int:
    if exponent < 0 or degree < 0 or degree > 2 * exponent:
        return 0
    return (powers[exponent] >> degree) & 1


@lru_cache(maxsize=None)
def ternary_coefficient_digit(exponent: int, degree: int) -> int:
    if exponent < 0 or degree < 0 or degree > 2 * exponent:
        return 0
    if exponent == 0:
        return int(degree == 0)
    half_exponent = exponent // 2
    half_degree = degree // 2
    if exponent % 2 == 0:
        if degree % 2:
            return 0
        return ternary_coefficient_digit(half_exponent, half_degree)
    if degree % 2:
        return ternary_coefficient_digit(half_exponent, half_degree)
    return ternary_coefficient_digit(
        half_exponent, half_degree
    ) ^ ternary_coefficient_digit(half_exponent, half_degree - 1)


def binomial_mod_two(top: int, bottom: int) -> int:
    if top < 0 or bottom < 0:
        return 0
    return int((bottom & ~top) == 0)


def hasse_value(polynomial: int, order: int) -> int:
    value = 0
    remaining = polynomial
    while remaining:
        bit = remaining & -remaining
        exponent = bit.bit_length() - 1
        value ^= binomial_mod_two(exponent, order)
        remaining ^= bit
    return value


def first_hasse_order(polynomial: int) -> int | None:
    if polynomial == 0:
        return None
    for order in range(polynomial.bit_length()):
        if hasse_value(polynomial, order):
            return order
    raise RuntimeError("CHECK FAILED: nonzero polynomial has no Hasse order")


def polynomial_exponents(polynomial: int) -> tuple[int, ...]:
    exponents = []
    remaining = polynomial
    while remaining:
        bit = remaining & -remaining
        exponents.append(bit.bit_length() - 1)
        remaining ^= bit
    return tuple(exponents)


def tail_target_data(
    rows: list[frozenset[int]], target: int
) -> tuple[int, tuple[int, int]]:
    """Return C_t(q) and the two reversed parity packets."""

    target_polynomial = 0
    channels = [0, 0]
    for source_depth in range(3, target):
        remainder = target - source_depth - 1
        parity = remainder & 1
        half_remainder = remainder // 2
        for distance in range(half_remainder + 1):
            slack = remainder - 2 * distance
            source = folded_source(rows, source_depth, distance)
            transport = ternary_coefficient_digit(distance + slack, slack)
            if source & transport:
                target_polynomial ^= 1 << slack
                channels[parity] ^= 1 << (half_remainder - distance)
    return target_polynomial, (channels[0], channels[1])


def parity_dilate(polynomial: int, parity: int) -> int:
    result = 0
    for exponent in polynomial_exponents(polynomial):
        result ^= 1 << (2 * exponent + parity)
    return result


def characteristic_root(truncation: int) -> int:
    """Small W satisfying W^2+(1+z)W+z^2=0 modulo z^truncation."""

    coefficients = [0] * truncation
    for degree in range(1, truncation):
        coefficients[degree] = (
            coefficients[degree - 1]
            ^ (coefficients[degree // 2] if degree % 2 == 0 else 0)
            ^ int(degree == 2)
        )
    return sum(coefficient << degree for degree, coefficient in enumerate(coefficients))


def series_multiply(left: int, right: int, truncation: int) -> int:
    result = 0
    remaining = right
    mask = (1 << truncation) - 1
    while remaining:
        bit = remaining & -remaining
        shift = bit.bit_length() - 1
        result ^= left << shift
        remaining ^= bit
    return result & mask


def expected_taylor_owner(order: int) -> tuple[int, tuple[int, int]]:
    if order == 0:
        return 3, (3, 0)
    residue = order % 64
    if residue % 4 == 2:
        return 3, (3, 0)
    if residue % 8 in (1, 5, 7):
        return 4, (4, 0)
    if residue % 8 == 4:
        return 5, (5, 0)
    if residue % 8 == 0:
        return 6, (4, 1)
    if residue in (3, 35):
        return 7, (7, 0)
    if residue in (19, 51):
        return 8, (8, 0)
    if residue in (43, 59):
        return 9, (9, 0)
    if residue == 11:
        return 11, (11, 0)
    if residue == 27:
        return 12, (4, 4)
    raise RuntimeError(f"CHECK FAILED: uncovered Taylor residue {residue}")


def audit_source_cycles(rows: list[frozenset[int]]) -> None:
    for source_depth, (head, cycle) in SOURCE_CYCLES.items():
        preperiod = len(head)
        period = len(cycle)
        count = max(preperiod + 2 * period + 1, 161)
        word, states = edge_states(source_depth, count)
        expected = tuple(int(bit) for bit in head + cycle + cycle)
        check(
            tuple(word[: len(expected)]) == expected,
            f"source word u={source_depth}",
        )
        check(
            states[preperiod] == states[preperiod + period],
            f"complete edge-state return u={source_depth}",
        )
        for distance in range(min(160, ROW_HORIZON - source_depth + 1)):
            source_time = source_depth + distance
            check(
                folded_source(rows, source_depth, distance)
                == source_from_cycle(source_depth, distance),
                f"centered/edge source u={source_depth} d={distance}",
            )
            right, left = states[distance]
            row = rows[source_time]
            for index, bit in enumerate(right):
                check(
                    bit == int(source_time - index in row),
                    f"right edge prefix u={source_depth} d={distance} r={index}",
                )
            for index, bit in enumerate(left):
                check(
                    bit == int(-source_time + index in row),
                    f"left edge prefix u={source_depth} d={distance} r={index}",
                )

    check(folded_source(rows, 3, 0) == 1, "alpha_3(0)=1")
    check(folded_source(rows, 4, 0) == 0, "alpha_4(0)=0")
    check(folded_source(rows, 5, 0) == 1, "alpha_5(0)=1")
    check(folded_source(rows, 3, 1) == 0, "alpha_3(1)=0")


def audit_green(powers: list[int]) -> None:
    for exponent in range(GREEN_CROSSCHECK_CAP + 1):
        for degree in range(2 * exponent + 1):
            check(
                ternary_coefficient_polynomial(powers, exponent, degree)
                == ternary_coefficient_digit(exponent, degree),
                f"Green engines n={exponent} r={degree}",
            )
    for exponent in range(GREEN_CENTRAL_CAP + 1):
        check(
            ternary_coefficient_polynomial(powers, exponent, exponent) == 1,
            f"central Green coefficient n={exponent}",
        )


def audit_targets(rows: list[frozenset[int]]) -> list[tuple[int, tuple[int, int]]]:
    data: list[tuple[int, tuple[int, int]]] = [(0, (0, 0))] * (TARGET_CAP + 1)
    for target in range(4, TARGET_CAP + 1):
        polynomial, channels = tail_target_data(rows, target)
        data[target] = (polynomial, channels)
        check(polynomial != 0, f"nonzero marked tail target={target}")
        check(
            polynomial.bit_length() - 1 == target - 4,
            f"monic degree target={target}",
        )
        check((polynomial >> (target - 4)) & 1, f"top face target={target}")
        if target >= 5:
            check(
                ((polynomial >> (target - 5)) & 1) == 0,
                f"empty weight-five face target={target}",
            )
        if target >= 6:
            check(
                (polynomial >> (target - 6)) & 1,
                f"weight-six face target={target}",
            )

    for order in range(GLOBAL_JET_CAP + 1):
        first_target = None
        for target in range(4, TARGET_CAP + 1):
            if hasse_value(data[target][0], order):
                first_target = target
                break
        check(first_target == order + 4, f"global jet onset order={order}")
        check(
            hasse_value(data[order + 4][0], order) == 1,
            f"global jet leading coefficient order={order}",
        )
    return data


def audit_cartier(data: list[tuple[int, tuple[int, int]]]) -> None:
    valuation_cases: set[str] = set()
    for target in range(4, CARTIER_TARGET_CAP + 1):
        polynomial, channels = data[target]
        reconstructed = parity_dilate(channels[0], 0) ^ parity_dilate(channels[1], 1)
        check(reconstructed == polynomial, f"Cartier reconstruction target={target}")

        max_channel_degree = max(channel.bit_length() for channel in channels)
        for order in range(max_channel_degree + 1):
            coefficient_zero = hasse_value(channels[0], order)
            coefficient_one = hasse_value(channels[1], order)
            check(
                hasse_value(polynomial, 2 * order)
                == (coefficient_zero ^ coefficient_one),
                f"Cartier even Hasse target={target} order={order}",
            )
            check(
                hasse_value(polynomial, 2 * order + 1) == coefficient_one,
                f"Cartier odd Hasse target={target} order={order}",
            )

        orders = (first_hasse_order(channels[0]), first_hasse_order(channels[1]))
        actual = first_hasse_order(polynomial)
        if orders == (None, None):
            predicted = None
            valuation_cases.add("zero")
        elif orders[0] is None:
            predicted = 2 * orders[1]
            valuation_cases.add("unequal")
        elif orders[1] is None:
            predicted = 2 * orders[0]
            valuation_cases.add("unequal")
        elif orders[0] != orders[1]:
            predicted = 2 * min(orders)
            valuation_cases.add("unequal")
        else:
            predicted = 2 * orders[0] + 1
            valuation_cases.add("equal")
        check(actual == predicted, f"Cartier valuation target={target}")

    check(data[6][0] == ((1 << 0) | (1 << 2)), "target-six hostile")
    check(data[9][1] == (1, 0b111), "target-nine channel hostile")
    check(
        data[9][0] == ((1 << 0) | (1 << 1) | (1 << 3) | (1 << 5)),
        "target-nine marked hostile",
    )
    check(valuation_cases >= {"equal", "unequal"}, "both Cartier valuation cases")


def audit_taylor(rows: list[frozenset[int]]) -> None:
    for order in range(TAYLOR_K_CAP + 1):
        best_weight = 10**9
        owners: list[tuple[int, int]] = []
        for source_depth in range(3, 13):
            for offset in range(5):
                weight = source_depth + 2 * offset
                if weight > 12:
                    continue
                distance = order + offset
                if (
                    source_from_cycle(source_depth, distance)
                    and binomial_mod_two(distance, order)
                ):
                    if weight < best_weight:
                        best_weight = weight
                        owners = [(source_depth, offset)]
                    elif weight == best_weight:
                        owners.append((source_depth, offset))
        expected_weight, expected_owner = expected_taylor_owner(order)
        check(best_weight == expected_weight, f"Taylor weight order={order}")
        check(owners == [expected_owner], f"Taylor unique owner order={order}")

    root = characteristic_root(Z_TRUNCATION)
    characteristic = (
        series_multiply(root, root, Z_TRUNCATION)
        ^ root
        ^ ((root << 1) & ((1 << Z_TRUNCATION) - 1))
        ^ (1 << 2)
    )
    check(characteristic == 0, "truncated small-root characteristic")
    root_powers = [1]
    for _ in range((Z_TRUNCATION - 3) // 2 + 1):
        root_powers.append(
            series_multiply(root_powers[-1], root, Z_TRUNCATION)
        )
    mask = (1 << Z_TRUNCATION) - 1
    for order in range(TAYLOR_K_CAP + 1):
        coordinate = 0
        for source_depth in range(3, Z_TRUNCATION):
            max_offset = (Z_TRUNCATION - 1 - source_depth) // 2
            for offset in range(max_offset + 1):
                distance = order + offset
                if (
                    folded_source(rows, source_depth, distance)
                    and binomial_mod_two(distance, order)
                ):
                    coordinate ^= (root_powers[offset] << source_depth) & mask
        first_term = (coordinate & -coordinate).bit_length() - 1 if coordinate else None
        expected_weight, _ = expected_taylor_owner(order)
        check(first_term == expected_weight, f"direct Taylor coordinate order={order}")


def main() -> None:
    rows = centered_rows(ROW_HORIZON)
    powers = ternary_powers(GREEN_CENTRAL_CAP)
    audit_source_cycles(rows)
    audit_green(powers)
    data = audit_targets(rows)
    audit_cartier(data)
    audit_taylor(rows)

    print("RULE 30 INWARD-SLACK CARTIER EXACT AUDIT")
    print(
        "ballistic_faces=omega4:z^4/(1+zq) omega5:0 "
        "omega6:z^6/(1+zq)"
    )
    print(
        f"tail_monic_targets=4..{TARGET_CAP} degree=t-4 "
        "top3_for_t>=6=[1,0,1]"
    )
    print(
        f"global_jet_staircase=orders_0..{GLOBAL_JET_CAP} "
        "ord_z(J_j)=j+4 leading=1"
    )
    print(
        f"parity_cartier_targets=4..{CARTIER_TARGET_CAP} "
        "J_2j=h0_j+h1_j J_(2j+1)=h1_j valuation=exact"
    )
    print(
        "source_cycles_u3..12="
        "preperiods_[1,0,1,0,0,1,0,0,0,0]_"
        "periods_[4,8,8,16,32,32,64,64,64,64] state_return=exact"
    )
    print(
        f"source_taylor_orders=0..{TAYLOR_K_CAP} residue_classes=64 "
        "ord_z_range=[3,12] unique_shallow_owner=exact"
    )
    print(
        f"green_crosscheck_exponents=0..{GREEN_CROSSCHECK_CAP} "
        f"central_coefficients=0..{GREEN_CENTRAL_CAP} exact"
    )
    print("hostiles=t6:1+q^2 t9:1+q+q^3+q^5 scalar_cancellation=exact")
    print(
        "scope=FINITE-EXACT audit of universal proofs; no Rule 30 prize, "
        "density theorem, or unrestricted complexity lower bound"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
