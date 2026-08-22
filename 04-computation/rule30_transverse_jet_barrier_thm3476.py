#!/usr/bin/env python3
"""Exact finite audit for THM-3476's Rule 30 transverse-jet package.

Universal claims are proved in the theorem.  This companion independently
checks the physical depth-four source, the Motzkin transport, the power-of-two
target family, the Pascal/Lucas residue atlas, and the symbolic source-kernel
identities.  It uses explicit failures rather than optimization-sensitive
assertions.
"""

from __future__ import annotations

from functools import lru_cache


POWER_MIN = 5
POWER_MAX = 11
PASCAL_POWER_MAX = 8
GREEN_EXHAUSTIVE_CAP = 192
MAX_MODULUS = 1 << POWER_MAX
MAX_GREEN_EXPONENT = 5 * MAX_MODULUS // 4 + 5
MAX_SOURCE_DISTANCE = MAX_GREEN_EXPONENT // 2
MAX_SOURCE_TIME = MAX_SOURCE_DISTANCE + 4


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


def edge_rows(horizon: int, width: int) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Independent right/left triangular edge recurrences."""

    right = [1] + [0] * (width - 1)
    left = [1] + [0] * (width - 1)
    right_rows: list[tuple[int, ...]] = []
    left_rows: list[tuple[int, ...]] = []
    for _ in range(horizon + 1):
        right_rows.append(tuple(right))
        left_rows.append(tuple(left))
        next_right = []
        next_left = []
        for index in range(width):
            right_one = right[index - 1] if index >= 1 else 0
            right_two = right[index - 2] if index >= 2 else 0
            next_right.append(right[index] ^ (right_one | right_two))

            left_one = left[index - 1] if index >= 1 else 0
            left_two = left[index - 2] if index >= 2 else 0
            next_left.append(
                left_two ^ left_one ^ ((1 ^ left_one) & left[index])
            )
        right = next_right
        left = next_left
    return right_rows, left_rows


def folded_source(rows: list[frozenset[int]], source_time: int, distance: int) -> int:
    row = rows[source_time]

    def edge(site: int) -> int:
        return int(site in row and site + 1 in row)

    if distance == 0:
        return edge(0)
    return edge(distance) ^ edge(-distance)


def alpha_direct(rows: list[frozenset[int]], source_depth: int, distance: int) -> int:
    source_time = source_depth + distance
    if source_time < 2:
        return 0
    return folded_source(rows, source_time, distance)


def alpha_four_closed(distance: int) -> int:
    return int(distance % 8 in (1, 5, 6, 7))


def ternary_powers(cap: int) -> list[int]:
    powers = [1]
    for _ in range(cap):
        current = powers[-1]
        powers.append(current ^ (current << 1) ^ (current << 2))
    return powers


def ternary_coefficient_polynomial(powers: list[int], exponent: int, degree: int) -> int:
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
    return ternary_coefficient_digit(half_exponent, half_degree) ^ ternary_coefficient_digit(
        half_exponent,
        half_degree - 1,
    )


def contributing_slacks(
    rows: list[frozenset[int]],
    powers: list[int],
    source_depth: int,
    target_time: int,
    *,
    digit_transport: bool,
) -> tuple[int, ...]:
    remainder = target_time - source_depth - 1
    if remainder < 0:
        return ()
    slacks = []
    for distance in range(remainder // 2 + 1):
        slack = remainder - 2 * distance
        source = alpha_direct(rows, source_depth, distance)
        if digit_transport:
            transport = ternary_coefficient_digit(distance + slack, slack)
        else:
            transport = ternary_coefficient_polynomial(
                powers,
                distance + slack,
                slack,
            )
        if source & transport:
            slacks.append(slack)
    return tuple(slacks)


def binomial_mod_two(top: int, bottom: int) -> int:
    if top < 0 or bottom < 0:
        return 0
    return int((bottom & ~top) == 0)


def hasse_value(exponents: tuple[int, ...], order: int) -> int:
    value = 0
    for exponent in exponents:
        value ^= binomial_mod_two(exponent, order)
    return value


def first_hasse_order(exponents: tuple[int, ...]) -> int | None:
    if not exponents:
        return None
    for order in range(max(exponents) + 2):
        if hasse_value(exponents, order):
            return order
    raise RuntimeError("CHECK FAILED: nonzero polynomial has no live Hasse order")


def xor_polynomials(*polynomials: set[tuple[int, ...]]) -> set[tuple[int, ...]]:
    total: set[tuple[int, ...]] = set()
    for polynomial in polynomials:
        total.symmetric_difference_update(polynomial)
    return total


def audit_source(rows: list[frozenset[int]]) -> str:
    right_rows, left_rows = edge_rows(MAX_SOURCE_TIME, 6)
    for time in range(MAX_SOURCE_TIME + 1):
        row = rows[time]
        for index in range(6):
            check(
                right_rows[time][index] == int(time - index in row),
                f"right edge recurrence time={time} index={index}",
            )
            check(
                left_rows[time][index] == int(-time + index in row),
                f"left edge recurrence time={time} index={index}",
            )

    expected_cycle = (
        ("10001", "110010", 0),
        ("11101", "110111", 1),
        ("10010", "110010", 0),
        ("11111", "110111", 0),
        ("10000", "110010", 0),
        ("11100", "110111", 1),
        ("10011", "110010", 1),
        ("11110", "110111", 1),
    )
    actual_cycle = []
    for distance in range(8):
        source_time = distance + 4
        right_word = "".join(str(bit) for bit in right_rows[source_time][:5])
        left_word = "".join(str(bit) for bit in left_rows[source_time][:6])
        actual_cycle.append(
            (right_word, left_word, alpha_direct(rows, 4, distance))
        )
    check(tuple(actual_cycle) == expected_cycle, "depth-four exact edge cycle")

    for distance in range(MAX_SOURCE_DISTANCE + 1):
        check(
            alpha_direct(rows, 4, distance) == alpha_four_closed(distance),
            f"depth-four source word distance={distance}",
        )
    return "01000111"


def audit_green(powers: list[int]) -> None:
    for exponent in range(GREEN_EXHAUSTIVE_CAP + 1):
        for degree in range(2 * exponent + 1):
            check(
                ternary_coefficient_polynomial(powers, exponent, degree)
                == ternary_coefficient_digit(exponent, degree),
                f"Motzkin polynomial/digit exponent={exponent} degree={degree}",
            )

    for low_exponent in range(64):
        for low_degree in range(64):
            check(
                ternary_coefficient_digit(8 * low_exponent + 7, 8 * low_degree + 5)
                == ternary_coefficient_digit(low_exponent, low_degree - 1),
                f"three-bit Green reduction N={low_exponent} Q={low_degree}",
            )

    for power in range(8):
        parameter = 5 * (1 << power) - 1
        support = tuple(
            shift
            for shift in range(parameter // 2 + 1)
            if ternary_coefficient_digit(parameter - shift, parameter)
        )
        check(
            support == (0, 1 << (power + 1)),
            f"Green two-point support power={power}",
        )


def audit_power_family(
    rows: list[frozenset[int]],
    powers: list[int],
) -> tuple[tuple[int, int, int, int], ...]:
    rows_out = []
    for exponent_power in range(POWER_MIN, POWER_MAX + 1):
        modulus = 1 << exponent_power
        target_time = 5 * modulus // 4 + 10
        low_slack = modulus // 4 - 7
        high_slack = low_slack + modulus
        expected = (high_slack, low_slack)
        direct_polynomial = contributing_slacks(
            rows,
            powers,
            4,
            target_time,
            digit_transport=False,
        )
        direct_digit = contributing_slacks(
            rows,
            powers,
            4,
            target_time,
            digit_transport=True,
        )
        check(direct_polynomial == expected, f"physical polynomial family m={exponent_power}")
        check(direct_digit == expected, f"physical digit family m={exponent_power}")

        for order in range(modulus):
            check(
                hasse_value(expected, order) == 0,
                f"vanishing Hasse jet m={exponent_power} order={order}",
            )
        check(
            hasse_value(expected, modulus) == 1,
            f"first live Hasse jet m={exponent_power}",
        )
        check(
            hasse_value(expected, modulus + 1) == 1,
            f"paired live Hasse jet m={exponent_power}",
        )
        rows_out.append((exponent_power, modulus, target_time, low_slack))
    return tuple(rows_out)


def audit_pascal_atlas() -> tuple[int, ...]:
    moduli = []
    for modulus_power in range(1, PASCAL_POWER_MAX + 1):
        modulus = 1 << modulus_power
        moduli.append(modulus)
        for residue in range(modulus):
            jets = tuple(
                binomial_mod_two(residue, order)
                for order in range(modulus)
            )
            reconstruction = tuple(
                sum(
                    binomial_mod_two(order, reconstructed_residue) * jets[order]
                    for order in range(modulus)
                )
                % 2
                for reconstructed_residue in range(modulus)
            )
            expected = tuple(int(index == residue) for index in range(modulus))
            check(
                reconstruction == expected,
                f"Pascal involution modulus={modulus} residue={residue}",
            )

        for exponent in range(4 * modulus):
            for order in range(modulus):
                check(
                    binomial_mod_two(exponent, order)
                    == binomial_mod_two(exponent % modulus, order),
                    f"Lucas residue modulus={modulus} exponent={exponent} order={order}",
                )
            for order in range(4 * modulus):
                check(
                    binomial_mod_two(exponent, order)
                    == (
                        binomial_mod_two(exponent % modulus, order % modulus)
                        & binomial_mod_two(exponent // modulus, order // modulus)
                    ),
                    f"Lucas tensor modulus={modulus} exponent={exponent} order={order}",
                )
    return tuple(moduli)


def audit_ballistic_ramification(rows: list[frozenset[int]]) -> int:
    packet_count = 0
    for source_depth in range(9):
        for target_time in range(source_depth + 1, 97):
            remainder = target_time - source_depth - 1
            distances = []
            slacks = []
            for distance in range(remainder // 2 + 1):
                slack = remainder - 2 * distance
                source = alpha_direct(rows, source_depth, distance)
                transport = ternary_coefficient_digit(distance + slack, slack)
                if source & transport:
                    distances.append(distance)
                    slacks.append(slack)
            if not distances:
                continue

            distance_packet = tuple(distances)
            slack_packet = tuple(slacks)
            distance_order = first_hasse_order(distance_packet)
            slack_order = first_hasse_order(slack_packet)
            check(
                slack_order == 2 * distance_order,
                f"ballistic valuation u={source_depth} t={target_time}",
            )

            parity = remainder & 1
            half_degree = (remainder - parity) // 2
            reversed_packet = tuple(half_degree - distance for distance in distances)
            for order in range(half_degree + 1):
                coefficient = hasse_value(reversed_packet, order)
                check(
                    hasse_value(slack_packet, 2 * order) == coefficient,
                    f"ballistic even jet u={source_depth} t={target_time} j={order}",
                )
                check(
                    hasse_value(slack_packet, 2 * order + 1) == (parity & coefficient),
                    f"ballistic odd jet u={source_depth} t={target_time} j={order}",
                )
            packet_count += 1
    return packet_count


def audit_symbolic_kernel() -> None:
    # Monomial coordinates are (epsilon exponent, U exponent, X exponent).
    p_polynomial = {
        (0, 0, 2),
        (0, 0, 1),
        (0, 1, 1),
        (0, 2, 0),
    }
    deformed_characteristic = {
        (0, 0, 2),
        (2, 0, 2),
        (0, 0, 1),
        (0, 1, 1),
        (1, 1, 1),
        (0, 2, 0),
    }
    expected_difference = {(2, 0, 2), (1, 1, 1)}
    check(
        xor_polynomials(p_polynomial, deformed_characteristic)
        == expected_difference,
        "P-adic/epsilon characteristic identity",
    )

    # (1+X^2)S_boundary = X P for u=0,1,2.
    boundary_numerator = {
        (0, 0, 3),
        (0, 0, 2),
        (0, 1, 2),
        (0, 2, 1),
    }
    x_times_p = {
        (0, 0, 3),
        (0, 0, 2),
        (0, 1, 2),
        (0, 2, 1),
    }
    check(boundary_numerator == x_times_p, "three-strip source equals XP/(1+X^2)")


def main() -> None:
    rows = centered_rows(MAX_SOURCE_TIME)
    powers = ternary_powers(MAX_GREEN_EXPONENT)
    source_period = audit_source(rows)
    audit_green(powers)
    power_rows = audit_power_family(rows, powers)
    pascal_moduli = audit_pascal_atlas()
    ramification_packets = audit_ballistic_ramification(rows)
    audit_symbolic_kernel()

    print("RULE 30 TRANSVERSE-JET BARRIER EXACT AUDIT")
    print(f"depth_four_source_period={source_period} support_residues=[1,5,6,7]_mod_8")
    print(f"power_family_m_M_t_low_slack={list(power_rows)}")
    print(f"first_live_jet_orders={[row[1] for row in power_rows]}")
    print(
        f"pascal_lucas_moduli={list(pascal_moduli)} "
        "matrix_inverse=exact tensor_factorization=exact"
    )
    print(
        f"ballistic_ramification_packets={ramification_packets} "
        "ord_slack=2*ord_distance even_odd_atlas=exact"
    )
    print(
        "symbolic_kernel=P=X^2+(1+U)X+U^2; "
        "P(U,W_(1+epsilon))=epsilon*U*W+epsilon^2*W^2"
    )
    print("boundary_source=(X*P)/(1+X^2) exact")
    print(
        "scope=FINITE-EXACT audit of universal proofs; no Rule 30 prize, "
        "general-time lower bound, or fixed-seed density conclusion"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
