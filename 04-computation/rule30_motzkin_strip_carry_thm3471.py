#!/usr/bin/env python3
"""Exact verifier for THM-3471's Rule 30 strip/carry package.

The universal statements belong to the theorem.  This companion keeps the
finite checks typed and separate:

* centered Rule 30 rows supply the first three folded source strips;
* polynomial multiplication and an independent binary-digit recurrence
  supply the ternary/Motzkin Green coefficients;
* the three strip identities and the transport-slack first jet are checked
  without substituting a precomputed center word;
* the innovation atlas is compiled into a cyclic F_2 cochain problem;
* innovative carries are checked for odd support, full ANF degree, and full
  Walsh support (outside the one-variable boundary case); and
* the distinguished center is independently recovered as a moving marked
  arc of the transition current and as a physical left-front ray.

All checks are finite exact certificates.  No Rule 30 prize conclusion is
extrapolated from them.
"""

from __future__ import annotations

from functools import lru_cache


STRIP_TARGET_CAP = 512
STRIP_SOURCE_CAP = (STRIP_TARGET_CAP + 3) // 2 + 2
GREEN_RECURRENCE_CAP = STRIP_TARGET_CAP
INNOVATION_DEPTH_CAP = 30
MARKED_ARC_CAP = 30

EXPECTED_LIFT_WORD = "101101101000001100000001101010"
EXPECTED_INNOVATIONS = (1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29)
EXPECTED_CARRY_WEIGHTS = (1, 1, 3, 5, 11, 21, 41, 93, 187, 381, 757, 1533)
EXPECTED_MARKED_INNOVATION_WORD = "111001111001"


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def centered_rows(horizon: int) -> list[frozenset[int]]:
    """Direct centered-coordinate Rule 30, independent of packed edges."""

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


def folded_source(rows: list[frozenset[int]], source_time: int, distance: int) -> int:
    """F_s(d), retaining the oriented half-integer bond convention."""

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


def alpha_closed(source_depth: int, distance: int) -> int:
    if source_depth == 0:
        return int(distance >= 2)
    if source_depth == 1:
        return int(distance >= 2 and distance % 2 == 0)
    if source_depth == 2:
        return int(distance >= 1 and distance % 2 == 1)
    raise ValueError("closed alpha is only defined for source depths 0,1,2")


def ternary_powers(cap: int) -> list[int]:
    """Bitsets for (1+x+x^2)^n over F_2."""

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
    """Independent Frobenius/digit recurrence for [x^degree](1+x+x^2)^exponent."""

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


def green_shell_polynomial(powers: list[int], distance: int, time: int) -> int:
    """H_distance(time), using the palindromic ternary chart."""

    if time < 0 or distance < 0 or distance > time:
        return 0
    return ternary_coefficient_polynomial(powers, time, time + distance)


def strip_value(
    rows: list[frozenset[int]],
    powers: list[int],
    source_depth: int,
    target_time: int,
    *,
    use_closed_source: bool,
    use_digit_green: bool,
) -> int:
    """R_u(t)=sum alpha_u(d) K(d+v,v), t=u+1+2d+v."""

    total = 0
    remainder = target_time - source_depth - 1
    if remainder < 0:
        return 0
    for distance in range(remainder // 2 + 1):
        slack = remainder - 2 * distance
        if use_closed_source:
            source = alpha_closed(source_depth, distance)
        else:
            source = alpha_direct(rows, source_depth, distance)
        if use_digit_green:
            transport = ternary_coefficient_digit(distance + slack, slack)
        else:
            transport = ternary_coefficient_polynomial(powers, distance + slack, slack)
        total ^= source & transport
    return total


def first_transport_jet(
    rows: list[frozenset[int]],
    powers: list[int],
    target_time: int,
) -> int:
    """The q=1 first jet: weight every strip term by v mod 2."""

    total = 0
    for source_depth in range(3):
        remainder = target_time - source_depth - 1
        if remainder < 0:
            continue
        for distance in range(remainder // 2 + 1):
            slack = remainder - 2 * distance
            total ^= (
                (slack & 1)
                & alpha_direct(rows, source_depth, distance)
                & ternary_coefficient_polynomial(powers, distance + slack, slack)
            )
    return total


def audit_motzkin_strips() -> tuple[str, str, str, str, str]:
    rows = centered_rows(STRIP_SOURCE_CAP)
    powers = ternary_powers(GREEN_RECURRENCE_CAP)

    for exponent in range(GREEN_RECURRENCE_CAP + 1):
        for degree in range(2 * exponent + 1):
            check(
                ternary_coefficient_polynomial(powers, exponent, degree)
                == ternary_coefficient_digit(exponent, degree),
                f"ternary polynomial/digit n={exponent} q={degree}",
            )

    for source_depth in range(3):
        for distance in range(STRIP_SOURCE_CAP - source_depth + 1):
            check(
                alpha_direct(rows, source_depth, distance)
                == alpha_closed(source_depth, distance),
                f"folded source alpha u={source_depth} d={distance}",
            )

    strip_words = [[], [], []]
    backbone_word = []
    jet_word = []
    for target_time in range(2, STRIP_TARGET_CAP + 1):
        strips = []
        for source_depth in range(3):
            direct_polynomial = strip_value(
                rows,
                powers,
                source_depth,
                target_time,
                use_closed_source=False,
                use_digit_green=False,
            )
            closed_digit = strip_value(
                rows,
                powers,
                source_depth,
                target_time,
                use_closed_source=True,
                use_digit_green=True,
            )
            check(
                direct_polynomial == closed_digit,
                f"strip independent paths u={source_depth} t={target_time}",
            )
            strips.append(direct_polynomial)
            strip_words[source_depth].append(str(direct_polynomial))

        backbone = green_shell_polynomial(powers, 1, target_time - 2)
        check(
            strips[0] ^ strips[1] ^ strips[2] == 0,
            f"three-strip cancellation t={target_time}",
        )
        check(
            (backbone ^ strips[2]) == int(target_time >= 3 and target_time % 2 == 1),
            f"backbone plus strip two alternating t={target_time}",
        )
        jet = first_transport_jet(rows, powers, target_time)
        check(jet == strips[1], f"transport first jet equals strip one t={target_time}")
        backbone_word.append(str(backbone))
        jet_word.append(str(jet))

    prefix_length = 64
    return (
        "".join(backbone_word[:prefix_length]),
        "".join(strip_words[0][:prefix_length]),
        "".join(strip_words[1][:prefix_length]),
        "".join(strip_words[2][:prefix_length]),
        "".join(jet_word[:prefix_length]),
    )


def packed_phi_width(row: int, width: int) -> int:
    return (row ^ ((row << 1) | (row << 2))) & ((1 << width) - 1)


def seed_period(width: int) -> int:
    row = 1
    period = 0
    while True:
        row = packed_phi_width(row, width)
        period += 1
        if row == 1:
            return period


def period_lift_bit(width: int, period: int) -> int:
    row = 1
    for _ in range(period):
        row = packed_phi_width(row, width + 1)
    return (row >> width) & 1


def packed_edge_cycle(width: int, period: int) -> list[int]:
    rows = []
    row = 1
    for _ in range(period):
        rows.append(row)
        row = packed_phi_width(row, width)
    check(row == 1, f"packed edge full cycle width={width}")
    return rows


def edge_bit(rows: list[int], periods: list[int], depth: int, time: int) -> int:
    return (rows[time % periods[depth]] >> depth) & 1


def arrival_readout(rows: list[int], periods: list[int], depth: int, phase: int) -> int:
    return edge_bit(rows, periods, depth, phase + depth)


def transition_current(rows: list[int], periods: list[int], depth: int, phase: int) -> int:
    check(depth >= 2, "transition current depth at least two")
    return arrival_readout(rows, periods, depth - 1, phase + 1) | arrival_readout(
        rows,
        periods,
        depth - 2,
        phase + 2,
    )


def mobius_anf(values: bytes) -> bytearray:
    size = len(values)
    variables = size.bit_length() - 1
    check(size == 1 << variables, "ANF table size is power of two")
    coefficients = bytearray(values)
    for bit in range(variables):
        pivot = 1 << bit
        for mask in range(size):
            if mask & pivot:
                coefficients[mask] ^= coefficients[mask ^ pivot]
    return coefficients


def walsh_spectrum(values: bytes) -> list[int]:
    spectrum = [1 - 2 * value for value in values]
    stride = 1
    while stride < len(spectrum):
        for base in range(0, len(spectrum), 2 * stride):
            for offset in range(stride):
                left = spectrum[base + offset]
                right = spectrum[base + stride + offset]
                spectrum[base + offset] = left + right
                spectrum[base + stride + offset] = left - right
        stride *= 2
    return spectrum


def phase_masks(
    rows: list[int],
    periods: list[int],
    innovation_depths: list[int],
    period: int,
) -> tuple[list[int], list[int]]:
    masks = []
    inverse = [-1] * period
    for phase in range(period):
        mask = sum(
            arrival_readout(rows, periods, depth, phase) << index
            for index, depth in enumerate(innovation_depths)
        )
        check(0 <= mask < period, f"innovation mask range phase={phase}")
        check(inverse[mask] == -1, f"innovation atlas collision phase={phase}")
        inverse[mask] = phase
        masks.append(mask)
    check(all(phase >= 0 for phase in inverse), "innovation atlas is surjective")
    return masks, inverse


def audit_innovation_cochains() -> tuple[str, list[tuple[int, int, int, int, int]], list[int]]:
    periods = [seed_period(width) for width in range(1, INNOVATION_DEPTH_CAP + 2)]
    lift_bits = [
        period_lift_bit(width, periods[width - 1])
        for width in range(1, INNOVATION_DEPTH_CAP + 1)
    ]
    lift_word = "".join(map(str, lift_bits))
    check(lift_word == EXPECTED_LIFT_WORD, "frozen innovation lift word")

    full_width = INNOVATION_DEPTH_CAP + 1
    full_period = periods[INNOVATION_DEPTH_CAP]
    rows = packed_edge_cycle(full_width, full_period)

    innovation_depths = [
        depth
        for depth in range(1, INNOVATION_DEPTH_CAP + 1)
        if lift_bits[depth - 1]
    ]
    check(tuple(innovation_depths) == EXPECTED_INNOVATIONS, "frozen innovation depths")

    carry_rows: list[tuple[int, int, int, int, int]] = []
    observed_carry_weights = [1]

    # The first carry has no lower variables and toggles every phase step.
    check(arrival_readout(rows, periods, 1, 1) ^ arrival_readout(rows, periods, 1, 0) == 1,
          "zero-variable first carry")
    carry_rows.append((1, 0, 1, 0, 0))

    for depth in range(2, INNOVATION_DEPTH_CAP + 1):
        lower_innovations = [item for item in innovation_depths if item < depth]
        lower_period = periods[depth - 1]
        check(lower_period == 1 << len(lower_innovations), f"lower cube size k={depth}")
        masks, inverse = phase_masks(rows, periods, lower_innovations, lower_period)

        current_values = bytearray(lower_period)
        readout_values = bytearray(lower_period)
        for phase in range(lower_period):
            mask = masks[phase]
            current = transition_current(rows, periods, depth, phase)
            current_values[mask] = current
            readout_values[mask] = arrival_readout(rows, periods, depth, phase)
            check(
                arrival_readout(rows, periods, depth, phase + 1)
                ^ arrival_readout(rows, periods, depth, phase)
                == current,
                f"transition current identity k={depth} h={phase}",
            )

        current_parity = sum(current_values) & 1
        epsilon = lift_bits[depth - 1]
        check(current_parity == epsilon, f"cochain holonomy k={depth}")

        tau = [masks[(inverse[mask] + 1) % lower_period] for mask in range(lower_period)]
        seen = set()
        cursor = 0
        for _ in range(lower_period):
            check(cursor not in seen, f"lower odometer early return k={depth}")
            seen.add(cursor)
            cursor = tau[cursor]
        check(cursor == 0 and len(seen) == lower_period, f"lower odometer one cycle k={depth}")

        if epsilon == 0:
            for mask in range(lower_period):
                check(
                    readout_values[tau[mask]] ^ readout_values[mask] == current_values[mask],
                    f"noninnovation primitive k={depth} mask={mask}",
                )
                check(
                    (readout_values[tau[mask]] ^ 1) ^ (readout_values[mask] ^ 1)
                    == current_values[mask],
                    f"complement primitive k={depth} mask={mask}",
                )
        else:
            variables = len(lower_innovations)
            values = bytes(current_values)
            weight = sum(values)
            observed_carry_weights.append(weight)
            coefficients = mobius_anf(values)
            full_mask = (1 << variables) - 1
            check(coefficients[full_mask] == 1, f"innovative carry full monomial k={depth}")
            degree = max(
                (mask.bit_count() for mask, value in enumerate(coefficients) if value),
                default=-1,
            )
            check(degree == variables, f"innovative carry full degree k={depth}")

            spectrum = walsh_spectrum(values)
            zero_count = sum(coefficient == 0 for coefficient in spectrum)
            if variables >= 2:
                check(zero_count == 0, f"innovative carry full Walsh support k={depth}")
                check(
                    all(coefficient % 4 == 2 for coefficient in spectrum),
                    f"innovative carry Walsh mod four k={depth}",
                )

            extended_period = periods[depth]
            check(extended_period == 2 * lower_period, f"innovative double cover k={depth}")
            extended_seen = set()
            for phase in range(extended_period):
                lower_mask = masks[phase % lower_period]
                new_bit = arrival_readout(rows, periods, depth, phase)
                state = lower_mask | (new_bit << variables)
                check(state not in extended_seen, f"double-cover atlas collision k={depth}")
                extended_seen.add(state)
                next_lower = masks[(phase + 1) % lower_period]
                next_bit = arrival_readout(rows, periods, depth, phase + 1)
                check(next_lower == tau[lower_mask], f"skew-product lower coordinate k={depth}")
                check(next_bit == (new_bit ^ current_values[lower_mask]),
                      f"skew-product carry coordinate k={depth}")
            check(len(extended_seen) == extended_period, f"double cover one cycle k={depth}")
            carry_rows.append((depth, variables, weight, degree, zero_count))

    check(tuple(observed_carry_weights) == EXPECTED_CARRY_WEIGHTS, "frozen carry weights")
    return lift_word, carry_rows, periods


def audit_marked_arcs(periods: list[int]) -> tuple[str, str]:
    full_width = INNOVATION_DEPTH_CAP + 1
    full_period = periods[INNOVATION_DEPTH_CAP]
    packed_rows = packed_edge_cycle(full_width, full_period)
    direct_rows = centered_rows(MARKED_ARC_CAP)

    lift_bits = [
        period_lift_bit(width, periods[width - 1])
        for width in range(1, INNOVATION_DEPTH_CAP + 1)
    ]
    innovation_depths = [
        depth
        for depth in range(1, INNOVATION_DEPTH_CAP + 1)
        if lift_bits[depth - 1]
    ]

    center_word = []
    for depth in range(2, MARKED_ARC_CAP + 1):
        check(arrival_readout(packed_rows, periods, depth, -depth) == 0,
              f"light-cone boundary anchor k={depth}")

        arc = 0
        physical_ray = 0
        for phase in range(-depth, 0):
            arc ^= transition_current(packed_rows, periods, depth, phase)
            time = phase + depth
            row = direct_rows[time]
            physical_ray ^= int(
                (time - depth + 1 in row) or (time - depth + 2 in row)
            )

        center = int(0 in direct_rows[depth])
        check(arrival_readout(packed_rows, periods, depth, 0) == center,
              f"packed/direct center agreement k={depth}")
        check(arc == center, f"marked current arc k={depth}")
        check(physical_ray == center, f"independent physical left-front ray k={depth}")

        period = periods[depth - 1]
        remainder = depth % period
        reduced_arc = 0
        for phase in range(-remainder, 0):
            reduced_arc ^= transition_current(packed_rows, periods, depth, phase)
        quotient_term = ((depth // period) & 1) & lift_bits[depth - 1]
        check(center == (quotient_term ^ reduced_arc), f"period-reduced marked arc k={depth}")
        center_word.append(str(center))

    marked_innovation_word = "".join(
        str(arrival_readout(packed_rows, periods, depth, 0))
        for depth in innovation_depths
    )
    check(
        marked_innovation_word == EXPECTED_MARKED_INNOVATION_WORD,
        "frozen marked innovation word",
    )
    return "".join(center_word), marked_innovation_word


def main() -> None:
    strip_prefixes = audit_motzkin_strips()
    lift_word, carry_rows, periods = audit_innovation_cochains()
    center_word, marked_innovation_word = audit_marked_arcs(periods)

    print("THM-3471 RULE 30 MOTZKIN STRIP/CARRY EXACT AUDIT")
    print(
        "strip_universe="
        f"targets_t=2..{STRIP_TARGET_CAP}; direct_sources_s=0..{STRIP_SOURCE_CAP}; "
        f"ternary_coefficients_n=0..{GREEN_RECURRENCE_CAP}"
    )
    print("strip_source_words=alpha0=1[d>=2],alpha1=1[d>=2_even],alpha2=1[d>=1_odd]")
    print("strip_identities=R0+R1+R2=0;B+R2=1[t>=3_odd];transport_q_first_jet=R1")
    print(f"backbone_B_prefix_t2..65={strip_prefixes[0]}")
    print(f"strip_R0_prefix_t2..65={strip_prefixes[1]}")
    print(f"strip_R1_prefix_t2..65={strip_prefixes[2]}")
    print(f"strip_R2_prefix_t2..65={strip_prefixes[3]}")
    print(f"transport_first_jet_prefix_t2..65={strip_prefixes[4]}")
    print(f"innovation_lift_bits_eps1..{INNOVATION_DEPTH_CAP}={lift_word}")
    print(f"innovative_carries_k_variables_weight_degree_walshzeros={carry_rows}")
    print(f"marked_arc_universe=k=2..{MARKED_ARC_CAP}")
    print(f"center_word_c2..c{MARKED_ARC_CAP}={center_word}")
    print(f"marked_innovation_word={marked_innovation_word}")
    print("scope=FINITE-EXACT; no Rule30 prize conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
