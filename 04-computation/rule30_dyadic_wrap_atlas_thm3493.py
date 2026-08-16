#!/usr/bin/env python3
"""Exact companion for THM-3493.

The long gate follows the packed Rule 30 orbit modulo 2^68 for exactly 2^27
steps.  All checks use explicit exceptions; optimized Python therefore runs
the same verification path.
"""

from __future__ import annotations


LOW_WIDTH = 68
MAX_DYADIC_EXPONENT = 27
PERIOD_WIDTH_CAP = 33
BOUNDARY_ROW_CAP = 256
FACE_EXPONENT_CAP = 7

EXPECTED_VALUATIONS = (
    1,
    3,
    4,
    6,
    7,
    9,
    15,
    16,
    24,
    25,
    27,
    29,
    34,
    36,
    37,
    39,
    41,
    43,
    48,
    49,
    51,
    54,
    55,
    58,
    60,
    63,
    64,
    66,
)


explicit_checks = 0
checked_hasse_basis_pairs: set[tuple[int, int]] = set()


def check(condition: bool, label: str) -> None:
    """Fail explicitly, including under ``python -O``."""

    global explicit_checks
    explicit_checks += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def packed_step(row: int, mask: int | None = None) -> int:
    """One inward-packed Rule 30 step."""

    successor = row ^ ((row << 1) | (row << 2))
    return successor if mask is None else successor & mask


def valuation_two(value: int) -> int:
    check(value != 0, "2-adic valuation has nonzero input")
    return (value & -value).bit_length() - 1


def is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def direct_boundary_pair_audit() -> int:
    """Check the two persistent left-boundary bits on full packed rows."""

    row = 1
    for time in range(1, BOUNDARY_ROW_CAP + 1):
        row = packed_step(row)
        check(row.bit_length() == 2 * time + 1, f"top support time={time}")
        check((row >> (2 * time)) & 1 == 1, f"leftmost bit time={time}")
        check((row >> (2 * time - 1)) & 1 == 1, f"left neighbor time={time}")
    return BOUNDARY_ROW_CAP


def seed_period(width: int) -> int:
    """Exact return period of the seed modulo 2^width."""

    mask = (1 << width) - 1
    row = 1
    period = 0
    while True:
        row = packed_step(row, mask)
        period += 1
        if row == 1:
            return period


def finite_period_floor_audit() -> tuple[list[int], int]:
    """Finite hostile control for the period floor and lift law."""

    periods = [0]
    for width in range(1, PERIOD_WIDTH_CAP + 2):
        periods.append(seed_period(width))

    lift_checks = 0
    for width in range(1, PERIOD_WIDTH_CAP + 1):
        period = periods[width]
        dyadic_floor = 1 << (width.bit_length() - 1)
        check(is_power_of_two(period), f"power-of-two period width={width}")
        check(period >= dyadic_floor, f"dyadic period floor width={width}")
        check(width != 2 * period, f"no two-cycle endpoint width={width}")

        row = 1
        mask = (1 << (width + 1)) - 1
        for _ in range(period):
            row = packed_step(row, mask)
        epsilon = (row >> width) & 1
        check(
            periods[width + 1] == period << epsilon,
            f"period lift width={width}",
        )
        lift_checks += 1
    return periods, lift_checks


def dyadic_valuation_scout() -> list[int]:
    """Run the declared 2^27-step finite-exact universe."""

    mask = (1 << LOW_WIDTH) - 1
    row = 1
    valuations: list[int] = []
    next_dyadic_time = 1
    final_time = 1 << MAX_DYADIC_EXPONENT

    for time in range(1, final_time + 1):
        row = packed_step(row, mask)
        if time == next_dyadic_time:
            difference = row - 1
            check(difference != 0, f"visible dyadic displacement time={time}")
            value = valuation_two(difference)
            check(value < LOW_WIDTH, f"valuation inside truncation time={time}")
            valuations.append(value)
            next_dyadic_time <<= 1

    check(
        len(valuations) == MAX_DYADIC_EXPONENT + 1,
        "one valuation at every declared dyadic time",
    )
    check(
        all(valuations[index] < valuations[index + 1] for index in range(len(valuations) - 1)),
        "innovation depths strictly increase",
    )
    check(tuple(valuations) == EXPECTED_VALUATIONS, "committed dyadic valuations")
    return valuations


def wrap_ledger(valuations: list[int]) -> tuple[list[int], int, int]:
    """Apply the exact block formula to the computed valuations."""

    wrap_depths: list[int] = []
    wrapped_ones = 0
    for exponent, value in enumerate(valuations):
        block_start = 1 << exponent
        block_end = 2 * block_start - 1
        check(value <= block_end, f"strict top-neighbor bound m={exponent}")
        if value < block_start:
            wrap_length = 0
        else:
            wrap_length = value - block_start + 1
            wrapped_ones += 1
            wrap_depths.extend(range(block_start, value + 1))
        check(0 <= wrap_length <= block_start, f"wrap length m={exponent}")

    check(wrap_depths == [1, 2, 3, 4], "finite exact wrap depths")
    check(wrapped_ones == 3, "finite exact wrapped-one count")
    hard_last = (1 << (MAX_DYADIC_EXPONENT + 1)) - 1
    check(hard_last == 268_435_455, "declared hard-range endpoint")
    return wrap_depths, wrapped_ones, hard_last


def submasks(mask: int) -> list[int]:
    values: list[int] = []
    submask = mask
    while True:
        values.append(submask)
        if submask == 0:
            return values
        submask = (submask - 1) & mask


def superset_functional(period: int, order: int) -> int:
    """Truth-table mask of the Hasse coordinate M_order."""

    result = 0
    for phase in range(period):
        if phase & order == order:
            result |= 1 << phase
    return result


def hasse_basis_profile(order: int) -> int:
    """X-coordinate profile of Y^order, where Y=X+1."""

    result = 0
    phase = order
    while True:
        result |= 1 << phase
        if phase == 0:
            return result
        phase = (phase - 1) & order


def moment(profile: int, period: int, order: int) -> int:
    return (profile & superset_functional(period, order)).bit_count() & 1


def check_hasse_basis_delta(period: int, order: int) -> None:
    """Independently check that Y^order has one named Hasse coordinate."""

    key = (period, order)
    if key in checked_hasse_basis_pairs:
        return
    profile = hasse_basis_profile(order)
    coordinates = [(profile >> phase) & 1 for phase in range(period)]
    bit = 1
    while bit < period:
        for index in range(period):
            if index & bit == 0:
                coordinates[index] ^= coordinates[index | bit]
        bit <<= 1
    check(
        coordinates[order] == 1 and sum(coordinates) == 1,
        f"one-coordinate Hasse basis p={period} j={order}",
    )
    checked_hasse_basis_pairs.add(key)


def xor_masks(values: list[int]) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def mersenne_face_audit() -> tuple[int, int, int]:
    """Check both Mersenne Hasse query interfaces on exact masks."""

    cases = 0
    profile_perturbations = 0
    current_perturbations = 0
    for exponent in range(1, FACE_EXPONENT_CAP + 1):
        block_size = 1 << exponent
        depth = 2 * block_size - 1
        for period_multiplier_exponent in range(1, 4):
            period = block_size << period_multiplier_exponent
            endpoint = period - depth

            profile_orders = sorted(
                endpoint + offset for offset in submasks(depth - 1)
            )
            check(len(profile_orders) == block_size, f"profile face size m={exponent} p={period}")
            check(all(order & 1 for order in profile_orders), f"odd profile face m={exponent} p={period}")
            profile_mask = xor_masks(
                [superset_functional(period, order) for order in profile_orders]
            )
            check(profile_mask == 1 << endpoint, f"profile point functional m={exponent} p={period}")

            current_orders = sorted(
                1 + ((period - depth - 1) | offset)
                for offset in submasks(depth)
                if offset != depth
            )
            check(len(current_orders) == depth, f"current face size m={exponent} p={period}")
            check(min(current_orders) >= 1, f"current face excludes holonomy m={exponent} p={period}")
            terminal_mask = ((1 << period) - 1) ^ ((1 << endpoint) - 1)
            current_mask = xor_masks(
                [superset_functional(period, order) for order in current_orders]
            )
            check(current_mask == terminal_mask, f"current terminal functional m={exponent} p={period}")

            for order in profile_orders:
                perturbation = hasse_basis_profile(order)
                check_hasse_basis_delta(period, order)
                check(moment(perturbation, period, order) == 1, f"profile named coordinate m={exponent} j={order}")
                check((perturbation >> endpoint) & 1 == 1, f"profile perturbation flips point m={exponent} j={order}")
                profile_perturbations += 1

            for order in current_orders:
                perturbation = hasse_basis_profile(order)
                check_hasse_basis_delta(period, order)
                check(moment(perturbation, period, 0) == 0, f"current perturbation fixes holonomy m={exponent} j={order}")
                check(
                    (perturbation & terminal_mask).bit_count() & 1 == 1,
                    f"current perturbation flips terminal value m={exponent} j={order}",
                )
                current_perturbations += 1
            cases += 1
    return cases, profile_perturbations, current_perturbations


def scalar_schedule_hostiles() -> tuple[int, int]:
    """Show that the scalar tower laws alone allow both density extremes."""

    cutoff = 4095

    def audit(epsilon_word: list[int], expect_all_wrap: bool) -> int:
        period = 1
        hard = 0
        for depth, epsilon in enumerate(epsilon_word, start=1):
            dyadic_floor = 1 << (depth.bit_length() - 1)
            check(period >= dyadic_floor, f"hostile dyadic floor depth={depth}")
            check(depth != 2 * period, f"hostile no endpoint depth={depth}")
            if depth >= 3:
                check(
                    epsilon_word[depth - 3 : depth] != [1, 1, 1],
                    f"hostile no-111 depth={depth}",
                )
            if period > depth:
                hard += 1
            if epsilon:
                period *= 2
        if expect_all_wrap:
            check(hard == 0, "wrap-dense hostile has no hard depths")
        else:
            check(hard > cutoff - 32, "hard-dense hostile is eventually hard")
        return hard

    wrap_dense = [
        int(is_power_of_two(depth + 1)) for depth in range(1, cutoff + 1)
    ]
    hard_dense = [depth & 1 for depth in range(1, cutoff + 1)]
    wrap_hard_count = audit(wrap_dense, True)
    hard_hard_count = audit(hard_dense, False)
    return wrap_hard_count, hard_hard_count


def main() -> None:
    boundary_rows = direct_boundary_pair_audit()
    periods, lift_checks = finite_period_floor_audit()
    valuations = dyadic_valuation_scout()
    wraps, wrapped_ones, hard_last = wrap_ledger(valuations)
    face_cases, profile_perturbations, current_perturbations = mersenne_face_audit()
    wrap_hostile_hard, hard_hostile_hard = scalar_schedule_hostiles()

    print("THM-3493 RULE 30 DYADIC WRAP ATLAS EXACT AUDIT")
    print(f"boundary_pair_universe=full packed rows t=1..{boundary_rows}")
    print(
        "period_floor_universe="
        f"exact seed cycles widths=1..{PERIOD_WIDTH_CAP + 1};"
        f"floor/no-endpoint/lift checks widths=1..{PERIOD_WIDTH_CAP}"
    )
    print(f"period_floor_lift_checks={lift_checks};P_1..P_{PERIOD_WIDTH_CAP}={periods[1:PERIOD_WIDTH_CAP + 1]}")
    print(
        "finite_exact_universe="
        f"t=0..2^{MAX_DYADIC_EXPONENT};packed recurrence modulo 2^{LOW_WIDTH};"
        f"record t=2^m,m=0..{MAX_DYADIC_EXPONENT}"
    )
    print(f"finite_exact_step_updates={1 << MAX_DYADIC_EXPONENT}")
    print("dyadic_valuations_v0..v27=" + ",".join(map(str, valuations)))
    print(f"wrap_depths_through_2^28-1={wraps};wrapped_ones={wrapped_ones}")
    print(f"finite_exact_hard_range=5..{hard_last}")
    print(
        "mersenne_face_universe="
        f"m=1..{FACE_EXPONENT_CAP};p/2^m in {{2,4,8}};exact functional masks"
    )
    print(
        f"mersenne_face_cases={face_cases};"
        f"profile_load_bearing={profile_perturbations};"
        f"current_fixed_holonomy_load_bearing={current_perturbations}"
    )
    print(
        "scalar_hostiles=depths_1..4095;"
        f"wrap_dense_hard_count={wrap_hostile_hard};"
        f"hard_dense_hard_count={hard_hostile_hard}"
    )
    print(f"explicit_nonassert_checks={explicit_checks}")
    print(
        "scope=UNIVERSAL_PROOFS_PLUS_DECLARED_FINITE_EXACT_SCOUT;"
        "no_balance,nonperiodicity,random_access,or_fixed_seed_lower_bound_prize"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
