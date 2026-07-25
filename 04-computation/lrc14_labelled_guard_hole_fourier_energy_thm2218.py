#!/usr/bin/env python3
"""Exact stress audit for the THM-2204 Fourier/signed-energy addendum.

This is the canonical finite stress artifact for THM-2218.  It verifies:

* the oriented spread-2460 lift and guard-hole vectors;
* the failure of a lower-fibre top-five margin and the common-lift regret;
* the exact signed-energy certificate on seven concentrated base families
  over every (1,1,3) and (1,2,3) two-blocker residual;
* the exact integer family-knapsack bound on the mixed hostile row (799,46).

All certificate decisions use integers.  Floats are printed only to make
the size of a proved integer margin readable.  ``require`` remains active
under ``python -O``.
"""

from collections import defaultdict
import math

import numpy as np


P = 13
N = P**4
Q = P**3
Q0 = P**2


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign_classes(modulus):
    values = np.arange(1, (modulus + 1) // 2, dtype=np.int64)
    return values[values % P != 0]


def norm_numerators(values, modulus):
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def row_to_bitset(row):
    packed = np.packbits(row, bitorder="little")
    return int.from_bytes(packed.tobytes(), "little")


def rows_to_bitsets(rows):
    packed = np.packbits(rows, axis=1, bitorder="little")
    return tuple(int.from_bytes(row.tobytes(), "little") for row in packed)


def ceil_sqrt(value):
    root = math.isqrt(value)
    return root if root * root == value else root + 1


def integer_signed_certificate(capacities, residual_size):
    """Return the exact k=5 signed-energy certificate data."""
    total = sum(capacities)
    centered = tuple(13 * value - total for value in capacities)
    positive = sum(value * value for value in centered if value > 0)
    negative = sum(value * value for value in centered if value < 0)
    delta = 13 * residual_size - 5 * total
    radicand = min(5 * positive, 8 * negative)
    gap = delta * delta - radicand
    exact = delta > 0 and gap > 0
    readable_bound = (5 * total + math.sqrt(radicand)) / 13
    return {
        "total": total,
        "centered": centered,
        "positive": positive,
        "negative": negative,
        "delta": delta,
        "radicand": radicand,
        "gap": gap,
        "exact": exact,
        "bound": readable_bound,
        "slack": residual_size - readable_bound,
    }


def run():
    phases = np.arange(1, Q, dtype=np.int64)
    phases = phases[phases % P != 0]
    sheets = np.arange(P, dtype=np.int64)
    roots = phases[:, None] + Q * sheets[None, :]
    guard = 7 * norm_numerators(roots, N) > N
    guard_counts = guard.sum(axis=1)
    guard_nine = row_to_bitset(guard_counts == 9)
    phase_full = (1 << len(phases)) - 1

    unit_labels = sign_classes(N)
    base_labels = sign_classes(Q)
    shallow_labels = sign_classes(P**2)

    require(len(phases) == 2028, "bad phase count")
    require(len(unit_labels) == 13182, "bad unit count")
    require(len(base_labels) == 1014, "bad base count")
    require(len(shallow_labels) == 78, "bad shallow count")
    require(
        tuple(
            (size, int(np.sum(guard_counts == size)))
            for size in (9, 10)
        )
        == ((9, 1450), (10, 578)),
        "bad guard histogram",
    )

    unit_blocker_rows = (
        14
        * norm_numerators(base_labels[:, None] * phases[None, :], Q)
        < Q
    )
    divisible_blocker_rows = (
        14
        * norm_numerators(
            (P * shallow_labels[:, None]) * phases[None, :],
            Q,
        )
        < Q
    )
    unit_blockers = rows_to_bitsets(unit_blocker_rows)
    divisible_blockers = rows_to_bitsets(divisible_blocker_rows)

    def residual_size(residual):
        return (
            10 * residual.bit_count()
            - (residual & guard_nine).bit_count()
        )

    def oriented_family(base):
        return tuple((base + lift * Q) % N for lift in range(P))

    def family_weight_bitsets(base):
        at_least_one = []
        at_least_two = []
        for coefficient in oriented_family(base):
            dangerous = (
                14
                * norm_numerators(coefficient * roots, N)
                < N
            )
            counts = np.sum(dangerous & guard, axis=1)
            require(int(counts.max()) <= 2, "unit window too large")
            at_least_one.append(row_to_bitset(counts >= 1))
            at_least_two.append(row_to_bitset(counts >= 2))
        return tuple(at_least_one), tuple(at_least_two)

    def family_capacities(residual, at_least_one, at_least_two):
        return tuple(
            (residual & one).bit_count()
            + (residual & two).bit_count()
            for one, two in zip(at_least_one, at_least_two)
        )

    # ------------------------------------------------------------------
    # Spread-2460 row and its exact local obstruction.
    # ------------------------------------------------------------------
    shallow_four_index = tuple(map(int, shallow_labels)).index(4)
    spread_residual_bits = (
        phase_full ^ divisible_blockers[shallow_four_index]
    )
    spread_residual_row = ~divisible_blocker_rows[shallow_four_index]
    spread_coefficients = oriented_family(1098)
    spread_capacities = []
    spread_holes = []
    spread_endpoint_mass = []
    spread_phase_weights = []

    for coefficient in spread_coefficients:
        dangerous = (
            14 * norm_numerators(coefficient * roots, N) < N
        )
        phase_weights = np.sum(dangerous & guard, axis=1)
        spread_phase_weights.append(phase_weights)
        spread_capacities.append(
            int(np.sum(dangerous & guard & spread_residual_row[:, None]))
        )
        spread_holes.append(
            int(
                np.sum(
                    dangerous
                    & (~guard)
                    & spread_residual_row[:, None]
                )
            )
        )
        spread_endpoint_mass.append(
            int(np.sum(dangerous & spread_residual_row[:, None]))
        )

    spread_capacities = tuple(spread_capacities)
    spread_holes = tuple(spread_holes)
    spread_endpoint_mass = tuple(spread_endpoint_mass)
    require(
        spread_capacities
        == (
            2460,
            2456,
            2454,
            2450,
            2456,
            2452,
            0,
            2452,
            2450,
            2456,
            2452,
            2452,
            2454,
        ),
        "spread capacity vector changed",
    )
    require(
        spread_holes
        == (
            730,
            734,
            736,
            740,
            734,
            738,
            3190,
            738,
            740,
            734,
            738,
            738,
            736,
        ),
        "spread hole vector changed",
    )
    require(set(spread_endpoint_mass) == {3190}, "endpoint mass changed")

    spread_size = residual_size(spread_residual_bits)
    spread_certificate = integer_signed_certificate(
        spread_capacities, spread_size
    )
    require(spread_size == 15932, "spread residual size changed")
    require(spread_certificate["exact"], "spread certificate failed")
    require(
        (
            spread_certificate["positive"],
            spread_certificate["negative"],
            spread_certificate["delta"],
            spread_certificate["radicand"],
            spread_certificate["gap"],
        )
        == (
            72261760,
            866949136,
            59896,
            361308800,
            3226222016,
        ),
        "spread signed-energy invoice changed",
    )

    phase_to_index = {
        int(phase): index for index, phase in enumerate(phases)
    }
    lower_phases = np.arange(1, Q0, dtype=np.int64)
    lower_phases = lower_phases[lower_phases % P != 0]
    lower_groups = tuple(
        tuple(
            phase_to_index[int(lower + digit * Q0)]
            for digit in range(P)
        )
        for lower in lower_phases
    )

    local_margins = []
    local_rows = []
    for lower, group in zip(map(int, lower_phases), lower_groups):
        selected = tuple(
            index for index in group if spread_residual_row[index]
        )
        local_residual = int(guard_counts[list(selected)].sum())
        local_capacities = tuple(
            int(weights[list(selected)].sum())
            for weights in spread_phase_weights
        )
        local_top = sum(sorted(local_capacities, reverse=True)[:5])
        margin = local_residual - local_top
        local_margins.append(margin)
        if margin == -5:
            local_rows.append(
                (lower, local_residual, local_top, local_capacities)
            )

    require(sum(margin < 0 for margin in local_margins) == 22,
            "negative local-fibre count changed")
    require(min(local_margins) == -5, "local minimum changed")
    require(sum(local_margins) == -18, "local margin sum changed")
    require(
        local_rows[:2]
        == [
            (
                4,
                120,
                125,
                (25, 7, 25, 13, 20, 25, 0, 25, 20, 13, 25, 7, 25),
            ),
            (
                6,
                120,
                125,
                (25, 25, 25, 20, 13, 7, 0, 7, 13, 20, 25, 25, 25),
            ),
        ],
        "named local hostile rows changed",
    )
    spread_top_five = sum(sorted(spread_capacities, reverse=True)[:5])
    local_envelope = spread_size - sum(local_margins)
    common_lift_regret = local_envelope - spread_top_five
    require(
        (spread_top_five, local_envelope, common_lift_regret)
        == (12282, 15950, 3668),
        "common-lift regret changed",
    )

    # ------------------------------------------------------------------
    # Exact signed-energy audit on seven concentrated base families.
    # ------------------------------------------------------------------
    stress_bases = (1, 1098, 799, 599, 1015, 597, 1013)
    stress_rows = []
    for base in stress_bases:
        at_least_one, at_least_two = family_weight_bitsets(base)
        for profile, left_masks, right_masks, symmetric in (
            ("113", unit_blockers, unit_blockers, True),
            ("123", unit_blockers, divisible_blockers, False),
        ):
            checked = 0
            failures = 0
            worst_slack = None
            worst_record = None
            for left_index, left_mask in enumerate(left_masks):
                right_indices = (
                    range(left_index, len(right_masks))
                    if symmetric
                    else range(len(right_masks))
                )
                for right_index in right_indices:
                    residual = (
                        phase_full
                        ^ (left_mask | right_masks[right_index])
                    )
                    capacities = family_capacities(
                        residual, at_least_one, at_least_two
                    )
                    size = residual_size(residual)
                    certificate = integer_signed_certificate(
                        capacities, size
                    )
                    checked += 1
                    if not certificate["exact"]:
                        failures += 1
                    if (
                        worst_slack is None
                        or certificate["slack"] < worst_slack
                    ):
                        worst_slack = certificate["slack"]
                        worst_record = (
                            left_index,
                            right_index,
                            size,
                            capacities,
                            certificate["gap"],
                        )
            expected = 514605 if profile == "113" else 79092
            require(checked == expected, "bad stress universe count")
            require(failures == 0, "signed-energy stress failure")
            stress_rows.append(
                (
                    base,
                    profile,
                    checked,
                    failures,
                    worst_slack,
                    worst_record,
                )
            )

    # ------------------------------------------------------------------
    # Exact global family-knapsack certificate on the mixed hostile row.
    # ------------------------------------------------------------------
    depth_one_799_index = tuple(map(int, base_labels)).index(799)
    depth_two_46_index = tuple(map(int, shallow_labels)).index(46)
    mixed_residual_row = ~(
        unit_blocker_rows[depth_one_799_index]
        | divisible_blocker_rows[depth_two_46_index]
    )
    mixed_residual_bits = (
        phase_full
        ^ (
            unit_blockers[depth_one_799_index]
            | divisible_blockers[depth_two_46_index]
        )
    )
    mixed_size = residual_size(mixed_residual_bits)

    all_capacities = []
    for start in range(0, len(unit_labels), 64):
        batch = unit_labels[start : start + 64]
        dangerous = (
            14
            * norm_numerators(batch[:, None, None] * roots[None, :, :], N)
            < N
        )
        capacities = np.sum(
            dangerous
            & guard[None, :, :]
            & mixed_residual_row[None, :, None],
            axis=(1, 2),
        )
        all_capacities.extend(map(int, capacities))

    actual_global_rows = sorted(
        zip(all_capacities, map(int, unit_labels)),
        key=lambda item: (-item[0], item[1]),
    )[:5]
    require(
        (mixed_size, tuple(actual_global_rows))
        == (
            13526,
            (
                (2604, 5193),
                (2472, 10386),
                (2292, 7773),
                (2288, 10388),
                (2262, 7775),
            ),
        ),
        "mixed hostile row changed",
    )

    families = defaultdict(list)
    for label, capacity in zip(map(int, unit_labels), all_capacities):
        residue = label % Q
        base = min(residue, (-residue) % Q)
        families[base].append(capacity)
    require(len(families) == 1014, "bad lift-family count")
    require(set(map(len, families.values())) == {13},
            "bad lift-family size")

    negative_infinity = -(10**100)
    dp = [0] + [negative_infinity] * 5
    paths = [[] for _ in range(6)]
    for base in map(int, base_labels):
        capacities = families[base]
        total = sum(capacities)
        centered = tuple(13 * value - total for value in capacities)
        positive = sum(value * value for value in centered if value > 0)
        negative = sum(value * value for value in centered if value < 0)
        upper_numerators = [0]
        for count in range(1, 6):
            radicand = min(
                count * positive,
                (13 - count) * negative,
            )
            upper_numerators.append(
                count * total + ceil_sqrt(radicand)
            )

        next_dp = [negative_infinity] * 6
        next_paths = [None] * 6
        for total_count in range(6):
            for count in range(total_count + 1):
                candidate = (
                    dp[total_count - count]
                    + upper_numerators[count]
                )
                if candidate > next_dp[total_count]:
                    next_dp[total_count] = candidate
                    next_paths[total_count] = (
                        paths[total_count - count]
                        + (
                            [(base, count, upper_numerators[count])]
                            if count
                            else []
                        )
                    )
        dp = next_dp
        paths = next_paths

    require(dp[5] == 159477, "mixed knapsack numerator changed")
    require(13 * mixed_size - dp[5] == 16361,
            "mixed knapsack margin changed")
    require(
        paths[5]
        == [
            (1, 1, 31769),
            (599, 1, 32136),
            (799, 1, 33852),
            (1015, 1, 29796),
            (1098, 1, 31924),
        ],
        "mixed knapsack optimizer changed",
    )

    print("THM-2204 labelled guard-hole Fourier stress audit")
    print("arithmetic: exact integers for every certificate decision")
    print()
    print("spread-2460 oriented family")
    print(f"  coefficients={spread_coefficients}")
    print(f"  capacities={spread_capacities}")
    print(f"  hole_vector={spread_holes}")
    print(f"  endpoint_mass={spread_endpoint_mass[0]}")
    print(
        f"  residual_size={spread_size} "
        f"actual_top5={spread_top_five} "
        f"actual_margin={spread_size-spread_top_five}"
    )
    print(
        "  signed_invoice="
        f"(Zplus={spread_certificate['positive']},"
        f"Zminus={spread_certificate['negative']},"
        f"D={spread_certificate['delta']},"
        f"radicand={spread_certificate['radicand']},"
        f"gap={spread_certificate['gap']})"
    )
    print(
        f"  readable_signed_bound={spread_certificate['bound']:.12f} "
        f"bound_margin={spread_certificate['slack']:.12f}"
    )
    print()
    print("lower-fibre obstruction")
    print(
        f"  negative_fibres={sum(margin < 0 for margin in local_margins)} "
        f"minimum_margin={min(local_margins)} "
        f"sum_local_margins={sum(local_margins)}"
    )
    print(f"  first_two_minimum_rows={tuple(local_rows[:2])}")
    print(
        f"  local_envelope={local_envelope} "
        f"global_common_lift_top5={spread_top_five} "
        f"common_lift_regret={common_lift_regret}"
    )
    print()
    print("seven-family exhaustive signed-energy stress")
    for (
        base,
        profile,
        checked,
        failures,
        worst_slack,
        worst_record,
    ) in stress_rows:
        print(
            f"  base={base} profile={profile} checked={checked} "
            f"failures={failures} "
            f"worst_readable_slack={worst_slack:.12f} "
            f"worst_indices={worst_record[0:2]} "
            f"worst_residual={worst_record[2]} "
            f"worst_exact_gap={worst_record[4]}"
        )
    print()
    print("mixed hostile family-knapsack certificate")
    print(f"  blocker_labels=(799,46) residual_size={mixed_size}")
    print(f"  actual_top_five={tuple(actual_global_rows)}")
    print(f"  actual_margin={mixed_size-sum(x[0] for x in actual_global_rows)}")
    print(
        f"  exact_ceil_knapsack_numerator={dp[5]} "
        f"bound={dp[5]}/13={dp[5]/13:.12f}"
    )
    print(
        f"  exact_margin_numerator={13*mixed_size-dp[5]} "
        f"margin={(13*mixed_size-dp[5])/13:.12f}"
    )
    print(f"  optimizing_bound_path={tuple(paths[5])}")
    print()
    print("CONCLUSION: all exact stress checks PASS")


if __name__ == "__main__":
    run()
