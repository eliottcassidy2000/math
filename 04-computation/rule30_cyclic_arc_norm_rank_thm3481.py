#!/usr/bin/env python3
"""Exact verifier for THM-3481's cyclic arc-norm and Rule 30 coupling.

The universal rank theorem is proved in THM-3481.  This companion checks two
finite, deliberately independent surfaces:

* direct circulant arc matrices on every phase cycle of size 1,2,...,128 and
  every arc length through two full cycles, including their ranks, complete
  Hasse constraints, and constructive odd-length inverses; and
* the physical Rule 30 edge cycles through depth 30, comparing the cyclic
  current integral with an independently evolved centered light cone.

No finite computation is extrapolated to a Rule 30 prize conclusion.
"""

from __future__ import annotations


OPERATOR_MAX_EXPONENT = 7
RULE30_DEPTH_CAP = 30


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def valuation_two(value: int) -> int:
    check(value > 0, "two-adic valuation has positive input")
    return (value & -value).bit_length() - 1


def gf2_rank(rows: list[int]) -> int:
    """Rank of a bit-packed binary matrix, independently of the theorem formula."""

    pivots: dict[int, int] = {}
    for row in rows:
        current = row
        while current:
            pivot = current.bit_length() - 1
            if pivot in pivots:
                current ^= pivots[pivot]
            else:
                pivots[pivot] = current
                break
    return len(pivots)


def direct_arc_rows(period: int, length: int) -> list[int]:
    """Rows of q |-> (h |-> xor_{0<=j<length} q(h+j))."""

    rows = []
    for phase in range(period):
        row = 0
        for offset in range(length):
            row ^= 1 << ((phase + offset) % period)
        rows.append(row)
    return rows


def hasse_pullback(rows: list[int], order: int) -> int:
    """Input mask of the output Hasse moment sum_h binom(h,order) A(h)."""

    total = 0
    for phase, row in enumerate(rows):
        if order & ~phase == 0:  # Lucas: binom(phase,order)=1 mod 2.
            total ^= row
    return total


def cyclic_norm_mask(period: int, count: int, stride: int) -> int:
    mask = 0
    for index in range(count):
        mask ^= 1 << ((index * stride) % period)
    return mask


def cyclic_convolution(left: int, right: int, period: int) -> int:
    product = 0
    for left_index in range(period):
        if (left >> left_index) & 1:
            for right_index in range(period):
                if (right >> right_index) & 1:
                    product ^= 1 << ((left_index + right_index) % period)
    return product


def audit_cyclic_operators() -> tuple[int, int, int]:
    matrix_count = 0
    hasse_count = 0
    inverse_count = 0

    for exponent in range(OPERATOR_MAX_EXPONENT + 1):
        period = 1 << exponent
        for length in range(1, 2 * period + 1):
            rows = direct_arc_rows(period, length)
            power = 1 << valuation_two(length)
            radical_order = power - 1
            predicted_rank = max(period - radical_order, 0)
            actual_rank = gf2_rank(rows)
            check(
                actual_rank == predicted_rank,
                f"arc rank p={period} m={length}",
            )
            matrix_count += 1

            constraint_count = min(radical_order, period)
            for order in range(constraint_count):
                check(
                    hasse_pullback(rows, order) == 0,
                    f"Hasse constraint p={period} m={length} j={order}",
                )
                hasse_count += 1

            if radical_order < period:
                check(
                    hasse_pullback(rows, radical_order) != 0,
                    f"first nonforced Hasse moment p={period} m={length}",
                )
            else:
                check(all(row == 0 for row in rows), f"zero arc map p={period} m={length}")

            if length & 1:
                inverse_length = pow(length, -1, 2 * period)
                norm = cyclic_norm_mask(period, length, 1)
                inverse_norm = cyclic_norm_mask(period, inverse_length, length)
                check(
                    cyclic_convolution(norm, inverse_norm, period) == 1,
                    f"constructive odd inverse p={period} m={length}",
                )
                inverse_count += 1

            if length == period:
                all_inputs = (1 << period) - 1
                check(all(row == all_inputs for row in rows), f"one-cycle norm p={period}")
            if length == 2 * period:
                check(all(row == 0 for row in rows), f"two-cycle norm p={period}")

    return matrix_count, hasse_count, inverse_count


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


def packed_edge_cycle(width: int, period: int) -> list[int]:
    rows = []
    row = 1
    for _ in range(period):
        rows.append(row)
        row = packed_phi_width(row, width)
    check(row == 1, f"packed full cycle width={width}")
    return rows


def arc_profile(current: list[int], length: int) -> list[int]:
    period = len(current)
    value = 0
    for offset in range(length):
        value ^= current[offset % period]
    profile = []
    for phase in range(period):
        profile.append(value)
        value ^= current[phase % period]
        value ^= current[(phase + length) % period]
    check(value == profile[0], f"sliding arc closes p={period} m={length}")
    return profile


def mobius_anf(values: list[int]) -> list[int]:
    coefficients = list(values)
    variables = len(values).bit_length() - 1
    check(len(values) == 1 << variables, "ANF table has power-of-two size")
    for bit in range(variables):
        pivot = 1 << bit
        for mask in range(len(values)):
            if mask & pivot:
                coefficients[mask] ^= coefficients[mask ^ pivot]
    return coefficients


def walsh_spectrum(values: list[int]) -> list[int]:
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


def audit_rule30_coupling() -> tuple[str, str, list[tuple[int, int, int, int, int]]]:
    periods = [0] + [seed_period(width) for width in range(1, RULE30_DEPTH_CAP + 2)]
    lift_bits = [0] + [
        int(periods[depth + 1] == 2 * periods[depth])
        for depth in range(1, RULE30_DEPTH_CAP + 1)
    ]
    innovations = [
        depth for depth in range(1, RULE30_DEPTH_CAP + 1) if lift_bits[depth]
    ]

    full_width = RULE30_DEPTH_CAP + 1
    full_period = periods[full_width]
    packed_rows = packed_edge_cycle(full_width, full_period)
    direct_rows = centered_rows(RULE30_DEPTH_CAP)

    def arrival(depth: int, phase: int) -> int:
        if depth < 0:
            return 0
        return (packed_rows[(phase + depth) % full_period] >> depth) & 1

    def current(depth: int, phase: int) -> int:
        check(depth >= 2, "transition current depth at least two")
        return arrival(depth - 1, phase + 1) | arrival(depth - 2, phase + 2)

    center_word = []
    innovative_rows: list[tuple[int, int, int, int, int]] = []

    for depth in range(2, RULE30_DEPTH_CAP + 1):
        period = periods[depth]
        q_values = [current(depth, phase) for phase in range(period)]
        profile = arc_profile(q_values, depth)
        epsilon = lift_bits[depth]

        check(period >= (depth + 1) // 2, f"edge-period lower bound k={depth}")
        check((sum(q_values) & 1) == epsilon, f"current holonomy k={depth}")
        for phase in range(period):
            check(
                current(depth, phase + period) == q_values[phase],
                f"current phase period k={depth} h={phase}",
            )

        for phase in range(period):
            direct_arc = 0
            for offset in range(depth):
                direct_arc ^= q_values[(phase + offset) % period]
            check(profile[phase] == direct_arc, f"direct/sliding arc k={depth} h={phase}")

        check(
            (sum(profile) & 1) == ((depth & 1) & epsilon),
            f"arc support parity k={depth}",
        )
        radical_order = (1 << valuation_two(depth)) - 1
        for order in range(min(radical_order, period)):
            moment = 0
            for phase, value in enumerate(profile):
                if value and order & ~phase == 0:
                    moment ^= 1
            check(moment == 0, f"Rule 30 terminal Hasse moment k={depth} j={order}")

        center = int(0 in direct_rows[depth])
        marked = profile[(-depth) % period]
        check(marked == center, f"marked terminal profile k={depth}")
        check(arrival(depth, 0) == center, f"packed/direct center k={depth}")
        check(arrival(depth, -depth) == 0, f"light-cone basepoint k={depth}")

        physical_arc = 0
        for phase in range(-depth, 0):
            time = phase + depth
            row = direct_rows[time]
            physical_arc ^= int((phase + 1 in row) or (phase + 2 in row))
        check(physical_arc == center, f"independent physical arc k={depth}")

        if depth == period:
            check(all(value == epsilon for value in profile), f"one-wrap profile k={depth}")
            check(center == epsilon, f"one-wrap center k={depth}")
        elif period < depth < 2 * period:
            remainder_profile = arc_profile(q_values, depth - period)
            check(
                all(profile[h] == (epsilon ^ remainder_profile[h]) for h in range(period)),
                f"one-cycle plus residual k={depth}",
            )
            check(depth - period < depth / 2, f"strictly shorter residual k={depth}")
        elif depth == 2 * period:
            check(not any(profile), f"two-wrap zero profile k={depth}")
            check(center == 0, f"two-wrap zero center k={depth}")
            check(
                ((packed_rows[period % full_period] >> depth) & 1) == 1,
                f"saturated extreme-left lift k={depth}",
            )
            check(epsilon == 1, f"saturated innovation k={depth}")
        else:
            check(depth < period, f"period lower bound trichotomy k={depth}")

        if epsilon:
            lower_innovations = [item for item in innovations if item < depth]
            check(period == 1 << len(lower_innovations), f"innovation cube size k={depth}")
            phase_masks = []
            values_by_mask = [-1] * period
            for phase in range(period):
                mask = sum(
                    arrival(item, phase) << index
                    for index, item in enumerate(lower_innovations)
                )
                check(values_by_mask[mask] == -1, f"innovation atlas injection k={depth}")
                values_by_mask[mask] = profile[phase]
                phase_masks.append(mask)
            check(all(value >= 0 for value in values_by_mask), f"innovation atlas onto k={depth}")

            coefficients = mobius_anf(values_by_mask)
            degree = max(
                (mask.bit_count() for mask, value in enumerate(coefficients) if value),
                default=-1,
            )
            spectrum = walsh_spectrum(values_by_mask)
            walsh_zeros = sum(value == 0 for value in spectrum)
            variables = len(lower_innovations)
            if depth & 1:
                check(degree == variables, f"odd innovative terminal full degree k={depth}")
                check(
                    all(
                        any(value and ((mask >> bit) & 1) for mask, value in enumerate(coefficients))
                        for bit in range(variables)
                    ),
                    f"odd innovative terminal all variables essential k={depth}",
                )
                if variables >= 2:
                    check(walsh_zeros == 0, f"odd innovative terminal full Walsh k={depth}")
                    check(
                        all(value % 4 == 2 for value in spectrum),
                        f"odd innovative terminal Walsh mod four k={depth}",
                    )
            innovative_rows.append(
                (depth, period, sum(values_by_mask), degree, walsh_zeros)
            )

        center_word.append(str(center))

    lift_word = "".join(str(lift_bits[depth]) for depth in range(1, RULE30_DEPTH_CAP + 1))
    return lift_word, "".join(center_word), innovative_rows


def main() -> None:
    matrix_count, hasse_count, inverse_count = audit_cyclic_operators()
    lift_word, center_word, innovative_rows = audit_rule30_coupling()

    print("THM-3481 RULE 30 CYCLIC ARC-NORM EXACT AUDIT")
    print(
        "operator_universe="
        f"p=2^d,d=0..{OPERATOR_MAX_EXPONENT};m=1..2p;"
        "all direct circulant matrices"
    )
    print(f"operator_matrices_checked={matrix_count}")
    print("operator_rank=rank(N_m)=max(p-2^nu2(m)+1,0)")
    print(f"complete_hasse_constraints_checked={hasse_count}")
    print(f"constructive_odd_inverses_checked={inverse_count}")
    print("boundary_controls=N_p has rank 1;N_2p=0")
    print(f"rule30_universe=depths_k=2..{RULE30_DEPTH_CAP};all phases_mod_Pk")
    print(f"innovation_lift_bits_eps1..{RULE30_DEPTH_CAP}={lift_word}")
    print(f"terminal_innovative_k_p_weight_degree_walshzeros={innovative_rows}")
    print(f"center_word_c2..c{RULE30_DEPTH_CAP}={center_word}")
    print("marked_controls=cyclic_profile=physical_left_cone_arc=center")
    print("scope=UNIVERSAL_OPERATOR_PROOF_PLUS_FINITE_EXACT_RULE30_AUDIT;no prize conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
