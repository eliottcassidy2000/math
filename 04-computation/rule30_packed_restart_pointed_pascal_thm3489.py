#!/usr/bin/env python3
"""Exact verifier for THM-3489's restart and pointed-Pascal laws.

The universal statements are proved in THM-3489.  This companion keeps four
finite, independently checkable surfaces:

* every pointed Pascal face for p=2^d, 1<=d<=8 and 0<s<p, represented as
  exact bit-packed linear functionals (so the check covers every input);
* direct cyclic arc matrices through p=64, including ranks and an explicit
  independently load-bearing image perturbation for every face coordinate;
* every Rule 30 restart phase through depth 30, with centered-light-cone,
  terminal-current, wrap, Hasse-face, and innovation-address controls; and
* exhaustive small odd-support Boolean tables plus translated-singleton
  hostiles for the ANF/Walsh point-evaluation boundary.

No finite computation is extrapolated to a Rule 30 prize conclusion.
"""

from __future__ import annotations


PASCAL_MAX_EXPONENT = 8
DIRECT_ARC_MAX_EXPONENT = 6
RULE30_DEPTH_CAP = 30
BOOLEAN_EXHAUSTIVE_MAX_DIMENSION = 3
SINGLETON_MAX_DIMENSION = 6


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def valuation_two(value: int) -> int:
    check(value > 0, "two-adic valuation has positive input")
    return (value & -value).bit_length() - 1


def submasks(mask: int) -> list[int]:
    values = []
    current = mask
    while True:
        values.append(current)
        if current == 0:
            return sorted(values)
        current = (current - 1) & mask


def superset_mask(period: int, order: int) -> int:
    """Bit mask of h with binom(h,order)=1 modulo two."""

    mask = 0
    for phase in range(period):
        if order & ~phase == 0:
            mask |= 1 << phase
    return mask


def profile_face(period: int, length: int) -> list[int]:
    check(0 < length < period, "profile face has a proper nonempty arc")
    endpoint = period - length
    return sorted(endpoint + item for item in submasks(length - 1))


def current_face(period: int, length: int) -> list[int]:
    check(0 < length < period, "current face has a proper nonempty arc")
    complement = period - length - 1
    return sorted(
        1 + (complement | item)
        for item in submasks(length)
        if item != length
    )


def xor_masks(masks: list[int]) -> int:
    total = 0
    for mask in masks:
        total ^= mask
    return total


def gf2_span_basis(vectors: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for vector in vectors:
        current = vector
        while current:
            pivot = current.bit_length() - 1
            if pivot in basis:
                current ^= basis[pivot]
            else:
                basis[pivot] = current
                break
    return basis


def gf2_in_span(vector: int, basis: dict[int, int]) -> bool:
    current = vector
    while current:
        pivot = current.bit_length() - 1
        if pivot not in basis:
            return False
        current ^= basis[pivot]
    return True


def direct_arc_columns(period: int, length: int) -> list[int]:
    """Columns of q -> (h -> xor_{0<=j<length}q(h+j))."""

    columns = []
    for input_phase in range(period):
        column = 0
        for offset in range(length):
            column ^= 1 << ((input_phase - offset) % period)
        columns.append(column)
    return columns


def y_power_profile(power: int) -> int:
    """X-coefficient mask of (X+1)^power over F_2."""

    return sum(1 << phase for phase in submasks(power))


def hasse_vector(profile: int, period: int) -> list[int]:
    return [
        (profile & superset_mask(period, order)).bit_count() & 1
        for order in range(period)
    ]


def audit_pascal_faces() -> tuple[int, int, int, int]:
    face_pairs = 0
    odd_face_pairs = 0
    direct_matrices = 0
    load_bearing_profiles = 0

    for exponent in range(1, PASCAL_MAX_EXPONENT + 1):
        period = 1 << exponent
        all_phases = (1 << period) - 1
        for length in range(1, period):
            endpoint = period - length
            pface = profile_face(period, length)
            cface = current_face(period, length)

            check(
                len(pface) == 1 << (length - 1).bit_count(),
                f"profile-face size p={period} s={length}",
            )
            check(
                len(cface) == (1 << length.bit_count()) - 1,
                f"current-face size p={period} s={length}",
            )
            check(
                xor_masks([superset_mask(period, order) for order in pface])
                == 1 << endpoint,
                f"profile point functional p={period} s={length}",
            )
            terminal_mask = all_phases ^ ((1 << endpoint) - 1)
            check(
                xor_masks([superset_mask(period, order) for order in cface])
                == terminal_mask,
                f"current terminal functional p={period} s={length}",
            )

            radical_order = (1 << valuation_two(length)) - 1
            check(
                min(pface) >= radical_order + 1,
                f"point face above rank-loss sector p={period} s={length}",
            )
            if length & 1:
                check(
                    all(order & 1 for order in pface),
                    f"odd arc has odd Hasse face p={period} s={length}",
                )
                odd_face_pairs += 1

            if exponent <= DIRECT_ARC_MAX_EXPONENT:
                columns = direct_arc_columns(period, length)
                basis = gf2_span_basis(columns)
                check(
                    len(basis) == period - radical_order,
                    f"direct arc rank p={period} s={length}",
                )
                for order in pface:
                    perturbation = y_power_profile(order)
                    check(
                        gf2_in_span(perturbation, basis),
                        f"load-bearing perturbation in image p={period} s={length} j={order}",
                    )
                    moments = hasse_vector(perturbation, period)
                    check(
                        moments == [int(index == order) for index in range(period)],
                        f"one-coordinate Hasse perturbation p={period} s={length} j={order}",
                    )
                    check(
                        (perturbation >> endpoint) & 1 == 1,
                        f"load-bearing perturbation flips point p={period} s={length} j={order}",
                    )
                    load_bearing_profiles += 1
                direct_matrices += 1

            face_pairs += 1

    return face_pairs, odd_face_pairs, direct_matrices, load_bearing_profiles


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def phi_width(row: int, width: int) -> int:
    return (row ^ ((row << 1) | (row << 2))) & ((1 << width) - 1)


def seed_cycle(width: int) -> list[int]:
    rows = []
    seen: set[int] = set()
    row = 1
    while row not in seen:
        seen.add(row)
        rows.append(row)
        row = phi_width(row, width)
    check(row == 1, f"seed returns without preperiod width={width}")
    return rows


def centered_packed_rows(horizon: int) -> list[int]:
    """Independent local-rule evolution, repacked inward from the right edge."""

    physical = frozenset({0})
    packed = []
    for time in range(horizon + 1):
        value = 0
        for inward in range(2 * time + 1):
            if time - inward in physical:
                value |= 1 << inward
        packed.append(value)
        physical = frozenset(
            site
            for site in range(-time - 1, time + 2)
            if rule30(
                int(site - 1 in physical),
                int(site in physical),
                int(site + 1 in physical),
            )
        )
    return packed


def orbit_prefix(width: int, horizon: int) -> list[int]:
    rows = []
    row = 1
    for _ in range(horizon + 1):
        rows.append(row)
        row = phi_width(row, width)
    return rows


def arc_profile(current: list[int], length: int) -> list[int]:
    period = len(current)
    value = 0
    for offset in range(length):
        value ^= current[offset % period]
    profile = []
    for phase in range(period):
        profile.append(value)
        value ^= current[phase]
        value ^= current[(phase + length) % period]
    check(value == profile[0], f"sliding arc closes p={period} m={length}")
    return profile


def hasse_transform(values: list[int]) -> list[int]:
    """M_j=sum_{h superset j} values[h], the Taylor coefficients at X=1."""

    result = list(values)
    dimension = len(values).bit_length() - 1
    check(len(values) == 1 << dimension, "Hasse table has power-of-two size")
    for bit in range(dimension):
        pivot = 1 << bit
        for mask in range(len(result)):
            if mask & pivot == 0:
                result[mask] ^= result[mask | pivot]
    return result


def mobius_anf(values: list[int]) -> list[int]:
    coefficients = list(values)
    dimension = len(values).bit_length() - 1
    check(len(values) == 1 << dimension, "ANF table has power-of-two size")
    for bit in range(dimension):
        pivot = 1 << bit
        for mask in range(len(coefficients)):
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


def anf_degree(coefficients: list[int]) -> int:
    return max(
        (mask.bit_count() for mask, value in enumerate(coefficients) if value),
        default=-1,
    )


def audit_rule30() -> tuple[
    int,
    str,
    list[tuple[int, int, int, int]],
    list[tuple[int, int, int]],
    list[tuple[int, int, int, tuple[tuple[int, int], ...]]],
    int,
    int,
]:
    cycles = {
        width: seed_cycle(width)
        for width in range(1, RULE30_DEPTH_CAP + 2)
    }
    periods = [0] + [len(cycles[width]) for width in range(1, RULE30_DEPTH_CAP + 2)]
    innovations = [0] + [
        int(periods[depth + 1] == 2 * periods[depth])
        for depth in range(1, RULE30_DEPTH_CAP + 1)
    ]
    innovative_depths = [
        depth for depth in range(1, RULE30_DEPTH_CAP + 1) if innovations[depth]
    ]

    direct_rows = centered_packed_rows(RULE30_DEPTH_CAP)
    unbounded_row = 1
    for time in range(RULE30_DEPTH_CAP + 1):
        check(unbounded_row == direct_rows[time], f"packed/direct full row t={time}")
        unbounded_row ^= (unbounded_row << 1) | (unbounded_row << 2)

    def arrival(depth: int, phase: int) -> int:
        if depth < 0:
            return 0
        cycle = cycles[depth + 1]
        return (cycle[(phase + depth) % len(cycle)] >> depth) & 1

    restart_checks = 0
    wrap_rows: list[tuple[int, int, int, int]] = []
    hard_odd_innovative: list[tuple[int, int, int]] = []
    odd_face_examples: list[tuple[int, int, int, tuple[tuple[int, int], ...]]] = []
    marked_face_checks = 0
    self_addressed_checks = 0
    full_profiles: dict[int, list[int]] = {}
    centers = {
        depth: (direct_rows[depth] >> depth) & 1
        for depth in range(1, RULE30_DEPTH_CAP + 1)
    }

    for depth in range(2, RULE30_DEPTH_CAP + 1):
        period = periods[depth]
        horizon = max(2 * period, period + depth)
        orbit = orbit_prefix(depth + 1, horizon)
        epsilon = innovations[depth]
        top_toggle = epsilon << depth

        check(period >= (depth + 1) // 2, f"edge-period lower bound k={depth}")
        check(
            orbit[period] == (1 ^ top_toggle),
            f"one-cycle seed lift k={depth}",
        )
        for time in range(period + 1):
            check(
                orbit[period + time] == (orbit[time] ^ top_toggle),
                f"packed restart k={depth} t={time}",
            )
            restart_checks += 1

        center = (orbit[depth] >> depth) & 1
        check(center == centers[depth], f"packed/direct center k={depth}")
        quotient, residual = divmod(depth, period)
        if depth >= period:
            check(
                center == ((quotient & 1) & epsilon),
                f"wrapped center formula k={depth}",
            )
            wrap_rows.append((depth, period, epsilon, center))

        current = [
            ((orbit[phase + depth] >> (depth - 1)) & 1)
            | ((orbit[phase + depth] >> (depth - 2)) & 1)
            for phase in range(period)
        ]
        check((sum(current) & 1) == epsilon, f"current holonomy k={depth}")
        profile = arc_profile(current, depth)
        full_profiles[depth] = profile
        marked_phase = (-depth) % period
        check(profile[marked_phase] == center, f"marked terminal profile k={depth}")

        physical_arc = 0
        for phase in range(-depth, 0):
            row = orbit[phase + depth]
            physical_arc ^= (
                ((row >> (depth - 1)) & 1)
                | ((row >> (depth - 2)) & 1)
            )
        check(physical_arc == center, f"physical current arc k={depth}")

        if depth < period:
            check(quotient == 0 and residual == depth, f"hard split k={depth}")
        elif depth == period:
            check(all(value == epsilon for value in profile), f"one-cycle profile k={depth}")
            check(center == epsilon, f"one-cycle center k={depth}")
        elif period < depth < 2 * period:
            check(quotient == 1 and residual > 0, f"strict-wrap split k={depth}")
            residual_profile = arc_profile(current, residual)
            endpoint = period - residual
            check(residual_profile[endpoint] == 0, f"strict-wrap residual zero k={depth}")
            check(center == epsilon, f"strict-wrap center k={depth}")
        elif depth == 2 * period:
            check(not any(profile), f"two-cycle profile k={depth}")
            check(center == 0, f"two-cycle center k={depth}")
            check(epsilon == 1, f"two-cycle endpoint innovative k={depth}")
        else:
            check(False, f"period lower bound exhausted k={depth}")

        if residual:
            endpoint = period - residual
            pface = profile_face(period, residual)
            profile_moments = hasse_transform(profile)
            check(
                (sum(profile_moments[order] for order in pface) & 1) == center,
                f"physical profile Pascal face k={depth}",
            )

            cface = current_face(period, residual)
            current_moments = hasse_transform(current)
            current_terminal = sum(current_moments[order] for order in cface) & 1
            check(
                center == (((quotient & 1) & epsilon) ^ current_terminal),
                f"physical current Pascal face k={depth}",
            )
            if period < depth < 2 * period:
                check(current_terminal == 0, f"strict-wrap high-face cancellation k={depth}")
            if depth & 1:
                check(all(order & 1 for order in pface), f"odd physical Hasse face k={depth}")
            if depth in (7, 9):
                odd_face_examples.append(
                    (
                        depth,
                        period,
                        center,
                        tuple((order, profile_moments[order]) for order in pface),
                    )
                )
            marked_face_checks += 1

        if epsilon and depth & 1 and depth < period:
            hard_odd_innovative.append((depth, period, center))

    for depth in innovative_depths:
        if depth < 2:
            continue
        period = periods[depth]
        lower_innovations = [item for item in innovative_depths if item < depth]
        check(period == 1 << len(lower_innovations), f"innovation cube size k={depth}")
        values_by_mask = [-1] * period
        phase_to_mask = []
        for phase in range(period):
            mask = sum(
                arrival(item, phase) << index
                for index, item in enumerate(lower_innovations)
            )
            check(values_by_mask[mask] == -1, f"innovation atlas injective k={depth}")
            values_by_mask[mask] = full_profiles[depth][phase]
            phase_to_mask.append(mask)
        check(all(value >= 0 for value in values_by_mask), f"innovation atlas onto k={depth}")

        origin_mask = sum(
            centers[item] << index
            for index, item in enumerate(lower_innovations)
        )
        check(origin_mask == phase_to_mask[0], f"origin is center prefix k={depth}")
        marked_mask = phase_to_mask[(-depth) % period]
        center = centers[depth]
        check(values_by_mask[marked_mask] == center, f"self-addressed recursion k={depth}")

        if depth & 1:
            coefficients = mobius_anf(values_by_mask)
            check(
                anf_degree(coefficients) == len(lower_innovations),
                f"terminal full ANF degree k={depth}",
            )
            check(
                (sum(coefficients[item] for item in submasks(marked_mask)) & 1)
                == center,
                f"ANF point face k={depth}",
            )
            spectrum = walsh_spectrum(values_by_mask)
            signed_sum = sum(
                coefficient
                * (-1 if (mode & marked_mask).bit_count() & 1 else 1)
                for mode, coefficient in enumerate(spectrum)
            )
            check(
                signed_sum == period * (1 - 2 * center),
                f"Walsh point inversion k={depth}",
            )
            dimension = len(lower_innovations)
            if dimension >= 2:
                check(
                    all(value and value % 4 == 2 for value in spectrum),
                    f"terminal full Walsh support k={depth}",
                )
                half_sum = sum(
                    (coefficient // 2)
                    * (-1 if (mode & marked_mask).bit_count() & 1 else 1)
                    for mode, coefficient in enumerate(spectrum)
                )
                check(
                    half_sum == (1 << (dimension - 1)) * (1 - 2 * center),
                    f"Walsh final signed carry k={depth}",
                )
                check(
                    half_sum % (1 << dimension) == 1 << (dimension - 1),
                    f"Walsh low two-adic layers point-blind k={depth}",
                )
        self_addressed_checks += 1

    lift_word = "".join(
        str(innovations[depth]) for depth in range(1, RULE30_DEPTH_CAP + 1)
    )
    return (
        restart_checks,
        lift_word,
        wrap_rows,
        hard_odd_innovative,
        odd_face_examples,
        marked_face_checks,
        self_addressed_checks,
    )


def audit_boolean_point_boundary() -> tuple[int, int]:
    exhaustive_tables = 0
    singleton_pairs = 0

    for dimension in range(2, BOOLEAN_EXHAUSTIVE_MAX_DIMENSION + 1):
        size = 1 << dimension
        for truth_mask in range(1 << size):
            if truth_mask.bit_count() & 1 == 0:
                continue
            values = [(truth_mask >> point) & 1 for point in range(size)]
            coefficients = mobius_anf(values)
            spectrum = walsh_spectrum(values)
            check(coefficients[-1] == 1, f"odd support top ANF d={dimension}")
            check(
                all(value and value % 4 == 2 for value in spectrum),
                f"odd support full Walsh d={dimension}",
            )
            for point in range(size):
                check(
                    (sum(coefficients[item] for item in submasks(point)) & 1)
                    == values[point],
                    f"ANF point inversion d={dimension} point={point}",
                )
                signed_sum = sum(
                    coefficient
                    * (-1 if (mode & point).bit_count() & 1 else 1)
                    for mode, coefficient in enumerate(spectrum)
                )
                check(
                    signed_sum == size * (1 - 2 * values[point]),
                    f"Walsh point inversion d={dimension} point={point}",
                )
                half_sum = signed_sum // 2
                check(
                    half_sum % size == size // 2,
                    f"Walsh sign invisible modulo cube size d={dimension} point={point}",
                )
            exhaustive_tables += 1

    for dimension in range(2, SINGLETON_MAX_DIMENSION + 1):
        size = 1 << dimension
        for marked_point in range(size):
            other_point = marked_point ^ 1
            marked_singleton = [int(point == marked_point) for point in range(size)]
            other_singleton = [int(point == other_point) for point in range(size)]
            marked_anf = mobius_anf(marked_singleton)
            other_anf = mobius_anf(other_singleton)
            marked_walsh = walsh_spectrum(marked_singleton)
            other_walsh = walsh_spectrum(other_singleton)
            check(marked_singleton[marked_point] == 1, f"singleton marked one d={dimension}")
            check(other_singleton[marked_point] == 0, f"singleton marked zero d={dimension}")
            check(anf_degree(marked_anf) == dimension, f"marked singleton degree d={dimension}")
            check(anf_degree(other_anf) == dimension, f"other singleton degree d={dimension}")
            check(
                [abs(value) for value in marked_walsh]
                == [abs(value) for value in other_walsh],
                f"singleton Walsh magnitudes agree d={dimension}",
            )
            check(
                all(value and value % 4 == 2 for value in marked_walsh + other_walsh),
                f"singleton full Walsh support d={dimension}",
            )
            singleton_pairs += 1

    return exhaustive_tables, singleton_pairs


def main() -> None:
    (
        face_pairs,
        odd_face_pairs,
        direct_matrices,
        load_bearing_profiles,
    ) = audit_pascal_faces()
    (
        restart_checks,
        lift_word,
        wrap_rows,
        hard_odd_innovative,
        odd_face_examples,
        marked_face_checks,
        self_addressed_checks,
    ) = audit_rule30()
    exhaustive_tables, singleton_pairs = audit_boolean_point_boundary()

    print("THM-3489 RULE 30 PACKED RESTART AND POINTED PASCAL EXACT AUDIT")
    print(
        "pascal_universe="
        f"p=2^d,d=1..{PASCAL_MAX_EXPONENT};s=1..p-1;"
        "exact bit-packed functionals"
    )
    print(f"pascal_face_pairs_checked={face_pairs}")
    print(f"odd_length_all_odd_face_pairs_checked={odd_face_pairs}")
    print(
        "direct_arc_universe="
        f"p=2^d,d=1..{DIRECT_ARC_MAX_EXPONENT};s=1..p-1"
    )
    print(f"direct_arc_matrices_checked={direct_matrices}")
    print(f"independently_load_bearing_image_profiles_checked={load_bearing_profiles}")
    print(
        "rule30_restart_universe="
        f"depths_k=2..{RULE30_DEPTH_CAP};all t=0..P_k"
    )
    print(f"rule30_restart_phase_checks={restart_checks}")
    print(f"innovation_lift_bits_eps1..{RULE30_DEPTH_CAP}={lift_word}")
    print(f"actual_wrap_rows_k_p_epsilon_center={wrap_rows}")
    print(f"hard_odd_innovative_rows_k_p_center={hard_odd_innovative}")
    print(f"odd_face_examples_k_p_center_order_value={odd_face_examples}")
    print(f"physical_marked_pascal_faces_checked={marked_face_checks}")
    print(f"self_addressed_innovation_recursions_checked={self_addressed_checks}")
    print(
        "boolean_universe="
        f"all odd-support tables d=2..{BOOLEAN_EXHAUSTIVE_MAX_DIMENSION};"
        f"all translated singleton pairs d=2..{SINGLETON_MAX_DIMENSION}"
    )
    print(f"odd_support_boolean_tables_checked={exhaustive_tables}")
    print(f"translated_singleton_hostile_pairs_checked={singleton_pairs}")
    print(
        "scope=UNIVERSAL_PROOFS_PLUS_FINITE_EXACT_CONTROLS;"
        "no wrap-frequency,balance,nonperiodicity,or random-access conclusion"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
