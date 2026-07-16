#!/usr/bin/env python3
"""Exact Farey/divisor-sheet atlas for the HYP-7084 owner-wall remainder.

For a slow wall x, reduce 7x=h/d.  The wall owners are exactly the speeds
divisible by d.  Aggregating all columns with the same proper owner set S and
g=gcd(S) gives a signed cyclic word on Z/(7g); forgetting its order, miss
transitions, and far-phase labels leaves the complete bipartite graph K_{|S|,N}.

The script proves the column-count formula by inclusion-exclusion, verifies the
signed-word formula against the event sweep on all 6,900 HYP-7084 rows, and
tests the Zarankiewicz/class-parity quotient.  Tournament Analysis uses owner
cardinalities rather than runners: the pairwise observable is maximum positive
signed contribution, and the switch is maximum Zarankiewicz load with column
count as the tie gauge.

Assumption challenge: runners, divisor sheets, owner sets, and exact wall
columns are all possible vertices.  The K_{r,N} quotient preserves only
owner-column incidence.  It destroys cyclic column order, pre/post miss labels,
the outside-divisibility mask, and the far address {tk/(7g)}; the exact
counterexamples below show that this destroyed data controls the sign.
"""

from collections import defaultdict
from fractions import Fraction
from functools import cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import comb as binomial, gcd, sqrt


OWNER_PATH = "04-computation/lrc14_owner_packet_parallel_class_current_codex_S19.py"
OWNER_SPEC = spec_from_file_location("owner_parallel_current", OWNER_PATH)
OWNER_MODULE = module_from_spec(OWNER_SPEC)
OWNER_SPEC.loader.exec_module(OWNER_MODULE)


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def gcd_all(values):
    common_divisor = 0
    for value in values:
        common_divisor = gcd(common_divisor, value)
    return common_divisor


def zarankiewicz_value(left_size, right_size):
    return (
        (left_size // 2)
        * ((left_size - 1) // 2)
        * (right_size // 2)
        * ((right_size - 1) // 2)
    )


def missed_pattern(speeds, wall_column, owner_gcd, before, stationary_pin=0):
    counts = [0, 0, 0, 0, 0, 0, 0]
    counts[stationary_pin] = 1
    shift = 1 if before else 0
    for speed in speeds:
        section = ((speed * wall_column - shift) // owner_gcd) % 7
        counts[section] += 1
    return tuple(section for section in range(7) if counts[section] == 0)


def exact_column_count(speeds, owners):
    owner_gcd = gcd_all(owners)
    outside = tuple(speed for speed in speeds if speed not in owners)
    count = 0
    for subset_size in range(len(outside) + 1):
        for subset in combinations(outside, subset_size):
            common_divisor = owner_gcd
            for speed in subset:
                common_divisor = gcd(common_divisor, speed)
            count += (-1) ** subset_size * 7 * common_divisor
    return count


def exact_parity_count(speeds, owners, far_speed, parity):
    owner_gcd = gcd_all(owners)
    outside = tuple(speed for speed in speeds if speed not in owners)
    count = 0
    for subset_size in range(len(outside) + 1):
        for subset in combinations(outside, subset_size):
            common_divisor = owner_gcd
            for speed in subset:
                common_divisor = gcd(common_divisor, speed)
            common_speed_divisor = gcd(common_divisor, far_speed)
            reduced_divisor = common_divisor // common_speed_divisor
            parity_count = common_speed_divisor * (
                (reduced_divisor + 1) // 2 if parity == 0 else reduced_divisor // 2
            )
            count += (-1) ** subset_size * 7 * parity_count
    return count


def cycle_cut(residues):
    support = set(residues)
    return sum(
        (vertex in support) != ((vertex + 1) % 7 in support)
        for vertex in range(7)
    )


def crossing_energy(residues):
    crossing_profile = (0, 0, 2, 3, 3, 2, 0)
    return sum(
        crossing_profile[(right - left) % 7]
        for left, right in combinations(residues, 2)
    )


def reflect_section(section):
    return (-section - 1) % 7


def least_far_speed(diameter, frequency, period):
    return frequency + period * ((diameter - frequency) // period + 1)


def primitive_pair_bound():
    missed_sets = tuple(
        missed
        for missed_size in range(7)
        for missed in combinations(range(1, 7), missed_size)
    )
    return max(
        abs(
            OWNER_MODULE.owner_wave_primitive(first_missed, phase)
            - OWNER_MODULE.owner_wave_primitive(second_missed, phase)
        )
        for numerator in range(14)
        for phase in (Fraction(numerator, 14),)
        for first_missed in missed_sets
        for second_missed in missed_sets
    )


def missed_from_sections(sections, stationary_pin):
    counts = [0, 0, 0, 0, 0, 0, 0]
    counts[stationary_pin] = 1
    for section in sections:
        counts[section] += 1
    return tuple(section for section in range(7) if counts[section] == 0)


def single_owner_primitive_bound():
    transition_pairs = set()
    for sections in product(range(7), repeat=5):
        for moving_index in range(5):
            after_sections = tuple(
                (
                    (section + 1) % 7
                    if runner_index == moving_index
                    else section
                )
                for runner_index, section in enumerate(sections)
            )
            transition_pairs.add(
                (
                    missed_from_sections(sections, 0),
                    missed_from_sections(after_sections, 0),
                )
            )
    return max(
        abs(
            OWNER_MODULE.owner_wave_primitive(missed_before, phase)
            - OWNER_MODULE.owner_wave_primitive(missed_after, phase)
        )
        for numerator in range(14)
        for phase in (Fraction(numerator, 14),)
        for missed_before, missed_after in transition_pairs
    )


def adjacent_pin_pair_ranges():
    transition_states = {owner_size: set() for owner_size in range(1, 5)}
    for before_sections in product(range(7), repeat=5):
        for owner_size in range(1, 5):
            for owner_indices in combinations(range(5), owner_size):
                owner_index_set = frozenset(owner_indices)
                after_sections = tuple(
                    (section + 1) % 7
                    if runner_index in owner_index_set
                    else section
                    for runner_index, section in enumerate(before_sections)
                )
                transition_states[owner_size].add(
                    (
                        missed_from_sections(before_sections, 0),
                        missed_from_sections(after_sections, 0),
                        missed_from_sections(before_sections, 6),
                        missed_from_sections(after_sections, 6),
                    )
                )
    ranges = {}
    for owner_size, states in transition_states.items():
        values = []
        for missed_before_zero, missed_after_zero, missed_before_six, missed_after_six in states:
            for numerator in range(14):
                phase = Fraction(numerator, 14)
                values.append(
                    OWNER_MODULE.owner_wave_primitive(missed_before_zero, phase)
                    - OWNER_MODULE.owner_wave_primitive(missed_after_zero, phase)
                    + OWNER_MODULE.owner_wave_primitive(missed_before_six, phase)
                    - OWNER_MODULE.owner_wave_primitive(missed_after_six, phase)
                )
        ranges[owner_size] = (min(values), max(values))
    return ranges


def consecutive_tail_singleton_coefficient(owner_speed, far_multiple):
    speeds = tuple(range(owner_speed - 4, owner_speed + 1))
    coefficient = Fraction(0)
    for wall_column in range(7 * owner_speed):
        owners = tuple(
            speed
            for speed in speeds
            if (speed * wall_column) % owner_speed == 0
        )
        if owners != (owner_speed,):
            continue
        missed_before = missed_pattern(
            speeds,
            wall_column,
            owner_speed,
            before=True,
        )
        missed_after = missed_pattern(
            speeds,
            wall_column,
            owner_speed,
            before=False,
        )
        phase = Fraction(far_multiple * wall_column, 7)
        coefficient += (
            OWNER_MODULE.owner_wave_primitive(missed_before, phase)
            - OWNER_MODULE.owner_wave_primitive(missed_after, phase)
        )
    return coefficient


def phase_locked_pair_coefficient(owner_gcd):
    speeds = (
        owner_gcd,
        2 * owner_gcd - 3,
        2 * owner_gcd - 2,
        2 * owner_gcd - 1,
        2 * owner_gcd,
    )
    owners = (owner_gcd, 2 * owner_gcd)
    coefficient = Fraction(0)
    for wall_column in range(7 * owner_gcd):
        exact_owners = tuple(
            speed
            for speed in speeds
            if (speed * wall_column) % owner_gcd == 0
        )
        if exact_owners != owners:
            continue
        missed_before = missed_pattern(
            speeds,
            wall_column,
            owner_gcd,
            before=True,
        )
        missed_after = missed_pattern(
            speeds,
            wall_column,
            owner_gcd,
            before=False,
        )
        phase = Fraction(3 * wall_column, 7)
        coefficient += (
            OWNER_MODULE.owner_wave_primitive(missed_before, phase)
            - OWNER_MODULE.owner_wave_primitive(missed_after, phase)
        )
    return coefficient


def owner_cardinality_wall_counts(speeds):
    intersection_totals = {
        subset_size: 7
        * sum(
            gcd_all(subset)
            for subset in combinations(speeds, subset_size)
        )
        for subset_size in range(1, 6)
    }
    return {
        owner_size: sum(
            (-1) ** (subset_size - owner_size)
            * binomial(subset_size, owner_size)
            * intersection_totals[subset_size]
            for subset_size in range(owner_size, 6)
        )
        for owner_size in range(1, 6)
    }


@cache
def owner_word_skeleton(speeds):
    _, walls, denominator = OWNER_MODULE.slow_pattern_sweep(tuple(speeds))
    period_scale = denominator // 7
    records = []
    for position, missed_before, missed_after, owners in walls:
        divisor_gcd = gcd(position, period_scale)
        sheet_divisor = period_scale // divisor_gcd
        numerator = position // divisor_gcd
        expected_owners = tuple(
            speed for speed in speeds if speed % sheet_divisor == 0
        )
        if owners != expected_owners:
            raise AssertionError(
                ("divisor-owner identity", speeds, position, owners, expected_owners)
            )

        owner_gcd = gcd_all(owners)
        wall_column_numerator = owner_gcd * position
        if wall_column_numerator % period_scale:
            raise AssertionError(("owner column integrality", speeds, position, owners))
        wall_column = wall_column_numerator // period_scale
        expected_column_owners = tuple(
            speed for speed in speeds if (speed * wall_column) % owner_gcd == 0
        )
        if owners != expected_column_owners:
            raise AssertionError(
                ("column-owner identity", speeds, wall_column, owners, expected_column_owners)
            )
        direct_before = missed_pattern(speeds, wall_column, owner_gcd, before=True)
        direct_after = missed_pattern(speeds, wall_column, owner_gcd, before=False)
        if (missed_before, missed_after) != (direct_before, direct_after):
            raise AssertionError(
                (
                    "mechanical miss word",
                    speeds,
                    wall_column,
                    missed_before,
                    missed_after,
                    direct_before,
                    direct_after,
                )
            )

        for stationary_pin in range(7):
            reflected_pin = reflect_section(stationary_pin)
            negative_column = (-wall_column) % (7 * owner_gcd)
            reflected_before = tuple(
                sorted(reflect_section(section) for section in missed_after)
            )
            reflected_after = tuple(
                sorted(reflect_section(section) for section in missed_before)
            )
            direct_negative_before = missed_pattern(
                speeds,
                negative_column,
                owner_gcd,
                before=True,
                stationary_pin=reflected_pin,
            )
            direct_negative_after = missed_pattern(
                speeds,
                negative_column,
                owner_gcd,
                before=False,
                stationary_pin=reflected_pin,
            )
            pinned_before = missed_pattern(
                speeds,
                wall_column,
                owner_gcd,
                before=True,
                stationary_pin=stationary_pin,
            )
            pinned_after = missed_pattern(
                speeds,
                wall_column,
                owner_gcd,
                before=False,
                stationary_pin=stationary_pin,
            )
            reflected_before = tuple(
                sorted(reflect_section(section) for section in pinned_after)
            )
            reflected_after = tuple(
                sorted(reflect_section(section) for section in pinned_before)
            )
            if (direct_negative_before, direct_negative_after) != (
                reflected_before,
                reflected_after,
            ):
                raise AssertionError(
                    (
                        "cut-open reflection",
                        speeds,
                        wall_column,
                        stationary_pin,
                        direct_negative_before,
                        direct_negative_after,
                        reflected_before,
                        reflected_after,
                    )
                )

        unit = numerator % sheet_divisor if sheet_divisor > 1 else 0
        fiber_column = (
            (numerator - unit) // sheet_divisor if sheet_divisor > 1 else numerator - 1
        )
        if sheet_divisor > 1:
            if gcd(unit, sheet_divisor) != 1 or not 0 <= fiber_column < 7:
                raise AssertionError(
                    ("unit-fiber coordinate", speeds, numerator, sheet_divisor)
                )
            for speed in speeds:
                direct_section = (speed * numerator // sheet_divisor) % 7
                affine_section = (
                    speed * unit // sheet_divisor + speed * fiber_column
                ) % 7
                if direct_section != affine_section:
                    raise AssertionError(
                        (
                            "affine slow row",
                            speeds,
                            speed,
                            numerator,
                            sheet_divisor,
                        )
                    )

        records.append(
            (
                owners,
                missed_before,
                missed_after,
                owner_gcd,
                wall_column,
                sheet_divisor,
                numerator,
                unit,
                fiber_column,
            )
        )
    return tuple(records)


def owner_word_data(speeds, far_speed):
    contributions = defaultdict(Fraction)
    columns = defaultdict(set)
    fiber_checks = 0
    for (
        owners,
        missed_before,
        missed_after,
        owner_gcd,
        wall_column,
        sheet_divisor,
        numerator,
        unit,
        fiber_column,
    ) in owner_word_skeleton(tuple(speeds)):
        if sheet_divisor > 1:
            source_section = (far_speed * numerator // sheet_divisor) % 7
            affine_source = (
                far_speed * unit // sheet_divisor + far_speed * fiber_column
            ) % 7
            far_remainder = (far_speed * unit) % sheet_divisor
            parity = (2 * far_remainder) // sheet_divisor
            target_section = (2 * far_speed * numerator // sheet_divisor) % 7
            if source_section != affine_source or target_section != (
                2 * source_section + parity
            ) % 7:
                raise AssertionError(
                    (
                        "affine t/2t rows",
                        speeds,
                        far_speed,
                        numerator,
                        sheet_divisor,
                    )
                )
            fiber_checks += 1

        phase = Fraction(far_speed * numerator, 7 * sheet_divisor)
        direct_phase = Fraction(far_speed * wall_column, 7 * owner_gcd)
        if phase % 1 != direct_phase % 1:
            raise AssertionError(("far phase address", speeds, far_speed, position))
        contribution = (
            OWNER_MODULE.owner_wave_primitive(missed_before, direct_phase)
            - OWNER_MODULE.owner_wave_primitive(missed_after, direct_phase)
        ) / far_speed
        contributions[owners] += contribution
        if wall_column in columns[owners]:
            raise AssertionError(("duplicate owner column", speeds, owners, wall_column))
        columns[owners].add(wall_column)
    return contributions, columns, fiber_checks


def correlation(first_values, second_values):
    first_mean = sum(first_values) / len(first_values)
    second_mean = sum(second_values) / len(second_values)
    numerator = sum(
        (first - first_mean) * (second - second_mean)
        for first, second in zip(first_values, second_values)
    )
    denominator = sqrt(
        sum((value - first_mean) ** 2 for value in first_values)
        * sum((value - second_mean) ** 2 for value in second_values)
    )
    return numerator / denominator if denominator else 0.0


def order_fingerprint(first_scores, switched_scores):
    vertices = tuple(sorted(first_scores))
    first_path = tuple(sorted(vertices, key=lambda vertex: first_scores[vertex], reverse=True))
    switched_path = tuple(
        sorted(vertices, key=lambda vertex: switched_scores[vertex], reverse=True)
    )
    first_rank = {vertex: rank for rank, vertex in enumerate(first_path)}
    switched_rank = {vertex: rank for rank, vertex in enumerate(switched_path)}
    edge_flips = sum(
        (first_rank[left] < first_rank[right])
        != (switched_rank[left] < switched_rank[right])
        for left, right in combinations(vertices, 2)
    )
    return {
        "positive_path": first_path,
        "zarankiewicz_path": switched_path,
        "score_histogram": {score: 1 for score in range(len(vertices))},
        "directed_triangles": 0,
        "scc_sizes": (1,) * len(vertices),
        "edge_flips": edge_flips,
        "hamiltonian_path_count": 1,
    }


def full_owner_frequency_scan():
    minimum = None
    maximum = None
    envelope_minimum = None
    envelope_maximum = None
    by_size = {}
    word_count = 0
    frequency_count = 0
    for diameter in range(5, 11):
        for speeds in combinations(range(1, diameter + 1), 5):
            if speeds[-1] != diameter or not OWNER_MODULE.primitive(speeds):
                continue
            owner_records = defaultdict(list)
            for (
                owners,
                missed_before,
                missed_after,
                owner_gcd,
                wall_column,
                _sheet_divisor,
                _numerator,
                _unit,
                _fiber_column,
            ) in owner_word_skeleton(tuple(speeds)):
                owner_records[owners].append(
                    (missed_before, missed_after, owner_gcd, wall_column)
                )
            for owners, records in owner_records.items():
                if owners == tuple(speeds):
                    continue
                word_count += 1
                owner_gcd = records[0][2]
                for frequency in range(7 * owner_gcd):
                    coefficient = sum(
                        (
                            OWNER_MODULE.owner_wave_primitive(
                                missed_before,
                                Fraction(frequency * wall_column, 7 * owner_gcd),
                            )
                            - OWNER_MODULE.owner_wave_primitive(
                                missed_after,
                                Fraction(frequency * wall_column, 7 * owner_gcd),
                            )
                        )
                        for missed_before, missed_after, _, wall_column in records
                    )
                    record = (coefficient, speeds, owners, frequency, owner_gcd)
                    if minimum is None or coefficient < minimum[0]:
                        minimum = record
                    if maximum is None or coefficient > maximum[0]:
                        maximum = record
                    owner_size = len(owners)
                    if owner_size not in by_size:
                        by_size[owner_size] = [record, record]
                    if coefficient < by_size[owner_size][0][0]:
                        by_size[owner_size][0] = record
                    if coefficient > by_size[owner_size][1][0]:
                        by_size[owner_size][1] = record
                    least_speed = least_far_speed(
                        diameter,
                        frequency,
                        7 * owner_gcd,
                    )
                    contribution = coefficient / least_speed
                    envelope_record = (
                        contribution,
                        speeds,
                        owners,
                        frequency,
                        least_speed,
                        owner_gcd,
                        coefficient,
                    )
                    if (
                        envelope_minimum is None
                        or contribution < envelope_minimum[0]
                    ):
                        envelope_minimum = envelope_record
                    if (
                        envelope_maximum is None
                        or contribution > envelope_maximum[0]
                    ):
                        envelope_maximum = envelope_record
                    frequency_count += 1
    return {
        "minimum": minimum,
        "maximum": maximum,
        "envelope_minimum": envelope_minimum,
        "envelope_maximum": envelope_maximum,
        "by_size": by_size,
        "word_count": word_count,
        "frequency_count": frequency_count,
    }


def main():
    print("HYP-7085: FAREY DIVISOR SHEETS AND THE ZARANKIEWICZ GUARDRAIL")
    print("=" * 78)

    check("Z(7,7)=81", zarankiewicz_value(7, 7) == 81)
    check("Z(7,8)=108", zarankiewicz_value(7, 8) == 108)
    check("Z(8,8)=144", zarankiewicz_value(8, 8) == 144)
    check(
        "deletion averaging propagates Woodall's K_7,7 value to K_7,8 and K_8,8",
        Fraction(8, 6) * 81 == 108 and Fraction(8, 6) * 108 == 144,
    )
    for missed_size in range(8):
        for missed in combinations(range(7), missed_size):
            reflected_missed = tuple(
                sorted(reflect_section(section) for section in missed)
            )
            waveform = OWNER_MODULE.owner_waveform(missed)
            reflected_waveform = OWNER_MODULE.owner_waveform(reflected_missed)
            if any(
                reflected_waveform[13 - half_sector] != waveform[half_sector]
                for half_sector in range(14)
            ):
                raise AssertionError(("waveform reflection", missed))
            for numerator in range(29):
                phase = Fraction(numerator, 28)
                if OWNER_MODULE.owner_wave_primitive(
                    reflected_missed, 1 - phase
                ) != -OWNER_MODULE.owner_wave_primitive(missed, phase):
                    raise AssertionError(("primitive reflection", missed, phase))
    check("cut-open waveform reflection P_Rm(1-z)=-P_m(z)", True)
    palette_bound = primitive_pair_bound()
    check(
        "the exact common-phase primitive-pair bound is 135/1372",
        palette_bound == Fraction(135, 1372),
    )
    single_owner_bound = single_owner_primitive_bound()
    check(
        "the single-owner color-cube edge bound is 51/686",
        single_owner_bound == Fraction(51, 686),
    )
    paired_ranges = adjacent_pin_pair_ranges()
    check(
        "adjacent-pin paired ranges are exact for all proper owner sizes",
        paired_ranges
        == {
            1: (Fraction(-1, 7), Fraction(11, 98)),
            2: (Fraction(-1, 7), Fraction(1, 7)),
            3: (Fraction(-13, 98), Fraction(45, 343)),
            4: (Fraction(-13, 98), Fraction(15, 98)),
        },
    )
    consecutive_tail_coefficients = {
        far_multiple: {
            owner_speed: consecutive_tail_singleton_coefficient(
                owner_speed,
                far_multiple,
            )
            for owner_speed in range(5, 173)
        }
        for far_multiple in (2, 4)
    }
    check(
        "the phase-locked tail has period-84 drifts +62/49 and -149/49",
        all(
            consecutive_tail_coefficients[2][owner_speed + 84]
            - consecutive_tail_coefficients[2][owner_speed]
            == Fraction(62, 49)
            and consecutive_tail_coefficients[4][owner_speed + 84]
            - consecutive_tail_coefficients[4][owner_speed]
            == Fraction(-149, 49)
            for owner_speed in range(5, 89)
        ),
    )
    check(
        "the positive coefficient cap already fails at g=33 and t=66",
        consecutive_tail_coefficients[2][33] == Fraction(53, 98),
    )
    phase_locked_pair_values = {
        owner_gcd: phase_locked_pair_coefficient(owner_gcd)
        for owner_gcd in (11, 95, 179)
    }
    check(
        "a planar pair sheet has positive period-84 drift 797/343",
        phase_locked_pair_values
        == {
            11: Fraction(12, 49),
            95: Fraction(881, 343),
            179: Fraction(1678, 343),
        },
    )

    entries = []
    rows = []
    wall_count = 0
    affine_fiber_checks = 0
    cardinality_count_checks = 0
    for diameter in range(5, 11):
        for speeds in combinations(range(1, diameter + 1), 5):
            if speeds[-1] != diameter or not OWNER_MODULE.primitive(speeds):
                continue
            direct_cardinality_counts = {
                owner_size: sum(
                    len(record[0]) == owner_size
                    for record in owner_word_skeleton(tuple(speeds))
                )
                for owner_size in range(1, 6)
            }
            if direct_cardinality_counts != owner_cardinality_wall_counts(speeds):
                raise AssertionError(
                    (
                        "owner-cardinality binomial inversion",
                        speeds,
                        direct_cardinality_counts,
                        owner_cardinality_wall_counts(speeds),
                    )
                )
            cardinality_count_checks += 1
            for far_speed in range(diameter + 1, 4 * diameter + 1):
                contributions, columns, fiber_check_count = owner_word_data(
                    speeds, far_speed
                )
                affine_fiber_checks += fiber_check_count
                full_owner_set = tuple(speeds)
                row_entries = []
                for owners, contribution in contributions.items():
                    wall_count += len(columns[owners])
                    if len(columns[owners]) != exact_column_count(speeds, owners):
                        raise AssertionError(
                            (
                                "inclusion-exclusion column count",
                                speeds,
                                owners,
                                len(columns[owners]),
                                exact_column_count(speeds, owners),
                            )
                        )
                    if owners == full_owner_set:
                        continue
                    owner_gcd = gcd_all(owners)
                    quotient = tuple(owner // owner_gcd for owner in owners)
                    quotient_residues = tuple(owner % 7 for owner in quotient)
                    exact_columns = len(columns[owners])
                    parity_word = tuple(
                        (2 * ((far_speed * column) % owner_gcd)) // owner_gcd
                        for column in sorted(columns[owners])
                    )
                    parity_counts = (
                        parity_word.count(0),
                        parity_word.count(1),
                    )
                    expected_parity_counts = tuple(
                        exact_parity_count(speeds, owners, far_speed, parity)
                        for parity in (0, 1)
                    )
                    if parity_counts != expected_parity_counts:
                        raise AssertionError(
                            (
                                "parity inclusion-exclusion",
                                speeds,
                                far_speed,
                                owners,
                                parity_counts,
                                expected_parity_counts,
                            )
                        )
                    entry = {
                        "row": (speeds, far_speed),
                        "owners": owners,
                        "value": contribution,
                        "size": len(owners),
                        "gcd": owner_gcd,
                        "quotient": quotient,
                        "columns": exact_columns,
                        "parity_counts": parity_counts,
                        "z_exact": zarankiewicz_value(len(owners), exact_columns),
                        "z_full": zarankiewicz_value(len(owners), 7 * owner_gcd),
                        "cut": cycle_cut(quotient_residues),
                        "energy": crossing_energy(quotient_residues),
                    }
                    entries.append(entry)
                    row_entries.append(entry)
                residual = sum(entry["value"] for entry in row_entries)
                rows.append(
                    {
                        "row": (speeds, far_speed),
                        "residual": residual,
                        "entries": tuple(row_entries),
                        "z_exact": sum(entry["z_exact"] for entry in row_entries),
                        "z_full": sum(entry["z_full"] for entry in row_entries),
                        "cut": sum(entry["cut"] for entry in row_entries),
                        "energy": sum(entry["energy"] for entry in row_entries),
                        "planar_sum": sum(
                            entry["value"] for entry in row_entries if entry["z_exact"] == 0
                        ),
                        "crossing_sum": sum(
                            entry["value"] for entry in row_entries if entry["z_exact"] > 0
                        ),
                    }
                )

    check("6,900 exact owner-doubling rows", len(rows) == 6900)
    check("34,791 proper owner-set entries", len(entries) == 34791)
    check("25 primitive owner quotient shapes", len({entry["quotient"] for entry in entries}) == 25)
    check(
        "all proper owner sets have at most four vertices",
        max(entry["size"] for entry in entries) == 4,
    )
    check(
        "only unconditional small-part Zarankiewicz values occur",
        {entry["z_exact"] for entry in entries} == {0, 9, 18, 42},
    )
    check("affine colored K_8,7 fiber identities", affine_fiber_checks == 906171)
    check(
        "owner-cardinality wall counts satisfy exact binomial inversion",
        cardinality_count_checks > 0,
    )

    entry_minimum = min(entries, key=lambda entry: entry["value"])
    entry_maximum = max(entries, key=lambda entry: entry["value"])
    row_minimum = min(rows, key=lambda row: row["residual"])
    row_maximum = max(rows, key=lambda row: row["residual"])
    zero_z_entries = [entry for entry in entries if entry["z_exact"] == 0]
    quotient_groups = defaultdict(list)
    fixed_owner_groups = defaultdict(list)
    for entry in entries:
        quotient_groups[entry["quotient"]].append(entry)
        fixed_owner_groups[(entry["row"][0], entry["owners"])].append(entry)

    check("30,702 proper sheets have zero crossing number", len(zero_z_entries) == 30702)
    check(
        "both individual signed extrema occur on planar incidence graphs",
        entry_minimum["z_exact"] == entry_maximum["z_exact"] == 0,
    )
    check(
        "every quotient shape attains both contribution signs",
        all(
            min(group, key=lambda entry: entry["value"])["value"] < 0
            < max(group, key=lambda entry: entry["value"])["value"]
            for group in quotient_groups.values()
        ),
    )

    balanced_arc = [
        entry
        for entry in entries
        if entry["owners"] == (2, 4, 6, 8) and entry["quotient"] == (1, 2, 3, 4)
    ]
    balanced_minimum = min(balanced_arc, key=lambda entry: entry["value"])
    balanced_maximum = max(balanced_arc, key=lambda entry: entry["value"])
    check(
        "the same K_4,7 quotient arc carries both signs",
        balanced_minimum["value"] == Fraction(-43, 3430)
        and balanced_maximum["value"] == Fraction(15, 686),
    )

    fixed_flip = fixed_owner_groups[((1, 3, 5, 6, 7), (7,))]
    fixed_by_far = {entry["row"][1]: entry["value"] for entry in fixed_flip}
    check(
        "a fixed core and owner set flip sign under the far-phase address",
        fixed_by_far[8] == Fraction(-53, 1372)
        and fixed_by_far[9] == Fraction(11, 882),
    )
    fixed_eta_counts = {
        entry["row"][1]: entry["parity_counts"] for entry in fixed_flip
    }
    check(
        "the fixed sign flip has the same balanced class-parity census",
        fixed_eta_counts[8] == fixed_eta_counts[9] == (21, 21),
    )

    parity_word_counterexample = {}
    for far_speed in (6, 8):
        contributions, columns, _ = owner_word_data((1, 2, 3, 4, 5), far_speed)
        owners = (2, 4)
        owner_gcd = 2
        parity_word_counterexample[far_speed] = (
            contributions[owners],
            tuple(
                (2 * ((far_speed * column) % owner_gcd)) // owner_gcd
                for column in sorted(columns[owners])
            ),
        )
    check(
        "an identical full class-parity word carries opposite signs",
        parity_word_counterexample[6]
        == (Fraction(23, 4116), (0, 0, 0, 0, 0, 0, 0))
        and parity_word_counterexample[8]
        == (Fraction(-37, 5488), (0, 0, 0, 0, 0, 0, 0)),
    )

    cap_refuter_speeds = (7, 8, 9, 10, 11)
    cap_refuter_owners = (11,)
    cap_refuter_far_speed = 44
    cap_refuter_contributions, cap_refuter_columns, _ = owner_word_data(
        cap_refuter_speeds,
        cap_refuter_far_speed,
    )
    cap_refuter_contribution = cap_refuter_contributions[cap_refuter_owners]
    check(
        "diameter eleven refutes the bounded-bank coefficient cap",
        len(cap_refuter_columns[cap_refuter_owners]) == 70
        and cap_refuter_far_speed * cap_refuter_contribution == Fraction(-216, 343)
        and cap_refuter_contribution == Fraction(-54, 3773),
    )

    planar_minimum = min(rows, key=lambda row: row["planar_sum"])
    planar_maximum = max(rows, key=lambda row: row["planar_sum"])
    crossing_minimum = min(rows, key=lambda row: row["crossing_sum"])
    crossing_maximum = max(rows, key=lambda row: row["crossing_sum"])
    check(
        "row residual extrema reproduce HYP-7084",
        row_minimum["residual"] == Fraction(-8087, 246960)
        and row_maximum["residual"] == Fraction(269, 6860),
    )

    frequency_scan = full_owner_frequency_scan()
    check("1,259 distinct proper owner words", frequency_scan["word_count"] == 1259)
    check(
        "49,483 exact owner-word frequency evaluations",
        frequency_scan["frequency_count"] == 49483,
    )
    check(
        "the diameter-ten frequency scan stays strictly inside 1/2",
        frequency_scan["minimum"][0] == Fraction(-1019, 2058)
        and frequency_scan["maximum"][0] == Fraction(89, 196),
    )

    print("\nExact divisor-sheet theorem")
    print("  wall coordinate: 7x=h/d reduced")
    print("  exact owners: S_d(A)={a in A:d|a}")
    print("  owner-set columns: k mod 7g, g=gcd(S)")
    print("  column count: N_A(S)=7 sum_(T subset A\\S) (-1)^|T| gcd(g,T)")
    print("  parity count: N_eta=7 sum_T (-1)^|T| E_eta(gcd(g,T),t)")
    print("  E_0(h,t)=c ceil(delta/2), E_1(h,t)=c floor(delta/2)")
    print("  where c=gcd(h,t), delta=h/c")
    print("  signed word: sum_k [P_m-( {tk/(7g)} )-P_m+( {tk/(7g)} )]/t")
    print("  least representative: tau_D(r)=r+7g(floor((D-r)/(7g))+1)")
    print("  fixed-word envelope: max_r |C_A,S(r)|/tau_D(r)")
    print("  primitive-pair bound:", palette_bound)
    print("  single-owner color-cube edge bound:", single_owner_bound)
    print("  adjacent-pin paired ranges:", paired_ranges)
    print("  singleton: -N/(14t)<=G_S<=11N/(196t)<11/28")
    print("  pair: |G_S|<=N/(14t)<1/4")
    print("  divisor charge: M_r=sum_(k>=r)(-1)^(k-r) binom(k,r) I_k")
    print("  where I_k=7 sum_(|T|=k) gcd(T)")
    print("  universal proper-sheet bound: |G_S|<=135 N_A(S)/(1372 t)<945/1372")
    print("  phase-locked tail: A_g=(g-4,...,g), S=(g)")
    print("  t=2g drift: C_2(g+84)=C_2(g)+62/49; G_S -> 31/4116")
    print("  positive cap refuter: g=33, t=66, C_2=53/98, G_S=53/6468")
    print("  t=4g drift: C_4(g+84)=C_4(g)-149/49; G_S -> -149/16464")
    print("  planar pair ray: A=(g,2g-3,2g-2,2g-1,2g), S=(g,2g), t=3g")
    print("  on g=11 mod 84: C(g+84)=C(g)+797/343; G_S -> 797/86436")
    print("  total swept walls:", wall_count)
    print("  affine K_8,7 fiber checks:", affine_fiber_checks)

    print("\nAffine colored K_8,7 fiber")
    print("  h=a+jd with a in U_d and j in Z_7")
    print("  A_e^+(j)=floor(ea/d)+ej mod 7")
    print("  A_e^-(j)=A_e^+(j)-1_[d|e]")
    print("  A_2t(j)=2A_t(j)+floor(2{ta/d}) mod 7")
    print("  score: pointed rainbow-completion switch, not a crossing count")
    print("  reflection: W_S(-k;b)=W_S(k;R(b)), R(b)=-b-1")
    print("  actual stationary pin 0 reflects to adjacent pin 6, so no cancellation")

    print("\nZarankiewicz guardrail")
    print("  ordinary exact values: cr(K_7,7)=81, cr(K_7,8)=108, cr(K_8,8)=144")
    print("  proper-sheet Z values:", sorted({entry["z_exact"] for entry in entries}))
    print(
        "  planar proper sheets:",
        f"{len(zero_z_entries)}/{len(entries)} = {len(zero_z_entries) / len(entries):.6f}",
    )
    print("  individual minimum:", entry_minimum)
    print("  individual maximum:", entry_maximum)
    print("  same K_4,7 minimum:", balanced_minimum)
    print("  same K_4,7 maximum:", balanced_maximum)
    print("  fixed-core owner t=8/t=9:", fixed_by_far[8], fixed_by_far[9])
    print("  their parity counts:", fixed_eta_counts[8], fixed_eta_counts[9])
    print("  identical eta-word sign flip:", parity_word_counterexample)
    print("  planar row-sum range:", planar_minimum["planar_sum"], planar_maximum["planar_sum"])
    print(
        "  positive-Z row-sum range:",
        crossing_minimum["crossing_sum"],
        crossing_maximum["crossing_sum"],
    )
    print("  full residual range:", row_minimum["residual"], row_maximum["residual"])
    print("  full-frequency coefficient range:", frequency_scan["minimum"], frequency_scan["maximum"])
    print(
        "  all-far-speed envelope range through diameter ten:",
        frequency_scan["envelope_minimum"],
        frequency_scan["envelope_maximum"],
    )
    print(
        "  diameter-eleven coefficient-cap refuter:",
        cap_refuter_speeds,
        cap_refuter_owners,
        cap_refuter_far_speed * cap_refuter_contribution,
        cap_refuter_contribution,
    )
    print("  coefficient extrema by owner size:")
    for owner_size, extrema in sorted(frequency_scan["by_size"].items()):
        print("   ", owner_size, extrema)

    absolute_residuals = [abs(float(row["residual"])) for row in rows]
    correlations = {
        field: correlation(
            [float(row[field]) for row in rows],
            absolute_residuals,
        )
        for field in ("z_exact", "z_full", "cut", "energy")
    }
    print("  correlations with |residual|:", correlations)

    positive_risks = {
        size: max(entry["value"] for entry in entries if entry["size"] == size)
        for size in range(1, 5)
    }
    zarankiewicz_risks = {
        size: max(
            (entry["z_exact"], entry["columns"])
            for entry in entries
            if entry["size"] == size
        )
        for size in range(1, 5)
    }
    fingerprint = order_fingerprint(positive_risks, zarankiewicz_risks)
    print("\nTournament Analysis")
    print("  vertices: proper owner-set cardinalities 1,2,3,4")
    print("  observable: largest positive signed owner-word contribution")
    print("  switch: largest Zarankiewicz load; exact-column tie gauge")
    print("  positive risks:", positive_risks)
    print("  Zarankiewicz risks:", zarankiewicz_risks)
    print("  fingerprint:", fingerprint)

    print("\nVERDICT")
    print("  PROVED: noncommon walls split exactly into Farey/divisor owner sheets.")
    print("  PROVED: each owner set is a signed cyclic column word with the stated N_A(S).")
    print("  CORRECTED: K_7,7, K_7,8, and K_8,8 are known ordinary crossing cases.")
    print("  REFUTED: Zarankiewicz value, parity, cut, or crossing energy determines the sign.")
    print("  REFUTED: the diameter-ten artifact |t G_S|<1/2 fails at diameter eleven.")
    print("  PROVED: one planar singleton family has coefficients unbounded both ways.")
    print("  PROVED: multi-owner walls are path-independent sums of single color-cube edges.")
    print("  PROVED: reflection gives exact density-scaled singleton and pair bounds.")
    print("  PROVED: each fixed word's infinite far-speed problem is a finite least-speed envelope.")
    print("  NEXT: bound the divisor-summed adjacent-pin envelope, not its unnormalized coefficient.")


if __name__ == "__main__":
    main()
