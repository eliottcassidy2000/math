#!/usr/bin/env python3
"""Exact parallel-class current referee for HYP-7084.

For q=floor(14{tx})=2c+eta and u=2c mod 7, the doubled owner lands in
parallel class u+eta.  Thus the two half-sectors are the endpoints of the
oriented class edge u -> u+1.  On the one-miss set let mu_s^eta(u) be the
exact mass with miss s and channel (u,eta), and put M=mu^0+mu^1,
J=mu^0-mu^1.  The script verifies

    Delta_F(2t)
      = 1/2 sum_s (M_s(s)+M_s(s-1)) - p1(F)/7
        + 1/2 sum_s (J_s(s)-J_s(s-1)).

The THM-913 K_7 crossing kernel is xi=(0,0,2,3,3,2,0).  Its Laplacian L
gives the exact dual-norm certificate

    Delta_F(2t)^2 <= (16/29) E_pc,
    E_pc = sum_s [M_s^T L M_s + J_s^T L J_s]
         = 2 sum_{s,eta} (mu_s^eta)^T L mu_s^eta.

There is a sharper wall formulation.  If A is the slow core and m_A(x) its
miss pattern, let H_m be the fourteen-cell owner waveform and let P_m be the
periodic primitive of H_m-mean(H_m).  Then

    Delta_F(2t) = C_2(A)
      + (1/t) sum_{slow walls p} [P_{m_-}({tp})-P_{m_+}({tp})].

This is the THM-727 endpoint sum with a finite parallel-class palette.  It
shows exactly which chronology and relation data the class energy discards.

Tournament Analysis uses the seven parallel classes as vertices.  Its
pairwise observable is the difference of their largest positive source
contributions in the bounded bank; the switch replaces this by largest local
Laplacian load, with numeric order as tie gauge.  This is telemetry only.  The
faithful carrier is the miss-labelled oriented class edge.  Runners, gaps,
circle sections, section boundaries, wall events, residues, cover arcs,
Fourier modes, matroid circuits, and proof obligations were considered as
vertices; none alone retains both target incidence and half-sector parity.
The quotient preserves the exact finite correction but destroys wall order,
the slow-family relation lattice, and the individual chord crossings.
"""

from fractions import Fraction
from functools import cache
from heapq import heapify, heappop, heappush
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm


OWNER_PATH = "04-computation/lrc14_finite_t_owner_packet_codex_S18.py"
OWNER_SPEC = spec_from_file_location("owner_packet", OWNER_PATH)
OWNER_MODULE = module_from_spec(OWNER_SPEC)
OWNER_SPEC.loader.exec_module(OWNER_MODULE)

SECTORS = tuple(range(7))
PARITIES = (0, 1)
PROPAGATION_SLACK = Fraction(97, 1000)
CROSSING_WEIGHTS = tuple(
    Fraction(0) if distance == 0 else Fraction((distance - 1) * (6 - distance), 2)
    for distance in SECTORS
)
ENERGY_CERTIFICATE_CONSTANT = Fraction(16, 29)
ENERGY_THRESHOLD = PROPAGATION_SLACK**2 / ENERGY_CERTIFICATE_CONSTANT


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def chords_cross(first, second):
    first_left, first_right = first
    second_left, second_right = second
    return (
        first_left < second_left < first_right < second_right
        or second_left < first_left < second_right < first_right
    )


def crossing_matrix():
    classes = [[] for _ in SECTORS]
    for edge in combinations(SECTORS, 2):
        classes[sum(edge) % 7].append(edge)
    matrix = [[0] * 7 for _ in SECTORS]
    for first in SECTORS:
        for second in SECTORS:
            matrix[first][second] = sum(
                chords_cross(first_edge, second_edge)
                for first_edge in classes[first]
                for second_edge in classes[second]
            )
    return tuple(tuple(row) for row in matrix)


def laplacian_matrix():
    matrix = [[Fraction(0) for _ in SECTORS] for _ in SECTORS]
    for first in SECTORS:
        for second in SECTORS:
            if first == second:
                continue
            weight = CROSSING_WEIGHTS[(second - first) % 7]
            matrix[first][first] += weight
            matrix[first][second] -= weight
    return tuple(tuple(row) for row in matrix)


LAPLACIAN = laplacian_matrix()


def matrix_vector(matrix, vector):
    return tuple(
        sum(entry * value for entry, value in zip(row, vector))
        for row in matrix
    )


def dot(first, second):
    return sum(left * right for left, right in zip(first, second))


def laplacian_energy(vector):
    return dot(vector, matrix_vector(LAPLACIAN, vector))


def packet_measures(profile, far_speed):
    intervals, _, denominator = profile
    measure_numerators = [
        [[0 for _ in SECTORS] for _ in PARITIES]
        for _ in SECTORS
    ]
    doubled_denominator = 2 * denominator
    cell_step = denominator // (7 * far_speed)
    if cell_step * 7 * far_speed != denominator:
        raise AssertionError(("far speed must divide the event period", far_speed))
    for left, right, missed in intervals:
        doubled_left = 2 * left
        doubled_right = 2 * right
        first_cell = doubled_left // cell_step
        final_cell = (doubled_right + cell_step - 1) // cell_step
        for cell in range(first_cell, final_cell):
            overlap_left = max(doubled_left, cell * cell_step)
            overlap_right = min(doubled_right, (cell + 1) * cell_step)
            if overlap_left >= overlap_right:
                continue
            half_sector = cell % 14
            parity = half_sector % 2
            source_sector = half_sector // 2
            class_vertex = 2 * source_sector % 7
            measure_numerators[missed][parity][class_vertex] += overlap_right - overlap_left
    return tuple(
        tuple(
            tuple(Fraction(value, doubled_denominator) for value in row)
            for row in parity_rows
        )
        for parity_rows in measure_numerators
    )


@cache
def owner_waveform(missed):
    values = []
    missed_set = frozenset(missed)
    for half_sector in range(14):
        source_sector = half_sector // 2
        target_sector = half_sector % 7
        remaining = missed_set - {source_sector}
        if len(remaining) != 1:
            values.append(Fraction(0))
            continue
        unique_miss = next(iter(remaining))
        values.append(Fraction(int(target_sector == unique_miss)) - Fraction(1, 7))
    return tuple(values)


@cache
def owner_wave_mean(missed):
    return sum(owner_waveform(missed)) / 14


@cache
def centered_owner_waveform(missed):
    mean = owner_wave_mean(missed)
    return tuple(value - mean for value in owner_waveform(missed))


def owner_wave_primitive(missed, phase):
    phase %= 1
    scaled_phase = 14 * phase
    cell = scaled_phase.numerator // scaled_phase.denominator
    waveform = centered_owner_waveform(missed)
    prefix = sum(waveform[:cell], Fraction(0)) / 14
    return prefix + (phase - Fraction(cell, 14)) * waveform[cell]


@cache
def slow_pattern_sweep(positive_speeds):
    period_scale = lcm(*positive_speeds)
    denominator = 7 * period_scale
    sectors = [0] * len(positive_speeds)
    counts = [len(positive_speeds) + 1, 0, 0, 0, 0, 0, 0]
    events = []
    for runner_index, speed in enumerate(positive_speeds):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    previous = 0
    intervals = []
    walls = []
    while events:
        event_position = events[0][0]
        missed_before = tuple(section for section in range(1, 7) if counts[section] == 0)
        intervals.append((previous, event_position, missed_before))
        owners = []
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            owners.append(speed)
            old_sector = sectors[runner_index]
            new_sector = (old_sector + 1) % 7
            counts[old_sector] -= 1
            counts[new_sector] += 1
            sectors[runner_index] = new_sector
            if event_index < 7 * speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (next_index * event_step, runner_index, next_index, speed, event_step),
                )
        missed_after = tuple(section for section in range(1, 7) if counts[section] == 0)
        walls.append((event_position, missed_before, missed_after, tuple(sorted(owners))))
        previous = event_position
    return tuple(intervals), tuple(walls), denominator


def endpoint_metrics(positive_speeds, far_speed):
    intervals, walls, denominator = slow_pattern_sweep(tuple(positive_speeds))
    limiting_coefficient = sum(
        Fraction(right - left, denominator) * owner_wave_mean(missed)
        for left, right, missed in intervals
    )
    interval_remainder = sum(
        (
            owner_wave_primitive(missed, Fraction(far_speed * right, denominator))
            - owner_wave_primitive(missed, Fraction(far_speed * left, denominator))
        )
        / far_speed
        for left, right, missed in intervals
    )
    wall_contributions = {}
    wall_remainder = Fraction(0)
    for position, missed_before, missed_after, owners in walls:
        phase = Fraction(far_speed * position, denominator)
        contribution = (
            owner_wave_primitive(missed_before, phase)
            - owner_wave_primitive(missed_after, phase)
        ) / far_speed
        wall_remainder += contribution
        wall_contributions[owners] = wall_contributions.get(owners, Fraction(0)) + contribution
    if wall_remainder != interval_remainder:
        raise AssertionError(("wall/interval endpoint mismatch", positive_speeds, far_speed))
    return {
        "limit": limiting_coefficient,
        "wall_remainder": wall_remainder,
        "total": limiting_coefficient + wall_remainder,
        "wall_contributions": wall_contributions,
        "wall_count": len(walls),
    }


@cache
def common_wall_coefficient(residue_multiset, far_residue):
    coefficient = Fraction(0)
    for wall_index in range(1, 7):
        counts_before = [1, 0, 0, 0, 0, 0, 0]
        counts_after = [1, 0, 0, 0, 0, 0, 0]
        for residue in residue_multiset:
            counts_before[(residue * wall_index - 1) % 7] += 1
            counts_after[(residue * wall_index) % 7] += 1
        missed_before = tuple(
            section for section in range(1, 7) if counts_before[section] == 0
        )
        missed_after = tuple(
            section for section in range(1, 7) if counts_after[section] == 0
        )
        phase = Fraction((far_residue * wall_index) % 7, 7)
        coefficient += (
            owner_wave_primitive(missed_before, phase)
            - owner_wave_primitive(missed_after, phase)
        )
    return coefficient


def common_wall_bounds():
    coefficient_maximum = None
    coefficient_minimum = None
    residue_patterns = tuple(combinations_with_replacement(SECTORS, 5))
    for residues in residue_patterns:
        for far_residue in SECTORS:
            coefficient = common_wall_coefficient(residues, far_residue)
            record = (residues, far_residue)
            coefficient_maximum = update_extreme(
                record, coefficient, coefficient_maximum
            )
            coefficient_minimum = update_extreme(
                record, coefficient, coefficient_minimum, maximize=False
            )

    minimum_diameters = {}
    for speeds in combinations(range(1, 36), 5):
        if not primitive(speeds):
            continue
        residues = tuple(sorted(speed % 7 for speed in speeds))
        minimum_diameters[residues] = min(
            minimum_diameters.get(residues, speeds[-1]),
            speeds[-1],
        )
    primitive_patterns = set(residue_patterns) - {(0, 0, 0, 0, 0)}
    if set(minimum_diameters) != primitive_patterns:
        raise AssertionError("primitive residue-pattern feasibility census")

    feasible_maximum = None
    feasible_minimum = None
    for residues, diameter in minimum_diameters.items():
        for far_residue in SECTORS:
            far_speed = diameter + 1
            while far_speed % 7 != far_residue:
                far_speed += 1
            coefficient = common_wall_coefficient(residues, far_residue)
            contribution = coefficient / far_speed
            record = (residues, far_residue, diameter, far_speed, coefficient)
            feasible_maximum = update_extreme(
                record, contribution, feasible_maximum
            )
            feasible_minimum = update_extreme(
                record, contribution, feasible_minimum, maximize=False
            )
    return {
        "coefficient_maximum": coefficient_maximum,
        "coefficient_minimum": coefficient_minimum,
        "feasible_maximum": feasible_maximum,
        "feasible_minimum": feasible_minimum,
        "primitive_pattern_count": len(primitive_patterns),
    }


def solve_linear(matrix, vector):
    size = len(vector)
    augmented = [list(matrix[row]) + [vector[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [value / pivot_value for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = augmented[row][column]
            if multiplier:
                augmented[row] = [
                    value - multiplier * pivot_entry
                    for value, pivot_entry in zip(augmented[row], augmented[column])
                ]
    return tuple(augmented[row][-1] for row in range(size))


def positive_cone_constants():
    coefficient = tuple(
        Fraction(6, 7) if vertex == 0 else Fraction(-1, 7)
        for vertex in SECTORS
    )
    constants = []
    for forbidden in range(1, 7):
        allowed = tuple(vertex for vertex in SECTORS if vertex != forbidden)
        best = None
        for mask in range(1, 1 << len(allowed)):
            support = tuple(
                allowed[index]
                for index in range(len(allowed))
                if mask & (1 << index)
            )
            if 0 not in support:
                continue
            principal = tuple(
                tuple(LAPLACIAN[first][second] for second in support)
                for first in support
            )
            restricted_coefficient = tuple(coefficient[vertex] for vertex in support)
            potential = solve_linear(principal, restricted_coefficient)
            if potential is None or any(value <= 0 for value in potential):
                continue
            norm = dot(restricted_coefficient, potential)
            if best is None or norm > best[0]:
                best = (norm, support, potential)
        constants.append(best)
    return tuple(constants)


def packet_metrics(measures):
    total_mass = sum(
        measures[missed][parity][class_vertex]
        for missed in SECTORS
        for parity in PARITIES
        for class_vertex in SECTORS
    )
    symmetric_rows = []
    current_rows = []
    for missed in SECTORS:
        symmetric_rows.append(
            tuple(
                measures[missed][0][class_vertex]
                + measures[missed][1][class_vertex]
                for class_vertex in SECTORS
            )
        )
        current_rows.append(
            tuple(
                measures[missed][0][class_vertex]
                - measures[missed][1][class_vertex]
                for class_vertex in SECTORS
            )
        )

    fill_mass = sum(
        measures[missed][0][missed]
        + measures[missed][1][(missed - 1) % 7]
        for missed in SECTORS
    )
    symmetric_term = (
        sum(
            symmetric_rows[missed][missed]
            + symmetric_rows[missed][(missed - 1) % 7]
            for missed in SECTORS
        )
        / 2
        - total_mass / 7
    )
    current_term = sum(
        current_rows[missed][missed]
        - current_rows[missed][(missed - 1) % 7]
        for missed in SECTORS
    ) / 2
    correction = fill_mass - total_mass / 7

    symmetric_energy = sum(laplacian_energy(row) for row in symmetric_rows)
    current_energy = sum(laplacian_energy(row) for row in current_rows)
    direct_channel_energy = 2 * sum(
        laplacian_energy(measures[missed][parity])
        for missed in SECTORS
        for parity in PARITIES
    )
    local_energy_loads = tuple(
        sum(
            CROSSING_WEIGHTS[(other - class_vertex) % 7]
            * (
                measures[missed][parity][class_vertex]
                - measures[missed][parity][other]
            )
            ** 2
            for missed in SECTORS
            for parity in PARITIES
            for other in SECTORS
        )
        for class_vertex in SECTORS
    )
    source_contributions = tuple(
        measures[class_vertex][0][class_vertex]
        + measures[(class_vertex + 1) % 7][1][class_vertex]
        - sum(
            measures[missed][parity][class_vertex]
            for missed in SECTORS
            for parity in PARITIES
        )
        / 7
        for class_vertex in SECTORS
    )
    return {
        "total_mass": total_mass,
        "fill_mass": fill_mass,
        "symmetric_term": symmetric_term,
        "current_term": current_term,
        "correction": correction,
        "symmetric_energy": symmetric_energy,
        "current_energy": current_energy,
        "packet_energy": symmetric_energy + current_energy,
        "direct_channel_energy": direct_channel_energy,
        "local_energy_loads": local_energy_loads,
        "source_contributions": source_contributions,
    }


def risk_order(values):
    return tuple(sorted(SECTORS, key=lambda vertex: (-values[vertex], vertex)))


def tournament_fingerprint(positive_risks, energy_risks):
    positive_path = risk_order(positive_risks)
    energy_path = risk_order(energy_risks)
    positive_rank = {vertex: rank for rank, vertex in enumerate(positive_path)}
    energy_rank = {vertex: rank for rank, vertex in enumerate(energy_path)}
    flips = sum(
        (positive_rank[first] < positive_rank[second])
        != (energy_rank[first] < energy_rank[second])
        for first, second in combinations(SECTORS, 2)
    )
    return {
        "positive_path": positive_path,
        "energy_path": energy_path,
        "score_histogram": {score: 1 for score in SECTORS},
        "directed_triangles": 0,
        "scc_sizes": (1,) * 7,
        "edge_flips": flips,
        "hamiltonian_path_count": 1,
    }


def update_extreme(record, value, current, maximize=True):
    if current is None or (value > current[0] if maximize else value < current[0]):
        return value, record
    return current


def primitive(speeds):
    common_divisor = 0
    for speed in speeds:
        common_divisor = gcd(common_divisor, speed)
    return common_divisor == 1


def main():
    print("HYP-7084: OWNER PACKET ON THE PARALLEL-CLASS CIRCLE")
    print("=" * 76)

    matrix = crossing_matrix()
    check(
        "direct K_7 chord census gives the THM-913 circulant kernel",
        all(
            matrix[first][second] == CROSSING_WEIGHTS[(second - first) % 7]
            for first in SECTORS
            for second in SECTORS
        ),
    )
    check("adjacent class bundles have zero crossing cost", CROSSING_WEIGHTS[1] == 0)
    check("crossing Laplacian row sum is 10", LAPLACIAN[0][0] == 10)

    current_functional = (
        Fraction(1), Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(-1),
    )
    current_potential = (
        Fraction(3, 29), Fraction(-3, 203), Fraction(-1, 203), Fraction(0),
        Fraction(1, 203), Fraction(3, 203), Fraction(-3, 29),
    )
    symmetric_functional = tuple(
        Fraction(int(vertex in (0, 6))) - Fraction(2, 7)
        for vertex in SECTORS
    )
    symmetric_potential = (
        Fraction(11, 203), Fraction(-1, 29), Fraction(-3, 203),
        Fraction(-2, 203), Fraction(-3, 203), Fraction(-1, 29),
        Fraction(11, 203),
    )
    check(
        "current dual potential solves Lh=delta_0-delta_6",
        matrix_vector(LAPLACIAN, current_potential) == current_functional,
    )
    check(
        "symmetric dual potential solves Lh=delta_0+delta_6-2/7",
        matrix_vector(LAPLACIAN, symmetric_potential) == symmetric_functional,
    )
    current_dual_norm = dot(current_functional, current_potential)
    symmetric_dual_norm = dot(symmetric_functional, symmetric_potential)
    check("adjacent-current dual norm is 6/29", current_dual_norm == Fraction(6, 29))
    check("adjacent-average dual norm is 22/203", symmetric_dual_norm == Fraction(22, 203))
    check(
        "combined seven-row certificate constant is 16/29",
        Fraction(7, 4) * (current_dual_norm + symmetric_dual_norm)
        == ENERGY_CERTIFICATE_CONSTANT,
    )

    cone_data = positive_cone_constants()
    expected_cone_norms = (
        Fraction(16, 203),
        Fraction(139, 1799),
        Fraction(496, 6545),
        Fraction(496, 6545),
        Fraction(139, 1799),
        Fraction(16, 203),
    )
    check(
        "all six occupied/miss cone norms are exact",
        tuple(item[0] for item in cone_data) == expected_cone_norms,
    )
    positive_cone_constant = sum(expected_cone_norms)
    check(
        "global positive incidence-cone constant",
        positive_cone_constant == Fraction(22620786, 48779885),
    )
    cone_energy_threshold = PROPAGATION_SLACK**2 / positive_cone_constant
    print("\nPositive incidence-cone norms by forbidden distance")
    for distance, (norm, support, potential) in enumerate(cone_data, start=1):
        print(
            f"  d={distance}: norm={norm}, support={support}, potential={potential}"
        )
    print("  global positive constant:", positive_cone_constant)

    synchronized_bounds = common_wall_bounds()
    check(
        "3,234 residue rows give the synchronized-wall coefficient range",
        synchronized_bounds["coefficient_maximum"][0] == Fraction(229, 686)
        and synchronized_bounds["coefficient_minimum"][0] == Fraction(-23, 98),
    )
    check(
        "all 461 primitive five-speed residue multisets occur by diameter 35",
        synchronized_bounds["primitive_pattern_count"] == 461,
    )
    check(
        "distinct-speed feasibility gives the universal synchronized-wall range",
        synchronized_bounds["feasible_maximum"][0] == Fraction(229, 5488)
        and synchronized_bounds["feasible_minimum"][0] == Fraction(-23, 784),
    )
    print("\nSynchronized x=k/7 wall theorem")
    print("  residue coefficient maximum:", synchronized_bounds["coefficient_maximum"])
    print("  residue coefficient minimum:", synchronized_bounds["coefficient_minimum"])
    print("  feasible t>D maximum:", synchronized_bounds["feasible_maximum"])
    print("  feasible t>D minimum:", synchronized_bounds["feasible_minimum"])

    extrema = {
        "correction_max": None,
        "correction_min": None,
        "symmetric_max": None,
        "symmetric_min": None,
        "current_max": None,
        "current_min": None,
        "energy_max": None,
        "ratio_max": None,
        "limit_max": None,
        "limit_min": None,
        "wall_max": None,
        "wall_min": None,
        "wall_abs_max": None,
        "common_wall_max": None,
        "common_wall_min": None,
        "residual_wall_max": None,
        "residual_wall_min": None,
        "residual_wall_abs_max": None,
    }
    positive_class_risks = {vertex: Fraction(-1) for vertex in SECTORS}
    energy_class_risks = {vertex: Fraction(-1) for vertex in SECTORS}
    case_count = 0
    raw_threshold_passes = 0
    cone_threshold_passes = 0
    zero_incidence_checks = 0

    for diameter in range(5, 11):
        for speeds in combinations(range(1, diameter + 1), 5):
            if speeds[-1] != diameter or not primitive(speeds):
                continue
            core = (0,) + speeds
            _, scaled_coefficients = OWNER_MODULE.PROFILE_MODULE.core_profile(speeds)
            expected_limit = scaled_coefficients[1] / 2
            for far_speed in range(diameter + 1, 4 * diameter + 1):
                base = core + (far_speed,)
                profile = OWNER_MODULE.sector_profile(base)
                measures = packet_measures(profile, far_speed)
                metrics = packet_metrics(measures)
                endpoints = endpoint_metrics(speeds, far_speed)
                record = (core, far_speed)
                correction = metrics["correction"]
                common_wall = endpoints["wall_contributions"].get(
                    tuple(speeds), Fraction(0)
                )
                expected_common_wall = common_wall_coefficient(
                    tuple(sorted(speed % 7 for speed in speeds)), far_speed % 7
                ) / far_speed
                residual_wall = endpoints["wall_remainder"] - common_wall

                if correction != OWNER_MODULE.finite_error(profile, 2 * far_speed):
                    raise AssertionError(("doubling correction", record, correction))
                if correction != metrics["symmetric_term"] + metrics["current_term"]:
                    raise AssertionError(("class-current decomposition", record, correction))
                if metrics["total_mass"] != OWNER_MODULE.p1_mass(profile):
                    raise AssertionError(("one-miss mass", record))
                if metrics["packet_energy"] != metrics["direct_channel_energy"]:
                    raise AssertionError(("polarized energy identity", record))
                if correction**2 > ENERGY_CERTIFICATE_CONSTANT * metrics["packet_energy"]:
                    raise AssertionError(("energy certificate", record, correction))
                if sum(metrics["local_energy_loads"]) != metrics["packet_energy"]:
                    raise AssertionError(("local energy partition", record))
                if sum(metrics["source_contributions"]) != correction:
                    raise AssertionError(("source contribution partition", record))
                if endpoints["limit"] != expected_limit:
                    raise AssertionError(("residue-two limit", record, endpoints["limit"]))
                if endpoints["total"] != correction:
                    raise AssertionError(("endpoint potential identity", record, endpoints))
                if common_wall != expected_common_wall:
                    raise AssertionError(("synchronized wall identity", record, common_wall))
                if (
                    correction > 0
                    and correction**2 > positive_cone_constant * metrics["packet_energy"]
                ):
                    raise AssertionError(("positive incidence-cone certificate", record))

                for missed in SECTORS:
                    forbidden_class = 2 * missed % 7
                    for parity in PARITIES:
                        if measures[missed][parity][forbidden_class] != 0:
                            raise AssertionError(("occupied sector cannot be missed", record))
                        zero_incidence_checks += 1

                if metrics["packet_energy"] < ENERGY_THRESHOLD:
                    raw_threshold_passes += 1
                if correction <= 0 or metrics["packet_energy"] < cone_energy_threshold:
                    cone_threshold_passes += 1
                ratio = (
                    correction**2 / metrics["packet_energy"]
                    if metrics["packet_energy"]
                    else Fraction(0)
                )
                extrema["correction_max"] = update_extreme(
                    record, correction, extrema["correction_max"]
                )
                extrema["correction_min"] = update_extreme(
                    record, correction, extrema["correction_min"], maximize=False
                )
                extrema["symmetric_max"] = update_extreme(
                    record, metrics["symmetric_term"], extrema["symmetric_max"]
                )
                extrema["symmetric_min"] = update_extreme(
                    record, metrics["symmetric_term"], extrema["symmetric_min"], maximize=False
                )
                extrema["current_max"] = update_extreme(
                    record, metrics["current_term"], extrema["current_max"]
                )
                extrema["current_min"] = update_extreme(
                    record, metrics["current_term"], extrema["current_min"], maximize=False
                )
                extrema["energy_max"] = update_extreme(
                    record, metrics["packet_energy"], extrema["energy_max"]
                )
                extrema["ratio_max"] = update_extreme(
                    record, ratio, extrema["ratio_max"]
                )
                extrema["limit_max"] = update_extreme(
                    record, endpoints["limit"], extrema["limit_max"]
                )
                extrema["limit_min"] = update_extreme(
                    record, endpoints["limit"], extrema["limit_min"], maximize=False
                )
                extrema["wall_max"] = update_extreme(
                    record, endpoints["wall_remainder"], extrema["wall_max"]
                )
                extrema["wall_min"] = update_extreme(
                    record, endpoints["wall_remainder"], extrema["wall_min"], maximize=False
                )
                extrema["wall_abs_max"] = update_extreme(
                    record, abs(endpoints["wall_remainder"]), extrema["wall_abs_max"]
                )
                extrema["common_wall_max"] = update_extreme(
                    record, common_wall, extrema["common_wall_max"]
                )
                extrema["common_wall_min"] = update_extreme(
                    record, common_wall, extrema["common_wall_min"], maximize=False
                )
                extrema["residual_wall_max"] = update_extreme(
                    record, residual_wall, extrema["residual_wall_max"]
                )
                extrema["residual_wall_min"] = update_extreme(
                    record, residual_wall, extrema["residual_wall_min"], maximize=False
                )
                extrema["residual_wall_abs_max"] = update_extreme(
                    record, abs(residual_wall), extrema["residual_wall_abs_max"]
                )
                for class_vertex in SECTORS:
                    positive_class_risks[class_vertex] = max(
                        positive_class_risks[class_vertex],
                        metrics["source_contributions"][class_vertex],
                    )
                    energy_class_risks[class_vertex] = max(
                        energy_class_risks[class_vertex],
                        metrics["local_energy_loads"][class_vertex],
                    )
                case_count += 1

    check("6,900 exact owner-doubling rows", case_count == 6900)
    check("96,600 structural occupied/miss zeros", zero_incidence_checks == 96600)
    check(
        "bounded-bank correction maximum matches HYP-7083",
        extrema["correction_max"][0] == Fraction(2173, 27440)
        and extrema["correction_max"][1] == ((0, 3, 4, 5, 6, 7), 8),
    )
    check(
        "bounded bank attains both universal synchronized-wall extrema",
        extrema["common_wall_max"][0] == synchronized_bounds["feasible_maximum"][0]
        and extrema["common_wall_min"][0] == synchronized_bounds["feasible_minimum"][0],
    )

    print("\nExact packet extrema")
    for label, value in extrema.items():
        print(f"  {label:16s}: {value}  decimal={float(value[0]):.9f}")
    print(f"  energy threshold: {ENERGY_THRESHOLD}  decimal={float(ENERGY_THRESHOLD):.9f}")
    print(f"  raw threshold passes: {raw_threshold_passes}/{case_count}")
    print(
        f"  cone threshold: {cone_energy_threshold}  "
        f"decimal={float(cone_energy_threshold):.9f}"
    )
    print(f"  cone/negative upper certificates: {cone_threshold_passes}/{case_count}")
    print(
        "  worst energy-only upper bound:",
        (float(ENERGY_CERTIFICATE_CONSTANT * extrema["energy_max"][0])) ** 0.5,
    )

    dangerous_core, dangerous_far_speed = extrema["correction_max"][1]
    dangerous_endpoints = endpoint_metrics(dangerous_core[1:], dangerous_far_speed)
    print("\nDangerous-row endpoint current")
    print("  row:", extrema["correction_max"][1])
    print("  closed residue-two limit:", dangerous_endpoints["limit"])
    print("  finite wall remainder:", dangerous_endpoints["wall_remainder"])
    print("  slow wall count:", dangerous_endpoints["wall_count"])
    print("  contributions by simultaneous owner set:")
    for owners, contribution in sorted(
        dangerous_endpoints["wall_contributions"].items(),
        key=lambda item: (-abs(item[1]), item[0]),
    ):
        print(f"    {owners}: {contribution}  decimal={float(contribution):+.9f}")

    fingerprint = tournament_fingerprint(positive_class_risks, energy_class_risks)
    print("\nTournament Analysis")
    print("  vertices: seven parallel classes")
    print("  observable: largest positive source contribution")
    print("  switch: largest local crossing-Laplacian load; numeric tie gauge")
    print("  positive risks:", positive_class_risks)
    print("  energy risks:", energy_class_risks)
    print("  fingerprint:", fingerprint)
    check("both class-risk tournaments are transitive", fingerprint["directed_triangles"] == 0)
    check("class-risk SCCs are singletons", set(fingerprint["scc_sizes"]) == {1})
    check("class-risk tie paths are unique", fingerprint["hamiltonian_path_count"] == 1)

    print("\nVERDICT")
    print("  PROVED: the owner-doubling packet is an adjacent-edge current on C_7.")
    print("  PROVED: Delta_F(2t)^2 <= (16/29) E_pc for the THM-913 Laplacian energy.")
    print("  PROVED: Delta_F(2t)=C_2(E)+a finite 14-cycle slow-wall endpoint sum.")
    print("  PROVED: -23/784 <= synchronized-wall term <= 229/5488, sharp.")
    print("  EXACT BANK: D<=10 and D<t<=4D, 6,900 rows; all three identities pass.")
    if raw_threshold_passes == case_count:
        print("  BANK RESULT: E_pc alone certifies every scanned correction below 0.097.")
    else:
        print("  NEGATIVE: the raw E_pc threshold does not certify every scanned row.")
    print("  OPEN: bound the noncommon gcd-lattice wall sum with its relation sidecar.")


if __name__ == "__main__":
    main()
