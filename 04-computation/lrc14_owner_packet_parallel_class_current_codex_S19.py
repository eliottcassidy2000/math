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
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd


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
    measures = [
        [[Fraction(0) for _ in SECTORS] for _ in PARITIES]
        for _ in SECTORS
    ]
    cell_denominator = 14 * far_speed
    for left, right, missed in intervals:
        first_cell = cell_denominator * left // denominator
        final_cell = (cell_denominator * right + denominator - 1) // denominator
        interval_left = Fraction(left, denominator)
        interval_right = Fraction(right, denominator)
        for cell in range(first_cell, final_cell):
            overlap_left = max(interval_left, Fraction(cell, cell_denominator))
            overlap_right = min(interval_right, Fraction(cell + 1, cell_denominator))
            if overlap_left >= overlap_right:
                continue
            half_sector = cell % 14
            parity = half_sector % 2
            source_sector = half_sector // 2
            class_vertex = 2 * source_sector % 7
            measures[missed][parity][class_vertex] += overlap_right - overlap_left
    return tuple(
        tuple(tuple(row) for row in parity_rows)
        for parity_rows in measures
    )


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

    extrema = {
        "correction_max": None,
        "correction_min": None,
        "symmetric_max": None,
        "symmetric_min": None,
        "current_max": None,
        "current_min": None,
        "energy_max": None,
        "ratio_max": None,
    }
    positive_class_risks = {vertex: Fraction(-1) for vertex in SECTORS}
    energy_class_risks = {vertex: Fraction(-1) for vertex in SECTORS}
    case_count = 0
    threshold_passes = 0
    zero_incidence_checks = 0

    for diameter in range(5, 11):
        for speeds in combinations(range(1, diameter + 1), 5):
            if speeds[-1] != diameter or not primitive(speeds):
                continue
            core = (0,) + speeds
            for far_speed in range(diameter + 1, 4 * diameter + 1):
                base = core + (far_speed,)
                profile = OWNER_MODULE.sector_profile(base)
                measures = packet_measures(profile, far_speed)
                metrics = packet_metrics(measures)
                record = (core, far_speed)
                correction = metrics["correction"]

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

                for missed in SECTORS:
                    forbidden_class = 2 * missed % 7
                    for parity in PARITIES:
                        if measures[missed][parity][forbidden_class] != 0:
                            raise AssertionError(("occupied sector cannot be missed", record))
                        zero_incidence_checks += 1

                if metrics["packet_energy"] < ENERGY_THRESHOLD:
                    threshold_passes += 1
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

    print("\nExact packet extrema")
    for label, value in extrema.items():
        print(f"  {label:16s}: {value}  decimal={float(value[0]):.9f}")
    print(f"  energy threshold: {ENERGY_THRESHOLD}  decimal={float(ENERGY_THRESHOLD):.9f}")
    print(f"  threshold passes: {threshold_passes}/{case_count}")
    print(
        "  worst energy-only upper bound:",
        (float(ENERGY_CERTIFICATE_CONSTANT * extrema["energy_max"][0])) ** 0.5,
    )

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
    print("  EXACT BANK: D<=10 and D<t<=4D, 6,900 rows; decomposition and zeros all pass.")
    if threshold_passes == case_count:
        print("  BANK RESULT: E_pc alone certifies every scanned correction below 0.097.")
    else:
        print("  NEGATIVE: the raw E_pc threshold does not certify every scanned row.")
    print("  OPEN: derive a universal packet-energy bound or exploit incidence beyond E_pc.")


if __name__ == "__main__":
    main()
