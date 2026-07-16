#!/usr/bin/env python3
"""
Exact referee for THM-891 and HYP-7024.

For a slow six-offset core E containing the stationary offset 0, append a far
offset t and evaluate the THM-727 two-scale error at w = a*t.  Conditioning on
the source section c of the fast offset gives the microcell integral

    J_a(c,s) = integral_[c/7,(c+1)/7] g_s(a*y) dy,

where g_s = 1_[s/7,(s+1)/7) - 1/7.  If N_a(c,s) counts the a microcells whose
image section is s, then

    a*J_a(c,s) = N_a(c,s)/7 - a/49.

It depends only on a modulo 7.  The limiting scaled coefficient a*C_a(E) is a
signed expectation over the slow core's one- and two-miss patterns.  This file
computes that coefficient exactly, compares it with finite far peels, and scans
all primitive five-speed cores in [1,20].

Tournament Analysis is diagnostic rather than the proof engine here.  Vertices
are the six nonzero owner-resonance residues.  On each core the pairwise
observable is |r*C_r|-|q*C_q|; the aggregate switch/gauge orients r -> q when a
strict majority of primitive cores prefer r, with lexicographic tie gauge.  The
tie Hamiltonian path is the lexicographically first directed Hamiltonian path.
This quotient preserves the resonant error ordering but destroys wall order,
finite-t aliasing, and the arithmetic identity of the core.  Alternative vertex
sets explicitly considered were runners, swing arcs, fixed sections, section
boundaries, wall-crossing events, miss patterns, Fourier modes, and proof
obligations.  Miss patterns are the faithful coefficient carrier; residues are
used as tournament vertices only because runners and arcs do not define a
canonical binary relation for this scalar extremal problem.
"""

from fractions import Fraction
from heapq import heapify, heappop, heappush
from itertools import combinations, permutations
from math import gcd, lcm


RESIDUES = tuple(range(1, 7))
INNER_SECTIONS = tuple(range(1, 7))
EXPECTED_CONSECUTIVE = (
    Fraction(-239, 5145),
    Fraction(209, 20580),
    Fraction(-39, 6860),
    Fraction(47, 20580),
    Fraction(-4, 1715),
    Fraction(83, 6860),
)
EXPECTED_EXTREMAL = (
    Fraction(-16, 343),
    Fraction(3, 343),
    Fraction(-1, 2058),
    Fraction(17, 1372),
    Fraction(13, 4116),
    Fraction(31, 4116),
)


def scaled_microcell_numerator(residue, source_section, target_section):
    extra = any(
        (source_section * residue + step) % 7 == target_section
        for step in range(residue)
    )
    return 7 * int(extra) - residue


def pattern_kernel_numerator(residue, missed):
    if len(missed) == 1:
        target = missed[0]
        return sum(
            scaled_microcell_numerator(residue, source, target)
            for source in range(7)
            if source != target
        )
    first, second = missed
    return scaled_microcell_numerator(
        residue, second, first
    ) + scaled_microcell_numerator(residue, first, second)


KERNEL_NUMERATORS = {
    (residue, missed): pattern_kernel_numerator(residue, missed)
    for residue in RESIDUES
    for size in (1, 2)
    for missed in combinations(INNER_SECTIONS, size)
}


def core_profile(positive_speeds):
    period_scale = lcm(*positive_speeds)
    sectors = [0] * len(positive_speeds)
    counts = [len(positive_speeds) + 1, 0, 0, 0, 0, 0, 0]
    events = []
    for runner_index, speed in enumerate(positive_speeds):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    previous = 0
    pattern_numerators = {}
    coefficient_numerators = [0] * 6
    while events:
        event_position = events[0][0]
        interval_length = event_position - previous
        missed = tuple(section for section in INNER_SECTIONS if counts[section] == 0)
        pattern_numerators[missed] = (
            pattern_numerators.get(missed, 0) + interval_length
        )
        if len(missed) in (1, 2):
            for residue in RESIDUES:
                coefficient_numerators[residue - 1] += (
                    interval_length * KERNEL_NUMERATORS[residue, missed]
                )
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            old_sector = sectors[runner_index]
            new_sector = (old_sector + 1) % 7
            counts[old_sector] -= 1
            counts[new_sector] += 1
            sectors[runner_index] = new_sector
            if event_index < 7 * speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (
                        next_index * event_step,
                        runner_index,
                        next_index,
                        speed,
                        event_step,
                    ),
                )
        previous = event_position
    measure_denominator = 7 * period_scale
    pattern_measures = {
        missed: Fraction(numerator, measure_denominator)
        for missed, numerator in pattern_numerators.items()
    }
    scaled_coefficients = tuple(
        Fraction(numerator, 49 * measure_denominator)
        for numerator in coefficient_numerators
    )
    return pattern_measures, scaled_coefficients


def missed_count_distribution(pattern_measures):
    distribution = {}
    for missed, measure in pattern_measures.items():
        distribution[len(missed)] = distribution.get(len(missed), Fraction(0)) + measure
    return distribution


def exact_swing_arcs(offsets, target_section):
    breakpoints = sorted(
        {
            Fraction(index, 7 * speed)
            for speed in offsets
            if speed > 0
            for index in range(7 * speed)
        }
        | {Fraction(0), Fraction(1)}
    )
    arcs = []
    active = False
    start = None
    for left, right in zip(breakpoints, breakpoints[1:]):
        midpoint = (left + right) / 2
        occupied = {
            int((speed * midpoint % 1) * 7) if speed > 0 else 0
            for speed in offsets
        }
        current = len(occupied) == 6 and target_section not in occupied
        if current and not active:
            start = left
            active = True
        if active and not current:
            arcs.append((start, left))
            active = False
    if active:
        arcs.append((start, Fraction(1)))
    return arcs


def centered_section_primitive(value, target_section):
    fractional = value % 1
    left = Fraction(target_section, 7)
    right = Fraction(target_section + 1, 7)
    intersection = max(Fraction(0), min(fractional, right) - left)
    return intersection - fractional / 7


def finite_error(offsets, peel):
    scaled_error = Fraction(0)
    for target_section in range(7):
        for left, right in exact_swing_arcs(offsets, target_section):
            scaled_error += centered_section_primitive(
                peel * right, target_section
            ) - centered_section_primitive(peel * left, target_section)
    return scaled_error / peel


def pair_collision_mass(first_speed, second_speed):
    positive = tuple(speed for speed in (first_speed, second_speed) if speed > 0)
    period_scale = lcm(*positive) if positive else 1
    states = [0, 0]
    events = []
    for runner_index, speed in enumerate((first_speed, second_speed)):
        if speed == 0:
            continue
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    previous = 0
    collision_numerator = 0
    while events:
        event_position = events[0][0]
        if states[0] == states[1]:
            collision_numerator += event_position - previous
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            states[runner_index] = (states[runner_index] + 1) % 7
            if event_index < 7 * speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (
                        next_index * event_step,
                        runner_index,
                        next_index,
                        speed,
                        event_step,
                    ),
                )
        previous = event_position
    return Fraction(collision_numerator, 7 * period_scale)


def primitive(positive_speeds):
    common_divisor = 0
    for speed in positive_speeds:
        common_divisor = gcd(common_divisor, speed)
    return common_divisor == 1


def empty_tournament_counts():
    return {
        (first, second): [0, 0, 0]
        for first in RESIDUES
        for second in RESIDUES
        if first < second
    }


def update_tournament_counts(counts, scaled_coefficients):
    magnitudes = tuple(abs(value) for value in scaled_coefficients)
    for first in RESIDUES:
        for second in RESIDUES:
            if first >= second:
                continue
            if magnitudes[first - 1] > magnitudes[second - 1]:
                counts[first, second][0] += 1
            elif magnitudes[second - 1] > magnitudes[first - 1]:
                counts[first, second][1] += 1
            else:
                counts[first, second][2] += 1


def majority_edges(counts):
    edges = set()
    for (first, second), (first_wins, second_wins, _) in counts.items():
        winner = first if first_wins >= second_wins else second
        loser = second if winner == first else first
        edges.add((winner, loser))
    return edges


def tournament_fingerprint(counts):
    edges = majority_edges(counts)
    scores = {vertex: 0 for vertex in RESIDUES}
    for winner, _ in edges:
        scores[winner] += 1
    directed_triangles = 0
    for triple in combinations(RESIDUES, 3):
        induced_scores = {
            vertex: sum(
                (vertex, other) in edges for other in triple if other != vertex
            )
            for vertex in triple
        }
        directed_triangles += sorted(induced_scores.values()) == [1, 1, 1]
    reachability = {vertex: {vertex} for vertex in RESIDUES}
    for winner, loser in edges:
        reachability[winner].add(loser)
    for pivot in RESIDUES:
        for source in RESIDUES:
            if pivot in reachability[source]:
                reachability[source] |= reachability[pivot]
    components = []
    unused = set(RESIDUES)
    while unused:
        vertex = min(unused)
        component = tuple(
            candidate
            for candidate in RESIDUES
            if candidate in reachability[vertex]
            and vertex in reachability[candidate]
        )
        components.append(component)
        unused -= set(component)
    paths = [
        path
        for path in permutations(RESIDUES)
        if all((path[index], path[index + 1]) in edges for index in range(5))
    ]
    return {
        "edges": edges,
        "scores": scores,
        "score_histogram": {
            score: tuple(sorted(vertex for vertex, value in scores.items() if value == score))
            for score in sorted(set(scores.values()))
        },
        "directed_triangles": directed_triangles,
        "components": tuple(components),
        "hamiltonian_path_count": len(paths),
        "tie_hamiltonian_path": min(paths),
    }


def edge_flip_count(first_edges, second_edges):
    return sum(
        ((first, second) in first_edges) != ((first, second) in second_edges)
        for first in RESIDUES
        for second in RESIDUES
        if first < second
    )


def scan_cores(diameter):
    records = []
    tournament_counts = empty_tournament_counts()
    nested_counts = empty_tournament_counts()
    winning_residues = {residue: 0 for residue in RESIDUES}
    nested_winning_residues = {residue: 0 for residue in RESIDUES}
    residue_maxima = {
        residue: (Fraction(-1), None, None) for residue in RESIDUES
    }
    for positive_speeds in combinations(range(1, diameter + 1), 5):
        if not primitive(positive_speeds):
            continue
        pattern_measures, scaled_coefficients = core_profile(positive_speeds)
        core = (0,) + positive_speeds
        magnitudes = tuple(abs(value) for value in scaled_coefficients)
        winning_residue = min(
            RESIDUES, key=lambda residue: (-magnitudes[residue - 1], residue)
        )
        winning_residues[winning_residue] += 1
        update_tournament_counts(tournament_counts, scaled_coefficients)
        if positive_speeds[-1] <= 14:
            nested_winning_residues[winning_residue] += 1
            update_tournament_counts(nested_counts, scaled_coefficients)
        for residue in RESIDUES:
            candidate = magnitudes[residue - 1]
            if candidate > residue_maxima[residue][0]:
                residue_maxima[residue] = (
                    candidate,
                    core,
                    scaled_coefficients[residue - 1],
                )
        distribution = missed_count_distribution(pattern_measures)
        records.append(
            (
                max(magnitudes),
                core,
                winning_residue,
                scaled_coefficients[winning_residue - 1],
                distribution.get(1, Fraction(0)),
                distribution.get(2, Fraction(0)),
            )
        )
    return {
        "records": records,
        "tournament_counts": tournament_counts,
        "nested_counts": nested_counts,
        "winning_residues": winning_residues,
        "nested_winning_residues": nested_winning_residues,
        "residue_maxima": residue_maxima,
    }


def format_fraction(value):
    return f"{value} = {float(value):.12f}"


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def main():
    print("THM-891 / HYP-7024: EXACT RESONANT CROSS-SECTION CANCELLATION")
    print("=" * 78)
    print("\n[1] Seven-microcell identity")
    for multiplier in range(1, 71):
        residue = multiplier % 7
        if residue == 0:
            residue = 7
        for source in range(7):
            for target in range(7):
                count = sum(
                    index % 7 == target
                    for index in range(source * multiplier, (source + 1) * multiplier)
                )
                direct = Fraction(count, 7) - Fraction(multiplier, 49)
                expected = (
                    Fraction(0)
                    if residue == 7
                    else Fraction(
                        scaled_microcell_numerator(residue, source, target), 49
                    )
                )
                check_label = "microcell formula" if multiplier == source == target == 1 else None
                if direct != expected:
                    raise AssertionError((multiplier, source, target, direct, expected))
                if check_label:
                    print("PASS  microcell formula over a<=70, all 49 section pairs")
    print("      Therefore a*J_a(c,s) depends only on a mod 7 and vanishes at 0 mod 7.")

    print("\n[2] Exact slow-core coefficient tables")
    consecutive_patterns, consecutive_coefficients = core_profile((1, 2, 3, 4, 5))
    extremal_patterns, extremal_coefficients = core_profile((1, 2, 3, 4, 6))
    check("consecutive residue table", consecutive_coefficients == EXPECTED_CONSECUTIVE)
    check("perforated extremal residue table", extremal_coefficients == EXPECTED_EXTREMAL)
    for name, patterns, coefficients in (
        ("consecutive {0,1,2,3,4,5}", consecutive_patterns, consecutive_coefficients),
        ("extremal {0,1,2,3,4,6}", extremal_patterns, extremal_coefficients),
    ):
        distribution = missed_count_distribution(patterns)
        print(f"  {name}")
        print(f"    p1={distribution[1]}, p2={distribution[2]}, 3p1+p2={3*distribution[1]+distribution[2]}")
        for residue, coefficient in zip(RESIDUES, coefficients):
            print(f"    r={residue}: r*C_r = {format_fraction(coefficient)}")
        check(
            f"binding r=1 identity for {name}",
            coefficients[0] == -2 * (3 * distribution[1] + distribution[2]) / 49,
        )
    check(
        "sharp candidate equality 3p1+p2=8/7",
        3 * missed_count_distribution(extremal_patterns)[1]
        + missed_count_distribution(extremal_patterns)[2]
        == Fraction(8, 7),
    )
    check("candidate constant is below 0.097", Fraction(16, 343) < Fraction(97, 1000))

    print("\n[3] Finite far-peel referee for the consecutive core, w=t")
    limiting_error = EXPECTED_CONSECUTIVE[0]
    finite_rows = []
    for far_speed in (25, 50, 100, 200, 400):
        error = finite_error((0, 1, 2, 3, 4, 5, far_speed), far_speed)
        finite_rows.append((far_speed, error, error - limiting_error))
        print(
            f"  t={far_speed:3d}: Error(t)={format_fraction(error)}, "
            f"Error-C1={format_fraction(error-limiting_error)}"
        )
    check(
        "finite values converge at the predicted O(1/t) scale on the audit ladder",
        max(far_speed * abs(difference) for far_speed, _, difference in finite_rows)
        < Fraction(1, 7),
    )

    print("\n[4] Refuted collision shortcut")
    collision_18 = pair_collision_mass(1, 8)
    collision_29 = pair_collision_mass(2, 9)
    print(f"  collision mass (1,8) = {collision_18}")
    print(f"  collision mass (2,9) = {collision_29}")
    check("pair-collision mass is not universally 1/7", collision_18 == Fraction(1, 4) and collision_29 == Fraction(2, 9))

    print("\n[5] Exhaustive primitive-core scan through diameter 20")
    scan = scan_cores(20)
    records = sorted(scan["records"], reverse=True)
    check("primitive core count", len(records) == 15246)
    check("sharp scanned maximum", records[0][0] == Fraction(16, 343))
    check("unique scanned maximizer", sum(record[0] == Fraction(16, 343) for record in records) == 1)
    check("scanned maximizer core", records[0][1] == (0, 1, 2, 3, 4, 6))
    print("  top ten (max |r*C_r|, core, winning residue, signed coefficient, p1, p2)")
    for record in records[:10]:
        print("   ", record)
    print("  per-residue exact maxima")
    for residue in RESIDUES:
        print(f"    r={residue}: {scan['residue_maxima'][residue]}")
    print("  winning residues through D=14:", scan["nested_winning_residues"])
    print("  winning residues through D=20:", scan["winning_residues"])
    check(
        "residue-1 pointwise dominance is false",
        scan["winning_residues"][2] > 0 and scan["winning_residues"][3] > 0,
    )

    print("\n[6] Tournament Analysis on resonance residues")
    nested_fingerprint = tournament_fingerprint(scan["nested_counts"])
    full_fingerprint = tournament_fingerprint(scan["tournament_counts"])
    print("  observable: |r*C_r|-|q*C_q| per core")
    print("  switch/gauge: strict primitive-core majority; lexicographic tie orientation")
    print("  D=14 fingerprint:", nested_fingerprint)
    print("  D=20 fingerprint:", full_fingerprint)
    flips = edge_flip_count(
        nested_fingerprint["edges"], full_fingerprint["edges"]
    )
    print("  edge flips D=14 -> D=20:", flips)
    check("majority tournament is stable", flips == 0)
    check("majority tournament has no directed 3-cycle", full_fingerprint["directed_triangles"] == 0)
    check("majority tournament SCCs are singletons", all(len(component) == 1 for component in full_fingerprint["components"]))
    check("tie Hamiltonian path is unique", full_fingerprint["hamiltonian_path_count"] == 1)
    check("stable tie Hamiltonian path", full_fingerprint["tie_hamiltonian_path"] == (1, 2, 3, 4, 6, 5))
    print("  Assumption challenge: residues preserve the resonant LRC error ordering, while")
    print("  destroying finite-t aliases, wall chronology, and core arithmetic. Miss patterns,")
    print("  not runners or arcs, are the faithful carrier of the exact coefficient.")

    print("\nVERDICT")
    print("  PROVED: exact cross-section limit kernel and seven-residue reduction.")
    print("  PROVED: consecutive-core max coefficient 239/5145 < 0.097.")
    print("  FINITE-EXACT: 15,246 primitive cores through D=20; unique max 16/343.")
    print("  REFUTED: fixed pair-collision first moment and universal residue-1 dominance.")
    print("  OPEN CRUX: prove all-core max_r |r*C_r| <= 16/343, including the")
    print("  separate residue-2 through residue-6 miss-pattern inequalities.")


if __name__ == "__main__":
    main()
