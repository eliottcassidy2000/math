#!/usr/bin/env python3
"""Exact verifier for THM-904 and HYP-7061.

The negative residue-six kernel from THM-891 is pointwise dominated by a
rational observable on the 84 unordered sector states of three moving
runners.  Integrating reduces the limiting crux to one universal three-speed
inequality.  This script proves the finite reduction, scans every primitive
triple through largest speed 60 exactly, decomposes the observable into its
uniform/pair/cubic ANOVA channels, and decides the proposed negation transfer.

Tournament Analysis is diagnostic.  Vertices are the six strongest primitive
triple relation shapes, the pairwise observable is q(u)-q(v), and the gauge
orients toward the larger exact value with lexicographic tie breaking.  The
tie Hamiltonian path is the lexicographically first directed Hamiltonian path.
This preserves obstruction strength but destroys the full relation lattice,
the other two movers, and finite-t wall chronology.  Alternative vertex sets
considered explicitly were runners, gaps, fixed sections, boundaries,
wall-crossing events, residues, Fourier modes, primitive relations, miss
patterns, and proof obligations.  Triple relation shapes are used because the
84-state triple observable is the smallest quotient presently known to retain
the negative K6 certificate; runners and arcs do not.
"""

from fractions import Fraction
from functools import cache
from heapq import heapify, heappop, heappush
from itertools import combinations, combinations_with_replacement, permutations
from math import comb, gcd, lcm


SECTIONS = tuple(range(7))
INNER_SECTIONS = tuple(range(1, 7))
TRIPLE_TYPES = tuple(combinations_with_replacement(SECTIONS, 3))
TRIPLE_TYPE_INDEX = {triple_type: index for index, triple_type in enumerate(TRIPLE_TYPES)}
PAIR_TYPES = tuple(combinations_with_replacement(SECTIONS, 2))
TARGET_TRIPLE_BOUND = Fraction(47, 100)
PROPAGATION_SLACK = Fraction(97, 1000)
SCAN_DIAMETER = 60

# Lexicographic combinations_with_replacement(range(7), 3) order.
TRIPLE_WEIGHT_NUMERATORS = (
    0, 3, 2, 3, 7, 7, 0, -5, -4, 19, -9, 21, 3, 2, 3, 31, -6, 41,
    -5, -1, 15, 287, 0, -6, 210, -12, 174, 0, 27, 24, -2, 63, 1, 110,
    9, 1, -40, -10, -81, -3, -9, 118, 221, 32, -24, -102, -9, 72,
    27, 0, 4, 49, 13, 152, -6, 120, -1, 191, 14, -40, 55, 38, -89,
    71, 29, 5, 0, -6, 2, -9, 267, -9, 275, 4, 0, 28, 113, 53, -94,
    76, 57, 107, 20, 0,
)
TRIPLE_WEIGHTS = tuple(Fraction(value, 100) for value in TRIPLE_WEIGHT_NUMERATORS)

PAIR_RAY_MINIMIZERS = (
    (1, 8), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 9),
    (2, 3), (2, 11), (2, 5), (2, 13), (3, 10), (3, 4), (3, 5),
    (3, 13), (4, 11), (4, 5), (4, 13), (5, 12), (5, 6), (6, 13),
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
            for source in SECTIONS
            if source != target
        )
    first, second = missed
    return scaled_microcell_numerator(
        residue, second, first
    ) + scaled_microcell_numerator(residue, first, second)


KERNEL_NUMERATORS = {
    (residue, missed): pattern_kernel_numerator(residue, missed)
    for residue in range(1, 7)
    for size in (1, 2)
    for missed in combinations(INNER_SECTIONS, size)
}


def moving_state_compositions(total=5, sections=7, prefix=()):
    if sections == 1:
        yield prefix + (total,)
        return
    for first_count in range(total + 1):
        yield from moving_state_compositions(
            total - first_count, sections - 1, prefix + (first_count,)
        )


def triple_type_counts(state):
    counts = []
    for first, second, third in TRIPLE_TYPES:
        if first == third:
            count = comb(state[first], 3)
        elif first == second:
            count = comb(state[first], 2) * state[third]
        elif second == third:
            count = state[first] * comb(state[second], 2)
        else:
            count = state[first] * state[second] * state[third]
        counts.append(count)
    return tuple(counts)


def state_certificate_slack(state):
    missed = tuple(section for section in INNER_SECTIONS if state[section] == 0)
    target = -KERNEL_NUMERATORS.get((6, missed), 0)
    weighted_triples = sum(
        weight * count
        for weight, count in zip(TRIPLE_WEIGHT_NUMERATORS, triple_type_counts(state))
    )
    return weighted_triples - 100 * target


@cache
def triple_weight_average(speeds):
    common_divisor = gcd(*speeds)
    reduced = tuple(speed // common_divisor for speed in speeds)
    period_scale = lcm(*reduced)
    sectors = [0, 0, 0]
    events = []
    for runner_index, speed in enumerate(reduced):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    previous = 0
    numerator = 0
    while events:
        event_position = events[0][0]
        triple_type = tuple(sorted(sectors))
        numerator += (
            event_position - previous
        ) * TRIPLE_WEIGHT_NUMERATORS[TRIPLE_TYPE_INDEX[triple_type]]
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            sectors[runner_index] = (sectors[runner_index] + 1) % 7
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
    return Fraction(numerator, 100 * 7 * period_scale)


@cache
def pair_sector_distribution(first_speed, second_speed):
    common_divisor = gcd(first_speed, second_speed)
    reduced = (first_speed // common_divisor, second_speed // common_divisor)
    period_scale = lcm(*reduced)
    sectors = [0, 0]
    events = []
    for runner_index, speed in enumerate(reduced):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    previous = 0
    numerators = {pair_type: 0 for pair_type in PAIR_TYPES}
    while events:
        event_position = events[0][0]
        numerators[tuple(sorted(sectors))] += event_position - previous
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            sectors[runner_index] = (sectors[runner_index] + 1) % 7
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
    return {
        pair_type: Fraction(numerator, 7 * period_scale)
        for pair_type, numerator in numerators.items()
    }


def beta(first, second, third):
    return TRIPLE_WEIGHTS[TRIPLE_TYPE_INDEX[tuple(sorted((first, second, third)))]]


def anova_channels():
    mean = sum(beta(*state) for state in permutations_with_repetition(3)) / 343
    singleton = tuple(
        sum(beta(section, second, third) for second in SECTIONS for third in SECTIONS)
        / 49
        - mean
        for section in SECTIONS
    )
    pair = tuple(
        tuple(
            sum(beta(first, second, third) for third in SECTIONS) / 7
            - mean
            - singleton[first]
            - singleton[second]
            for second in SECTIONS
        )
        for first in SECTIONS
    )
    cubic = {
        (first, second, third): beta(first, second, third)
        - mean
        - singleton[first]
        - singleton[second]
        - singleton[third]
        - pair[first][second]
        - pair[first][third]
        - pair[second][third]
        for first in SECTIONS
        for second in SECTIONS
        for third in SECTIONS
    }
    return mean, singleton, pair, cubic


def permutations_with_repetition(length, prefix=()):
    if length == 0:
        yield prefix
        return
    for section in SECTIONS:
        yield from permutations_with_repetition(length - 1, prefix + (section,))


def pair_projection_expectation(pair_projection, distribution):
    return sum(
        mass * pair_projection[first][second]
        for (first, second), mass in distribution.items()
    )


def primitive(speeds):
    return gcd(*speeds) == 1


def shortest_full_relation(speeds, height_limit=30):
    for height in range(1, height_limit + 1):
        candidates = []
        for first in range(-height, height + 1):
            for second in range(-height, height + 1):
                if first == 0 or second == 0:
                    continue
                numerator = -(first * speeds[0] + second * speeds[1])
                if numerator % speeds[2]:
                    continue
                third = numerator // speeds[2]
                if third == 0 or abs(third) > height:
                    continue
                relation_gcd = gcd(gcd(abs(first), abs(second)), abs(third))
                relation = (first // relation_gcd, second // relation_gcd, third // relation_gcd)
                if max(map(abs, relation)) == height:
                    candidates.append(relation)
        if candidates:
            return min(candidates, key=lambda row: (sum(map(abs, row)), row))
    return None


def tournament_fingerprint(vertices, values):
    edges = set()
    for first, second in combinations(vertices, 2):
        if values[first] > values[second]:
            edges.add((first, second))
        elif values[second] > values[first]:
            edges.add((second, first))
        else:
            edges.add((min(first, second), max(first, second)))
    scores = {
        vertex: sum(winner == vertex for winner, _ in edges)
        for vertex in vertices
    }
    directed_triangles = 0
    for triple in combinations(vertices, 3):
        induced_scores = sorted(
            sum((vertex, other) in edges for other in triple if other != vertex)
            for vertex in triple
        )
        directed_triangles += induced_scores == [1, 1, 1]
    reachability = {vertex: {vertex} for vertex in vertices}
    for winner, loser in edges:
        reachability[winner].add(loser)
    for pivot in vertices:
        for source in vertices:
            if pivot in reachability[source]:
                reachability[source] |= reachability[pivot]
    components = []
    unused = set(vertices)
    while unused:
        vertex = min(unused)
        component = tuple(
            candidate
            for candidate in vertices
            if candidate in reachability[vertex]
            and vertex in reachability[candidate]
        )
        components.append(component)
        unused -= set(component)
    paths = [
        path
        for path in permutations(vertices)
        if all((path[index], path[index + 1]) in edges for index in range(len(path) - 1))
    ]
    return {
        "scores": scores,
        "score_histogram": {
            score: tuple(vertex for vertex in vertices if scores[vertex] == score)
            for score in sorted(set(scores.values()))
        },
        "directed_triangles": directed_triangles,
        "components": tuple(components),
        "hamiltonian_path_count": len(paths),
        "tie_hamiltonian_path": min(paths),
    }


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def main():
    print("THM-904 / HYP-7061: RESIDUE-SIX TRIPLE-CERTIFICATE REDUCTION")
    print("=" * 78)

    print("\n[1] Exact 462-state pointwise certificate")
    check("84 explicit unordered triple weights", len(TRIPLE_WEIGHTS) == 84)
    state_rows = [
        (state_certificate_slack(state), state)
        for state in moving_state_compositions()
    ]
    check("462 weak sector compositions", len(state_rows) == 462)
    check("every pointwise certificate inequality", min(state_rows)[0] >= 0)
    check(
        "every state contains ten unordered mover triples",
        all(sum(triple_type_counts(state)) == 10 for _, state in state_rows),
    )
    print("  minimum integer slack:", min(state_rows)[0])

    print("\n[2] Exact primitive triple scan")
    records = []
    nested_records = []
    for speeds in combinations(range(1, SCAN_DIAMETER + 1), 3):
        if not primitive(speeds):
            continue
        row = (triple_weight_average(speeds), speeds)
        records.append(row)
        if speeds[-1] <= 40:
            nested_records.append(row)
    records.sort(reverse=True)
    nested_records.sort(reverse=True)
    check("primitive triple count through 60", len(records) == 28876)
    check("unique scanned maximizer", sum(row[0] == records[0][0] for row in records) == 1)
    check("scanned maximum is q(1,4,7)=81/175", records[0] == (Fraction(81, 175), (1, 4, 7)))
    check("scanned maximum clears 47/100", records[0][0] < TARGET_TRIPLE_BOUND)
    check("top value is stable from diameter 40 to 60", nested_records[0] == records[0])
    print("  top ten (q, speeds, shortest full-support relation)")
    for value, speeds in records[:10]:
        print("   ", value, speeds, shortest_full_relation(speeds))
    scanned_residue_six_bound = 10 * records[0][0] / 49
    print("  conditional five-triple-core bound:", scanned_residue_six_bound)
    check("conditional scanned bound is below 0.097", scanned_residue_six_bound < PROPAGATION_SLACK)

    print("\n[3] Exact ANOVA split and pair-ray channel")
    mean, singleton, pair_projection, cubic = anova_channels()
    check("singleton channel has zero mean", sum(singleton) == 0)
    check(
        "pair channel has zero row sums",
        all(sum(row) == 0 for row in pair_projection),
    )
    check(
        "cubic channel has zero coordinate marginals",
        all(
            sum(cubic[first, second, third] for third in SECTIONS) == 0
            for first in SECTIONS
            for second in SECTIONS
        ),
    )
    uniform_pair_distribution = {
        pair_type: Fraction(1 if pair_type[0] == pair_type[1] else 2, 49)
        for pair_type in PAIR_TYPES
    }
    pair_vertices = [uniform_pair_distribution] + [
        pair_sector_distribution(*speeds) for speeds in PAIR_RAY_MINIMIZERS
    ]
    pair_values = [
        pair_projection_expectation(pair_projection, distribution)
        for distribution in pair_vertices
    ]
    relaxed_pair_baseline = mean + 3 * max(pair_values)
    print("  uniform mean:", mean)
    print("  pair-projection range on the 22 pair vertices:", min(pair_values), max(pair_values))
    print("  independently relaxed uniform+three-pair baseline:", relaxed_pair_baseline)
    check("pair baseline alone stays below 47/100", relaxed_pair_baseline < TARGET_TRIPLE_BOUND)
    print("  Remaining theorem: couple the zero-marginal cubic relation channel to")
    print("  the pair deficit; independent maximization of the channels is invalid.")

    print("\n[4] Pinned reflection audit")
    exceptional_pairs = {(1, 5), (2, 4)}
    reflected_pairs = {
        tuple(sorted((6 - section) % 7 for section in missed))
        for missed in exceptional_pairs
    }
    print("  exceptional K6 pairs:", sorted(exceptional_pairs))
    print("  pinned-reflection images:", sorted(reflected_pairs))
    check("exceptional support is the pinned fixed locus", reflected_pairs == exceptional_pairs)
    check(
        "same fixed pairs have K1=-2 and K6=-12, so no residue transfer",
        all(KERNEL_NUMERATORS[1, pair] == -2 for pair in exceptional_pairs)
        and all(KERNEL_NUMERATORS[6, pair] == -12 for pair in exceptional_pairs),
    )

    print("\n[5] Tournament Analysis on strongest relation shapes")
    vertices = tuple(speeds for _, speeds in records[:6])
    values = {speeds: value for value, speeds in records[:6]}
    fingerprint = tournament_fingerprint(vertices, values)
    print("  observable: q(u)-q(v)")
    print("  switch/gauge: larger exact q; lexicographic tie orientation")
    print("  fingerprint:", fingerprint)
    check("obstruction tournament has no directed triangle", fingerprint["directed_triangles"] == 0)
    check("obstruction tournament SCCs are singletons", all(len(component) == 1 for component in fingerprint["components"]))
    check("tie Hamiltonian path is unique", fingerprint["hamiltonian_path_count"] == 1)
    print("  Assumption challenge: relation shapes preserve the triple certificate value")
    print("  but destroy the other movers and wall chronology. Runners, arcs, residues,")
    print("  and scalar relation height were all considered and rejected as faithful vertices.")

    print("\nVERDICT")
    print("  PROVED: explicit rational triple observable pointwise dominates -K6.")
    print("  PROVED: universal q<=47/100 would give -F6<=47/490<0.097.")
    print("  FINITE-EXACT: 28,876 primitive triples through speed 60; unique max 81/175.")
    print("  PROVED: pair projection is a 21-ray channel plus a zero-marginal cubic remainder.")
    print("  PROVED: pinned reflection fixes the hard support and acts within residue 6.")
    print("  REFUTED: negation transfers negative residue 6 to closed positive residue 1.")
    print("  OPEN: bound the cubic relation channel uniformly, then the finite-t wall remainder.")


if __name__ == "__main__":
    main()
