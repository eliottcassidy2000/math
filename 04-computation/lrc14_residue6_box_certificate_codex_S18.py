#!/usr/bin/env python3
"""Exact verifier for THM-905 and HYP-7062.

The residue-six negative kernel is pointwise dominated by three occupation
monomials.  Expanding them reduces the limiting sign to two four-runner and
one five-runner sector-box inequalities.  This script proves the finite
algebra, verifies the occupation expansion on all labelled sector states,
and scans primitive speed tuples exactly on nested diameter boxes.

Tournament Analysis is diagnostic.  Vertices are the strongest scanned
four-speed tuples.  The pairwise observable is the difference of exact box
masses; switching the target from the exceptional-pair complements to the
singleton-six complement changes the gauge.  Ties point lexicographically,
and the tie Hamiltonian path is the lexicographically first directed path.
This quotient preserves obstruction strength for one box obligation but
destroys the other mover, wall chronology, and most of the relation lattice.
Alternative vertices considered were runners, gaps, fixed sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, primitive relations, and proof obligations.  Target boxes
are retained because they are the exact monomials in the pointwise theorem;
relation fingerprints are carried as a sidecar because box vertices alone
forget the nondecaying plateaus proved in THM-898 and THM-899.
"""

from fractions import Fraction
from heapq import heapify, heappop, heappush
from itertools import combinations, permutations, product
from math import gcd, lcm


SECTIONS = tuple(range(7))
INNER_SECTIONS = tuple(range(1, 7))
FOUR_SCAN_DIAMETER = 32
FIVE_SCAN_DIAMETER = 24

TARGET_A = frozenset((1, 3, 5, 6))
TARGET_B = frozenset((2, 3, 4, 6))
TARGET_C = frozenset((1, 2, 4, 5))
TARGET_ZERO_A = frozenset((0, 1, 3, 5, 6))
TARGET_ZERO_B = frozenset((0, 2, 3, 4, 6))

FOUR_UNION_CAP = Fraction(1, 12)
FOUR_C_CAP = Fraction(5, 42)
FIVE_UNION_CAP = Fraction(40, 441)
PROPAGATION_SLACK = Fraction(97, 1000)


def moving_state_compositions(total=5, sections=7, prefix=()):
    if sections == 1:
        yield prefix + (total,)
        return
    for first_count in range(total + 1):
        yield from moving_state_compositions(
            total - first_count, sections - 1, prefix + (first_count,)
        )


def residue_six_negative_kernel(state):
    missed = tuple(section for section in INNER_SECTIONS if state[section] == 0)
    if len(missed) == 1:
        return -6 if missed == (3,) else 1
    if len(missed) == 2:
        return 12 if missed in ((1, 5), (2, 4)) else -2
    return 0


def occupation_terms(state):
    first_box = state[1] * state[3] * state[5] * state[6]
    second_box = state[2] * state[3] * state[4] * state[6]
    exceptional_complements = first_box + second_box
    pinned_exceptional_complements = state[0] * exceptional_complements
    singleton_six_complement = state[1] * state[2] * state[4] * state[5]
    return (
        exceptional_complements,
        pinned_exceptional_complements,
        singleton_six_complement,
    )


def labelled_occupation_expansion(sectors):
    exceptional_complements = 0
    singleton_six_complement = 0
    for mover_indices in combinations(range(5), 4):
        sector_set = frozenset(sectors[index] for index in mover_indices)
        if len(sector_set) != 4:
            continue
        exceptional_complements += sector_set in (TARGET_A, TARGET_B)
        singleton_six_complement += sector_set == TARGET_C
    full_sector_set = frozenset(sectors)
    pinned_exceptional_complements = int(
        len(full_sector_set) == 5
        and full_sector_set in (TARGET_ZERO_A, TARGET_ZERO_B)
    )
    return (
        exceptional_complements,
        pinned_exceptional_complements,
        singleton_six_complement,
    )


def box_masses(speeds, targets):
    common_divisor = gcd(*speeds)
    reduced = tuple(speed // common_divisor for speed in speeds)
    period_scale = lcm(*reduced)
    sectors = [0] * len(reduced)
    events = []
    for runner_index, speed in enumerate(reduced):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    target_indices = {target: index for index, target in enumerate(targets)}
    numerators = [0] * len(targets)
    previous = 0
    while events:
        event_position = events[0][0]
        if len(set(sectors)) == len(sectors):
            target_index = target_indices.get(frozenset(sectors))
            if target_index is not None:
                numerators[target_index] += event_position - previous
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
    denominator = 7 * period_scale
    return tuple(Fraction(numerator, denominator) for numerator in numerators)


def primitive(speeds):
    return gcd(*speeds) == 1


def canonical_relation(relation):
    relation_divisor = gcd(*(abs(coefficient) for coefficient in relation))
    reduced = tuple(coefficient // relation_divisor for coefficient in relation)
    first_nonzero = next(coefficient for coefficient in reduced if coefficient)
    return tuple(-coefficient for coefficient in reduced) if first_nonzero < 0 else reduced


def short_relation_fingerprint(speeds, height=2):
    relations = set()
    for relation in product(range(-height, height + 1), repeat=len(speeds)):
        support = sum(coefficient != 0 for coefficient in relation)
        if support < 3 or sum(coefficient * speed for coefficient, speed in zip(relation, speeds)):
            continue
        relations.add(canonical_relation(relation))
    support_histogram = {
        support: sum(sum(coefficient != 0 for coefficient in relation) == support for relation in relations)
        for support in range(3, len(speeds) + 1)
    }
    return tuple(sorted(support_histogram.items())), len(relations)


def oriented_edges(vertices, values):
    edges = set()
    for first, second in combinations(vertices, 2):
        if values[first] > values[second]:
            edges.add((first, second))
        elif values[second] > values[first]:
            edges.add((second, first))
        else:
            edges.add((min(first, second), max(first, second)))
    return edges


def tournament_fingerprint(vertices, values):
    edges = oriented_edges(vertices, values)
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
    score_histogram = {
        score: sum(value == score for value in scores.values())
        for score in sorted(set(scores.values()))
    }
    return {
        "score_histogram": score_histogram,
        "directed_triangles": directed_triangles,
        "components": tuple(components),
        "hamiltonian_path_count": len(paths),
        "tie_hamiltonian_path": min(paths),
    }


def nested_maxima(records, diameters):
    return tuple(
        (diameter, max(record for record in records if record[1][-1] <= diameter))
        for diameter in diameters
    )


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def main():
    print("THM-905 / HYP-7062: RESIDUE-SIX BOX-CERTIFICATE REDUCTION")
    print("=" * 78)

    print("\n[1] Exact pointwise and occupation-expansion algebra")
    states = list(moving_state_compositions())
    check("462 weak sector compositions", len(states) == 462)
    slacks = []
    for state in states:
        exceptional, pinned, singleton_six = occupation_terms(state)
        slacks.append(
            6 * exceptional
            + 6 * pinned
            + singleton_six
            - residue_six_negative_kernel(state)
        )
    check("6H+6J+H6 dominates the negative K6 kernel", min(slacks) >= 0)
    check("pointwise certificate is sharp", min(slacks) == 0)
    print("  minimum integer slack:", min(slacks))

    labelled_states = product(SECTIONS, repeat=5)
    check(
        "occupation monomials equal their labelled box expansions",
        all(
            occupation_terms(tuple(sectors.count(section) for section in SECTIONS))
            == labelled_occupation_expansion(sectors)
            for sectors in labelled_states
        ),
    )

    print("\n[2] Exact primitive four-speed scan")
    four_records = []
    four_values = {}
    for speeds in combinations(range(1, FOUR_SCAN_DIAMETER + 1), 4):
        if not primitive(speeds):
            continue
        first_mass, second_mass, third_mass = box_masses(
            speeds, (TARGET_A, TARGET_B, TARGET_C)
        )
        union_mass = first_mass + second_mass
        four_records.append((union_mass, speeds, third_mass))
        four_values[speeds] = (union_mass, third_mass)
    union_records = [(union_mass, speeds) for union_mass, speeds, _ in four_records]
    third_records = [(third_mass, speeds) for _, speeds, third_mass in four_records]
    union_maximum = max(union_records)
    third_maximum = max(third_records)
    check("scanned four-speed union cap", union_maximum[0] <= FOUR_UNION_CAP)
    check("union cap is attained at (2,3,4,6)", four_values[(2, 3, 4, 6)][0] == FOUR_UNION_CAP)
    check("scanned singleton-six box cap", third_maximum[0] <= FOUR_C_CAP)
    check("singleton-six cap is attained at (1,2,3,4)", four_values[(1, 2, 3, 4)][1] == FOUR_C_CAP)
    print("  primitive tuples:", len(four_records))
    print("  nested union maxima:", nested_maxima(union_records, (12, 20, FOUR_SCAN_DIAMETER)))
    print("  nested singleton-six maxima:", nested_maxima(third_records, (12, 20, FOUR_SCAN_DIAMETER)))
    print("  top six union rows:")
    for value, speeds in sorted(union_records, reverse=True)[:6]:
        print("   ", value, speeds, short_relation_fingerprint(speeds))
    print("  top six singleton-six rows:")
    for value, speeds in sorted(third_records, reverse=True)[:6]:
        print("   ", value, speeds, short_relation_fingerprint(speeds))

    print("\n[3] Exact primitive five-speed scan")
    five_records = []
    for speeds in combinations(range(1, FIVE_SCAN_DIAMETER + 1), 5):
        if not primitive(speeds):
            continue
        first_mass, second_mass = box_masses(
            speeds, (TARGET_ZERO_A, TARGET_ZERO_B)
        )
        five_records.append((first_mass + second_mass, speeds))
    five_maximum = max(five_records)
    check("scanned five-speed union cap", five_maximum[0] <= FIVE_UNION_CAP)
    check(
        "five-speed cap is attained at (3,5,7,9,11)",
        dict((speeds, value) for value, speeds in five_records)[(3, 5, 7, 9, 11)]
        == FIVE_UNION_CAP,
    )
    print("  primitive tuples:", len(five_records))
    print("  nested maxima:", nested_maxima(five_records, (12, 18, FIVE_SCAN_DIAMETER)))
    print("  top six rows:")
    for value, speeds in sorted(five_records, reverse=True)[:6]:
        print("   ", value, speeds, short_relation_fingerprint(speeds))

    print("\n[4] Conditional residue-six closure")
    conditional_kernel_bound = (
        30 * FOUR_UNION_CAP + 6 * FIVE_UNION_CAP + 5 * FOUR_C_CAP
    )
    conditional_residue_bound = conditional_kernel_bound / 49
    loose_five_cap = (
        49 * PROPAGATION_SLACK - 30 * FOUR_UNION_CAP - 5 * FOUR_C_CAP
    ) / 6
    print("  conditional 49(-F6) bound:", conditional_kernel_bound)
    print("  conditional -F6 bound:", conditional_residue_bound)
    print("  five-box cap actually sufficient after four-box caps:", loose_five_cap)
    check("candidate caps imply -F6<0.097", conditional_residue_bound < PROPAGATION_SLACK)
    check("scanned five-box cap has ample slack", FIVE_UNION_CAP < loose_five_cap)

    print("\n[5] Tournament Analysis under target switch")
    selected_vertices = tuple(
        sorted(
            set(speeds for _, speeds in sorted(union_records, reverse=True)[:4])
            | set(speeds for _, speeds in sorted(third_records, reverse=True)[:4])
        )
    )
    union_values = {speeds: four_values[speeds][0] for speeds in selected_vertices}
    third_values = {speeds: four_values[speeds][1] for speeds in selected_vertices}
    union_fingerprint = tournament_fingerprint(selected_vertices, union_values)
    third_fingerprint = tournament_fingerprint(selected_vertices, third_values)
    union_edges = oriented_edges(selected_vertices, union_values)
    third_edges = oriented_edges(selected_vertices, third_values)
    edge_flips = sum(
        ((first, second) in union_edges) != ((first, second) in third_edges)
        for first, second in combinations(selected_vertices, 2)
    )
    print("  vertices:", selected_vertices)
    print("  observable: exact target-box mass difference")
    print("  switch/gauge: exceptional complements -> singleton-six complement")
    print("  exceptional-complement fingerprint:", union_fingerprint)
    print("  singleton-six-complement fingerprint:", third_fingerprint)
    print("  target-switch edge flips:", edge_flips)
    check("both scalar tournaments are transitive", union_fingerprint["directed_triangles"] == 0 and third_fingerprint["directed_triangles"] == 0)
    check("both have singleton SCCs", all(len(component) == 1 for component in union_fingerprint["components"] + third_fingerprint["components"]))
    print("  Assumption challenge: target boxes preserve (1) but erase the relation lattice;")
    print("  the short-relation histogram is therefore retained as a non-tournament sidecar.")

    print("\nVERDICT")
    print("  PROVED: -K6 <= 6H+6J+H6 on every five-mover sector state.")
    print("  PROVED: occupation expansion gives the exact three-box reduction.")
    print("  FINITE-EXACT: four-speed caps through 32 and five-speed cap through 24.")
    print("  CONDITIONAL: the three universal caps give -F6<=535/7203<0.097.")
    print("  OPEN: prove the caps uniformly by Bernoulli relation-lattice stratification.")


if __name__ == "__main__":
    main()
