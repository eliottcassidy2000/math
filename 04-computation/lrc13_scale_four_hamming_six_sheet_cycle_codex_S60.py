#!/usr/bin/env python3
"""Exact scale-four AP-centred Hamming-six sheet reduction (THM-952).

The program performs five independent exact checks.

1. It scans every c=4 effective-order/unit state allowed by the hereditary
   leave-one-out lcm law and applies the literal four-sheet masks.
2. It identifies the surviving label sets with the c=2 signed-doubling
   Hamiltonian-cycle bank and reconstructs the edge-coloured K6 carrier.
3. It proves, by direct CRT/parity equivalence, that the unit fibre is an
   affine rank-four code supported on two disjoint signed triangles.
4. It audits multiplier/reflection orbits and a deliberately completed
   tournament, including its declared tie Hamiltonian path.
5. It counts the complete THM-815 root-cap first layer.  This is only a
   workload/reduction certificate: no terminal metric verdict is claimed.

Only integer arithmetic and fractions.Fraction are used.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


P = 13
SCALE = 4
LABELS = tuple(range(1, P))
FULL_SHEET_MASK = (1 << SCALE) - 1

# (effective order, unit modulo that order).  D=1 has its trivial unit and
# D=2 has its unique unit.  The two D=4 states are the signs +1 and -1.
STATES = ((1, 0), (2, 1), (4, 1), (4, 3))
D4_STATES = frozenset((2, 3))

# In an owner-normalized four-sheet gauge, these two ratio classes supply the
# two unit-sensitive sheets.  The signs record which class is reversed.
VARIABLE_SIGN = {3: 1, 4: 1, 9: -1, 10: -1}
FIXED_PROVIDER_RATIOS = frozenset((6, 7))
DOUBLING_OWNER_RATIOS = frozenset((2, 11))

LINES: list[str] = []


def emit(line: str = "") -> None:
    LINES.append(line)
    print(line)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    """Least nonnegative solution u mod 13D to u=Dr mod 13, u=e mod D."""
    return next(
        u
        for u in range(P * order)
        if u % P == order * label % P and u % order == unit % order
    )


def sheet_mask(label: int, state: int, owner: int) -> int:
    order, unit = STATES[state]
    u = crt_base(label, order, unit)
    inverse_owner = pow(owner, -1, P)
    answer = 0
    for sheet in range(SCALE):
        z = centered(u * (inverse_owner + P * sheet), P * order)
        if -order < z <= order:
            answer |= 1 << sheet
    return answer


MASK = {
    (label, state, owner): sheet_mask(label, state, owner)
    for label in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}


def fnv64(data: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in data:
        value ^= byte
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def direct_context_cover(
    labels: tuple[int, ...], states: tuple[int, ...]
) -> bool:
    return all(
        reduce(
            int.__or__,
            (
                MASK[label, state, owner]
                for label, state in zip(labels, states)
            ),
            0,
        )
        == FULL_SHEET_MASK
        for owner in labels
    )


def context_row(
    labels: tuple[int, ...], states: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    return (
        labels,
        tuple(STATES[state][0] for state in states),
        tuple(STATES[state][1] for state in states),
    )


def full_sheet_bank():
    """Scan all 924*3648 hereditary-lcm-compatible state words."""
    contexts = set()
    presentations = set()
    tested = 0
    state_words = tuple(
        states
        for states in product(range(len(STATES)), repeat=6)
        if sum(state in D4_STATES for state in states) >= 2
    )
    require(len(state_words) == 3_648, "unexpected hereditary-lcm state count")

    for labels in combinations(LABELS, 6):
        # Pack all four owner sheets at six owners into one 24-bit vector.
        packed = {
            (label, state): sum(
                MASK[label, state, owner] << (SCALE * index)
                for index, owner in enumerate(labels)
            )
            for label in labels
            for state in range(len(STATES))
        }
        full = (1 << (SCALE * len(labels))) - 1
        for states in state_words:
            tested += 1
            union = 0
            for label, state in zip(labels, states):
                union |= packed[label, state]
            if union != full:
                continue
            row = context_row(labels, states)
            contexts.add(row)
            presentations.add(row[:2])

    require(tested == 3_370_752, "unexpected full-bank candidate count")
    require(len(contexts) == 256, "unexpected c=4 context count")
    require(len(presentations) == 64, "unexpected c=4 presentation count")
    require(
        {orders for _, orders, _ in contexts} == {(4, 4, 4, 4, 4, 4)},
        "a mixed-order context survived",
    )
    require(
        Counter(units for _, _, units in contexts).total() == 256,
        "context unit census mismatch",
    )
    return tuple(sorted(contexts)), tuple(sorted(presentations)), tested


CONTEXTS, PRESENTATIONS, TESTED_CONTEXTS = full_sheet_bank()
SUPPORTS = tuple(sorted(labels for labels, _ in PRESENTATIONS))


def doubling_edges(labels: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    """Provider -> owner law from the c=2 signed-doubling bank."""
    return frozenset(
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner
        and owner * pow(provider, -1, P) % P in DOUBLING_OWNER_RATIOS
    )


def c2_cycle_support(labels: tuple[int, ...]) -> bool:
    edges = doubling_edges(labels)
    return (
        len(edges) == 6
        and all(sum((provider, owner) in edges for provider in labels) == 1 for owner in labels)
        and all(sum((provider, owner) in edges for owner in labels) == 1 for provider in labels)
    )


C2_SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if c2_cycle_support(labels)
)
require(len(C2_SUPPORTS) == 64, "unexpected c=2 cycle-support count")
require(SUPPORTS == C2_SUPPORTS, "c=4 supports differ from the c=2 cycle bank")


def undirected_edges(
    labels: tuple[int, ...], allowed_ratios: frozenset[int]
) -> frozenset[tuple[int, int]]:
    return frozenset(
        (a, b)
        for a, b in combinations(labels, 2)
        if b * pow(a, -1, P) % P in allowed_ratios
    )


def component_sets(
    labels: tuple[int, ...], edges: frozenset[tuple[int, int]]
) -> tuple[frozenset[int], ...]:
    adjacency = {label: set() for label in labels}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    unseen = set(labels)
    answer = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component = {seed}
        unseen.remove(seed)
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        answer.append(frozenset(component))
    return tuple(sorted(answer, key=lambda component: tuple(sorted(component))))


def variable_graph(labels: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    return undirected_edges(labels, frozenset(VARIABLE_SIGN))


def zero_matching(
    labels: tuple[int, ...],
    cycle_edges: frozenset[tuple[int, int]],
    triangle_edges: frozenset[tuple[int, int]],
) -> frozenset[tuple[int, int]]:
    cycle_pairs = {tuple(sorted(edge)) for edge in cycle_edges}
    return frozenset(
        pair
        for pair in combinations(labels, 2)
        if pair not in cycle_pairs and pair not in triangle_edges
    )


def unit_sign(unit: int) -> int:
    require(unit in (1, 3), "non-D4 unit in parity code")
    return 1 if unit == 1 else -1


def parity_constraints(
    labels: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], int, int], ...]:
    """(opposite edge, rhs, owner) for the six owner obligations."""
    answer = []
    for owner in labels:
        neighbors = tuple(
            provider
            for provider in labels
            if provider != owner and ratio(provider, owner) in VARIABLE_SIGN
        )
        require(len(neighbors) == 2, "owner lacks its two variable providers")
        a, b = neighbors
        rhs = -VARIABLE_SIGN[ratio(a, owner)] * VARIABLE_SIGN[ratio(b, owner)]
        answer.append((tuple(sorted((a, b))), rhs, owner))
    return tuple(answer)


def parity_cover(labels: tuple[int, ...], units: tuple[int, ...]) -> bool:
    signs = {label: unit_sign(unit) for label, unit in zip(labels, units)}
    return all(
        signs[a] * signs[b] == rhs
        for (a, b), rhs, _ in parity_constraints(labels)
    )


def gf2_rank(rows: tuple[int, ...]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


PARITY_PAYLOAD = []
for labels in SUPPORTS:
    cycle = doubling_edges(labels)
    triangle_edges = variable_graph(labels)
    components = component_sets(labels, triangle_edges)
    matching = zero_matching(labels, cycle, triangle_edges)

    require(len(cycle) == 6, "cycle edge count mismatch")
    require(len({tuple(sorted(edge)) for edge in cycle}) == 6, "cycle has a two-cycle")
    require(len(triangle_edges) == 6, "variable edge count mismatch")
    require(tuple(sorted(map(len, components))) == (3, 3), "variable graph is not 2K3")
    require(
        all(
            sum(vertex in edge for edge in triangle_edges) == 2
            for vertex in labels
        ),
        "variable graph is not two triangles",
    )
    require(len(matching) == 3, "zero-pair complement is not a matching")
    require(
        all(sum(vertex in edge for edge in matching) == 1 for vertex in labels),
        "zero-pair complement is not 3K2",
    )

    constraints = parity_constraints(labels)
    require(
        {edge for edge, _, _ in constraints} == set(triangle_edges),
        "owner constraints do not index the opposite triangle edges",
    )
    index = {label: position for position, label in enumerate(labels)}
    coefficient_rows = tuple(
        (1 << index[a]) | (1 << index[b]) for (a, b), _, _ in constraints
    )
    require(gf2_rank(coefficient_rows) == 4, "unit parity rank is not four")

    rhs_by_edge = {edge: rhs for edge, rhs, _ in constraints}
    for component in components:
        component_edges = tuple(
            edge for edge in triangle_edges if edge[0] in component
        )
        require(len(component_edges) == 3, "triangle edge count mismatch")
        require(
            reduce(int.__mul__, (rhs_by_edge[edge] for edge in component_edges), 1)
            == 1,
            "signed triangle is inconsistent",
        )

    direct_words = {
        units
        for row_labels, _, units in CONTEXTS
        if row_labels == labels
    }
    parity_words = {
        units for units in product((1, 3), repeat=6) if parity_cover(labels, units)
    }
    require(len(direct_words) == 4, "direct unit fibre does not have size four")
    require(direct_words == parity_words, "direct masks and parity code disagree")
    PARITY_PAYLOAD.append(
        f"{','.join(map(str, labels))}:"
        + ";".join(
            f"{a}-{b}:{rhs:+d}@{owner}"
            for (a, b), rhs, owner in sorted(constraints)
        )
    )


def multiply_context(row, multiplier: int):
    labels, orders, units = row
    triples = sorted(
        (multiplier * label % P, order, unit)
        for label, order, unit in zip(labels, orders, units)
    )
    return (
        tuple(value[0] for value in triples),
        tuple(value[1] for value in triples),
        tuple(value[2] for value in triples),
    )


def reflect_context(row):
    labels, orders, units = row
    return labels, orders, tuple(4 - unit for unit in units)


def orbit_sizes(space, group, action) -> tuple[int, ...]:
    unseen = set(space)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {action(seed, element) for element in group}
        require(orbit <= set(space), "claimed action leaves the context bank")
        sizes.append(len(orbit))
        unseen -= orbit
    return tuple(sorted(sizes))


MULTIPLIER_ORBITS = orbit_sizes(
    CONTEXTS,
    tuple(range(1, P)),
    multiply_context,
)


def sheet_action(row, group_element):
    multiplier, reflection = group_element
    image = multiply_context(row, multiplier)
    return reflect_context(image) if reflection else image


SHEET_ORBITS = orbit_sizes(
    CONTEXTS,
    tuple((multiplier, reflection) for multiplier in range(1, P) for reflection in (0, 1)),
    sheet_action,
)
require(Counter(MULTIPLIER_ORBITS) == {4: 4, 12: 20}, "multiplier orbit mismatch")
require(Counter(SHEET_ORBITS) == {8: 2, 24: 10}, "sheet orbit mismatch")


def circle_norm(numerator: int, denominator: int) -> F:
    residue = numerator % denominator
    return F(min(residue, denominator - residue), denominator)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    """Exact piecewise-linear maximum of min_v ||vt||."""
    denominators = {2 * speed for speed in speeds}
    denominators |= {
        denominator
        for u, v in combinations(speeds, 2)
        for denominator in (u + v, abs(u - v))
        if denominator
    }
    best = F(0)
    witnesses: set[F] = set()
    for denominator in sorted(denominators):
        for numerator in range(denominator):
            value = min(
                circle_norm(speed * numerator, denominator) for speed in speeds
            )
            witness = F(numerator, denominator)
            if value > best:
                best = value
                witnesses = {witness}
            elif value == best:
                witnesses.add(witness)
    return best, tuple(sorted(witnesses))


def least_packet(labels: tuple[int, ...], units: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sorted(
            [4 * label for label in LABELS if label not in labels]
            + [crt_base(label, 4, unit) for label, unit in zip(labels, units)]
        )
    )


METRIC_GUARD_LABELS = (1, 2, 3, 4, 5, 6)
METRIC_GUARD_UNITS = (1, 1, 1, 3, 1, 1)
METRIC_GUARD_FLIP = tuple(4 - unit for unit in METRIC_GUARD_UNITS)
METRIC_GUARD_PACKET = least_packet(METRIC_GUARD_LABELS, METRIC_GUARD_UNITS)
METRIC_GUARD_FLIP_PACKET = least_packet(METRIC_GUARD_LABELS, METRIC_GUARD_FLIP)
METRIC_GUARD_VALUE = exact_maximin(METRIC_GUARD_PACKET)
METRIC_GUARD_FLIP_VALUE = exact_maximin(METRIC_GUARD_FLIP_PACKET)
require(METRIC_GUARD_VALUE == (F(7, 31), (F(8, 31), F(23, 31))), "metric guard value mismatch")
require(
    METRIC_GUARD_FLIP_VALUE == (F(14, 79), (F(39, 79), F(40, 79))),
    "reflected metric guard value mismatch",
)
require(
    METRIC_GUARD_VALUE[0] != METRIC_GUARD_FLIP_VALUE[0],
    "global sheet flip accidentally became metric",
)


def directed_reachable(
    source: int,
    target: int,
    labels: tuple[int, ...],
    edges: set[tuple[int, int]],
) -> bool:
    stack = [source]
    seen = {source}
    while stack:
        vertex = stack.pop()
        for neighbor in labels:
            if (vertex, neighbor) not in edges or neighbor in seen:
                continue
            if neighbor == target:
                return True
            seen.add(neighbor)
            stack.append(neighbor)
    return False


def completed_tournament_fingerprint(labels: tuple[int, ...]):
    """Complete the sparse cycle using a declared tie Hamiltonian path.

    The tie graph is K6 minus the undirected cycle, hence a triangular prism.
    Its lexicographically first Hamiltonian path gives the total tie gauge.
    """
    sparse = doubling_edges(labels)

    def is_tie(a: int, b: int) -> bool:
        return (a, b) not in sparse and (b, a) not in sparse

    tie_path = min(
        path
        for path in permutations(labels)
        if all(is_tie(path[index], path[index + 1]) for index in range(5))
    )
    position = {vertex: index for index, vertex in enumerate(tie_path)}
    tournament: set[tuple[int, int]] = set()
    sparse_flips = 0
    for a, b in combinations(labels, 2):
        if (a, b) in sparse:
            edge = (a, b)
        elif (b, a) in sparse:
            edge = (b, a)
        elif position[a] < position[b]:
            edge = (a, b)
        else:
            edge = (b, a)
        tournament.add(edge)
        if edge in sparse and position[edge[0]] > position[edge[1]]:
            sparse_flips += 1

    score_histogram = tuple(
        sorted(sum((vertex, other) in tournament for other in labels) for vertex in labels)
    )
    directed_triangles = sum(
        all(
            sum((vertex, other) in tournament for other in triple) == 1
            for vertex in triple
        )
        for triple in combinations(labels, 3)
    )
    strong = all(
        directed_reachable(a, b, labels, tournament)
        for a in labels
        for b in labels
        if a != b
    )
    hamiltonian_paths = sum(
        all((path[index], path[index + 1]) in tournament for index in range(5))
        for path in permutations(labels)
    )
    return (
        score_histogram,
        directed_triangles,
        (6,) if strong else (),
        sparse_flips,
        hamiltonian_paths,
    )


TOURNAMENT_FINGERPRINTS = Counter(
    completed_tournament_fingerprint(labels) for labels in SUPPORTS
)
require(len(TOURNAMENT_FINGERPRINTS) == 5, "unexpected tournament fingerprint count")
require(
    {fingerprint[0] for fingerprint in TOURNAMENT_FINGERPRINTS}
    == {(1, 2, 2, 3, 3, 4)},
    "completed tournament score histogram mismatch",
)
require(
    {fingerprint[1] for fingerprint in TOURNAMENT_FINGERPRINTS} == {6},
    "completed tournament triangle count mismatch",
)
require(
    {fingerprint[2] for fingerprint in TOURNAMENT_FINGERPRINTS} == {(6,)},
    "completed tournament SCC mismatch",
)


def danger_intervals(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, P * speed)
    answer = []
    for integer in range(speed):
        low = F(integer, speed) - radius
        high = F(integer, speed) + radius
        if low < 0:
            answer.extend(((F(0), high), (low + 1, F(1))))
        elif high > 1:
            answer.extend(((low, F(1)), (F(0), high - 1)))
        else:
            answer.append((low, high))
    return tuple(answer)


DANGER = {speed: danger_intervals(speed) for speed in LABELS}


def longest_strict_safe_component(speeds: tuple[int, ...]) -> F:
    intervals = sorted(interval for speed in speeds for interval in DANGER[speed])
    merged: list[list[F]] = []
    for low, high in intervals:
        if not merged or low > merged[-1][1]:
            merged.append([low, high])
        elif high > merged[-1][1]:
            merged[-1][1] = high
    last = F(0)
    longest = F(0)
    for low, high in merged:
        longest = max(longest, low - last)
        last = max(last, high)
    return max(longest, F(1) - last)


ROOT_LENGTH = {
    labels: longest_strict_safe_component(
        tuple(label for label in LABELS if label not in labels)
    )
    for labels in SUPPORTS
}

LOGICAL_FIRST_EDGES = 0
GEOMETRY_MULTIPLICITY: Counter[tuple[tuple[int, ...], int]] = Counter()
PER_CONTEXT_EDGES = Counter()
ROOT_REAL_CAPS = {}
for labels, orders, units in CONTEXTS:
    require(orders == (4, 4, 4, 4, 4, 4), "non-D4 first-layer row")
    # L(4P)=L(P)/4, and at s=6 the THM-815 real cap is 132/(13L).
    real_cap = F(528, P) / ROOT_LENGTH[labels]
    ROOT_REAL_CAPS[labels] = real_cap
    cap = real_cap.numerator // real_cap.denominator
    edge_count = 0
    for label, unit in zip(labels, units):
        base = crt_base(label, 4, unit)
        for speed in range(base, cap + 1, 52):
            LOGICAL_FIRST_EDGES += 1
            edge_count += 1
            GEOMETRY_MULTIPLICITY[labels, speed] += 1
    PER_CONTEXT_EDGES[edge_count] += 1

require(LOGICAL_FIRST_EDGES == 25_132, "unexpected first-edge count")
require(len(GEOMETRY_MULTIPLICITY) == 12_566, "unexpected geometry-key count")
require(
    Counter(GEOMETRY_MULTIPLICITY.values()) == {2: 12_566},
    "first geometry fibres are not uniformly two-to-one",
)
require(min(PER_CONTEXT_EDGES) == 60, "unexpected minimum context edge count")
require(max(PER_CONTEXT_EDGES) == 168, "unexpected maximum context edge count")
require(min(ROOT_REAL_CAPS.values()) == F(528), "unexpected minimum root cap")
require(max(ROOT_REAL_CAPS.values()) == F(1440), "unexpected maximum root cap")


CONTEXT_PAYLOAD = "\n".join(
    f"{','.join(map(str, labels))}|{','.join(map(str, orders))}|{','.join(map(str, units))}"
    for labels, orders, units in CONTEXTS
) + "\n"
PARITY_TEXT = "\n".join(PARITY_PAYLOAD) + "\n"


emit("THM-952 SCALE-FOUR HAMMING-SIX SHEET-CYCLE EXACT REDUCTION")
emit("SECTION exhaustive common-sheet bank")
emit(
    "label_sets=924 state_words_per_label_set=3648 "
    f"tested_contexts={TESTED_CONTEXTS}"
)
emit(
    "surviving_presentations=64 surviving_unit_contexts=256 "
    "order_strata={(D4^6):256} mixed_order_contexts=0"
)
emit("units_per_presentation={4:64}")
emit(f"context_payload_sha256={sha256(CONTEXT_PAYLOAD.encode()).hexdigest()}")
emit("SECTION edge-coloured K6 and affine unit code")
emit(
    "supports=64 equal_c2_signed_doubling_bank=yes "
    "fixed_provider_graph=directed_C6 variable_graph=2K3 zero_graph=3K2"
)
emit(
    "tie_graph=C6_complement=triangular_prism "
    "coverage_union=C6_union_2K3=octahedral_graph"
)
emit("parity_constraints=6 gf2_rank=4 signed_triangle_products=(+1,+1) words=4")
emit(f"parity_payload_fnv64={fnv64(PARITY_TEXT.encode()):016x}")
emit("SECTION sheet actions")
emit(f"multiplier_context_orbit_sizes={dict(sorted(Counter(MULTIPLIER_ORBITS).items()))}")
emit(f"multiplier_x_global_flip_orbit_sizes={dict(sorted(Counter(SHEET_ORBITS).items()))}")
emit(
    f"metric_guard_packet={METRIC_GUARD_PACKET} exact_M={METRIC_GUARD_VALUE[0]} "
    f"witnesses={METRIC_GUARD_VALUE[1]}"
)
emit(
    f"global_flip_packet={METRIC_GUARD_FLIP_PACKET} "
    f"exact_M={METRIC_GUARD_FLIP_VALUE[0]} "
    f"witnesses={METRIC_GUARD_FLIP_VALUE[1]} sheet_symmetry_is_metric=no"
)
emit("SECTION tournament audit")
emit(
    "observable=signed_doubling_provider_arc switch=complete_missing_pairs "
    "tie_gauge=lexicographically_first_Hamiltonian_path_of_triangular_prism"
)
emit(
    "completed_score_histogram={(1,2,2,3,3,4):64} "
    "directed_triangles={6:64} SCCs={(6,):64}"
)
emit(
    "sparse_edge_flip_histogram="
    f"{dict(sorted(Counter(fp[3] for fp in map(completed_tournament_fingerprint, SUPPORTS)).items()))}"
)
emit(
    "Hamiltonian_path_histogram="
    f"{dict(sorted(Counter(fp[4] for fp in map(completed_tournament_fingerprint, SUPPORTS)).items()))}"
)
emit(f"joint_tournament_fingerprints={len(TOURNAMENT_FINGERPRINTS)}")
emit("SECTION exact THM-815 root-cap first layer")
emit(
    f"root_contexts={len(CONTEXTS)} distinct_roots={len(ROOT_LENGTH)} "
    f"real_cap_range=[{min(ROOT_REAL_CAPS.values())},{max(ROOT_REAL_CAPS.values())}]"
)
emit(
    f"logical_first_edges={LOGICAL_FIRST_EDGES} "
    f"per_context_range=[{min(PER_CONTEXT_EDGES)},{max(PER_CONTEXT_EDGES)}]"
)
emit(
    f"geometric_keys={len(GEOMETRY_MULTIPLICITY)} "
    "fibre_multiplicity={2:12566} exact_cache_factor=2"
)
emit("terminal_metric_recursion=NOT_RUN verdict=OPEN")
emit("SECTION integrity")
emit(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
