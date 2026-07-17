#!/usr/bin/env python3
"""Exact c=5 AP-centred Hamming-six common-sheet obstruction (THM-958).

The verifier has two independent layers:

1. a literal 30-bit scan of all 924 label sets and all 15,600
   leave-one-out-lcm-compatible effective-order/unit words; and
2. a structural reduction through sheet capacity, the 64 signed-cycle
   supports, and their six owner all-different obligations.

The structural layer also audits the sparse carrier, the forbidden-owner
triple hypergraph, and a declared tournament completion.  No metric height or
floating-point arithmetic occurs because the common-sheet bank is empty.
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


P = 13
SCALE = 5
LABELS = tuple(range(1, P))
FULL_SHEET_MASK = (1 << SCALE) - 1
STATES = ((1, 0), (5, 1), (5, 2), (5, 3), (5, 4))
D5_STATES = frozenset(range(1, 5))
ZERO_RATIOS = frozenset((4, 9, 12))
SYMMETRIC_NONZERO_RATIOS = frozenset((2, 5, 6, 7, 8, 11))
DOUBLING_OWNER_RATIOS = frozenset((2, 11))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    return next(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )


def sheet_mask(label: int, state: int, owner: int) -> int:
    order, unit = STATES[state]
    base = crt_base(label, order, unit)
    inverse_owner = pow(owner, -1, P)
    answer = 0
    for sheet in range(SCALE):
        value = centered(base * (inverse_owner + P * sheet), P * order)
        if -order < value <= order:
            answer |= 1 << sheet
    return answer


MASK = {
    (label, state, owner): sheet_mask(label, state, owner)
    for label in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def fnv64(data: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in data:
        value ^= byte
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


# Freeze the literal owner-one mask table.  D=1 is the automatic owner mask.
# A D=5 provider supplies one sheet except at ratios 4,9,12; its self sheet is
# unit-independent, and every nonzero off-owner unit map is a permutation of
# the other four sheets.
LOCAL_MASK_PAYLOAD = []
for state, (order, unit) in enumerate(STATES):
    grouped: dict[int, list[int]] = {}
    for label in LABELS:
        grouped.setdefault(MASK[label, state, 1], []).append(label)
    LOCAL_MASK_PAYLOAD.append(
        f"D={order},e={unit}:"
        + ";".join(
            f"{mask:02x}={','.join(map(str, grouped[mask]))}"
            for mask in sorted(grouped)
        )
    )

for owner in LABELS:
    self_masks = {MASK[owner, state, owner] for state in D5_STATES}
    require(len(self_masks) == 1, "D=5 self sheet depends on the unit")
    self_mask = next(iter(self_masks))
    require(self_mask.bit_count() == 1, "D=5 self mask is not a singleton")
    for provider in LABELS:
        if provider == owner:
            continue
        masks = tuple(MASK[provider, state, owner] for state in D5_STATES)
        if ratio(provider, owner) in ZERO_RATIOS:
            require(set(masks) == {0}, "zero ratio supplies a sheet")
        else:
            require(
                all(mask.bit_count() == 1 for mask in masks),
                "nonzero D=5 mask is not a singleton",
            )
            require(len(set(masks)) == 4, "unit map is not a four-sheet permutation")
            require(self_mask not in masks, "off-owner mask hits the self sheet")


# Literal exhaustive scan.  Hereditary lcm=5 is exactly the requirement that
# at least two coordinates have effective order five.
STATE_WORDS_WITH_COUNT = tuple(
    (word, sum(state in D5_STATES for state in word))
    for word in product(range(len(STATES)), repeat=6)
    if sum(state in D5_STATES for state in word) >= 2
)
STATE_WORDS = tuple(word for word, _ in STATE_WORDS_WITH_COUNT)
require(len(STATE_WORDS) == 15_600, "unexpected c=5 state-word count")
EXPECTED_TESTED_BY_D5 = {
    count: 924 * sum(d5_count == count for _, d5_count in STATE_WORDS_WITH_COUNT)
    for count in range(2, 7)
}
require(
    EXPECTED_TESTED_BY_D5
    == {2: 221_760, 3: 1_182_720, 4: 3_548_160, 5: 5_677_056, 6: 3_784_704},
    "unexpected c=5 stratum sizes",
)

TESTED_BY_D5 = Counter()
SURVIVORS_BY_D5 = Counter()
for labels in combinations(LABELS, 6):
    packed = {
        (label, state): sum(
            MASK[label, state, owner] << (SCALE * index)
            for index, owner in enumerate(labels)
        )
        for label in labels
        for state in range(len(STATES))
    }
    full = (1 << (SCALE * 6)) - 1
    for word, d5_count in STATE_WORDS_WITH_COUNT:
        TESTED_BY_D5[d5_count] += 1
        union = 0
        for label, state in zip(labels, word):
            union |= packed[label, state]
        if union == full:
            SURVIVORS_BY_D5[d5_count] += 1

require(dict(TESTED_BY_D5) == EXPECTED_TESTED_BY_D5, "literal stratum scan mismatch")
require(sum(TESTED_BY_D5.values()) == 14_414_400, "literal candidate count mismatch")
require(not SURVIVORS_BY_D5, "a literal c=5 common-sheet context survived")


# Capacity reduction.  With k D=5 providers, a D=5 owner sees at most k
# singleton masks, so k<=4 is impossible.  At k=5 every ordered pair of D=5
# labels must be nonzero.  Both directions are nonzero exactly at the six
# odd powers of generator 2, giving K6,6 and clique number two.
generator_powers = []
value = 1
for exponent in range(12):
    generator_powers.append(value)
    value = 2 * value % P
require(len(set(generator_powers)) == 12 and value == 1, "2 is not a generator")
exponent = {value: index for index, value in enumerate(generator_powers)}
require(
    {exponent[value] for value in SYMMETRIC_NONZERO_RATIOS}
    == {1, 3, 5, 7, 9, 11},
    "symmetric nonzero graph is not the odd-exponent graph",
)
SYMMETRIC_GRAPH_EDGES = {
    (a, b)
    for a, b in combinations(LABELS, 2)
    if b * pow(a, -1, P) % P in SYMMETRIC_NONZERO_RATIOS
}
require(len(SYMMETRIC_GRAPH_EDGES) == 36, "symmetric graph is not K6,6")
require(
    all(
        sum(vertex in edge for edge in SYMMETRIC_GRAPH_EDGES) == 6
        for vertex in LABELS
    ),
    "symmetric graph degree mismatch",
)
require(
    not any(
        all(tuple(sorted(pair)) in SYMMETRIC_GRAPH_EDGES for pair in combinations(triple, 2))
        for triple in combinations(LABELS, 3)
    ),
    "symmetric graph contains a triangle",
)


def zero_edges(labels: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    return frozenset(
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner and ratio(provider, owner) in ZERO_RATIOS
    )


def capacity_support(labels: tuple[int, ...]) -> bool:
    edges = zero_edges(labels)
    return all(sum((provider, owner) in edges for provider in labels) <= 1 for owner in labels)


CAPACITY_SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if capacity_support(labels)
)
require(len(CAPACITY_SUPPORTS) == 64, "unexpected all-D5 capacity-support count")
for labels in CAPACITY_SUPPORTS:
    edges = zero_edges(labels)
    require(len(edges) == 6, "capacity support does not have six zero arcs")
    require(
        all(sum((provider, owner) in edges for provider in labels) == 1 for owner in labels),
        "capacity support does not have zero indegree one",
    )
    require(
        all(sum((provider, owner) in edges for owner in labels) == 1 for provider in labels),
        "capacity support does not have zero outdegree one",
    )
    require(
        {ratio(provider, owner) for provider, owner in edges} <= {4, 9},
        "capacity support contains an antipodal zero arc",
    )


def doubling_edges(labels: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    return frozenset(
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner
        and owner * pow(provider, -1, P) % P in DOUBLING_OWNER_RATIOS
    )


def cycle_support(labels: tuple[int, ...]) -> bool:
    edges = doubling_edges(labels)
    return (
        len(edges) == 6
        and all(sum((provider, owner) in edges for provider in labels) == 1 for owner in labels)
        and all(sum((provider, owner) in edges for owner in labels) == 1 for provider in labels)
    )


CYCLE_SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if cycle_support(labels)
)
require(CAPACITY_SUPPORTS == CYCLE_SUPPORTS, "c=5 capacity bank differs from c=2 cycle bank")


def directed_components(
    labels: tuple[int, ...], edges: frozenset[tuple[int, int]] | set[tuple[int, int]]
) -> tuple[int, ...]:
    reach = {(a, b): a == b or (a, b) in edges for a in labels for b in labels}
    for middle in labels:
        for a in labels:
            for b in labels:
                reach[a, b] = reach[a, b] or (reach[a, middle] and reach[middle, b])
    unused = set(labels)
    sizes = []
    while unused:
        seed = min(unused)
        component = {b for b in labels if reach[seed, b] and reach[b, seed]}
        unused -= component
        sizes.append(len(component))
    return tuple(sorted(sizes, reverse=True))


def directed_hamiltonian_paths(
    labels: tuple[int, ...], edges: frozenset[tuple[int, int]] | set[tuple[int, int]]
) -> int:
    return sum(
        all((path[index], path[index + 1]) in edges for index in range(5))
        for path in permutations(labels)
    )


def directed_cycle_order(labels: tuple[int, ...]) -> tuple[int, ...]:
    edges = doubling_edges(labels)
    successor = {provider: owner for provider, owner in edges}
    start = min(labels)
    answer = [start]
    while successor[answer[-1]] != start:
        answer.append(successor[answer[-1]])
    require(len(answer) == 6, "doubling relation is not a Hamiltonian cycle")
    return tuple(answer)


ZERO_SPARSE_FINGERPRINTS = Counter()
DOUBLING_SPARSE_FINGERPRINTS = Counter()
for labels in CAPACITY_SUPPORTS:
    zero = zero_edges(labels)
    doubling = doubling_edges(labels)
    ZERO_SPARSE_FINGERPRINTS[
        (
            len(zero),
            directed_components(labels, zero),
            directed_hamiltonian_paths(labels, zero),
        )
    ] += 1
    DOUBLING_SPARSE_FINGERPRINTS[
        (
            len(doubling),
            directed_components(labels, doubling),
            directed_hamiltonian_paths(labels, doubling),
        )
    ] += 1
require(ZERO_SPARSE_FINGERPRINTS == {(6, (3, 3), 0): 64}, "zero graph mismatch")
require(DOUBLING_SPARSE_FINGERPRINTS == {(6, (6,), 6): 64}, "cycle graph mismatch")


def owner_cover_mask(labels: tuple[int, ...], units: tuple[int, ...]) -> int:
    answer = 0
    for index, owner in enumerate(labels):
        union = reduce(
            int.__or__,
            (
                MASK[provider, unit, owner]
                for provider, unit in zip(labels, units)
            ),
            0,
        )
        if union == FULL_SHEET_MASK:
            answer |= 1 << index
    return answer


SUPPORT_PAYLOAD = []
SATISFACTION_PROFILE = Counter()
FORBIDDEN_TRIPLE_COUNT = Counter()
MAX_FAILURE_PAIR_COUNT = Counter()
for labels in CAPACITY_SUPPORTS:
    masks = tuple(
        owner_cover_mask(labels, units)
        for units in product(range(1, 5), repeat=6)
    )
    profile = Counter(mask.bit_count() for mask in masks)
    require(
        profile == {0: 2_332, 1: 1_248, 2: 504, 4: 12},
        "owner-satisfaction profile mismatch",
    )
    SATISFACTION_PROFILE.update(profile)
    require(max(profile) == 4, "more than four owner obligations are compatible")
    require(
        all(
            any(mask & pair_mask == pair_mask for mask in masks)
            for pair in combinations(range(6), 2)
            for pair_mask in (sum(1 << index for index in pair),)
        ),
        "an owner pair is already incompatible",
    )

    forbidden = set()
    for triple in combinations(range(6), 3):
        triple_mask = sum(1 << index for index in triple)
        if not any(mask & triple_mask == triple_mask for mask in masks):
            forbidden.add(frozenset(labels[index] for index in triple))

    cycle = directed_cycle_order(labels)
    expected_forbidden = {
        frozenset((cycle[index], cycle[(index + 1) % 6], cycle[(index + 2) % 6]))
        for index in range(6)
    }
    expected_forbidden |= {
        frozenset((cycle[0], cycle[2], cycle[4])),
        frozenset((cycle[1], cycle[3], cycle[5])),
    }
    require(forbidden == expected_forbidden, "forbidden triple hypergraph mismatch")
    FORBIDDEN_TRIPLE_COUNT[len(forbidden)] += 1

    opposite_pairs = {
        frozenset((cycle[index], cycle[(index + 3) % 6])) for index in range(3)
    }
    failure_pairs = Counter()
    for mask in masks:
        if mask.bit_count() != 4:
            continue
        failed = frozenset(
            labels[index] for index in range(6) if not ((mask >> index) & 1)
        )
        failure_pairs[failed] += 1
    require(
        set(failure_pairs) == opposite_pairs and set(failure_pairs.values()) == {4},
        "maximal near-solutions do not fail opposite owner pairs",
    )
    MAX_FAILURE_PAIR_COUNT[len(failure_pairs)] += 1

    SUPPORT_PAYLOAD.append(
        f"{','.join(map(str, labels))}|cycle={','.join(map(str, cycle))}|"
        f"zero={','.join(f'{a}>{b}' for a,b in sorted(zero_edges(labels)))}|"
        + "forbidden="
        + ";".join(
            ",".join(map(str, sorted(triple))) for triple in sorted(forbidden, key=lambda x: tuple(sorted(x)))
        )
    )

require(SATISFACTION_PROFILE == {0: 149_248, 1: 79_872, 2: 32_256, 4: 768}, "global profile mismatch")
require(FORBIDDEN_TRIPLE_COUNT == {8: 64}, "forbidden triple count mismatch")
require(MAX_FAILURE_PAIR_COUNT == {3: 64}, "near-solution pair count mismatch")


def support_orbit_sizes() -> tuple[int, ...]:
    unseen = set(CAPACITY_SUPPORTS)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {
            tuple(sorted(multiplier * label % P for label in seed))
            for multiplier in LABELS
        }
        require(orbit <= set(CAPACITY_SUPPORTS), "multiplier leaves capacity bank")
        sizes.append(len(orbit))
        unseen -= orbit
    return tuple(sizes)


ORBIT_SIZES = support_orbit_sizes()
require(ORBIT_SIZES == (12, 12, 12, 4, 12, 12), "support orbit mismatch")


def completed_zero_tournament(labels: tuple[int, ...]):
    sparse = zero_edges(labels)

    def tie(a: int, b: int) -> bool:
        return (a, b) not in sparse and (b, a) not in sparse

    # The tie graph is the complement of two triangles, hence K3,3.
    tie_path = min(
        path
        for path in permutations(labels)
        if all(tie(path[index], path[index + 1]) for index in range(5))
    )
    position = {vertex: index for index, vertex in enumerate(tie_path)}
    tournament = set()
    flips = 0
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
            flips += 1

    scores = tuple(
        sorted(sum((vertex, other) in tournament for other in labels) for vertex in labels)
    )
    triangles = sum(
        all(
            sum((vertex, other) in tournament for other in triple) == 1
            for vertex in triple
        )
        for triple in combinations(labels, 3)
    )
    return (
        scores,
        triangles,
        directed_components(labels, tournament),
        flips,
        directed_hamiltonian_paths(labels, tournament),
    )


TOURNAMENT_FINGERPRINTS = Counter(
    completed_zero_tournament(labels) for labels in CAPACITY_SUPPORTS
)
require(len(TOURNAMENT_FINGERPRINTS) == 3, "completed tournament fingerprint mismatch")
require(
    Counter(fingerprint[0] for fingerprint in map(completed_zero_tournament, CAPACITY_SUPPORTS))
    == {(1, 2, 2, 3, 3, 4): 64},
    "completed score histogram mismatch",
)
require(
    Counter(fingerprint[1] for fingerprint in map(completed_zero_tournament, CAPACITY_SUPPORTS))
    == {6: 64},
    "completed triangle count mismatch",
)
require(
    Counter(fingerprint[2] for fingerprint in map(completed_zero_tournament, CAPACITY_SUPPORTS))
    == {(6,): 64},
    "completed SCC mismatch",
)


LOCAL_TEXT = "\n".join(LOCAL_MASK_PAYLOAD) + "\n"
SUPPORT_TEXT = "\n".join(SUPPORT_PAYLOAD) + "\n"
flip_histogram = Counter(
    fingerprint[3] for fingerprint in map(completed_zero_tournament, CAPACITY_SUPPORTS)
)
hp_histogram = Counter(
    fingerprint[4] for fingerprint in map(completed_zero_tournament, CAPACITY_SUPPORTS)
)

print("THM-958 SCALE-FIVE HAMMING-SIX COMMON-SHEET OBSTRUCTION")
print("arithmetic=integer+bitmask floating_point=none metric_heights=none")
print("SECTION literal exhaustive bank")
print(f"label_sets=924 state_words_per_label_set={len(STATE_WORDS)}")
print(f"tested_by_D5_count={dict(sorted(TESTED_BY_D5.items()))}")
print(f"tested_contexts={sum(TESTED_BY_D5.values())}")
print("surviving_presentations=0 surviving_unit_contexts=0")
print("SECTION capacity reduction")
print("D5_count_2_to_4=impossible_by_singleton_sheet_capacity")
print("D5_count_5=impossible_symmetric_nonzero_graph=K6,6_clique_number=2")
print("D5_count_6_capacity_supports=64 equal_c2_signed_cycle_bank=yes")
print(f"multiplicative_support_orbit_sizes={ORBIT_SIZES}")
print("SECTION owner-obligation carrier")
print("zero_provider_graph=(edges=6,SCCs=(3,3),Hamiltonian_paths=0) count=64")
print("doubling_graph=(edges=6,SCCs=(6,),Hamiltonian_paths=6) count=64")
print("unit_words_per_support=4096 full_owner_covers=0 maximum_satisfied_owners=4")
print("per_support_satisfaction_profile={0:2332,1:1248,2:504,3:0,4:12,5:0,6:0}")
print("forbidden_owner_triples=8=(six_C6_windows)+(two_alternating_zero_SCCs)")
print("maximal_near_solutions=12 failed_pairs=three_C6_opposites multiplicity_each=4")
print("SECTION tournament audit")
print("observable=directed_zero_provider_arc switch=complete_missing_pairs")
print("tie_graph=K3,3 tie_gauge=lexicographically_first_tie_Hamiltonian_path")
print("completed_score_histogram={(1,2,2,3,3,4):64} directed_triangles={6:64} SCCs={(6,):64}")
print(f"sparse_edge_flip_histogram={dict(sorted(flip_histogram.items()))}")
print(f"Hamiltonian_path_histogram={dict(sorted(hp_histogram.items()))}")
print(f"joint_tournament_fingerprints={len(TOURNAMENT_FINGERPRINTS)}")
print("SECTION workload and verdict")
print("metric_root_contexts=0 metric_first_edges=0 terminal_recursion=unnecessary")
print("verdict=common_sheet_bank_empty_scale_five_face_closed")
print("SECTION integrity")
print(f"local_mask_payload_fnv64={fnv64(LOCAL_TEXT.encode()):016x}")
print(f"support_payload_sha256={sha256(SUPPORT_TEXT.encode()).hexdigest()}")
print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
