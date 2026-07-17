#!/usr/bin/env python3
"""Independent structural audit of the THM-960 c=6 AP-centred H6 sheet bank.

This scratch verifier does not perform the 37,710,288-row literal scan in the
companion C++ scout.  It derives the local D=6 grammar, identifies the 64
owner-locally-feasible supports, and computes their exact binary-unit nerve.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations, product


P = 13
C = 6
LABELS = tuple(range(1, P))
D6_UNITS = (1, 5)
FULL = (1 << C) - 1
DOUBLING_RATIOS = frozenset((2, 11))
OPPOSITE_RATIOS = frozenset((5, 8))
CENTRAL_RATIOS = frozenset((6, 7))
VARIABLE_GROUPS = (
    (frozenset((2, 3)), frozenset((10, 11))),
    (frozenset((4, 5)), frozenset((8, 9))),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, unit: int) -> int:
    return next(
        value
        for value in range(P * C)
        if value % P == C * label % P and value % C == unit
    )


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


def sheet_mask(provider: int, bit: int, owner: int) -> int:
    base = crt_base(provider, D6_UNITS[bit])
    inverse = pow(owner, -1, P)
    return sum(
        1 << sheet
        for sheet in range(C)
        if -C
        < centered(base * (inverse + P * sheet), P * C)
        <= C
    )


# Local owner-one table, independently frozen from the literal scout.
LOCAL_TABLE = {}
for bit in range(2):
    grouped: dict[int, list[int]] = {}
    for provider in LABELS:
        grouped.setdefault(sheet_mask(provider, bit, 1), []).append(provider)
    LOCAL_TABLE[bit] = tuple(
        (mask, tuple(grouped[mask])) for mask in sorted(grouped)
    )
require(
    LOCAL_TABLE
    == {
        0: ((0, (12,)), (1, (10, 11)), (2, (8, 9)), (4, (6, 7)),
            (8, (4, 5)), (16, (2, 3)), (32, (1,))),
        1: ((0, (12,)), (1, (2, 3)), (2, (4, 5)), (4, (6, 7)),
            (8, (8, 9)), (16, (10, 11)), (32, (1,))),
    },
    "unexpected D=6 local mask grammar",
)


def structural_support(labels: tuple[int, ...]) -> bool:
    """The unit-free local exact-cover criterion at every owner.

    Besides the self sheet, the five providers must consist of one central
    fixed-sheet provider and two providers in each reflected sheet-pair.
    """
    for owner in labels:
        off_ratios = [ratio(provider, owner) for provider in labels if provider != owner]
        if 12 in off_ratios:
            return False
        if sum(value in CENTRAL_RATIOS for value in off_ratios) != 1:
            return False
        for left, right in VARIABLE_GROUPS:
            if sum(value in left | right for value in off_ratios) != 2:
                return False
    return True


def doubling_edges(labels: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    return frozenset(
        (provider, owner)
        for provider in labels
        for owner in labels
        if provider != owner
        and owner * pow(provider, -1, P) % P in DOUBLING_RATIOS
    )


def cycle_support(labels: tuple[int, ...]) -> bool:
    edges = doubling_edges(labels)
    return (
        len(edges) == 6
        and all(sum((p, o) in edges for p in labels) == 1 for o in labels)
        and all(sum((p, o) in edges for o in labels) == 1 for p in labels)
    )


SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if structural_support(labels)
)
CYCLE_SUPPORTS = tuple(
    labels for labels in combinations(LABELS, 6) if cycle_support(labels)
)
require(len(SUPPORTS) == 64, "unexpected structural support count")
require(SUPPORTS == CYCLE_SUPPORTS, "structural supports differ from signed C6 bank")


def cycle_order(labels: tuple[int, ...]) -> tuple[int, ...]:
    successor = dict(doubling_edges(labels))
    start = min(labels)
    answer = [start]
    while successor[answer[-1]] != start:
        answer.append(successor[answer[-1]])
    require(len(answer) == 6, "doubling support is not a Hamiltonian cycle")
    return tuple(answer)


def owner_constraints(labels: tuple[int, ...], owner: int):
    """Return the two affine bit-pair equations for one owner.

    Bits 0/1 mean units 1/5.  Providers in the same reflected ratio bin must
    use opposite bits; providers in opposite bins must use equal bits.
    """
    answer = []
    for left, right in VARIABLE_GROUPS:
        left_providers = [
            p for p in labels if p != owner and ratio(p, owner) in left
        ]
        right_providers = [
            p for p in labels if p != owner and ratio(p, owner) in right
        ]
        providers = tuple(left_providers + right_providers)
        require(len(providers) == 2, "variable sheet-pair does not have two providers")
        rhs = int(len(left_providers) in (0, 2))
        answer.append((tuple(sorted(providers)), rhs))
    return tuple(answer)


def affine_owner_cover(labels: tuple[int, ...], bits: tuple[int, ...], owner: int) -> bool:
    bit_of = dict(zip(labels, bits))
    return all(bit_of[a] ^ bit_of[b] == rhs for (a, b), rhs in owner_constraints(labels, owner))


def literal_owner_cover(labels: tuple[int, ...], bits: tuple[int, ...], owner: int) -> bool:
    union = 0
    for provider, bit in zip(labels, bits):
        union |= sheet_mask(provider, bit, owner)
    return union == FULL


PROFILE_VARIANTS = Counter()
PAIR_NERVE_VARIANTS = Counter()
MINIMAL_OBSTRUCTION_VARIANTS = Counter()
SUPPORT_ORBIT_REPRESENTATIVES = []
for labels in SUPPORTS:
    owner_index = {owner: index for index, owner in enumerate(labels)}
    masks = []
    for bits in product(range(2), repeat=6):
        literal_mask = sum(
            1 << owner_index[owner]
            for owner in labels
            if literal_owner_cover(labels, bits, owner)
        )
        affine_mask = sum(
            1 << owner_index[owner]
            for owner in labels
            if affine_owner_cover(labels, bits, owner)
        )
        require(literal_mask == affine_mask, "affine grammar differs from literal masks")
        masks.append(literal_mask)

    profile = Counter(mask.bit_count() for mask in masks)
    require(profile == {0: 16, 2: 48}, "unexpected owner-satisfaction profile")
    PROFILE_VARIANTS[tuple(sorted(profile.items()))] += 1

    compatible_pairs = {
        tuple(labels[index] for index in pair)
        for pair in combinations(range(6), 2)
        if any(mask & sum(1 << i for i in pair) == sum(1 << i for i in pair) for mask in masks)
    }
    opposite_pairs = {
        pair
        for pair in combinations(labels, 2)
        if ratio(pair[0], pair[1]) in OPPOSITE_RATIOS
    }
    require(len(compatible_pairs) == 12, "owner nerve does not have 12 edges")
    require(
        set(combinations(labels, 2)) - compatible_pairs == opposite_pairs,
        "owner nerve nonedges are not the ratio-{5,8} matching",
    )
    require(len(opposite_pairs) == 3, "opposite relation is not a matching")
    require(
        all(
            sum(mask == sum(1 << owner_index[o] for o in pair) for mask in masks) == 4
            for pair in compatible_pairs
        ),
        "compatible owner edges do not have four unit witnesses",
    )
    PAIR_NERVE_VARIANTS[(len(compatible_pairs), len(opposite_pairs))] += 1

    cycle = cycle_order(labels)
    cycle_opposites = {
        tuple(sorted((cycle[i], cycle[(i + 3) % 6]))) for i in range(3)
    }
    require(opposite_pairs == cycle_opposites, "ratio matching is not cycle antipodality")
    consecutive_triples = {
        frozenset(cycle[(i + offset) % 6] for offset in range(3))
        for i in range(6)
    }
    alternating_triples = {
        frozenset(cycle[0::2]),
        frozenset(cycle[1::2]),
    }
    transversal_triples = {
        frozenset(triple)
        for triple in combinations(labels, 3)
        if all(not set(pair) <= set(triple) for pair in opposite_pairs)
    }
    require(
        transversal_triples == consecutive_triples | alternating_triples,
        "octahedral triangles do not have 6+2 cycle description",
    )
    require(
        all(
            not any(
                mask & sum(1 << owner_index[o] for o in triple)
                == sum(1 << owner_index[o] for o in triple)
                for mask in masks
            )
            for triple in transversal_triples
        ),
        "a transversal owner triple survived",
    )
    MINIMAL_OBSTRUCTION_VARIANTS[(len(opposite_pairs), len(transversal_triples))] += 1


def support_orbit_sizes() -> tuple[int, ...]:
    unseen = set(SUPPORTS)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {
            tuple(sorted(multiplier * label % P for label in seed))
            for multiplier in LABELS
        }
        require(orbit <= set(SUPPORTS), "multiplication leaves support bank")
        sizes.append(len(orbit))
        unseen -= orbit
    return tuple(sorted(sizes))


def directed_reachable(source, target, labels, edges) -> bool:
    stack, seen = [source], {source}
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


def tournament_fingerprint(labels: tuple[int, ...]):
    sparse = doubling_edges(labels)

    def tie(a, b):
        return (a, b) not in sparse and (b, a) not in sparse

    tie_path = min(
        path
        for path in permutations(labels)
        if all(tie(path[i], path[i + 1]) for i in range(5))
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
        flips += edge in sparse and position[edge[0]] > position[edge[1]]
    scores = tuple(
        sorted(sum((v, w) in tournament for w in labels) for v in labels)
    )
    triangles = sum(
        all(sum((v, w) in tournament for w in triple) == 1 for v in triple)
        for triple in combinations(labels, 3)
    )
    strong = all(
        directed_reachable(a, b, labels, tournament)
        for a in labels for b in labels if a != b
    )
    paths = sum(
        all((path[i], path[i + 1]) in tournament for i in range(5))
        for path in permutations(labels)
    )
    return scores, triangles, (6,) if strong else (), flips, paths


TOURNAMENTS = Counter(tournament_fingerprint(labels) for labels in SUPPORTS)
require({fp[0] for fp in TOURNAMENTS} == {(1, 2, 2, 3, 3, 4)}, "score mismatch")
require({fp[1] for fp in TOURNAMENTS} == {6}, "triangle mismatch")
require({fp[2] for fp in TOURNAMENTS} == {(6,)}, "SCC mismatch")


print("THM-960 SCALE-SIX HAMMING-SIX AFFINE-NERVE AUDIT")
print(f"local_table={LOCAL_TABLE}")
print(f"owner_locally_feasible_supports={len(SUPPORTS)}")
print(f"supports_equal_signed_doubling_C6={SUPPORTS == CYCLE_SUPPORTS}")
print(f"support_multiplier_orbits={support_orbit_sizes()}")
print(f"unit_contexts={len(SUPPORTS) * 64}")
print(f"per_support_owner_profile={{0:16,2:48}}")
print(f"owner_nerve=octahedral_graph_K6_minus_3_cycle_antipodes")
print(f"compatible_edge_witnesses=4_each_across_12_edges")
print(f"minimal_obstructions=3_opposite_pairs_plus_8_octahedral_triangles")
print(f"triangle_cycle_split=6_consecutive_plus_2_alternating")
print(f"PROFILE_VARIANTS={dict(PROFILE_VARIANTS)}")
print(f"PAIR_NERVE_VARIANTS={dict(PAIR_NERVE_VARIANTS)}")
print(f"MINIMAL_OBSTRUCTION_VARIANTS={dict(MINIMAL_OBSTRUCTION_VARIANTS)}")
print(f"tournament_joint_fingerprints={len(TOURNAMENTS)}")
print(f"tournament_scores={(1, 2, 2, 3, 3, 4)}")
print(f"tournament_directed_triangles=6")
print(f"tournament_SCCs={(6,)}")
print(f"tournament_sparse_flip_hist={dict(sorted(Counter(fp[3] for fp in map(tournament_fingerprint, SUPPORTS)).items()))}")
print(f"tournament_Hamilton_path_hist={dict(sorted(Counter(fp[4] for fp in map(tournament_fingerprint, SUPPORTS)).items()))}")
