#!/usr/bin/env python3
"""Exact forest-lift boundary obstruction for the THM-3110 Ewens current.

The rank-four zeta current is lifted to canonically oriented four-edge
forests.  A literal lift assigns every forest the sign of the macro
partition under it and a nonnegative magnitude.  The script proves, by
exact one-sign face propagation, that no nonzero literal lift can be a
cycle in either product-Gamma bank.

This is a stopping theorem for the literal lift, not a no-go for signed,
rooted, or local-system lifts.
"""

import importlib.util
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
SPEC = importlib.util.spec_from_file_location("thm3110_schur", UPSTREAM)
THM3110 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(THM3110)


@lru_cache(maxsize=None)
def set_partitions(items):
    """All canonically ordered set partitions of the tuple ``items``."""

    answer = []

    def recurse(position, blocks):
        if position == len(items):
            answer.append(tuple(tuple(block) for block in blocks))
            return
        item = items[position]
        for index in range(len(blocks)):
            blocks[index].append(item)
            recurse(position + 1, blocks)
            blocks[index].pop()
        blocks.append([item])
        recurse(position + 1, blocks)
        blocks.pop()

    recurse(0, [])
    return tuple(answer)


def canonical_partition(blocks):
    return tuple(
        sorted(
            (tuple(sorted(block)) for block in blocks),
            key=lambda block: (block[0], len(block), block),
        )
    )


def colour_type(partition, a_count=5):
    return tuple(
        sorted(
            (
                (
                    sum(item < a_count for item in block),
                    sum(item >= a_count for item in block),
                )
                for block in partition
            ),
            reverse=True,
        )
    )


def ewens_zeta(bank, vertex_count):
    """Reconstruct THM-3110's labelled macro-partition zeta current."""

    universe = tuple(range(vertex_count))
    partitions = set_partitions(universe)
    multiplicities = Counter(colour_type(partition) for partition in partitions)
    coefficient_by_type = {
        tuple(sorted(row, reverse=True)): coefficient for coefficient, row in bank
    }

    zeta = defaultdict(Fraction)
    for partition in partitions:
        row_type = colour_type(partition)
        coefficient = coefficient_by_type.get(row_type, 0)
        if not coefficient:
            continue
        weight = Fraction(coefficient, multiplicities[row_type])
        refinements = tuple(set_partitions(block) for block in partition)
        for choices in product(*refinements):
            refined = canonical_partition(
                sum((list(piece) for piece in choices), [])
            )
            zeta[refined] += weight

    nonzero = {
        partition: zeta[partition] for partition in partitions if zeta[partition]
    }
    require(
        all(vertex_count - len(partition) == 4 for partition in nonzero),
        "zeta current escaped rank four",
    )
    return partitions, nonzero


def component_partition(vertex_count, edges):
    """Return the component partition, or None when ``edges`` has a cycle."""

    parent = list(range(vertex_count))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right in edges:
        root_left = find(left)
        root_right = find(right)
        if root_left == root_right:
            return None
        parent[root_right] = root_left

    blocks = defaultdict(list)
    for vertex in range(vertex_count):
        blocks[find(vertex)].append(vertex)
    return canonical_partition(blocks.values())


def forest_incidence(vertex_count, zeta):
    """Build all literal variables and their oriented three-face incidences."""

    complete_edges = tuple(combinations(range(vertex_count), 2))
    forest_count = 0
    variables = []
    for forest in combinations(complete_edges, 4):
        partition = component_partition(vertex_count, forest)
        if partition is None:
            continue
        forest_count += 1
        if partition in zeta:
            sign = 1 if zeta[partition] > 0 else -1
            variables.append((forest, partition, sign))

    faces = defaultdict(list)
    for variable, (forest, _, sign) in enumerate(variables):
        for position in range(4):
            face = forest[:position] + forest[position + 1 :]
            faces[face].append((variable, sign * (-1) ** position))
    return forest_count, tuple(variables), dict(faces)


def simultaneous_one_sign_pruning(variables, faces):
    """Propagate exact zeros from faces whose live incidences have one sign."""

    active = set(range(len(variables)))
    rounds = []
    while active:
        forcing_faces = []
        for face, incidences in faces.items():
            live = tuple(
                (variable, sign)
                for variable, sign in incidences
                if variable in active
            )
            if live and len({sign for _, sign in live}) == 1:
                forcing_faces.append((face, live))
        if not forcing_faces:
            break
        killed = {
            variable
            for _, incidences in forcing_faces
            for variable, _ in incidences
        }
        rounds.append((len(forcing_faces), len(killed)))
        active.difference_update(killed)
    return tuple(rounds), len(active)


def example_face_profile(vertex_count, zeta, faces):
    """Check the three one-sign faces killing one explicit three-tree fibre."""

    singleton_blocks = [(vertex,) for vertex in range(7, vertex_count)]
    target = canonical_partition(
        [(0, 1, 2), (3, 4), (5, 6)] + singleton_blocks
    )
    witness_faces = (
        ((0, 1), (0, 2), (5, 6)),
        ((0, 1), (1, 2), (5, 6)),
        ((0, 2), (1, 2), (5, 6)),
    )
    profiles = []
    for face in witness_faces:
        incidences = faces[face]
        signs = {sign for _, sign in incidences}
        require(len(signs) == 1, "illustrative face is not one-signed")
        profiles.append((len(incidences), next(iter(signs))))
    return target, zeta[target], tuple(profiles)


EXPECTED = {
    8: {
        "partition_count": 4140,
        "nonzero": 480,
        "signs": (285, 195),
        "forest_count": 18865,
        "variables": 1440,
        "faces": 1860,
        "rounds": ((1203, 1368), (120, 72)),
        "example_weight": Fraction(1, 30),
        "example_profiles": ((3, 1), (3, 1), (3, 1)),
    },
    9: {
        "partition_count": 21147,
        "nonzero": 1620,
        "signs": (720, 900),
        "forest_count": 55755,
        "variables": 4860,
        "faces": 3960,
        "rounds": ((2106, 4038), (496, 551), (228, 244), (73, 27)),
        "example_weight": Fraction(-1, 60),
        "example_profiles": ((6, -1), (6, -1), (6, -1)),
    },
}


records = []
for bank_index, bank in enumerate(THM3110.BANKS, 1):
    vertex_count = 7 + bank_index
    expected = EXPECTED[vertex_count]
    partitions, zeta = ewens_zeta(bank, vertex_count)
    sign_count = (
        sum(weight > 0 for weight in zeta.values()),
        sum(weight < 0 for weight in zeta.values()),
    )
    forest_count, variables, faces = forest_incidence(vertex_count, zeta)
    rounds, remaining = simultaneous_one_sign_pruning(variables, faces)
    target, target_weight, profiles = example_face_profile(
        vertex_count, zeta, faces
    )

    require(len(partitions) == expected["partition_count"], "Bell count drift")
    require(len(zeta) == expected["nonzero"], "nonzero zeta count drift")
    require(sign_count == expected["signs"], "zeta sign count drift")
    require(forest_count == expected["forest_count"], "forest count drift")
    require(len(variables) == expected["variables"], "literal variable count drift")
    require(len(faces) == expected["faces"], "incident face count drift")
    require(rounds == expected["rounds"], "one-sign propagation drift")
    require(remaining == 0, "one-sign propagation left a live variable")
    require(target_weight == expected["example_weight"], "example weight drift")
    require(profiles == expected["example_profiles"], "example face drift")
    records.append(
        (
            bank_index,
            vertex_count,
            len(zeta),
            sign_count,
            forest_count,
            len(variables),
            len(faces),
            rounds,
            target,
            target_weight,
            profiles,
        )
    )


print("rank4_forest_holotopy_boundary=exact")
for (
    bank_index,
    vertex_count,
    nonzero,
    signs,
    forest_count,
    variable_count,
    face_count,
    rounds,
    target,
    target_weight,
    profiles,
) in records:
    round_text = ",".join(f"{faces}/{killed}" for faces, killed in rounds)
    profile_text = ",".join(f"{count}/{sign:+d}" for count, sign in profiles)
    print(
        f"I{bank_index}=K{vertex_count}:zeta_nonzero={nonzero}:"
        f"positive={signs[0]}:negative={signs[1]}:"
        f"all_4forests={forest_count}:literal_variables={variable_count}:"
        f"incident_3faces={face_count}:pruning_faces/killed={round_text}:"
        "remaining=0"
    )
    print(
        f"I{bank_index}_example=partition={target}:weight={target_weight}:"
        f"three_face_profiles=count/sign={profile_text}"
    )
print("literal_same_sign_forest_cycle=IMPOSSIBLE")
print("survivor=rooted_or_signed_local_system_lift_required")
print("all_exact_checks=PASS")
