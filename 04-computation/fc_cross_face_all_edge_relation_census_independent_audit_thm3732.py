#!/usr/bin/env python3
"""Independent replay of the all-edge F12/F13 relation census.

This path imports the primary THM-3732 implementation (and hence the audited
THM-3286 atlas), not the THM-3238 reconstruction used by the scratch probe.
It independently rebuilds every edge relation and length-two intersection,
then compares stable digests and the exceptional full/metric edges.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / (
    "04-computation/fc_cross_face_shared_reset_edge_thm3732.py"
)

EXPECTED_GRAPH_DIGEST = (
    "8cffe16bf5b2301cf077948da8588edf61edbe837a9015a916b3bea3a963c0e3"
)
EXPECTED_EDGE_DIGEST = (
    "0d924494e3a689db3f818418d6696eb349ed20eb1de6cc1876448f464a6a6327"
)
EXPECTED_PATH_DIGEST = (
    "2af1486303b7d5332d7e4911aedbd6a16157086f751d9b894be0985ab11a49aa"
)
EXPECTED_ID_SIZES = (
    (0, 82), (1, 105), (2, 40), (3, 61), (4, 36), (5, 28), (6, 57),
    (7, 31), (8, 39), (9, 36), (10, 20), (11, 27), (12, 6),
)
EXPECTED_PATH_ID_SIZES = (
    (0, 733), (1, 120), (2, 66), (3, 65), (4, 27), (5, 27),
    (6, 27), (7, 14), (8, 13), (9, 10), (10, 16), (11, 12), (12, 2),
)
EXPECTED_RAW_MAG_EDGE = (
    ((1, 1, 1, 3, 3, 4, 5), (1, 1, 3, 3, 4, 5)),
    ((15, 1, 1440),),
)


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


primary = load_module("thm3732_primary_source", PRIMARY_PATH)


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


@lru_cache(maxsize=None)
def values(face, state):
    return tuple(primary.row_values(face, state))


def edge_record(edge):
    source, target = edge
    face_deltas = {}
    positive = {}
    for face in ("F12", "F13"):
        face_deltas[face] = tuple(
            after - before
            for before, after in zip(values(face, source), values(face, target))
        )
        positive[face] = tuple(
            index + 1 for index, value in enumerate(face_deltas[face])
            if value > 0
        )
    relation_id = tuple(
        (row, row) for row in positive["F12"]
        if row in set(positive["F13"])
    )
    relation_mag_tagged = tuple(
        (i + 1, j + 1, left)
        for i, left in enumerate(face_deltas["F12"]) if left > 0
        for j, right in enumerate(face_deltas["F13"])
        if right > 0 and left == right
    )
    relation_mag = tuple((i, j) for i, j, _value in relation_mag_tagged)
    return (
        edge,
        positive["F12"],
        positive["F13"],
        relation_id,
        relation_mag,
        relation_mag_tagged,
    )


def full_bijection(record, relation_index):
    relation = record[relation_index]
    left, right = record[1], record[2]
    return (
        len(left) == len(right) == len(relation)
        and {i for i, _j in relation} == set(left)
        and {j for _i, j in relation} == set(right)
    )


def main():
    poles12, reset12 = primary.face_profile("F12")
    poles13, reset13 = primary.face_profile("F13")
    states = tuple(sorted(
        set(primary.atlas.physical_states_product(poles12))
        & set(primary.atlas.physical_states_product(poles13)),
        key=lambda state: (len(state), state),
    ))
    graph = []
    for state in states:
        right = set(primary.atlas.toward_reset_neighbours(
            state, poles13, reset13
        ))
        graph.extend(
            (state, target)
            for target in primary.atlas.toward_reset_neighbours(
                state, poles12, reset12
            )
            if target in right
        )
    graph = tuple(graph)
    if (len(states), len(graph), digest(graph)) != (
        239, 568, EXPECTED_GRAPH_DIGEST
    ):
        raise RuntimeError("primary-path graph mismatch")

    records = tuple(edge_record(edge) for edge in graph)
    id_sizes = tuple(sorted(Counter(len(record[3]) for record in records).items()))
    if id_sizes != EXPECTED_ID_SIZES:
        raise RuntimeError(("R_id distribution mismatch", id_sizes))
    raw_mag = tuple(
        (record[0], record[5]) for record in records if record[4]
    )
    normalized_raw_mag = tuple(
        (edge, tuple((i, j, int(value)) for i, j, value in relation))
        for edge, relation in raw_mag
    )
    if normalized_raw_mag != (EXPECTED_RAW_MAG_EDGE,):
        raise RuntimeError(("R_mag exceptional edge mismatch", normalized_raw_mag))

    outgoing = defaultdict(list)
    record_by_edge = {}
    for record in records:
        source, target = record[0]
        outgoing[source].append(target)
        record_by_edge[record[0]] = record
    paths = tuple(
        (source, middle, target)
        for source, middle in graph
        for target in sorted(outgoing[middle], key=lambda v: (len(v), v))
    )
    path_records = []
    id_path_sizes = Counter()
    raw_mag_path_nonempty = 0
    identical_nonempty_id_paths = 0
    constant_full_id_paths = 0
    for source, middle, target in paths:
        first = record_by_edge[(source, middle)]
        second = record_by_edge[(middle, target)]
        identity_intersection = tuple(sorted(set(first[3]) & set(second[3])))
        magnitude_intersection = tuple(sorted(set(first[4]) & set(second[4])))
        tagged_intersection = tuple(sorted(
            set(first[5]) & set(second[5]),
            key=lambda triple: (triple[0], triple[1], str(triple[2])),
        ))
        id_path_sizes[len(identity_intersection)] += 1
        raw_mag_path_nonempty += bool(magnitude_intersection or tagged_intersection)
        identical_nonempty_id_paths += bool(first[3]) and first[3] == second[3]
        constant_full_id_paths += (
            first[3] == second[3]
            and full_bijection(first, 3)
            and full_bijection(second, 3)
        )
        path_records.append((
            (source, middle, target),
            identity_intersection,
            magnitude_intersection,
            tagged_intersection,
        ))

    if tuple(sorted(id_path_sizes.items())) != EXPECTED_PATH_ID_SIZES:
        raise RuntimeError(("length-two R_id mismatch", id_path_sizes))
    if raw_mag_path_nonempty:
        raise RuntimeError("raw-magnitude relation unexpectedly composes")
    if identical_nonempty_id_paths != 57:
        raise RuntimeError((
            "identical nonempty identity path mismatch",
            identical_nonempty_id_paths,
        ))
    if constant_full_id_paths != 2:
        raise RuntimeError(("constant full identity path mismatch", constant_full_id_paths))
    if digest(records) != EXPECTED_EDGE_DIGEST:
        raise RuntimeError(("edge digest mismatch", digest(records)))
    if digest(tuple(path_records)) != EXPECTED_PATH_DIGEST:
        raise RuntimeError(("path digest mismatch", digest(tuple(path_records))))

    full_id = tuple(record for record in records if full_bijection(record, 3))
    full_mag = tuple(record for record in records if full_bijection(record, 4))
    print("primary_source_sha256=" + sha256(PRIMARY_PATH.read_bytes()).hexdigest())
    print("universe=" + repr((len(states), len(graph), len(paths))))
    print("full_bijection_counts=" + repr((len(full_id), len(full_mag))))
    print("raw_magnitude_exception=" + repr(normalized_raw_mag))
    print("length2_identical_nonempty_identity_paths="
          + repr(identical_nonempty_id_paths))
    print("length2_constant_full_identity_paths=" + repr(constant_full_id_paths))
    print("edge_relation_digest=" + digest(records))
    print("length2_relation_digest=" + digest(tuple(path_records)))
    print("independent_primary_path_replay=PASS")


if __name__ == "__main__":
    main()
