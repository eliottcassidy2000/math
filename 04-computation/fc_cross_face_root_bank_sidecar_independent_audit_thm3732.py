#!/usr/bin/env python3
"""Independent primary-atlas replay of the root/bank sidecar probe.

Unlike fc_cross_face_root_bank_sidecar_census_thm3732.py, this audit imports the primary
THM-3732/THM-3286 atlas path rather than the THM-3238 reconstruction.
It recomputes the graph, all response relations, native certificate metadata,
and the targeted root/core compositions.
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
CORE_S = frozenset((2, 11, 16, 18, 22))
CORE_F = frozenset((3, 7, 10, 13, 17, 19, 21))
CORE = CORE_S | CORE_F


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


primary = load_module("fc_hfc_sidecar_primary_source", PRIMARY_PATH)


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


@lru_cache(maxsize=None)
def values(face, state):
    return tuple(primary.row_values(face, state))


def moved_root(edge):
    source, target = edge
    left = Counter(source)
    right = Counter(target)
    changed = tuple(
        (root, left[root], right[root])
        for root in range(1, 9) if left[root] != right[root]
    )
    if len(changed) != 1:
        raise RuntimeError(("non-unit edge support", edge, changed))
    root, before, after = changed[0]
    if abs(after - before) != 1:
        raise RuntimeError(("non-unit edge magnitude", edge, changed))
    return root, after - before, before, after, "I2"


def edge_record(edge):
    source, target = edge
    delta = {
        face: tuple(
            after - before
            for before, after in zip(values(face, source), values(face, target))
        )
        for face in ("F12", "F13")
    }
    positive = {
        face: tuple(index + 1 for index, value in enumerate(delta[face]) if value > 0)
        for face in delta
    }
    identity = tuple(
        (row, row) for row in positive["F12"]
        if row in set(positive["F13"])
    )
    magnitude = tuple(
        (i + 1, j + 1, left)
        for i, left in enumerate(delta["F12"]) if left > 0
        for j, right in enumerate(delta["F13"])
        if right > 0 and left == right
    )
    return edge, positive["F12"], positive["F13"], identity, magnitude, moved_root(edge)


def full_identity(record):
    rows = tuple(i for i, _j in record[3])
    return bool(rows) and rows == record[1] == record[2]


def core_relation(record):
    return tuple(pair for pair in record[3] if pair[0] in CORE)


def full_core(record):
    left = tuple(row for row in record[1] if row in CORE)
    right = tuple(row for row in record[2] if row in CORE)
    rows = tuple(i for i, _j in core_relation(record))
    return bool(rows) and rows == left == right


def root_coarse(record):
    root, direction, _before, _after, bank = record[5]
    return root, direction, bank


def root_composes(first, second):
    r1, d1, _s1, t1, b1 = first[5]
    r2, d2, s2, _t2, b2 = second[5]
    return (r1, d1, t1, b1) == (r2, d2, s2, b2)


def fixed_length_chains(edges, edge_length):
    outgoing = defaultdict(list)
    for source, target in edges:
        outgoing[source].append(target)
    answer = []

    def extend(path):
        if len(path) == edge_length + 1:
            answer.append(tuple(path))
            return
        for target in sorted(outgoing.get(path[-1], ())):
            extend(path + [target])

    for source in sorted({source for source, _target in edges}):
        extend([source])
    return tuple(answer)


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
        targets13 = set(primary.atlas.toward_reset_neighbours(
            state, poles13, reset13
        ))
        graph.extend(
            (state, target)
            for target in primary.atlas.toward_reset_neighbours(
                state, poles12, reset12
            )
            if target in targets13
        )
    graph = tuple(graph)
    if (len(states), len(graph), digest(graph)) != (
        239, 568, EXPECTED_GRAPH_DIGEST
    ):
        raise RuntimeError(("graph drift", len(states), len(graph), digest(graph)))

    records = tuple(edge_record(edge) for edge in graph)
    by_edge = {record[0]: record for record in records}
    outgoing = defaultdict(list)
    for source, target in graph:
        outgoing[source].append(target)
    paths = tuple(
        (source, middle, target)
        for source, middle in graph
        for target in outgoing[middle]
    )
    if len(paths) != 1132:
        raise RuntimeError(("path drift", len(paths)))

    certificate = primary.atlas.cross.T.CERTIFICATE
    degree = {row: entry[0] for row, entry in enumerate(certificate, 1)}
    upset_size = {row: entry[1] for row, entry in enumerate(certificate, 1)}
    bank_cell = {row: (degree[row], upset_size[row]) for row in degree}
    if len(set(upset_size.values())) != 22 or len(set(bank_cell.values())) != 22:
        raise RuntimeError("native bank coordinate lost injectivity")

    full_records = tuple(record for record in records if full_identity(record))
    groups = defaultdict(list)
    for record in full_records:
        groups[record[3]].append(record[0])
    if tuple(sorted(map(len, groups.values()))) != (1, 1, 3):
        raise RuntimeError(("full group drift", groups))
    corridor_relation, corridor_edges = next(
        (relation, edges) for relation, edges in groups.items() if len(edges) == 3
    )
    corridor_rows = tuple(i for i, _j in corridor_relation)
    corridor_chains = fixed_length_chains(corridor_edges, 3)
    if len(corridor_chains) != 1:
        raise RuntimeError(("corridor chain drift", corridor_chains))

    magnitude_records = tuple(record for record in records if record[4])
    if len(magnitude_records) != 1:
        raise RuntimeError(("magnitude count drift", magnitude_records))
    singleton = magnitude_records[0]
    normalized_mag = tuple((i, j, int(value)) for i, j, value in singleton[4])
    if normalized_mag != ((15, 1, 1440),):
        raise RuntimeError(("singleton drift", normalized_mag))
    if degree[15] != degree[1] or upset_size[15] == upset_size[1]:
        raise RuntimeError("singleton bank-sidecar separation drift")

    successor_records = tuple(
        record for record in records if record[0][0] == singleton[0][1]
    )
    if len(successor_records) != 1:
        raise RuntimeError(("singleton successor drift", successor_records))
    singleton_next = successor_records[0]
    if not root_composes(singleton, singleton_next) or singleton_next[4]:
        raise RuntimeError(("singleton chronology drift", singleton_next))

    constant_full_dynamic = []
    magnitude_dynamic_nonempty = 0
    core_constant_full = []
    core_dynamic_full = []
    for path in paths:
        first = by_edge[(path[0], path[1])]
        second = by_edge[(path[1], path[2])]
        dynamic = root_composes(first, second)
        if dynamic and first[3] == second[3] and full_identity(first) and full_identity(second):
            constant_full_dynamic.append(path)
        first_mag = {(i, j) for i, j, _value in first[4]}
        second_mag = {(i, j) for i, j, _value in second[4]}
        magnitude_dynamic_nonempty += bool(dynamic and first_mag & second_mag)
        if core_relation(first) == core_relation(second) and full_core(first) and full_core(second):
            core_constant_full.append(path)
            if dynamic:
                core_dynamic_full.append(path)
    if (len(constant_full_dynamic), magnitude_dynamic_nonempty) != (2, 0):
        raise RuntimeError(("dynamic relation drift", constant_full_dynamic))

    core_full_records = tuple(record for record in records if full_core(record))
    if (len(core_full_records), len(core_constant_full), len(core_dynamic_full)) != (
        24, 10, 4
    ):
        raise RuntimeError((
            "rooted-core census drift",
            len(core_full_records), len(core_constant_full), len(core_dynamic_full),
        ))
    core_dynamic_groups = defaultdict(list)
    for record in core_full_records:
        core_dynamic_groups[(core_relation(record), root_coarse(record))].append(record[0])
    s_key = (
        tuple((row, row) for row in sorted(CORE_S)),
        (1, -1, "I2"),
    )
    s_three_chains = fixed_length_chains(core_dynamic_groups[s_key], 3)
    if len(s_three_chains) != 2:
        raise RuntimeError(("rooted S-corridor drift", s_three_chains))

    summary = (
        corridor_chains,
        corridor_rows,
        tuple(by_edge[edge][5] for edge in corridor_edges),
        singleton[0], normalized_mag, singleton[5], singleton_next[0],
        degree[15], degree[1], upset_size[15], upset_size[1],
        tuple(sorted(constant_full_dynamic)),
        len(core_full_records), len(core_constant_full), len(core_dynamic_full),
        s_three_chains,
    )
    print("primary_source_sha256=" + sha256(PRIMARY_PATH.read_bytes()).hexdigest())
    print("universe=" + repr((len(states), len(graph), len(paths))))
    print("native_upset_size_injective=" + repr((len(set(upset_size.values())), 22)))
    print("corridor=" + repr(corridor_chains[0]))
    print("corridor_rows=" + repr(corridor_rows))
    print("singleton=" + repr((singleton[0], normalized_mag, singleton[5])))
    print("singleton_bank_fields=" + repr((
        (15, degree[15], upset_size[15]),
        (1, degree[1], upset_size[1]),
    )))
    print("singleton_successor=" + repr((singleton_next[0], singleton_next[5], singleton_next[4])))
    print("dynamic_full_identity_length2=" + repr(len(constant_full_dynamic)))
    print("dynamic_magnitude_length2_nonempty=" + repr(magnitude_dynamic_nonempty))
    print("rooted_core_counts=" + repr((
        len(core_full_records), len(core_constant_full), len(core_dynamic_full),
    )))
    print("rooted_core_S_dynamic_three_edge_chains=" + repr(s_three_chains))
    print("audit_summary_digest=" + digest(summary))
    print("RESULT=PASS;ARTIFACT=THM-3732_ROOT_BANK_SIDECAR_AUDIT")


if __name__ == "__main__":
    main()
