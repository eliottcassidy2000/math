#!/usr/bin/env python3
"""Exact audit of native root/bank sidecars on the THM-3732 seam.

This program does not import the promoted all-edge census.  It reconstructs
the common reset graph and exact response rows from the older independent
THM-3732 source, which itself rebuilds the THM-3238 certificate bank.

There are three deliberately separate operations below.

* Retention attaches a sidecar while keeping the named row.  It cannot create
  a relation pair.
* Compatibility filters a raw-magnitude pair by equality of a native row
  field (degree, bank cell, or provenance).
* Quotienting replaces the row by a native field.  This is audited separately
  because it can identify distinct rows and is not an augmentation.

For the physical root coordinate, a frozen edge tag is also distinguished
from the lawful chronological update c_r -> c_r +/- 1.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / (
    "04-computation/"
    "fc_cross_face_shared_reset_edge_independent_audit_thm3732.py"
)
BASE_PATH = ROOT / (
    "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
ROOTED_CORE_PATH = ROOT / (
    "01-canon/theorems/"
    "THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-"
    "character-boundary.md"
)

EXPECTED_STATES = 239
EXPECTED_EDGES = 568
EXPECTED_PATHS = 1132
EXPECTED_GRAPH_DIGEST = (
    "8cffe16bf5b2301cf077948da8588edf61edbe837a9015a916b3bea3a963c0e3"
)
AMBIENT_BANK = "I2"
CORE_S = frozenset((2, 11, 16, 18, 22))
CORE_F = frozenset((3, 7, 10, 13, 17, 19, 21))
CORE = CORE_S | CORE_F
CANONICAL_CORE_ROOT = 17


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = load_module("fc_hfc_sidecar_independent_source", AUDIT_PATH)


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def provenance(row):
    if row <= 11:
        return "inherited"
    if row <= 14:
        return "principal"
    return "nonprincipal"


ROW_METADATA = tuple(
    (
        row,
        entry[0],  # partition-bank degree
        entry[1],  # number of shapes in the selected upset
        provenance(row),
    )
    for row, entry in enumerate(audit.base.CERTIFICATE, 1)
)
DEGREE = {row: degree for row, degree, _size, _kind in ROW_METADATA}
UPSET_SIZE = {row: size for row, _degree, size, _kind in ROW_METADATA}
CELL = {row: (degree, size) for row, degree, size, _kind in ROW_METADATA}
PROVENANCE = {row: kind for row, _degree, _size, kind in ROW_METADATA}


@lru_cache(maxsize=None)
def values(face, state):
    a_value, b_value = audit.FACES[face]["parameters"]
    return tuple(audit.row_values(a_value, b_value, state))


@lru_cache(maxsize=None)
def deltas(face, source, target):
    before = values(face, source)
    after = values(face, target)
    return tuple(right - left for left, right in zip(before, after))


def positive_rows(delta):
    return tuple(index + 1 for index, value in enumerate(delta) if value > 0)


def moved_root(edge):
    source, target = edge
    before = audit.as_counts(source)
    after = audit.as_counts(target)
    changed = tuple(
        (root, left, right)
        for root, (left, right) in enumerate(zip(before, after), 1)
        if left != right
    )
    if len(changed) != 1:
        raise RuntimeError(("edge is not a one-root update", edge, changed))
    root, left, right = changed[0]
    if abs(right - left) != 1:
        raise RuntimeError(("edge is not unit root update", edge, changed))
    return root, right - left, left, right, AMBIENT_BANK


def edge_record(edge):
    source, target = edge
    delta12 = deltas("F12", source, target)
    delta13 = deltas("F13", source, target)
    left = positive_rows(delta12)
    right = positive_rows(delta13)
    right_set = set(right)
    relation_id = tuple((row, row) for row in left if row in right_set)
    relation_mag = tuple(
        (i + 1, j + 1, left_value)
        for i, left_value in enumerate(delta12)
        if left_value > 0
        for j, right_value in enumerate(delta13)
        if right_value > 0 and left_value == right_value
    )
    return {
        "edge": edge,
        "left": left,
        "right": right,
        "id": relation_id,
        "mag": relation_mag,
        "root": moved_root(edge),
    }


def full_identity(record):
    rows = tuple(i for i, _j in record["id"])
    return bool(rows) and rows == record["left"] == record["right"]


def paths_from(graph):
    outgoing = defaultdict(list)
    for source, target in graph:
        outgoing[source].append(target)
    return tuple(
        (source, middle, target)
        for source, middle in graph
        for target in outgoing.get(middle, ())
    )


def same_pairs(first, second, key):
    return tuple(sorted(set(first[key]) & set(second[key])))


def root_coarse(record):
    root, direction, _before, _after, bank = record["root"]
    return root, direction, bank


def root_dynamic_composable(first, second):
    r1, d1, _s1, t1, b1 = first["root"]
    r2, d2, s2, _t2, b2 = second["root"]
    return (r1, d1, t1, b1) == (r2, d2, s2, b2)


def core_relation(record):
    return tuple(pair for pair in record["id"] if pair[0] in CORE)


def full_core_identity(record):
    left = tuple(row for row in record["left"] if row in CORE)
    right = tuple(row for row in record["right"] if row in CORE)
    relation = tuple(i for i, _j in core_relation(record))
    return bool(relation) and relation == left == right


def state_height(state):
    counts = audit.as_counts(state)
    return sum(
        audit.l1(counts, audit.FACES[face]["reset"])
        for face in ("F12", "F13")
    )


def longest_edge_chain(edges):
    if not edges:
        return ()
    outgoing = defaultdict(list)
    vertices = set()
    for source, target in edges:
        outgoing[source].append(target)
        vertices.update((source, target))
    best = {vertex: (vertex,) for vertex in vertices}
    for source in sorted(vertices, key=lambda item: -state_height(item)):
        for target in sorted(outgoing[source], key=lambda item: (len(item), item)):
            candidate = best[source] + (target,)
            if len(candidate) > len(best[target]):
                best[target] = candidate
            elif len(candidate) == len(best[target]) and candidate < best[target]:
                best[target] = candidate
    return max(best.values(), key=lambda path: (len(path), path))


def fixed_length_chains(edges, edge_length):
    edge_set = set(edges)
    outgoing = defaultdict(list)
    for source, target in edges:
        outgoing[source].append(target)

    answer = []

    def extend(path):
        if len(path) == edge_length + 1:
            answer.append(tuple(path))
            return
        for target in sorted(outgoing.get(path[-1], ())):
            if (path[-1], target) in edge_set:
                extend(path + [target])

    for source in sorted({source for source, _target in edges}):
        extend([source])
    return tuple(answer)


def quotient_collision(rows, key):
    fibres = defaultdict(list)
    for row in rows:
        fibres[key(row)].append(row)
    return tuple(
        (label, tuple(fibre))
        for label, fibre in sorted(fibres.items(), key=lambda item: repr(item[0]))
        if len(fibre) > 1
    )


def relation_size_distribution(records, relation):
    return tuple(sorted(Counter(len(relation(record)) for record in records).items()))


def main():
    if len(ROW_METADATA) != 22:
        raise RuntimeError(("row count drift", len(ROW_METADATA)))
    if len(set(CELL.values())) != 22:
        raise RuntimeError(("native bank cell is not injective", CELL))
    if len(set(UPSET_SIZE.values())) != 22:
        raise RuntimeError(("native upset size is not injective", UPSET_SIZE))

    states, graph, sinks, outdegree = audit.graph_audit()
    if (len(states), len(graph)) != (EXPECTED_STATES, EXPECTED_EDGES):
        raise RuntimeError(("universe drift", len(states), len(graph)))
    if digest(graph) != EXPECTED_GRAPH_DIGEST:
        raise RuntimeError(("graph digest drift", digest(graph)))

    records = tuple(edge_record(edge) for edge in graph)
    by_edge = {record["edge"]: record for record in records}
    paths = paths_from(graph)
    if len(paths) != EXPECTED_PATHS:
        raise RuntimeError(("path count drift", len(paths)))

    full_edges = tuple(record for record in records if full_identity(record))
    full_groups = defaultdict(list)
    for record in full_edges:
        full_groups[record["id"]].append(record["edge"])
    group_sizes = tuple(sorted(len(edges) for edges in full_groups.values()))
    if group_sizes != (1, 1, 3):
        raise RuntimeError(("full identity grouping drift", group_sizes))
    corridor_relation, corridor_edges = next(
        (relation, tuple(edges))
        for relation, edges in full_groups.items()
        if len(edges) == 3
    )
    corridor = longest_edge_chain(corridor_edges)
    corridor_rows = tuple(i for i, _j in corridor_relation)
    expected_corridor = (
        (1, 1, 1, 1, 2, 2, 2, 3, 4, 5),
        (1, 1, 1, 2, 2, 2, 3, 4, 5),
        (1, 1, 2, 2, 2, 3, 4, 5),
        (1, 2, 2, 2, 3, 4, 5),
    )
    expected_rows = (1, 2, 5, 8, 9, 11, 12, 14, 15, 16, 18, 22)
    if corridor != expected_corridor or corridor_rows != expected_rows:
        raise RuntimeError(("corridor drift", corridor, corridor_rows))

    nonempty_mag = tuple(record for record in records if record["mag"])
    if len(nonempty_mag) != 1:
        raise RuntimeError(("magnitude singleton count drift", len(nonempty_mag)))
    singleton = nonempty_mag[0]
    expected_singleton_edge = (
        (1, 1, 1, 3, 3, 4, 5),
        (1, 1, 3, 3, 4, 5),
    )
    if singleton["edge"] != expected_singleton_edge:
        raise RuntimeError(("magnitude edge drift", singleton["edge"]))
    if singleton["mag"] != ((15, 1, 1440),):
        raise RuntimeError(("magnitude pair drift", singleton["mag"]))

    # The singleton has one chronological successor.  The physical root
    # update composes, but the equal-magnitude row relation does not.
    singleton_successors = tuple(
        record for record in records if record["edge"][0] == singleton["edge"][1]
    )
    if len(singleton_successors) != 1:
        raise RuntimeError(("singleton successor drift", singleton_successors))
    singleton_next = singleton_successors[0]
    if not root_dynamic_composable(singleton, singleton_next):
        raise RuntimeError("singleton root update unexpectedly fails to compose")
    if singleton_next["mag"]:
        raise RuntimeError(("singleton magnitude unexpectedly persists", singleton_next))

    path_records = tuple(
        (path, by_edge[(path[0], path[1])], by_edge[(path[1], path[2])])
        for path in paths
    )

    frozen_id_sizes = Counter()
    coarse_id_sizes = Counter()
    dynamic_id_sizes = Counter()
    coarse_mag_sizes = Counter()
    dynamic_mag_sizes = Counter()
    dynamic_constant_full_id = []
    for path, first, second in path_records:
        id_intersection = same_pairs(first, second, "id")
        mag_first = tuple((i, j) for i, j, _value in first["mag"])
        mag_second = tuple((i, j) for i, j, _value in second["mag"])
        mag_intersection = tuple(sorted(set(mag_first) & set(mag_second)))

        frozen_id_sizes[
            len(id_intersection) if first["root"] == second["root"] else 0
        ] += 1
        coarse_ok = root_coarse(first) == root_coarse(second)
        dynamic_ok = root_dynamic_composable(first, second)
        coarse_id_sizes[len(id_intersection) if coarse_ok else 0] += 1
        dynamic_id_sizes[len(id_intersection) if dynamic_ok else 0] += 1
        coarse_mag_sizes[len(mag_intersection) if coarse_ok else 0] += 1
        dynamic_mag_sizes[len(mag_intersection) if dynamic_ok else 0] += 1
        if (
            dynamic_ok
            and first["id"] == second["id"]
            and full_identity(first)
            and full_identity(second)
        ):
            dynamic_constant_full_id.append(path)

    dynamic_constant_full_id = sorted(dynamic_constant_full_id)
    if tuple(dynamic_constant_full_id) != (
        corridor[:3],
        corridor[1:],
    ):
        raise RuntimeError((
            "dynamic full-window drift", tuple(dynamic_constant_full_id)
        ))
    if any(size for size, count in frozen_id_sizes.items() if size and count):
        raise RuntimeError(("frozen fine root tag unexpectedly persists", frozen_id_sizes))
    if any(size for size, count in dynamic_mag_sizes.items() if size and count):
        raise RuntimeError(("magnitude unexpectedly persists", dynamic_mag_sizes))

    # Compatibility filters on raw equal magnitude.  Degree alone retains the
    # singleton.  The next native certificate field (upset size) separates it.
    def mag_degree(record):
        return tuple(
            (i, j, value) for i, j, value in record["mag"]
            if DEGREE[i] == DEGREE[j]
        )

    def mag_cell(record):
        return tuple(
            (i, j, value) for i, j, value in record["mag"]
            if CELL[i] == CELL[j]
        )

    def mag_upset_size(record):
        return tuple(
            (i, j, value) for i, j, value in record["mag"]
            if UPSET_SIZE[i] == UPSET_SIZE[j]
        )

    def mag_provenance(record):
        return tuple(
            (i, j, value) for i, j, value in record["mag"]
            if PROVENANCE[i] == PROVENANCE[j]
        )

    if sum(bool(mag_degree(record)) for record in records) != 1:
        raise RuntimeError("degree-compatible magnitude singleton drift")
    if any(mag_cell(record) for record in records):
        raise RuntimeError("bank-cell-compatible magnitude should be empty")
    if any(mag_upset_size(record) for record in records):
        raise RuntimeError("upset-size-compatible magnitude should be empty")
    if any(mag_provenance(record) for record in records):
        raise RuntimeError("provenance-compatible magnitude should be empty")

    # Rooted-core restriction is not an augmentation, but it tests the least
    # used canonical FC sidecar without pretending it is defined on all rows.
    core_full_edges = tuple(record for record in records if full_core_identity(record))
    core_full_groups = defaultdict(list)
    for record in core_full_edges:
        core_full_groups[core_relation(record)].append(record["edge"])
    core_best_relation, core_best_edges = max(
        core_full_groups.items(),
        key=lambda item: (len(longest_edge_chain(item[1])), len(item[1]), item[0]),
    )
    core_best_chain = longest_edge_chain(core_best_edges)
    core_constant_full_paths = tuple(
        path
        for path, first, second in path_records
        if (
            core_relation(first) == core_relation(second)
            and full_core_identity(first)
            and full_core_identity(second)
        )
    )
    core_dynamic_constant_full_paths = tuple(
        path
        for path, first, second in path_records
        if (
            root_dynamic_composable(first, second)
            and core_relation(first) == core_relation(second)
            and full_core_identity(first)
            and full_core_identity(second)
        )
    )
    core_dynamic_groups = defaultdict(list)
    for record in core_full_edges:
        core_dynamic_groups[(core_relation(record), root_coarse(record))].append(
            record["edge"]
        )
    core_dynamic_best_key, core_dynamic_best_edges = max(
        core_dynamic_groups.items(),
        key=lambda item: (len(longest_edge_chain(item[1])), len(item[1]), item[0]),
    )
    core_dynamic_best_chain = longest_edge_chain(core_dynamic_best_edges)
    core_dynamic_three_edge_chains = fixed_length_chains(
        core_dynamic_best_edges, 3
    )
    if len(core_dynamic_three_edge_chains) != 2:
        raise RuntimeError((
            "rooted-core dynamic three-edge chain drift",
            core_dynamic_three_edge_chains,
        ))

    corridor_degrees = tuple(sorted({DEGREE[row] for row in corridor_rows}))
    corridor_cells = tuple(CELL[row] for row in corridor_rows)
    corridor_provenance_counts = tuple(sorted(Counter(
        PROVENANCE[row] for row in corridor_rows
    ).items()))
    corridor_core = tuple(row for row in corridor_rows if row in CORE)
    if corridor_core != tuple(sorted(CORE_S)):
        raise RuntimeError(("corridor rooted-core intersection drift", corridor_core))
    if CANONICAL_CORE_ROOT in corridor_rows:
        raise RuntimeError("canonical rooted-core root unexpectedly in corridor")
    if any(row in CORE for row in (15, 1)):
        raise RuntimeError("magnitude singleton unexpectedly in rooted core")

    degree_collisions_all = quotient_collision(range(1, 23), lambda row: DEGREE[row])
    degree_collisions_corridor = quotient_collision(
        corridor_rows, lambda row: DEGREE[row]
    )
    cell_collisions_all = quotient_collision(range(1, 23), lambda row: CELL[row])
    size_collisions_all = quotient_collision(
        range(1, 23), lambda row: UPSET_SIZE[row]
    )
    if cell_collisions_all:
        raise RuntimeError(("bank-cell quotient collision", cell_collisions_all))
    if size_collisions_all:
        raise RuntimeError(("upset-size quotient collision", size_collisions_all))

    corridor_root_updates = tuple(
        by_edge[(corridor[index], corridor[index + 1])]["root"]
        for index in range(3)
    )
    if corridor_root_updates != (
        (1, -1, 4, 3, AMBIENT_BANK),
        (1, -1, 3, 2, AMBIENT_BANK),
        (1, -1, 2, 1, AMBIENT_BANK),
    ):
        raise RuntimeError(("corridor root updates drift", corridor_root_updates))

    print("source_independent_audit_sha256=" + file_digest(AUDIT_PATH))
    print("source_thm3238_sha256=" + file_digest(BASE_PATH))
    print("source_rooted_core_thm3278_sha256=" + file_digest(ROOTED_CORE_PATH))
    print("universe=" + repr((
        len(states), len(graph), len(paths), len(sinks), outdegree,
    )))
    print("graph_digest=" + digest(graph))
    print("row_metadata=" + repr(ROW_METADATA))
    print("row_metadata_digest=" + digest(ROW_METADATA))
    print("native_bank_cell_injective=" + repr((len(set(CELL.values())), 22)))
    print("native_upset_size_injective="
          + repr((len(set(UPSET_SIZE.values())), 22)))
    print("degree_collisions_all=" + repr(degree_collisions_all))
    print("bank_cell_collisions_all=" + repr(cell_collisions_all))
    print("upset_size_collisions_all=" + repr(size_collisions_all))
    print("full_identity_group_sizes=" + repr(group_sizes))
    print("corridor=" + repr(corridor))
    print("corridor_rows=" + repr(corridor_rows))
    print("corridor_root_updates=" + repr(corridor_root_updates))
    print("corridor_row_cells=" + repr(tuple(zip(corridor_rows, corridor_cells))))
    print("corridor_degree_quotient=" + repr((
        len(corridor_rows), len(corridor_degrees), corridor_degrees,
        degree_collisions_corridor,
    )))
    print("corridor_bank_cell_quotient=" + repr((
        len(corridor_rows), len(set(corridor_cells)),
    )))
    print("corridor_provenance_counts=" + repr(corridor_provenance_counts))
    print("corridor_rooted_core=" + repr((
        corridor_core, tuple(row for row in corridor_rows if row in CORE_F),
        CANONICAL_CORE_ROOT in corridor_rows,
    )))
    print("singleton=" + repr((
        singleton["edge"], singleton["left"], singleton["right"],
        singleton["mag"], singleton["root"],
    )))
    print("singleton_row_sidecars=" + repr((
        (15, DEGREE[15], CELL[15], PROVENANCE[15], 15 in CORE),
        (1, DEGREE[1], CELL[1], PROVENANCE[1], 1 in CORE),
    )))
    print("singleton_successor=" + repr((
        singleton_next["edge"], singleton_next["root"],
        singleton_next["id"], singleton_next["mag"],
    )))
    print("R_mag_degree_distribution="
          + repr(relation_size_distribution(records, mag_degree)))
    print("R_mag_bank_cell_distribution="
          + repr(relation_size_distribution(records, mag_cell)))
    print("R_mag_upset_size_distribution="
          + repr(relation_size_distribution(records, mag_upset_size)))
    print("R_mag_provenance_distribution="
          + repr(relation_size_distribution(records, mag_provenance)))
    print("length2_frozen_fine_root_R_id_sizes="
          + repr(tuple(sorted(frozen_id_sizes.items()))))
    print("length2_coarse_root_bank_R_id_sizes="
          + repr(tuple(sorted(coarse_id_sizes.items()))))
    print("length2_dynamic_root_bank_R_id_sizes="
          + repr(tuple(sorted(dynamic_id_sizes.items()))))
    print("length2_coarse_root_bank_R_mag_sizes="
          + repr(tuple(sorted(coarse_mag_sizes.items()))))
    print("length2_dynamic_root_bank_R_mag_sizes="
          + repr(tuple(sorted(dynamic_mag_sizes.items()))))
    print("dynamic_root_bank_constant_full_R_id_windows="
          + repr(tuple(dynamic_constant_full_id)))
    print("rooted_core_R_id_distribution="
          + repr(relation_size_distribution(records, core_relation)))
    print("rooted_core_full_edges=" + repr(len(core_full_edges)))
    print("rooted_core_constant_full_length2=" + repr(len(core_constant_full_paths)))
    print("rooted_core_best_constant_full=" + repr((
        core_best_relation, len(core_best_edges), core_best_chain,
    )))
    print("rooted_core_dynamic_root_constant_full_length2="
          + repr(len(core_dynamic_constant_full_paths)))
    print("rooted_core_dynamic_root_best_constant_full=" + repr((
        core_dynamic_best_key,
        len(core_dynamic_best_edges),
        core_dynamic_best_chain,
    )))
    print("rooted_core_dynamic_root_three_edge_chains="
          + repr(core_dynamic_three_edge_chains))
    print("face_origin_scope=F12_and_F13_have_opposite_static_origin_labels;"
          "ordered_(0,1)_retention_is_inert;equal-label_filter_is_empty")
    print("typed_scope=named_F12/F13_bank-I2;strict_positive_exact_integer_rows;"
          "physical_unit-root_updates;native_certificate_degree/upset-size/"
          "provenance;THM3278_rooted_core_restriction")
    print("FIRST_OBSTRUCTION=coarse_root+degree passes singleton edge and root "
          "update, but next R_mag is empty; native upset-size equality fails "
          "already on singleton because 5!=1 (equivalently bank cells "
          "(5,5)!=(5,1)); canonical rooted-core coordinate is undefined on "
          "rows 15 and 1")
    print("HFC_SCOPE=no maintained map to boundary edge, orientation, basepoint, "
          "primitive constant, sheet, or Krylov operator; identity bijections "
          "and a dynamic root filtration do not establish HFC nilpotence")
    print("NO_CLAIM=FC3;HFC3;other_faces;other_banks;arbitrary_rescaling;"
          "row-changing_history;groupoid;torsor;boundary-current_bridge")
    print("semantic_digest=" + digest((
        corridor,
        corridor_rows,
        corridor_root_updates,
        tuple((record["edge"], record["root"], record["id"], record["mag"])
              for record in records),
        tuple(dynamic_constant_full_id),
        tuple(sorted(dynamic_mag_sizes.items())),
        ROW_METADATA,
        core_best_relation,
        core_best_chain,
    )))
    print("RESULT=PASS;ARTIFACT=THM-3732_ROOT_BANK_SIDECAR")


if __name__ == "__main__":
    main()
