#!/usr/bin/env python3
"""Independent exact all-edge relation audit for the THM-3732 F12/F13 seam.

The primary THM-3732 probe is deliberately not imported.  Exact response
rows and the reset graph come from its independent audit, which itself
rebuilds the 22 product-Gamma rows from the older THM-3238 coefficient source.

For every common reset-directed edge e, this probe records

  R_id(e)  = {(i,i): Delta_12,i(e)>0 and Delta_13,i(e)>0},
  R_mag(e) = {(i,j): Delta_12,i(e)=Delta_13,j(e)>0}.

On a directed length-two path, an inherited fixed row-pair sidecar is
compatible exactly when the same pair lies in both edge relations.  This is
an explicitly scoped identity update within each face; the computation does
not manufacture an arbitrary row update, relabelling, groupoid, or torsor.
"""

from collections import Counter, defaultdict, deque
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / (
    "04-computation/"
    "fc_cross_face_shared_reset_edge_independent_audit_thm3732.py"
)
BASE_PATH = ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"

EXPECTED_STATES = 239
EXPECTED_EDGES = 568
EXPECTED_GRAPH_DIGEST = (
    "8cffe16bf5b2301cf077948da8588edf61edbe837a9015a916b3bea3a963c0e3"
)


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = load_module("thm3732_independent_source", AUDIT_PATH)


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


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


def edge_record(edge):
    source, target = edge
    left_delta = deltas("F12", source, target)
    right_delta = deltas("F13", source, target)
    left = positive_rows(left_delta)
    right = positive_rows(right_delta)
    relation_id = tuple((row, row) for row in left if row in set(right))
    relation_mag_tagged = tuple(
        (left_index + 1, right_index + 1, left_value)
        for left_index, left_value in enumerate(left_delta)
        if left_value > 0
        for right_index, right_value in enumerate(right_delta)
        if right_value > 0 and left_value == right_value
    )
    relation_mag = tuple((i, j) for i, j, _value in relation_mag_tagged)
    return (
        edge,
        left,
        right,
        relation_id,
        relation_mag,
        relation_mag_tagged,
    )


def is_full_bijection(relation, left, right):
    if len(left) != len(right) or len(relation) != len(left):
        return False
    return (
        {i for i, _j in relation} == set(left)
        and {j for _i, j in relation} == set(right)
        and len({i for i, _j in relation}) == len(relation)
        and len({j for _i, j in relation}) == len(relation)
    )


def is_partial_bijection(relation):
    return (
        len({i for i, _j in relation}) == len(relation)
        and len({j for _i, j in relation}) == len(relation)
    )


def state_height(state):
    counts = audit.as_counts(state)
    return sum(
        audit.l1(counts, audit.FACES[face]["reset"])
        for face in ("F12", "F13")
    )


def longest_path(edges):
    """Return a longest directed path; common edges lower state_height by 2."""
    if not edges:
        return ()
    outgoing = defaultdict(list)
    vertices = set()
    for source, target in edges:
        if state_height(target) != state_height(source) - 2:
            raise RuntimeError(("non-descending common edge", source, target))
        outgoing[source].append(target)
        vertices.update((source, target))
    best = {vertex: (vertex,) for vertex in vertices}
    for source in sorted(vertices, key=lambda v: (-state_height(v), len(v), v)):
        for target in sorted(outgoing[source], key=lambda v: (len(v), v)):
            candidate = best[source] + (target,)
            if len(candidate) > len(best[target]):
                best[target] = candidate
            elif len(candidate) == len(best[target]) and candidate < best[target]:
                best[target] = candidate
    return max(best.values(), key=lambda path: (len(path), tuple(path)))


def weak_components(edges):
    if not edges:
        return ()
    adjacency = defaultdict(set)
    vertices = set()
    for source, target in edges:
        adjacency[source].add(target)
        adjacency[target].add(source)
        vertices.update((source, target))
    unseen = set(vertices)
    answer = []
    while unseen:
        seed = min(unseen, key=lambda v: (len(v), v))
        queue = deque((seed,))
        component = set()
        while queue:
            vertex = queue.popleft()
            if vertex in component:
                continue
            component.add(vertex)
            unseen.discard(vertex)
            queue.extend(adjacency[vertex] - component)
        component_edges = tuple(
            edge for edge in edges
            if edge[0] in component and edge[1] in component
        )
        answer.append((
            tuple(sorted(component, key=lambda v: (len(v), v))),
            component_edges,
        ))
    return tuple(answer)


def edge_family_summary(edges):
    if not edges:
        return None
    components = weak_components(edges)
    best_component = max(
        components,
        key=lambda item: (len(item[1]), len(item[0]), digest(item[1])),
    )
    path = longest_path(edges)
    return (
        len(edges),
        len({vertex for edge in edges for vertex in edge}),
        len(best_component[1]),
        len(best_component[0]),
        len(path) - 1,
        path,
        digest(tuple(edges)),
        digest(best_component[1]),
    )


def best_fixed_pair(records, relation_index):
    by_pair = defaultdict(list)
    for record in records:
        edge = record[0]
        for pair in record[relation_index]:
            by_pair[pair].append(edge)
    if not by_pair:
        return None
    candidates = tuple(
        (pair, edge_family_summary(tuple(edges)))
        for pair, edges in sorted(by_pair.items())
    )
    return max(
        candidates,
        key=lambda item: (
            item[1][4],  # longest chronological path
            item[1][2],  # edges in largest weak component
            item[1][0],  # total edges
            item[0],
        ),
    )


def best_constant_relation(records, relation_index, require_full=False):
    groups = defaultdict(list)
    for record in records:
        edge, left, right = record[:3]
        relation = record[relation_index]
        if not relation:
            continue
        if require_full and not is_full_bijection(relation, left, right):
            continue
        groups[relation].append(edge)
    if not groups:
        return None
    candidates = tuple(
        (relation, edge_family_summary(tuple(edges)))
        for relation, edges in groups.items()
    )
    return max(
        candidates,
        key=lambda item: (
            item[1][4],
            item[1][2],
            item[1][0],
            len(item[0]),
            item[0],
        ),
    )


def relation_distributions(records, relation_index):
    sizes = Counter()
    projection_sizes = Counter()
    nonempty = 0
    partial_bijections = 0
    full_bijections = 0
    for record in records:
        _edge, left, right = record[:3]
        relation = record[relation_index]
        sizes[len(relation)] += 1
        projection_sizes[(
            len({i for i, _j in relation}),
            len({j for _i, j in relation}),
        )] += 1
        nonempty += bool(relation)
        partial_bijections += is_partial_bijection(relation)
        full_bijections += is_full_bijection(relation, left, right)
    return (
        tuple(sorted(sizes.items())),
        tuple(sorted(projection_sizes.items())),
        nonempty,
        partial_bijections,
        full_bijections,
    )


def main():
    states, graph, sinks, outdegree_histogram = audit.graph_audit()
    if len(states) != EXPECTED_STATES or len(graph) != EXPECTED_EDGES:
        raise RuntimeError(("universe mismatch", len(states), len(graph)))
    if digest(graph) != EXPECTED_GRAPH_DIGEST:
        raise RuntimeError(("graph digest mismatch", digest(graph)))

    records = tuple(edge_record(edge) for edge in graph)
    for record in records:
        _edge, left, right, relation_id, relation_mag, _tagged = record
        if relation_id != tuple((i, i) for i in sorted(set(left) & set(right))):
            raise RuntimeError(("identity relation mismatch", record[0]))
        if not is_partial_bijection(relation_id):
            raise RuntimeError(("identity relation is not partial bijection", record[0]))

    outgoing = defaultdict(tuple)
    mutable_outgoing = defaultdict(list)
    record_by_edge = {}
    for record in records:
        edge = record[0]
        record_by_edge[edge] = record
        mutable_outgoing[edge[0]].append(edge[1])
    outgoing.update({
        source: tuple(sorted(targets, key=lambda v: (len(v), v)))
        for source, targets in mutable_outgoing.items()
    })
    paths = tuple(
        (source, middle, target)
        for source, middle in graph
        for target in outgoing.get(middle, ())
    )

    path_id_sizes = Counter()
    path_mag_sizes = Counter()
    path_mag_tagged_sizes = Counter()
    path_id_identical_nonempty = 0
    path_id_constant_full = 0
    path_mag_constant_full = 0
    path_records = []
    for source, middle, target in paths:
        first = record_by_edge[(source, middle)]
        second = record_by_edge[(middle, target)]
        identity_intersection = tuple(sorted(
            set(first[3]) & set(second[3])
        ))
        magnitude_intersection = tuple(sorted(
            set(first[4]) & set(second[4])
        ))
        tagged_intersection = tuple(sorted(
            set(first[5]) & set(second[5]),
            key=lambda triple: (triple[0], triple[1], str(triple[2])),
        ))
        path_id_sizes[len(identity_intersection)] += 1
        path_mag_sizes[len(magnitude_intersection)] += 1
        path_mag_tagged_sizes[len(tagged_intersection)] += 1
        path_id_identical_nonempty += bool(first[3]) and first[3] == second[3]
        path_id_constant_full += (
            first[3] == second[3]
            and is_full_bijection(first[3], first[1], first[2])
            and is_full_bijection(second[3], second[1], second[2])
        )
        path_mag_constant_full += (
            first[4] == second[4]
            and is_full_bijection(first[4], first[1], first[2])
            and is_full_bijection(second[4], second[1], second[2])
        )
        path_records.append((
            (source, middle, target),
            identity_intersection,
            magnitude_intersection,
            tagged_intersection,
        ))

    availability_distribution = Counter(
        (len(record[1]), len(record[2])) for record in records
    )
    equal_availability = sum(
        len(record[1]) == len(record[2]) for record in records
    )

    identity_distribution = relation_distributions(records, 3)
    magnitude_distribution = relation_distributions(records, 4)
    full_identity_edges = tuple(
        (record[0], record[1], record[3])
        for record in records
        if is_full_bijection(record[3], record[1], record[2])
    )
    nonempty_magnitude_edges = tuple(
        (record[0], record[1], record[2], record[5])
        for record in records if record[4]
    )

    print("source_independent_audit_sha256=" + file_digest(AUDIT_PATH))
    print("source_thm3238_sha256=" + file_digest(BASE_PATH))
    print("universe=" + repr((
        len(states), len(graph), len(paths), len(sinks), outdegree_histogram,
    )))
    print("graph_digest=" + digest(graph))
    print("availability_size_pair_distribution="
          + repr(tuple(sorted(availability_distribution.items()))))
    print("equal_availability_edges=" + repr((equal_availability, len(graph))))
    print("R_id_distribution=" + repr(identity_distribution))
    print("R_mag_distribution=" + repr(magnitude_distribution))
    print("full_bijection_edges_R_id=" + repr(full_identity_edges))
    print("nonempty_edges_R_mag=" + repr(nonempty_magnitude_edges))
    print("length2_R_id_intersection_sizes="
          + repr(tuple(sorted(path_id_sizes.items()))))
    print("length2_R_mag_pair_intersection_sizes="
          + repr(tuple(sorted(path_mag_sizes.items()))))
    print("length2_R_mag_tagged_intersection_sizes="
          + repr(tuple(sorted(path_mag_tagged_sizes.items()))))
    print("length2_identical_nonempty_R_id=" + repr(path_id_identical_nonempty))
    print("length2_constant_full_bijection_counts="
          + repr((path_id_constant_full, path_mag_constant_full, len(paths))))
    print("best_fixed_pair_R_id=" + repr(best_fixed_pair(records, 3)))
    print("best_fixed_pair_R_mag=" + repr(best_fixed_pair(records, 4)))
    print("best_constant_relation_R_id="
          + repr(best_constant_relation(records, 3)))
    print("best_constant_relation_R_mag="
          + repr(best_constant_relation(records, 4)))
    print("best_constant_full_bijection_R_id="
          + repr(best_constant_relation(records, 3, require_full=True)))
    print("best_constant_full_bijection_R_mag="
          + repr(best_constant_relation(records, 4, require_full=True)))
    print("edge_relation_digest=" + digest(records))
    print("length2_relation_digest=" + digest(tuple(path_records)))
    print("typed_scope=F12/F13;239_shared_states;568_common_edges;"
          "strict_positive_response;named_identity_or_raw_equal_magnitude;"
          "identity_within_face_row_update_on_length2_paths")
    print("NO_CLAIM=arbitrary_relabelling;row_rescaling;variable_groupoid;"
          "torsor;other_faces;FC3;HFC3")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
