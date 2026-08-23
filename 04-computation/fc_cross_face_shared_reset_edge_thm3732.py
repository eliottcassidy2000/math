#!/usr/bin/env python3
"""Exact local probe of the first F12/F13 reset-selector seam.

This deliberately imports the audited THM-3286 computation rather than
reimplementing the product-Gamma response bank.  It computes only the two
hostile states and their reset-directed neighbours.
"""

from collections import Counter
import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = ROOT / (
    "04-computation/"
    "fc3_three_support_face_availability_hypergraph_scout_20260803.py"
)
AUDIT_PATH = ROOT / (
    "04-computation/"
    "fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py"
)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


atlas = load_module("thm3286_source", ATLAS_PATH)
cross = atlas.cross
audit = load_module("thm3286_independent_source", AUDIT_PATH)

FACE_PARAMETERS = {
    "F12": (1, 2),
    "F13": (1, 3),
}
N0 = (3, 4, 5)
N1 = (1, 3, 4, 5)


def face_profile(face):
    a_value, b_value = FACE_PARAMETERS[face]
    poles = tuple(cross.T.reduced_poles(
        1, cross.BANK, a_value, b_value
    )[0])
    reset = tuple(sorted(cross.T.residual_roots(
        1, cross.T.dominant_row(1), a_value, b_value
    )))
    return poles, reset


def row_values(face, state):
    a_value, b_value = FACE_PARAMETERS[face]
    vector = cross.T.coefficient_vectors(
        1, cross.BANK, a_value, b_value, state
    )
    return atlas.row_coordinates(vector)


def independent_row_values(face, state):
    face_data = audit.faces[face]
    vector = audit.base.coefficient_vectors(
        1, audit.base.BANK, face_data["a"], face_data["b"], state
    )
    return audit.row_coordinates(vector)


def local_record(face, state):
    poles, reset = face_profile(face)
    neighbours = atlas.toward_reset_neighbours(state, poles, reset)
    source = row_values(face, state)
    witness_fibres = []
    for row_index in range(22):
        witnesses = []
        for target in neighbours:
            delta = row_values(face, target)[row_index] - source[row_index]
            if delta > 0:
                witnesses.append((target, delta))
        if witnesses:
            witness_fibres.append((row_index + 1, tuple(witnesses)))
    return {
        "face": face,
        "state": state,
        "neighbours": neighbours,
        "availability": tuple(row for row, _ in witness_fibres),
        "witness_fibres": tuple(witness_fibres),
    }


def edge_delta(face, source, target):
    before = row_values(face, source)
    after = row_values(face, target)
    return tuple(after_value - before_value
                 for before_value, after_value in zip(before, after))


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def predicted_common_neighbours(state):
    """Coordinate description of N12(state) intersect N13(state)."""
    counts = Counter(state)
    answer = []
    for root in (1, 2, 3, 4, 5):
        direction = 0
        if root == 1 and counts[root] != 1:
            direction = 1 if counts[root] == 0 else -1
        elif root == 2 and counts[root] == 3:
            direction = -1
        elif root == 3 and counts[root] == 0:
            direction = 1
        elif root in (4, 5) and counts[root] == 0:
            direction = 1
        if direction:
            changed = counts.copy()
            changed[root] += direction
            if changed[root] == 0:
                del changed[root]
            answer.append(tuple(sorted(changed.elements())))
    return tuple(answer)


def main():
    records = {
        (face, state): local_record(face, state)
        for face in ("F12", "F13")
        for state in (N0, N1)
    }
    for face in ("F12", "F13"):
        poles, reset = face_profile(face)
        independent_profile = audit.faces[face]
        if (tuple(poles), reset) != (
            tuple(independent_profile["poles"]), independent_profile["reset"]
        ):
            raise RuntimeError(("independent profile mismatch", face))
        print("profile=" + repr((
            face, tuple(sorted(Counter(poles).items())), reset
        )))
        states_to_audit = {N0, N1}
        for state in (N0, N1):
            primary_neighbours = atlas.toward_reset_neighbours(
                state, poles, reset
            )
            independent_neighbours = audit.toward_reset_neighbours(
                state, poles, reset
            )
            if primary_neighbours != independent_neighbours:
                raise RuntimeError(("independent neighbour mismatch", face, state))
            states_to_audit.update(primary_neighbours)
        for state in states_to_audit:
            if row_values(face, state) != independent_row_values(face, state):
                raise RuntimeError(("independent response mismatch", face, state))
    for key in (("F12", N0), ("F13", N0),
                ("F12", N1), ("F13", N1)):
        record = records[key]
        print("local_record=" + repr((
            record["face"], record["state"], record["neighbours"],
            record["availability"], record["witness_fibres"]
        )))

    shared_n0 = tuple(sorted(
        set(records[("F12", N0)]["neighbours"])
        & set(records[("F13", N0)]["neighbours"])
    ))
    shared_n1 = tuple(sorted(
        set(records[("F12", N1)]["neighbours"])
        & set(records[("F13", N1)]["neighbours"])
    ))
    if shared_n0 != (N1,):
        raise RuntimeError(("unexpected source overlap", shared_n0))
    if shared_n1:
        raise RuntimeError(("unexpected target overlap", shared_n1))

    deltas = {
        face: edge_delta(face, N0, N1) for face in ("F12", "F13")
    }
    positive_rows = {
        face: tuple(index + 1 for index, value in enumerate(deltas[face])
                    if value > 0)
        for face in ("F12", "F13")
    }
    identity_relation = tuple(sorted(
        set(positive_rows["F12"]) & set(positive_rows["F13"])
    ))
    equal_delta_relation = tuple(
        (left_index + 1, right_index + 1, left_delta)
        for left_index, left_delta in enumerate(deltas["F12"])
        if left_delta > 0
        for right_index, right_delta in enumerate(deltas["F13"])
        if right_delta > 0 and left_delta == right_delta
    )
    equal_delta_left_projection = tuple(sorted({
        left for left, _, _ in equal_delta_relation
    }))
    equal_delta_right_projection = tuple(sorted({
        right for _, right, _ in equal_delta_relation
    }))

    print("shared_reset_neighbours_at_N0=" + repr(shared_n0))
    print("shared_reset_neighbours_at_N1=" + repr(shared_n1))
    print("shared_edge_deltas=" + repr(tuple(
        (face, deltas[face]) for face in ("F12", "F13")
    )))
    print("shared_edge_positive_rows=" + repr(tuple(
        (face, positive_rows[face]) for face in ("F12", "F13")
    )))
    print("identity_row_transition_fibre=" + repr(identity_relation))
    print("equal_positive_delta_transition_relation="
          + repr(equal_delta_relation))
    print("equal_delta_projection=" + repr((
        equal_delta_left_projection, equal_delta_right_projection
    )))

    availability_n0 = {
        face: records[(face, N0)]["availability"]
        for face in ("F12", "F13")
    }
    availability_n1 = {
        face: records[(face, N1)]["availability"]
        for face in ("F12", "F13")
    }
    same_row_n0 = tuple(sorted(
        set(availability_n0["F12"]) & set(availability_n0["F13"])
    ))
    same_row_n1 = tuple(sorted(
        set(availability_n1["F12"]) & set(availability_n1["F13"])
    ))
    print("same_row_selector_fibres=" + repr((same_row_n0, same_row_n1)))
    print("availability_cardinality_jump=" + repr((
        tuple((face, len(availability_n0[face]))
              for face in ("F12", "F13")),
        tuple((face, len(availability_n1[face]))
              for face in ("F12", "F13")),
    )))

    poles12, reset12 = face_profile("F12")
    poles13, reset13 = face_profile("F13")
    shared_states = tuple(sorted(
        set(atlas.physical_states_product(poles12))
        & set(atlas.physical_states_product(poles13)),
        key=lambda state: (len(state), state),
    ))
    shared_decisions = tuple(
        state for state in shared_states if state not in (reset12, reset13)
    )
    common_graph = []
    common_sinks = []
    outdegree_histogram = Counter()
    for state in shared_states:
        common_neighbours = tuple(
            target for target in atlas.toward_reset_neighbours(
                state, poles12, reset12
            )
            if target in set(atlas.toward_reset_neighbours(
                state, poles13, reset13
            ))
        )
        independent_common_neighbours = tuple(
            target for target in audit.toward_reset_neighbours(
                state, poles12, reset12
            )
            if target in set(audit.toward_reset_neighbours(
                state, poles13, reset13
            ))
        )
        if common_neighbours != independent_common_neighbours:
            raise RuntimeError((
                "independent common graph mismatch", state,
                common_neighbours, independent_common_neighbours,
            ))
        if common_neighbours != predicted_common_neighbours(state):
            raise RuntimeError((
                "coordinate common-neighbour formula failed",
                state, common_neighbours, predicted_common_neighbours(state),
            ))
        outdegree_histogram[len(common_neighbours)] += 1
        if common_neighbours:
            common_graph.extend((state, target) for target in common_neighbours)
        else:
            common_sinks.append(state)
    expected_sinks = tuple(sorted((
        tuple(sorted((1,) + (2,) * count2 + (3,) * count3 + (4, 5)))
        for count2 in (0, 1, 2)
        for count3 in (1, 2)
    ), key=lambda state: (len(state), state)))
    if tuple(common_sinks) != expected_sinks:
        raise RuntimeError(("common sink formula failed", common_sinks))
    print("common_reset_graph=" + repr((
        len(shared_states), len(shared_decisions), len(common_graph),
        tuple(sorted(outdegree_histogram.items())), tuple(common_sinks),
    )))
    print("common_reset_graph_digest=" + digest(tuple(common_graph)))
    print("record_digest=" + digest(tuple(
        (key, records[key]) for key in sorted(records)
    )))
    print("independent_profile_neighbour_response_audit=PASS")
    print("typed_scope=two_named_bank_I2_faces;full_239_state_common_graph;"
          "two_hostile_state_response_fibres;strict_positive_row_response;"
          "reset_directed_one_pole_edges")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
