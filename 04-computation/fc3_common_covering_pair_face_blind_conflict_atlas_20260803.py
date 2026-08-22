#!/usr/bin/env python3
"""Exact face-blind conflict atlas for all THM-3249 common row pairs."""

import ast
import contextlib
import hashlib
import importlib.util
import io
import runpy
from collections import Counter
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
CROSS_SCRIPT = ROOT / (
    "04-computation/gmc_cross_support_upset_atlas_holonomy_thm3249.py"
)
CROSS_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_cross_support_upset_atlas_holonomy_thm3249.out"
)
FULL_SCRIPT = ROOT / "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py"
FULL_OUTPUT = ROOT / "05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out"
JOINT_SCRIPT = ROOT / (
    "04-computation/fc3_joint_two_face_policy_sidecar_beach_20260803.py"
)
JOINT_OUTPUT = ROOT / (
    "05-knowledge/results/fc3_joint_two_face_policy_sidecar_beach_20260803.out"
)
DEPENDENCIES = {
    CROSS_SCRIPT:
        "bb901e92687544c69d67a55d057f8293ecbf516b80d491a636c4a62af19eebef",
    CROSS_OUTPUT:
        "76d037d4f0737b37ad48531e23be1a0a37509ce0a3d029d3cefdd42662ec04f2",
    FULL_SCRIPT:
        "3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba",
    FULL_OUTPUT:
        "27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9",
    JOINT_SCRIPT:
        "1d69c98d226648d1289c36c72250126e9dbf561c73aa6c3e846ddf3e1a6b47e8",
    JOINT_OUTPUT:
        "67aa3d2e7e6bc836d5a1d2e79caf16f9d1d9cfb5059cb1666fe2e8b74da792f0",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected,
            ("dependency drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
        "optimization-sensitive assert")
require(not any(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
), "floating literal")
require("dependency_hash_checks=8" in JOINT_OUTPUT.read_text(encoding="utf-8")
        and "all_exact_checks=PASS" in JOINT_OUTPUT.read_text(encoding="utf-8"),
        "repaired joint sidecar transcript drift")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cross = load_module("thm3249_cross_support_conflict_atlas", CROSS_SCRIPT)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    full = runpy.run_path(str(FULL_SCRIPT))
require("all_exact_checks=PASS" in captured.getvalue(), "THM-3244 replay failed")


def padded_counts(state):
    counts = Counter(state)
    return tuple(counts[value] for value in range(1, 9))


small_states = cross.SMALL_STATES
small_state_index = {state: index for index, state in enumerate(small_states)}
small_rows = tuple(
    cross.row_coordinates(cross.T.coefficient_vectors(
        1, cross.BANK, 1, 2, state
    ))
    for state in small_states
)
small_nonreset = frozenset(small_states) - {cross.SMALL_RESET}
small_cover_sets = []
for chart in range(22):
    covered = set()
    for state in small_nonreset:
        source = small_rows[small_state_index[state]][chart]
        if any(
            small_rows[small_state_index[target]][chart] > source
            for target in cross.toward_reset_neighbours(
                state, cross.SMALL_POLES, cross.SMALL_RESET
            )
        ):
            covered.add(state)
    small_cover_sets.append(frozenset(covered))
small_cover_sets = tuple(small_cover_sets)
require(tuple(map(len, small_cover_sets)) == (
    192, 219, 184, 186, 210, 186, 184, 202, 224, 186, 224,
    213, 190, 218, 171, 230, 186, 225, 193, 200, 186, 224,
), "small row-cover census drift")

small_pairs = tuple(
    (left + 1, right + 1)
    for left in range(22) for right in range(left + 1, 22)
    if small_cover_sets[left] | small_cover_sets[right] == small_nonreset
)
require(small_pairs == cross.EXPECTED_SMALL_COVERING_PAIRS,
        "small covering-pair drift")
full_cover_sets = full["row_cover_sets"]
full_pairs = full["covering_pairs"]
common_pairs = tuple(pair for pair in small_pairs if pair in frozenset(full_pairs))
require(common_pairs == cross.EXPECTED_COMMON_COVERING_PAIRS,
        "common covering-pair drift")

small_count_vectors = tuple(padded_counts(state) for state in small_states)
full_count_vectors = full["count_vectors"]
full_index_by_counts = {
    counts: index for index, counts in enumerate(full_count_vectors)
}
require(len(full_index_by_counts) == 4319
        and set(small_count_vectors) <= set(full_count_vectors),
        "common raw count-vector universe drift")


def available_rows_at_small(state):
    return tuple(
        row for row, cover in enumerate(small_cover_sets, 1) if state in cover
    )


def available_rows_at_full(counts):
    index = full_index_by_counts[counts]
    return tuple(
        row for row, cover in enumerate(full_cover_sets, 1) if index in cover
    )


universal_witness = (5,)
universal_counts = padded_counts(universal_witness)
small_witness_rows = available_rows_at_small(universal_witness)
full_witness_rows = available_rows_at_full(universal_counts)
require(small_witness_rows == (2, 8, 9, 11, 12, 14, 16, 18, 22),
        "small singleton availability drift")
require(full_witness_rows == (3, 4, 5, 6, 7, 10, 12, 13, 17, 19, 20, 21),
        "full singleton availability drift")

small_witness_row_set = frozenset(small_witness_rows)
full_witness_row_set = frozenset(full_witness_rows)
small_exclusive_witness_rows = small_witness_row_set - full_witness_row_set
full_exclusive_witness_rows = full_witness_row_set - small_witness_row_set
common_witness_rows = small_witness_row_set & full_witness_row_set
require(small_exclusive_witness_rows == frozenset((2, 8, 9, 11, 14, 16, 18, 22))
        and full_exclusive_witness_rows
        == frozenset((3, 4, 5, 6, 7, 10, 13, 17, 19, 20, 21))
        and common_witness_rows == frozenset((12,)),
        "singleton row partition drift")
require(all(
    (left in small_exclusive_witness_rows
     and right in full_exclusive_witness_rows)
    or (right in small_exclusive_witness_rows
        and left in full_exclusive_witness_rows)
    for left, right in common_pairs
), "common-pair graph lost singleton bipartition")
pair_vertices = frozenset(row for pair in common_pairs for row in pair)
pair_small_side = tuple(sorted(pair_vertices & small_exclusive_witness_rows))
pair_full_side = tuple(sorted(pair_vertices & full_exclusive_witness_rows))
require(pair_small_side == (2, 9, 11, 14, 16, 18, 22)
        and pair_full_side == (3, 7, 10, 13, 17, 19, 21)
        and pair_vertices
        == frozenset(pair_small_side) | frozenset(pair_full_side),
        "common-pair singleton bipartition sides drift")

expected_conflict_counts = (
    19, 13, 16, 16, 14, 23, 11, 11, 8, 19, 8, 24,
    13, 14, 8, 16, 16, 12, 19, 18, 10, 17, 16, 8,
)
pair_records = []
conflict_sets = []
for pair, expected_conflicts in zip(common_pairs, expected_conflict_counts):
    left, right = pair
    small_left_only = small_cover_sets[left - 1] - small_cover_sets[right - 1]
    small_right_only = small_cover_sets[right - 1] - small_cover_sets[left - 1]
    full_left_only = full_cover_sets[left - 1] - full_cover_sets[right - 1]
    full_right_only = full_cover_sets[right - 1] - full_cover_sets[left - 1]
    left_right = []
    right_left = []
    for state in small_nonreset:
        full_index = full_index_by_counts[padded_counts(state)]
        if state in small_left_only and full_index in full_right_only:
            left_right.append(state)
        if state in small_right_only and full_index in full_left_only:
            right_left.append(state)
    left_right = tuple(sorted(left_right, key=lambda state: (len(state), state)))
    right_left = tuple(sorted(right_left, key=lambda state: (len(state), state)))
    conflicts = left_right + right_left
    require(len(conflicts) == expected_conflicts and conflicts,
            ("pair conflict census drift", pair))
    require(not (left_right and right_left),
            ("pair gained both conflict orientations", pair))
    minimum_depth = min(map(len, conflicts))
    minimum_witnesses = tuple(
        state for state in conflicts if len(state) == minimum_depth
    )
    require(minimum_depth == 1 and minimum_witnesses == (universal_witness,),
            ("pair lost unique minimum singleton witness", pair))
    conflict_digest = hashlib.sha256(
        repr(conflicts).encode("ascii")
    ).hexdigest()
    record = (
        pair,
        (len(small_left_only), len(small_right_only)),
        (len(full_left_only), len(full_right_only)),
        (len(left_right), len(right_left)),
        minimum_depth,
        minimum_witnesses,
        conflict_digest,
    )
    pair_records.append(record)
    conflict_sets.append(frozenset(conflicts))
pair_records = tuple(pair_records)
conflict_sets = tuple(conflict_sets)

orientation_census = Counter(
    "small_left/full_right" if record[3][0]
    else "small_right/full_left"
    for record in pair_records
)
require(orientation_census == Counter({
    "small_left/full_right": 9,
    "small_right/full_left": 15,
}), "conflict orientation census drift")
face_blind_pairs = tuple(
    record[0] for record in pair_records if sum(record[3]) == 0
)
require(not face_blind_pairs, "a common pair unexpectedly became face-blind")

conflict_union = frozenset().union(*conflict_sets)
conflict_intersection = frozenset.intersection(*conflict_sets)
conflict_multiplicity = Counter(
    sum(state in locus for locus in conflict_sets)
    for state in conflict_union
)
atlas_digest = hashlib.sha256(repr(pair_records).encode("ascii")).hexdigest()
require(len(conflict_union) == 47, "conflict-union census drift")
require(tuple(sorted(
    conflict_intersection, key=lambda state: (len(state), state)
)) == (
    (5,), (4, 5), (1, 3, 5), (1, 4, 5), (3, 4, 5), (1, 3, 4, 5),
), "universal conflict locus drift")
require(tuple(sorted(conflict_multiplicity.items())) == (
    (1, 6), (2, 5), (3, 8), (4, 3), (5, 5), (6, 3),
    (7, 6), (8, 1), (9, 1), (11, 1), (20, 2), (24, 6),
), "conflict multiplicity histogram drift")
require(atlas_digest
        == "5faf9b8fc55cc49c83caa436d78ed98f4253e1c369e7a625ddce91cd69e61b1b",
        "conflict-atlas digest drift")

require(all(
    state in small_cover_sets[
        (left if state in small_cover_sets[left - 1] else right) - 1
    ]
    for left, right in common_pairs for state in small_nonreset
), "small face-tagged selector construction failed")
require(all(
    index in full_cover_sets[
        (left if index in full_cover_sets[left - 1] else right) - 1
    ]
    for left, right in common_pairs for index in full["nonreset_indices"]
), "full face-tagged selector construction failed")

print("FC3 all-common-pair face-blind conflict atlas")
print("dependency_hash_checks=6")
print("faces=(support_(1,2)_I2:239,support_(1,3)_I2:4319)")
print("common_nonreset_raw_vectors=238")
print("common_covering_pair_count=" + repr(len(common_pairs)))
print("common_covering_pairs=" + repr(common_pairs))
print("singleton_witness=(5):counts=" + repr(universal_counts))
print("singleton_small_available_rows=" + repr(small_witness_rows))
print("singleton_full_available_rows=" + repr(full_witness_rows))
print("singleton_common_available_rows=" + repr(tuple(sorted(common_witness_rows))))
print("singleton_common_pair_small_side=" + repr(pair_small_side))
print("singleton_common_pair_full_side=" + repr(pair_full_side))
print("singleton_common_pair_graph=bipartite_between_face-exclusive_row_sets:PASS")
for record in pair_records:
    print("pair_conflict=" + repr(record))
print("per_pair_conflict_counts=" + repr(tuple(sum(record[3]) for record in pair_records)))
print("orientation_census=" + repr(tuple(sorted(orientation_census.items()))))
print("conflict_union_count=" + repr(len(conflict_union)))
print("conflict_intersection_count=" + repr(len(conflict_intersection)))
print("conflict_intersection=" + repr(tuple(sorted(
    conflict_intersection, key=lambda state: (len(state), state)
))))
print("conflict_pair_multiplicity_histogram="
      + repr(tuple(sorted(conflict_multiplicity.items()))))
print("conflict_atlas_sha256=" + atlas_digest)
print("face_blind_count_vector_selector_pairs=()")
print("static_face_bit=NECESSARY_AND_SUFFICIENT_for_each_of_24_pairs")
print("sidecar_scope=pair_identity_is_given; arbitrary_finite_lookup_selector; no_tree_or_dynamic_memory_optimum")
print("typed_scope=two_promoted_I2_faces_and_their_24_common_covering_pairs; no_FC3_closure")
print("all_exact_checks=PASS")
