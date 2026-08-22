#!/usr/bin/env python3
"""Exact unrestricted-22-row face-blind selector obstruction on two FC3 faces."""

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
PAIR_SCRIPT = ROOT / (
    "04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py"
)
PAIR_OUTPUT = ROOT / (
    "05-knowledge/results/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.out"
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
    PAIR_SCRIPT:
        "c03a837f1ed2fbbadc4c9aaef8609a79b1411a1898a56a57546e5460f1fdca56",
    PAIR_OUTPUT:
        "40d074ea26d838bd43edea52b2cdaffeaf27a6a2cdb0dc9abedcd6156eed0e82",
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
require("conflict_intersection_count=6" in PAIR_OUTPUT.read_text(encoding="utf-8")
        and "all_exact_checks=PASS" in PAIR_OUTPUT.read_text(encoding="utf-8"),
        "THM-3266 pair-atlas transcript drift")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cross = load_module("thm3249_cross_support_unrestricted22", CROSS_SCRIPT)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    full = runpy.run_path(str(FULL_SCRIPT))
require("all_exact_checks=PASS" in captured.getvalue(), "THM-3244 replay failed")


def padded_counts(state):
    counts = Counter(state)
    return tuple(counts[value] for value in range(1, 9))


def multiset_leq(left, right):
    left_counts = padded_counts(left)
    right_counts = padded_counts(right)
    return all(a <= b for a, b in zip(left_counts, right_counts))


def l1_distance(left, right):
    return sum(abs(a - b) for a, b in zip(left, right))


# Reconstruct the support-(1,2), bank-I2 availability bank directly from the
# pinned THM-3249 product-Gamma response formula.
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

# Reconstruct the support-(1,3), bank-I2 bank by replaying the pinned
# THM-3244 certificate and retaining its exact response cover sets.
full_cover_sets = tuple(full["row_cover_sets"])
full_count_vectors = tuple(full["count_vectors"])
full_nonreset_indices = frozenset(full["nonreset_indices"])
full_index_by_counts = {
    counts: index for index, counts in enumerate(full_count_vectors)
}
small_count_vectors = tuple(padded_counts(state) for state in small_states)
small_state_by_counts = {
    padded_counts(state): state for state in small_states
}
require(len(small_states) == 239 and len(full_count_vectors) == 4319,
        "physical-bank census drift")
require(len(full_index_by_counts) == 4319
        and set(small_count_vectors) <= set(full_count_vectors),
        "shared padded-count universe drift")
require(frozenset().union(*small_cover_sets) == small_nonreset,
        "small face lost all-row coverage")
require(frozenset().union(*full_cover_sets) == full_nonreset_indices,
        "full face lost all-row coverage")


def available_rows_at_small(state):
    return tuple(
        row for row, cover in enumerate(small_cover_sets, 1) if state in cover
    )


def available_rows_at_full(counts):
    index = full_index_by_counts[counts]
    return tuple(
        row for row, cover in enumerate(full_cover_sets, 1) if index in cover
    )


# All 239 small vectors occur on the full face.  The small reset has no small
# decision obligation, leaving 238 vectors at which both faces must act.
shared_records = []
for state in sorted(small_nonreset, key=lambda item: (len(item), item)):
    counts = padded_counts(state)
    small_available = available_rows_at_small(state)
    full_available = available_rows_at_full(counts)
    common_available = tuple(sorted(
        frozenset(small_available) & frozenset(full_available)
    ))
    require(small_available and full_available,
            ("empty facewise availability", state))
    shared_records.append((
        state, small_available, full_available, common_available,
    ))
shared_records = tuple(shared_records)
require(len(shared_records) == 238, "joint nonreset universe drift")

availability_digest = hashlib.sha256(
    repr(shared_records).encode("ascii")
).hexdigest()
require(availability_digest
        == "9ee92c1bc68f02cb023031ec3f0cf6eb154408f49ff3e67f39db54e9c706bcd8",
        "unrestricted availability-atlas digest drift")

intersection_histogram = tuple(sorted(Counter(
    len(record[3]) for record in shared_records
).items()))
require(intersection_histogram == (
    (0, 2), (1, 4), (2, 2), (7, 1), (9, 8), (10, 8),
    (11, 11), (12, 11), (13, 15), (14, 10), (15, 7),
    (16, 11), (17, 10), (18, 16), (19, 13), (20, 29),
    (21, 52), (22, 28),
), "common-availability histogram drift")

conflicts = tuple(record for record in shared_records if not record[3])
expected_conflicts = (
    (
        (3, 4, 5),
        (2, 5, 8, 9, 11, 12, 14, 16, 18, 22),
        (3, 4, 6, 7, 10, 13, 17, 19, 20, 21),
        (),
    ),
    (
        (1, 3, 4, 5),
        (2, 5, 9, 11, 12, 14, 16, 18, 22),
        (3, 4, 6, 7, 10, 13, 17, 19, 20, 21),
        (),
    ),
)
require(conflicts == expected_conflicts, "unrestricted conflict locus drift")

conflict_states = tuple(record[0] for record in conflicts)
inclusion_minimal_conflicts = tuple(
    state for state in conflict_states
    if not any(
        other != state and multiset_leq(other, state)
        for other in conflict_states
    )
)
require(inclusion_minimal_conflicts == ((3, 4, 5),),
        "inclusion-minimal conflict drift")
require(multiset_leq((3, 4, 5), (1, 3, 4, 5))
        and (1, 3, 4, 5) in cross.toward_reset_neighbours(
            (3, 4, 5), cross.SMALL_POLES, cross.SMALL_RESET
        ), "conflict chain lost its Q-directed edge")

small_reset_counts = padded_counts(cross.SMALL_RESET)
full_reset_counts = padded_counts(cross.FULL_RESET)
conflict_distance_records = tuple(
    (
        state,
        l1_distance(padded_counts(state), small_reset_counts),
        l1_distance(padded_counts(state), full_reset_counts),
    )
    for state in conflict_states
)
require(conflict_distance_records == (
    ((3, 4, 5), 3, 5),
    ((1, 3, 4, 5), 2, 4),
), "conflict distance classification drift")

# THM-3266 found six states hostile to all 24 common pairs.  Four are repaired
# by an outside row; the remaining two are precisely the unrestricted locus.
pair_universal_states = (
    (5,), (4, 5), (1, 3, 5),
    (1, 4, 5), (3, 4, 5), (1, 3, 4, 5),
)
record_by_state = {record[0]: record for record in shared_records}
pair_universal_intersections = tuple(
    (state, record_by_state[state][3]) for state in pair_universal_states
)
require(pair_universal_intersections == (
    ((5,), (12,)),
    ((4, 5), (5,)),
    ((1, 3, 5), (5,)),
    ((1, 4, 5), (5,)),
    ((3, 4, 5), ()),
    ((1, 3, 4, 5), ()),
), "pair-to-unrestricted hostile control drift")

# A canonical partial face-blind policy chooses the least common row on all
# 236 compatible shared inputs.  A face bit is consulted only at the two
# conflicts; elsewhere both tagged copies deliberately make the same choice.
common_selector_census = Counter()
small_selector = {}
full_selector = {}
for state, small_available, full_available, common_available in shared_records:
    counts = padded_counts(state)
    if common_available:
        choice = min(common_available)
        small_selector[counts] = choice
        full_selector[counts] = choice
        common_selector_census[choice] += 1
    else:
        small_selector[counts] = min(small_available)
        full_selector[counts] = min(full_available)
require(tuple(sorted(common_selector_census.items())) == (
    (1, 184), (2, 25), (3, 21), (5, 5), (12, 1),
), "canonical common-selector census drift")

require(all(
    state in small_cover_sets[small_selector[padded_counts(state)] - 1]
    for state in small_nonreset
), "localized-bit selector failed on small face")

for index in full_nonreset_indices:
    counts = full_count_vectors[index]
    if counts not in full_selector:
        full_selector[counts] = min(available_rows_at_full(counts))
require(all(
    index in full_cover_sets[full_selector[full_count_vectors[index]] - 1]
    for index in full_nonreset_indices
), "localized-bit selector failed on full face")

face_dependent_vectors = tuple(sorted(
    (
        counts for counts in small_selector
        if small_selector[counts] != full_selector[counts]
    ),
    key=lambda counts: (sum(counts), counts),
))
require(face_dependent_vectors == tuple(
    padded_counts(state) for state in conflict_states
), "face-bit dependency locus drift")

print("FC3 unrestricted-22-row two-face selector obstruction")
print("dependency_hash_checks=6")
print("faces=(support_(1,2)_I2:239,support_(1,3)_I2:4319)")
print("shared_raw_count_vectors=239")
print("shared_joint_nonreset_vectors=238")
print("small_reset_shared_but_small_unconstrained=" + repr(small_reset_counts))
print("facewise_all_22_row_coverage=PASS")
print("common_availability_size_histogram=" + repr(intersection_histogram))
print("unrestricted_conflict_count=" + repr(len(conflicts)))
for record in conflicts:
    print("unrestricted_conflict=" + repr(record))
print("conflict_distance_records=(state,small_l1,full_l1):"
      + repr(conflict_distance_records))
print("unique_minimum_pole_cardinality_conflict=(3,4,5)")
print("unique_multiset_inclusion_minimal_conflict=(3,4,5)")
print("unique_closest_to_small_reset_conflict=(1,3,4,5):distance=2")
print("conflicts_form_one_small_Q_directed_edge=PASS")
print("pair_universal_intersections=" + repr(pair_universal_intersections))
print("positive_control=(5):common_available=(12,)")
print("hostile_controls=((3,4,5),(1,3,4,5))")
print("unrestricted_availability_atlas_sha256=" + availability_digest)
print("compatible_face_blind_shared_vectors=236")
print("canonical_common_selector_row_census="
      + repr(tuple(sorted(common_selector_census.items()))))
print("face_blind_unrestricted_22_row_selector=IMPOSSIBLE")
print("minimum_static_face_origin_sidecar_bits=1")
print("minimum_face_bit_dependency_locus_raw_vectors=2")
print("localized_one_bit_selector_cases_checked="
      + repr(len(small_nonreset) + len(full_nonreset_indices)))
print("localized_one_bit_selector=PASS")
print("typed_scope=two_promoted_bank_I2_faces; arbitrary_deterministic_count_vector_selector; no_other_face_or_FC3_closure")
print("all_exact_checks=PASS")
