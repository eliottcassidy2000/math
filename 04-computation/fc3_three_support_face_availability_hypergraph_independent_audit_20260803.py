#!/usr/bin/env python3
"""Independent hostile audit of the THM-3286 three-face Helly candidate.

This script does not import the candidate companion, THM-3249's face
constructor, or any frozen availability table.  It starts from the pinned
THM-3238 product--Gamma coefficient formula, rederives the three pole
multisets and resets for supports (1,2), (1,3), and (2,3), reconstructs all
8,397 physical response vectors, and forms the twenty-two row-availability
sets by direct strict one-pole comparison.

It then independently enumerates every pair and triple fibre, proper Helly
defects, static origin colourings, and the raw-vector dependency locus.  The
THM-3286 and THM-3278 artifacts are pinned for comparison only; no constants
or constructors are imported from them.

The result is a finite exact response-bank audit, not FC(3), SFC(3), GMC,
positivity, or a sequential physical controller.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from functools import lru_cache, reduce
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from multiprocessing import Pool, freeze_support
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"

BASE_SCRIPT = COMP / "gmc_complete_physical_bank_unique_reset_thm3238.py"
BASE_OUTPUT = RESULTS / "gmc_complete_physical_bank_unique_reset_thm3238.out"
CANDIDATE_SCRIPT = (
    COMP / "fc3_three_support_face_availability_hypergraph_scout_20260803.py"
)
CANDIDATE_OUTPUT = (
    RESULTS / "fc3_three_support_face_availability_hypergraph_scout_20260803.out"
)
THM3278_SCRIPT = COMP / "fc3_selector_origin_bipartition_phase_bridge_thm3278.py"
THM3278_OUTPUT = RESULTS / "fc3_selector_origin_bipartition_phase_bridge_thm3278.out"

PINNED = {
    BASE_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    BASE_OUTPUT:
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    CANDIDATE_SCRIPT:
        "c31eff8a10c6a6e1ab7e3ea6759388b02cbd9591695182f2d3756e08da38c8c3",
    CANDIDATE_OUTPUT:
        "a3e44f6d2eb0e26386e399e7d582e6eb67512125233eb1b49d2b62c4d0869e05",
    THM3278_SCRIPT:
        "07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f",
    THM3278_OUTPUT:
        "5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(path)) for path in PINNED)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "a pinned base, candidate, or THM-3278 artifact changed",
)

source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(
    isinstance(node, ast.Assert) for node in ast.walk(source_tree)
)
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(source_tree)
)
require(
    assert_nodes == 0 and float_literals == 0,
    "source validity gate found assert or floating-point literals",
)


base_spec = spec_from_file_location("thm3238_for_thm3286_independent_audit", BASE_SCRIPT)
base = module_from_spec(base_spec)
base_spec.loader.exec_module(base)


@lru_cache(maxsize=None)
def integer_partitions(total, maximum=None):
    """Independent partition constructor, in decreasing-part order."""
    if total == 0:
        return ((),)
    bound = total if maximum is None else min(total, maximum)
    return tuple(
        (first,) + tail
        for first in range(bound, 0, -1)
        for tail in integer_partitions(total - first, first)
    )


def audit_coarsens(fine, coarse):
    """Whether ``coarse`` is obtained by merging the parts of ``fine``."""
    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    capacities = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, remaining):
        if index == len(pieces):
            return not any(remaining)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(remaining):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            changed = list(remaining)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, capacities)


def build_audit_upsets():
    answer = []
    for degree, expected_size, generators, _factor in base.CERTIFICATE:
        shapes = integer_partitions(degree)
        upset = frozenset(
            shape for shape in shapes
            if any(audit_coarsens(generator, shape) for generator in generators)
        )
        require(
            len(upset) == expected_size,
            ("independent upset size drift", degree, generators),
        )
        answer.append(upset)
    return tuple(answer)


AUDIT_UPSETS = build_audit_upsets()
require(
    AUDIT_UPSETS == base.UPSETS,
    "independently reconstructed upset atlas differs from THM-3238",
)


def row_coordinates(vectors):
    return tuple(
        sum(vectors[entry[0]][shape] for shape in upset)
        for entry, upset in zip(base.CERTIFICATE, AUDIT_UPSETS)
    )


FACE_PARAMETERS = (
    ("F12", 1, 2),
    ("F13", 1, 3),
    ("F23", 2, 3),
)
FACE_ORDER = tuple(name for name, _a, _b in FACE_PARAMETERS)


def physical_states(poles):
    multiplicities = Counter(poles)
    roots = tuple(sorted(multiplicities))
    return tuple(sorted((
        tuple(
            root
            for root, count in zip(roots, choices)
            for _ in range(count)
        )
        for choices in product(
            *(range(multiplicities[root] + 1) for root in roots)
        )
        if any(choices)
    ), key=lambda state: (len(state), state)))


def padded_counts(state):
    counts = Counter(state)
    return tuple(counts[root] for root in range(1, 9))


def toward_reset_neighbours(state, poles, reset):
    """All nonempty one-pole moves decreasing padded-count l1 distance."""
    current = Counter(state)
    pole_counts = Counter(poles)
    reset_counts = Counter(reset)
    answer = []
    for root in sorted(set(pole_counts) | set(reset_counts)):
        changed = current.copy()
        if current[root] > reset_counts[root]:
            changed[root] -= 1
            if changed[root] == 0:
                del changed[root]
        elif current[root] < reset_counts[root]:
            require(
                current[root] < pole_counts[root],
                "reset-directed insertion left the pole bank",
            )
            changed[root] += 1
        else:
            continue
        target = tuple(sorted(changed.elements()))
        if target:
            require(
                sum(abs(left - right) for left, right in zip(
                    padded_counts(target), padded_counts(reset)
                ))
                + 1
                == sum(abs(left - right) for left, right in zip(
                    padded_counts(state), padded_counts(reset)
                )),
                "one-pole move did not lower reset distance exactly",
            )
            answer.append(target)
    return tuple(answer)


faces = {}
for name, a_value, b_value in FACE_PARAMETERS:
    poles, numerator = base.reduced_poles(
        1, base.BANK, a_value, b_value
    )
    reset = tuple(sorted(base.residual_roots(
        1, base.dominant_row(1), a_value, b_value
    )))
    states = physical_states(poles)
    require(
        Counter(reset) <= Counter(poles)
        and len(states) == len({padded_counts(state) for state in states}),
        f"{name} reset or padded-count universe is invalid",
    )
    faces[name] = {
        "a": a_value,
        "b": b_value,
        "poles": poles,
        "pole_profile": tuple(sorted(Counter(poles).items())),
        "numerator": numerator,
        "reset": reset,
        "states": states,
    }

require(
    faces["F12"]["poles"]
    == (5, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1)
    and faces["F12"]["reset"] == (1, 2, 2, 3, 4, 5)
    and len(faces["F12"]["states"]) == 239,
    "F12 product-Gamma universe changed",
)
require(
    faces["F13"]["poles"]
    == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1)
    and faces["F13"]["reset"] == (1, 3, 3, 4, 5, 6, 7, 8)
    and len(faces["F13"]["states"]) == 4319,
    "F13 product-Gamma universe changed",
)
require(
    faces["F23"]["poles"]
    == (8, 7, 6, 5, 5, 5, 4, 4, 4, 3, 2, 2, 2, 2, 1, 1)
    and faces["F23"]["pole_profile"]
    == ((1, 2), (2, 4), (3, 1), (4, 3),
        (5, 3), (6, 1), (7, 1), (8, 1))
    and faces["F23"]["reset"] == (2, 3, 4, 4, 5, 6, 7, 8)
    and len(faces["F23"]["states"]) == 3839,
    "F23 product-Gamma universe changed",
)


def response_task(task):
    face_name, a_value, b_value, state = task
    vectors = base.coefficient_vectors(
        1, base.BANK, a_value, b_value, state
    )
    return (face_name, state), row_coordinates(vectors)


def intersection_all(sets):
    sets = tuple(sets)
    require(sets, "empty family passed to intersection_all")
    return frozenset.intersection(*(frozenset(values) for values in sets))


def stable_digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def bits(values):
    return "".join(map(str, values))


def main():
    tasks = tuple(
        (name, face["a"], face["b"], state)
        for name, face in faces.items()
        for state in face["states"]
    )
    require(len(tasks) == 239 + 4319 + 3839 == 8397,
            "response task universe drift")
    with Pool(processes=4) as pool:
        responses = dict(pool.imap_unordered(
            response_task, tasks, chunksize=4
        ))
    require(len(responses) == len(tasks), "a response task was lost")

    availability = {}
    availability_records = {}
    for name in FACE_ORDER:
        face = faces[name]
        states = face["states"]
        state_set = frozenset(states)
        reset = face["reset"]
        require(reset in state_set, f"{name} reset absent")
        face_availability = {}
        for state in states:
            if state == reset:
                require(
                    responses[name, state] == (0,) * 22,
                    f"{name} reset response is nonzero",
                )
                continue
            neighbours = toward_reset_neighbours(
                state, face["poles"], reset
            )
            require(
                neighbours and all(target in state_set for target in neighbours),
                ("invalid reset-directed neighbourhood", name, state),
            )
            source = responses[name, state]
            available = tuple(
                row
                for row in range(1, 23)
                if any(
                    responses[name, target][row - 1] > source[row - 1]
                    for target in neighbours
                )
            )
            require(available, ("empty facewise availability", name, state))
            face_availability[state] = available
        availability[name] = face_availability
        availability_records[name] = tuple(
            (state, face_availability[state])
            for state in sorted(face_availability, key=lambda item: (len(item), item))
        )

    require(
        tuple(len(availability[name]) for name in FACE_ORDER)
        == (238, 4318, 3838),
        "facewise decision census drift",
    )

    pair_names = (("F12", "F13"), ("F12", "F23"), ("F13", "F23"))
    pair_records = {}
    pair_empty_loci = {}
    for left, right in pair_names:
        joint = tuple(sorted(
            set(availability[left]) & set(availability[right]),
            key=lambda item: (len(item), item),
        ))
        records = tuple(
            (
                state,
                availability[left][state],
                availability[right][state],
                tuple(sorted(
                    set(availability[left][state])
                    & set(availability[right][state])
                )),
            )
            for state in joint
        )
        pair_records[left, right] = records
        pair_empty_loci[left, right] = tuple(
            record[0] for record in records if not record[3]
        )

    require(
        tuple(len(pair_records[pair]) for pair in pair_names)
        == (238, 94, 1726),
        "pairwise joint-nonreset universe drift",
    )
    require(
        pair_empty_loci["F12", "F13"]
        == ((3, 4, 5), (1, 3, 4, 5))
        and pair_empty_loci["F12", "F23"] == ()
        and pair_empty_loci["F13", "F23"] == (),
        "pairwise empty-intersection locus drift",
    )

    triple_states = tuple(sorted(
        set(availability["F12"])
        & set(availability["F13"])
        & set(availability["F23"]),
        key=lambda item: (len(item), item),
    ))
    triple_records = []
    proper_helly = []
    for state in triple_states:
        face_sets = tuple(
            frozenset(availability[name][state]) for name in FACE_ORDER
        )
        pair_intersections = (
            tuple(sorted(face_sets[0] & face_sets[1])),
            tuple(sorted(face_sets[0] & face_sets[2])),
            tuple(sorted(face_sets[1] & face_sets[2])),
        )
        triple_intersection = tuple(sorted(intersection_all(face_sets)))
        triple_records.append((
            state,
            tuple(tuple(sorted(values)) for values in face_sets),
            pair_intersections,
            triple_intersection,
        ))
        if not triple_intersection and all(pair_intersections):
            proper_helly.append(state)
    triple_records = tuple(triple_records)
    triple_empty_locus = tuple(
        record[0] for record in triple_records if not record[3]
    )
    proper_helly = tuple(proper_helly)

    require(len(triple_states) == 94, "triple joint universe drift")
    require(
        triple_empty_locus == ((5,), (3, 4, 5), (1, 3, 4, 5))
        and proper_helly == ((5,),),
        "triple/proper-Helly locus drift",
    )

    expected_triple_rows = {
        (5,): (
            (2, 8, 9, 11, 12, 14, 16, 18, 22),
            (3, 4, 5, 6, 7, 10, 12, 13, 17, 19, 20, 21),
            (3, 4, 5, 6, 7, 8, 10, 13, 16, 17, 18, 19, 20, 21, 22),
        ),
        (3, 4, 5): (
            (2, 5, 8, 9, 11, 12, 14, 16, 18, 22),
            (3, 4, 6, 7, 10, 13, 17, 19, 20, 21),
            (3, 4, 5, 6, 7, 10, 12, 13, 17, 19, 21),
        ),
        (1, 3, 4, 5): (
            (2, 5, 9, 11, 12, 14, 16, 18, 22),
            (3, 4, 6, 7, 10, 13, 17, 19, 20, 21),
            tuple(range(1, 23)),
        ),
    }
    triple_by_state = {record[0]: record for record in triple_records}
    require(
        all(
            triple_by_state[state][1] == expected
            for state, expected in expected_triple_rows.items()
        ),
        "one of the three empty-triple row sets changed",
    )
    require(
        triple_by_state[(5,)][2]
        == ((12,), (8, 16, 18, 22),
            (3, 4, 5, 6, 7, 10, 13, 17, 19, 20, 21)),
        "proper Helly pair intersections changed",
    )

    # Every raw vector at which at least two faces have a decision obligation.
    all_vectors = tuple(sorted(
        set().union(*(set(availability[name]) for name in FACE_ORDER)),
        key=lambda item: (len(item), item),
    ))
    overlap_vectors = []
    dependency_locus = []
    local_branch_histogram = Counter()
    for state in all_vectors:
        present = tuple(name for name in FACE_ORDER if state in availability[name])
        if len(present) < 2:
            continue
        overlap_vectors.append(state)
        common = intersection_all(
            availability[name][state] for name in present
        )
        if common:
            minimum_local_labels = 1
        else:
            dependency_locus.append(state)
            minimum_local_labels = None
            for number in range(2, len(present) + 1):
                for colouring in product(range(number), repeat=len(present)):
                    if set(colouring) != set(range(number)):
                        continue
                    valid = True
                    for colour in range(number):
                        group = tuple(
                            name for name, assigned in zip(present, colouring)
                            if assigned == colour
                        )
                        if not intersection_all(
                            availability[name][state] for name in group
                        ):
                            valid = False
                            break
                    if valid:
                        minimum_local_labels = number
                        break
                if minimum_local_labels is not None:
                    break
            require(
                minimum_local_labels is not None,
                ("local face colouring failed", state),
            )
        local_branch_histogram[len(present), minimum_local_labels] += 1

    dependency_locus = tuple(dependency_locus)
    require(
        tuple(sorted(local_branch_histogram.items()))
        == (((2, 1), 1776), ((3, 1), 91), ((3, 2), 3))
        and dependency_locus == triple_empty_locus,
        "local branching or necessary dependency locus drift",
    )

    def colouring_valid(colouring, allowed_rows=None):
        for state in all_vectors:
            for colour in set(colouring):
                group = tuple(
                    name
                    for name, assigned in zip(FACE_ORDER, colouring)
                    if assigned == colour and state in availability[name]
                )
                if not group:
                    continue
                common = intersection_all(
                    availability[name][state] for name in group
                )
                if allowed_rows is not None:
                    common &= frozenset(allowed_rows)
                if not common:
                    return False
        return True

    one_label_colourings = tuple(
        colouring for colouring in ((0, 0, 0),)
        if colouring_valid(colouring)
    )
    binary_colourings = tuple(
        colouring
        for colouring in product((0, 1), repeat=3)
        if colouring_valid(colouring)
    )
    require(
        one_label_colourings == ()
        and binary_colourings
        == ((0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1)),
        "global static origin-colouring classification changed",
    )

    canonical_pair_extensions = tuple(
        colouring for colouring in binary_colourings
        if colouring_valid(colouring, allowed_rows=(16, 17))
    )

    # THM-3278 sees only F12/F13.  Its old conflict set is reproduced, but
    # the new proper defect (5,) is not one of those conflicts: the first two
    # faces share row 12 there, and only F23 deletes the simultaneous choice.
    require(
        pair_empty_loci["F12", "F13"]
        == ((3, 4, 5), (1, 3, 4, 5))
        and triple_by_state[(5,)][2][0] == (12,)
        and not triple_by_state[(5,)][3]
        and b"conflict_states=((3,4,5),(1,3,4,5))"
        in THM3278_OUTPUT.read_bytes().replace(b" ", b""),
        "THM-3278 non-subsumption witness changed",
    )

    face_digests = tuple(
        (name, stable_digest(availability_records[name]))
        for name in FACE_ORDER
    )
    pair_digests = tuple(
        (pair, stable_digest(pair_records[pair])) for pair in pair_names
    )
    triple_digest = stable_digest(triple_records)
    dependency_digest = stable_digest(dependency_locus)

    print("THM-3286 THREE-FACE AVAILABILITY INDEPENDENT HOSTILE AUDIT")
    print("status=FINITE-EXACT independent audit;candidate not promoted here")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path.relative_to(ROOT)}")
    print("dependency_hashes=END")
    print(
        "implementation=direct_THM3238_product_Gamma_coefficients;"
        "independent_partitions_and_upsets;independent_physical_submultisets;"
        "independent_reset_neighbours;no_candidate_import_or_availability_table"
    )
    print(
        "face_universes="
        + repr(tuple(
            (
                name,
                faces[name]["pole_profile"],
                faces[name]["reset"],
                len(faces[name]["states"]),
                len(availability[name]),
            )
            for name in FACE_ORDER
        ))
    )
    print(f"response_vectors_rebuilt={len(responses)}/8397;all_rows=22")
    print("face_availability_digests=" + repr(face_digests))
    print(
        "pair_joint_universes="
        + repr(tuple((pair, len(pair_records[pair])) for pair in pair_names))
    )
    print("pair_empty_loci=" + repr(tuple(
        (pair, pair_empty_loci[pair]) for pair in pair_names
    )))
    print("pair_record_digests=" + repr(pair_digests))
    print(
        f"triple_joint_universe={len(triple_states)};"
        f"triple_empty_locus={triple_empty_locus};"
        f"proper_Helly_locus={proper_helly}"
    )
    print(
        "proper_Helly_(5)_row_sets="
        + repr(triple_by_state[(5,)][1])
        + ";pair_intersections="
        + repr(triple_by_state[(5,)][2])
    )
    print("other_empty_triple_row_sets=" + repr(tuple(
        (state, triple_by_state[state][1])
        for state in triple_empty_locus if state != (5,)
    )))
    print(
        f"local_branch_histogram={tuple(sorted(local_branch_histogram.items()))};"
        "three_branch_required=0"
    )
    print(
        f"minimum_static_origin_labels=2;binary_width=1;"
        f"optimal_binary_assignments={binary_colourings};"
        "F23_may_share_either_F12_or_F13=yes"
    )
    print(
        f"necessary_origin_dependency_locus={dependency_locus};"
        f"dependency_digest={dependency_digest}"
    )
    print(
        "THM3278_non_subsumption=new_proper_state_(5)_has_"
        "F12_intersect_F13=(12)_but_F12_intersect_F13_intersect_F23=empty;"
        "old_two_face_conflicts_reproduced="
        + repr(pair_empty_loci["F12", "F13"])
    )
    print(
        "THM3278_canonical_pair_(16,17)_global_three_face_extensions="
        + repr(canonical_pair_extensions)
    )
    print(
        f"audit_digests=(triple:{triple_digest},dependency:{dependency_digest})"
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})"
    )
    print(
        "scope=three_named_bank_I2_response_faces_and_22_lawful_rows_only;"
        "no_FC3_SFC3_GMC_positivity_or_sequential_physical_controller"
    )
    print("VERDICT=ACCEPT_CANDIDATE_MATHEMATICS_IF_TRANSCRIPTS_MATCH")
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    freeze_support()
    main()
