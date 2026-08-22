#!/usr/bin/env python3
"""Exact three-face response-availability hypergraph scout for FC(3)."""

import ast
import hashlib
import importlib.util
from collections import Counter
from itertools import combinations_with_replacement, product
from multiprocessing import Pool, freeze_support
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
THM3249_SCRIPT = ROOT / (
    "04-computation/gmc_cross_support_upset_atlas_holonomy_thm3249.py"
)
THM3249_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_cross_support_upset_atlas_holonomy_thm3249.out"
)
THM3275_SCRIPT = ROOT / (
    "04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py"
)
THM3275_OUTPUT = ROOT / (
    "05-knowledge/results/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.out"
)
DEPENDENCIES = {
    THM3249_SCRIPT:
        "bb901e92687544c69d67a55d057f8293ecbf516b80d491a636c4a62af19eebef",
    THM3249_OUTPUT:
        "76d037d4f0737b37ad48531e23be1a0a37509ce0a3d029d3cefdd42662ec04f2",
    THM3275_SCRIPT:
        "2ad9eacf8e893f881b5616672d5fde872a80612b780372dd2318d57f75b6ea30",
    THM3275_OUTPUT:
        "664b7677e89873dda26cd004878c134e621086095d4155d307b52213471baf63",
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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cross = load_module("thm3249_three_face_source", THM3249_SCRIPT)

FACE_PARAMETERS = (
    ("F12", 1, 2),
    ("F13", 1, 3),
    ("F23", 2, 3),
)


def physical_states_product(poles):
    multiplicities = Counter(poles)
    roots = tuple(sorted(multiplicities))
    return tuple(sorted((
        tuple(root for root, count in zip(roots, choices)
              for _ in range(count))
        for choices in product(
            *(range(multiplicities[root] + 1) for root in roots)
        )
        if any(choices)
    ), key=lambda state: (len(state), state)))


def physical_states_filtered(poles):
    multiplicities = Counter(poles)
    roots = tuple(sorted(multiplicities))
    return tuple(
        state
        for depth in range(1, len(poles) + 1)
        for state in combinations_with_replacement(roots, depth)
        if Counter(state) <= multiplicities
    )


def row_coordinates(vector):
    return tuple(
        sum(vector[entry[0]][shape] for shape in upset)
        for entry, upset in zip(cross.T.CERTIFICATE, cross.T.UPSETS)
    )


def compute_task(task):
    face, a_value, b_value, state = task
    vector = cross.T.coefficient_vectors(
        1, cross.BANK, a_value, b_value, state
    )
    return face, state, row_coordinates(vector)


def toward_reset_neighbours(state, poles, reset):
    current = Counter(state)
    capacities = Counter(poles)
    target = Counter(reset)
    answer = []
    for root in sorted(capacities):
        changed = current.copy()
        if current[root] > target[root]:
            changed[root] -= 1
            if not changed[root]:
                del changed[root]
        elif current[root] < target[root]:
            changed[root] += 1
        else:
            continue
        candidate = tuple(sorted(changed.elements()))
        require(Counter(candidate) <= capacities,
                ("directed neighbour left physical bank", state, candidate))
        if candidate:
            answer.append(candidate)
    return tuple(answer)


def toward_reset_neighbours_counts(state, poles, reset):
    maximum = max(poles)
    counts = tuple(Counter(state)[root] for root in range(1, maximum + 1))
    capacities = tuple(Counter(poles)[root] for root in range(1, maximum + 1))
    target = tuple(Counter(reset)[root] for root in range(1, maximum + 1))
    answer = []
    for index, (value, target_value) in enumerate(zip(counts, target)):
        if value == target_value:
            continue
        changed = list(counts)
        changed[index] += 1 if value < target_value else -1
        require(0 <= changed[index] <= capacities[index],
                ("count neighbour left capacity", state, index))
        candidate = tuple(
            root
            for root, multiplicity in enumerate(changed, 1)
            for _ in range(multiplicity)
        )
        if candidate:
            answer.append(candidate)
    return tuple(answer)


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def common_rows(availability, faces, state):
    banks = [frozenset(availability[face][state]) for face in faces]
    return tuple(sorted(set.intersection(*map(set, banks))))


def minimum_local_labels(availability, faces, state):
    for label_count in range(1, len(faces) + 1):
        for labels in product(range(label_count), repeat=len(faces)):
            if set(labels) != set(range(label_count)):
                continue
            lawful = True
            for label in range(label_count):
                block = tuple(
                    face for face, face_label in zip(faces, labels)
                    if face_label == label
                )
                if not common_rows(availability, block, state):
                    lawful = False
                    break
            if lawful:
                return label_count
    raise RuntimeError(("no local origin labelling", faces, state))


def main():
    face_data = {}
    tasks = []
    for face, a_value, b_value in FACE_PARAMETERS:
        poles = tuple(cross.T.reduced_poles(
            1, cross.BANK, a_value, b_value
        )[0])
        reset = tuple(sorted(cross.T.residual_roots(
            1, cross.T.dominant_row(1), a_value, b_value
        )))
        states_product = physical_states_product(poles)
        states_filtered = physical_states_filtered(poles)
        require(states_product == states_filtered,
                ("independent physical-state constructors disagree", face))
        require(Counter(reset) <= Counter(poles) and reset in states_product,
                ("reset left physical bank", face))
        face_data[face] = {
            "parameters": (a_value, b_value),
            "poles": poles,
            "reset": reset,
            "states": states_product,
            "decisions": tuple(state for state in states_product if state != reset),
        }
        tasks.extend(
            (face, a_value, b_value, state) for state in states_product
        )

    with Pool(processes=4) as pool:
        computed = tuple(pool.imap_unordered(compute_task, tasks, chunksize=2))
    require(len(computed) == len(tasks), "response reconstruction lost tasks")
    rows = {face: {} for face, _, _ in FACE_PARAMETERS}
    for face, state, row in computed:
        require(state not in rows[face], ("duplicate response task", face, state))
        rows[face][state] = row

    availability = {face: {} for face, _, _ in FACE_PARAMETERS}
    cover_sets = {face: [set() for _ in range(22)]
                  for face, _, _ in FACE_PARAMETERS}
    for face, _, _ in FACE_PARAMETERS:
        data = face_data[face]
        require(len(rows[face]) == len(data["states"]),
                ("response bank incomplete", face))
        for state in data["states"]:
            require(
                toward_reset_neighbours(state, data["poles"], data["reset"])
                == toward_reset_neighbours_counts(
                    state, data["poles"], data["reset"]
                ),
                ("independent directed-neighbour constructors disagree", face, state),
            )
        for state in data["decisions"]:
            neighbours = toward_reset_neighbours(
                state, data["poles"], data["reset"]
            )
            available = tuple(
                chart + 1
                for chart in range(22)
                if any(rows[face][target][chart] > rows[face][state][chart]
                       for target in neighbours)
            )
            require(available, ("all-row coverage failed", face, state))
            availability[face][state] = available
            for chart in available:
                cover_sets[face][chart - 1].add(state)
        require(frozenset().union(*map(frozenset, cover_sets[face]))
                == frozenset(data["decisions"]),
                ("transposed cover reconstruction failed", face))

    # Positive control against the already promoted two-face theorem: direct
    # product-Gamma reconstruction must reproduce THM-3275's full atlas digest.
    joint_12_13 = tuple(sorted(
        set(face_data["F12"]["decisions"])
        & set(face_data["F13"]["decisions"]),
        key=lambda state: (len(state), state),
    ))
    canonical_two_face_records = tuple(
        (
            state,
            availability["F12"][state],
            availability["F13"][state],
            common_rows(availability, ("F12", "F13"), state),
        )
        for state in joint_12_13
    )
    require(digest(canonical_two_face_records)
            == "9ee92c1bc68f02cb023031ec3f0cf6eb154408f49ff3e67f39db54e9c706bcd8",
            "direct reconstruction failed THM-3275 atlas control")

    pair_records = {}
    pair_empty = {}
    for left, right in (("F12", "F13"), ("F12", "F23"), ("F13", "F23")):
        joint = tuple(sorted(
            set(face_data[left]["decisions"])
            & set(face_data[right]["decisions"]),
            key=lambda state: (len(state), state),
        ))
        records = tuple(
            (
                state,
                availability[left][state],
                availability[right][state],
                common_rows(availability, (left, right), state),
            )
            for state in joint
        )
        pair_records[(left, right)] = records
        pair_empty[(left, right)] = tuple(
            record[0] for record in records if not record[3]
        )

    triple_joint = tuple(sorted(
        set(face_data["F12"]["decisions"])
        & set(face_data["F13"]["decisions"])
        & set(face_data["F23"]["decisions"]),
        key=lambda state: (len(state), state),
    ))
    triple_records = tuple(
        (
            state,
            availability["F12"][state],
            availability["F13"][state],
            availability["F23"][state],
            common_rows(availability, ("F12", "F13"), state),
            common_rows(availability, ("F12", "F23"), state),
            common_rows(availability, ("F13", "F23"), state),
            common_rows(availability, ("F12", "F13", "F23"), state),
        )
        for state in triple_joint
    )
    triple_empty_records = tuple(record for record in triple_records if not record[7])
    proper_helly_records = tuple(
        record for record in triple_empty_records
        if record[4] and record[5] and record[6]
    )

    decision_face_sets = {}
    for face, _, _ in FACE_PARAMETERS:
        for state in face_data[face]["decisions"]:
            decision_face_sets.setdefault(state, []).append(face)
    shared_face_sets = tuple(
        (state, tuple(faces))
        for state, faces in sorted(
            decision_face_sets.items(), key=lambda item: (len(item[0]), item[0])
        )
        if len(faces) >= 2
    )
    local_branch_histogram = tuple(sorted(Counter(
        (len(faces), minimum_local_labels(availability, faces, state))
        for state, faces in shared_face_sets
    ).items()))

    dependency_locus = tuple(
        state for state, faces in shared_face_sets
        if not common_rows(availability, faces, state)
    )
    static_label_assignments = []
    for labels in product(range(2), repeat=3):
        if set(labels) != {0, 1}:
            continue
        face_labels = dict(zip(("F12", "F13", "F23"), labels))
        lawful = True
        for state, faces in shared_face_sets:
            for label in (0, 1):
                block = tuple(
                    face for face in faces if face_labels[face] == label
                )
                if block and not common_rows(availability, block, state):
                    lawful = False
                    break
            if not lawful:
                break
        if lawful:
            static_label_assignments.append(labels)
    static_label_assignments = tuple(static_label_assignments)

    # Every successful binary assignment yields a literal state+label lookup
    # on all three full decision banks.  This checks sufficiency, not just the
    # three hostile fibres.
    policy_case_count = 0
    for labels in static_label_assignments:
        face_labels = dict(zip(("F12", "F13", "F23"), labels))
        for state, faces in decision_face_sets.items():
            for label in (0, 1):
                block = tuple(
                    face for face in faces if face_labels[face] == label
                )
                if not block:
                    continue
                choices = common_rows(availability, block, state)
                require(choices, ("binary selector block failed", labels, state, block))
                choice = min(choices)
                for face in block:
                    require(choice in availability[face][state],
                            ("binary selector chose unavailable row", face, state))
                    policy_case_count += 1

    face_records = {
        face: tuple(
            (state, availability[face][state])
            for state in face_data[face]["decisions"]
        )
        for face, _, _ in FACE_PARAMETERS
    }
    face_digests = tuple(
        (face, digest(face_records[face])) for face, _, _ in FACE_PARAMETERS
    )
    triple_digest = digest(triple_records)
    dependency_digest = digest(dependency_locus)

    face_summary = tuple(
        (
            face,
            data["parameters"],
            tuple(sorted(Counter(data["poles"]).items())),
            data["reset"],
            len(data["states"]),
            len(data["decisions"]),
        )
        for face, _, _ in FACE_PARAMETERS
        for data in (face_data[face],)
    )
    require(face_summary == (
        ("F12", (1, 2), ((1, 4), (2, 3), (3, 2), (4, 1), (5, 1)),
         (1, 2, 2, 3, 4, 5), 239, 238),
        ("F13", (1, 3), ((1, 4), (2, 3), (3, 2), (4, 2), (5, 2),
                           (6, 1), (7, 1), (8, 1)),
         (1, 3, 3, 4, 5, 6, 7, 8), 4319, 4318),
        ("F23", (2, 3), ((1, 2), (2, 4), (3, 1), (4, 3), (5, 3),
                           (6, 1), (7, 1), (8, 1)),
         (2, 3, 4, 4, 5, 6, 7, 8), 3839, 3838),
    ), "three-face physical profile drift")
    require(tuple(
        (pair, len(pair_records[pair]), pair_empty[pair])
        for pair in (("F12", "F13"), ("F12", "F23"), ("F13", "F23"))
    ) == (
        (("F12", "F13"), 238, ((3, 4, 5), (1, 3, 4, 5))),
        (("F12", "F23"), 94, ()),
        (("F13", "F23"), 1726, ()),
    ), "pairwise Helly loci drift")
    require(len(triple_joint) == 94
            and tuple(record[0] for record in triple_empty_records)
            == ((5,), (3, 4, 5), (1, 3, 4, 5)),
            "triple obstruction locus drift")
    require(tuple(record[0] for record in proper_helly_records) == ((5,),),
            "proper three-way Helly witness drift")
    require(local_branch_histogram == (
        ((2, 1), 1776), ((3, 1), 91), ((3, 2), 3),
    ), "local branch-width histogram drift")
    require(static_label_assignments == (
        (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1),
    ), "global binary-origin assignment drift")
    require(dependency_locus == ((5,), (3, 4, 5), (1, 3, 4, 5)),
            "necessary origin-bit dependency locus drift")
    require(face_digests == (
        ("F12", "31d22b0d203f365b7dfc586a4acab31e591380fd24755dc32d7dc7bd946218c3"),
        ("F13", "8bf051c08c9f6f4ce04ce715faf5c92738c51717dc7c59091f9dfadc04bcc739"),
        ("F23", "3ee999fe0e397b8938daa5ca7b49b08be27fab2a07f08366b274a12a43790fe9"),
    ), "face availability digest drift")
    require(triple_digest
            == "d0ed0a5ace30b771637daaf388bdd794d787de5d2df566210fc2a413bc78ce2c",
            "triple availability atlas digest drift")
    require(dependency_digest
            == "d0427276be9fdf1018ee3cb29035aa596457f098b4358b156d640f3c6521edc4",
            "dependency-locus digest drift")
    require(policy_case_count == 4 * sum(
        len(face_data[face]["decisions"]) for face, _, _ in FACE_PARAMETERS
    ), "binary policy case census drift")

    print("FC3 three-support-face availability hypergraph exact scout")
    print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
    print("faces=" + repr(face_summary))
    print("physical_state_constructors=product_vs_filtered:PASS")
    print("directed_neighbour_constructors=counter_vs_padded_counts:PASS")
    print("all_nonreset_states_have_at_least_one_of_22_rows=PASS")
    print("face_availability_digests=" + repr(face_digests))
    for pair in (("F12", "F13"), ("F12", "F23"), ("F13", "F23")):
        print("pair=%s:joint_nonreset=%d:empty=%s" % (
            pair, len(pair_records[pair]), repr(pair_empty[pair])
        ))
    print("THM3275_two_face_atlas_control=PASS")
    print("triple_joint_nonreset=" + repr(len(triple_joint)))
    print("triple_empty_count=" + repr(len(triple_empty_records)))
    for record in triple_empty_records:
        print("triple_empty_record=" + repr(record))
    print("proper_three_way_helly_defect_count=" + repr(len(proper_helly_records)))
    for record in proper_helly_records:
        print("proper_three_way_helly_defect=" + repr(record))
    print("shared_raw_state_count=" + repr(len(shared_face_sets)))
    print("local_(active_faces,minimum_origin_labels)_histogram="
          + repr(local_branch_histogram))
    print("maximum_local_origin_labels="
          + repr(max(key[1] for key, _ in local_branch_histogram)))
    print("minimum_fixed_origin_labels=2")
    print("minimum_fixed_origin_bits=1")
    print("optimal_binary_face_assignments_F12_F13_F23="
          + repr(static_label_assignments))
    print("binary_policy_face_decision_checks=" + repr(policy_case_count))
    print("necessary_bit_dependency_locus=" + repr(dependency_locus))
    print("necessary_bit_dependency_locus_digest=" + dependency_digest)
    print("triple_atlas_digest=" + triple_digest)
    print("mechanism=proper_Cech_Helly_defect_at_(5);F23_shares_either_binary_origin_class")
    print("typed_scope=three_named_bank_I2_response_faces_and_22_lawful_rows;local_strict_reset_directed_availability_only;no_FC3_SFC3_GMC_or_positivity_conclusion")
    print("source_ast=(assert_nodes=0,float_literals=0)")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    freeze_support()
    main()
