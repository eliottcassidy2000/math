#!/usr/bin/env python3
"""Finite-exact partial scout for move-decorated THM-3286 availability.

For the three named bank-I2 faces and exactly the inherited 94-state triple
decision universe, replace each available row by every pair

    (row, actual reset-directed one-pole successor)

that realizes its strict response ascent.  Intersections now require both
the same row and the same physical successor.  The script recomputes the
decorated pair/triple defects, local and fixed supplied-origin widths, static
sections, and finite-graph termination.  Projection to the row coordinate is
checked against THM-3286's frozen triple atlas.

This is FINITE-EXACT PARTIAL evidence only.  It proves no FC(3), SFC(3), GMC,
positivity statement, or controller outside the inherited finite graph.
"""

from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from multiprocessing import Pool, freeze_support
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
THM3286_THEOREM = (
    ROOT
    / "01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md"
)
BASE_SCRIPT = COMP / "gmc_complete_physical_bank_unique_reset_thm3238.py"
BASE_OUTPUT = RESULTS / "gmc_complete_physical_bank_unique_reset_thm3238.out"
THM3286_SCRIPT = (
    COMP / "fc3_three_support_face_availability_hypergraph_scout_20260803.py"
)
THM3286_OUTPUT = (
    RESULTS / "fc3_three_support_face_availability_hypergraph_scout_20260803.out"
)
THM3278_SCRIPT = COMP / "fc3_selector_origin_bipartition_phase_bridge_thm3278.py"
THM3278_OUTPUT = RESULTS / "fc3_selector_origin_bipartition_phase_bridge_thm3278.out"
PINNED = {
    THM3286_THEOREM:
        "fd48e7e142d44b5c60baec7e7bc1b59e2640f672f8de019d3ed2a2e7d969163c",
    BASE_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    BASE_OUTPUT:
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    THM3286_SCRIPT:
        "c31eff8a10c6a6e1ab7e3ea6759388b02cbd9591695182f2d3756e08da38c8c3",
    THM3286_OUTPUT:
        "a3e44f6d2eb0e26386e399e7d582e6eb67512125233eb1b49d2b62c4d0869e05",
    THM3278_SCRIPT:
        "07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f",
    THM3278_OUTPUT:
        "5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a",
}

FACE_PARAMETERS = (
    ("F12", 1, 2),
    ("F13", 1, 3),
    ("F23", 2, 3),
)
FACE_ORDER = tuple(face for face, _a, _b in FACE_PARAMETERS)
PAIR_ORDER = (("F12", "F13"), ("F12", "F23"), ("F13", "F23"))
EXPECTED_HOSTILE_ROWS = {
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


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


def stable_digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(path)) for path in PINNED)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "a pinned inherited artifact changed",
)
require(
    "none is globally lawful on the three faces"
    in " ".join(THM3286_THEOREM.read_text(encoding="utf-8").split()),
    "pinned THM-3286 full-domain (16,17) failure wording changed",
)

SOURCE = Path(__file__).read_text(encoding="utf-8")
SOURCE_TREE = ast.parse(SOURCE)
ASSERT_NODES = sum(isinstance(node, ast.Assert) for node in ast.walk(SOURCE_TREE))
FLOAT_LITERALS = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(SOURCE_TREE)
)
require(
    ASSERT_NODES == 0 and FLOAT_LITERALS == 0,
    "source validity gate found assert or floating-point literal",
)


BASE_SPEC = spec_from_file_location("thm3238_for_move_decorated_partial", BASE_SCRIPT)
base = module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)


def row_coordinates(vectors):
    return tuple(
        sum(vectors[entry[0]][shape] for shape in upset)
        for entry, upset in zip(base.CERTIFICATE, base.UPSETS)
    )


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


def reset_distance(state, reset):
    return sum(
        abs(left - right)
        for left, right in zip(padded_counts(state), padded_counts(reset))
    )


def toward_reset_neighbours(state, poles, reset):
    """All nonempty physical one-pole moves lowering reset distance by one."""

    current = Counter(state)
    capacities = Counter(poles)
    target = Counter(reset)
    answer = []
    for root in sorted(capacities):
        changed = current.copy()
        if current[root] > target[root]:
            changed[root] -= 1
            if changed[root] == 0:
                del changed[root]
        elif current[root] < target[root]:
            changed[root] += 1
        else:
            continue
        successor = tuple(sorted(changed.elements()))
        if not successor:
            continue
        require(
            Counter(successor) <= capacities,
            ("reset-directed move left nonempty physical bank", state, successor),
        )
        require(
            reset_distance(successor, reset) + 1 == reset_distance(state, reset),
            ("reset-directed move failed distance test", state, successor),
        )
        answer.append(successor)
    return tuple(answer)


FACES = {}
for face, a_value, b_value in FACE_PARAMETERS:
    poles = tuple(base.reduced_poles(1, base.BANK, a_value, b_value)[0])
    reset = tuple(sorted(base.residual_roots(
        1, base.dominant_row(1), a_value, b_value
    )))
    states = physical_states(poles)
    require(reset in states and Counter(reset) <= Counter(poles), f"{face} reset")
    FACES[face] = {
        "a": a_value,
        "b": b_value,
        "poles": poles,
        "reset": reset,
        "states": states,
        "decisions": frozenset(state for state in states if state != reset),
    }


def response_task(task):
    face, a_value, b_value, state = task
    vectors = base.coefficient_vectors(1, base.BANK, a_value, b_value, state)
    return (face, state), row_coordinates(vectors)


def intersection_all(values):
    values = tuple(values)
    require(values, "intersection_all received no sets")
    return frozenset.intersection(*(frozenset(value) for value in values))


def common_decorated(decorated, faces, state, allowed_rows=None):
    common = intersection_all(decorated[face][state] for face in faces)
    if allowed_rows is not None:
        common = frozenset(
            pair for pair in common if pair[0] in frozenset(allowed_rows)
        )
    return tuple(sorted(common))


def common_rows(projected, faces, state):
    return tuple(sorted(intersection_all(projected[face][state] for face in faces)))


def common_successors(successors, faces, state):
    return tuple(sorted(intersection_all(successors[face][state] for face in faces)))


def surjective_colourings(number):
    return tuple(
        colouring
        for colouring in product(range(number), repeat=len(FACE_ORDER))
        if set(colouring) == set(range(number))
    )


def canonical_colouring(colouring):
    relabel = {}
    next_label = 0
    answer = []
    for label in colouring:
        if label not in relabel:
            relabel[label] = next_label
            next_label += 1
        answer.append(relabel[label])
    return tuple(answer)


def colouring_valid(decorated, triple_states, colouring, allowed_rows=None):
    face_labels = dict(zip(FACE_ORDER, colouring))
    for state in triple_states:
        for label in range(max(colouring) + 1):
            block = tuple(face for face in FACE_ORDER if face_labels[face] == label)
            if not common_decorated(
                decorated, block, state, allowed_rows=allowed_rows
            ):
                return False
    return True


def minimum_local_labels(decorated, state):
    for number in range(1, len(FACE_ORDER) + 1):
        for colouring in surjective_colourings(number):
            face_labels = dict(zip(FACE_ORDER, colouring))
            if all(
                common_decorated(
                    decorated,
                    tuple(
                        face for face in FACE_ORDER if face_labels[face] == label
                    ),
                    state,
                )
                for label in range(number)
            ):
                return number
    raise RuntimeError(("no local decorated colouring", state))


def section_iteration_record(decorated, triple_states, colouring):
    """Audit every allowed section edge and one lexicographic section."""

    triple_set = frozenset(triple_states)
    face_labels = dict(zip(FACE_ORDER, colouring))
    block_records = []
    for label in range(max(colouring) + 1):
        block = tuple(face for face in FACE_ORDER if face_labels[face] == label)
        adjacency = {}
        exit_choices = {}
        lexicographic = {}
        decorated_edge_count = 0
        for state in triple_states:
            choices = common_decorated(decorated, block, state)
            require(choices, ("section block unexpectedly empty", colouring, state))
            decorated_edge_count += len(choices)
            internal = tuple(sorted({target for _row, target in choices if target in triple_set}))
            exits = tuple(pair for pair in choices if pair[1] not in triple_set)
            adjacency[state] = internal
            exit_choices[state] = exits
            lexicographic[state] = min(choices)
            for _row, target in choices:
                require(
                    all(
                        reset_distance(target, FACES[face]["reset"]) + 1
                        == reset_distance(state, FACES[face]["reset"])
                        for face in block
                    ),
                    ("common edge is not reset-directed on its block", block, state),
                )

        visiting = set()
        longest_cache = {}

        def longest_internal(state):
            if state in longest_cache:
                return longest_cache[state]
            require(state not in visiting, ("decorated internal cycle", block, state))
            visiting.add(state)
            value = 0 if not adjacency[state] else 1 + max(
                longest_internal(target) for target in adjacency[state]
            )
            visiting.remove(state)
            longest_cache[state] = value
            return value

        longest_all_choice = max(longest_internal(state) for state in triple_states)
        lex_lengths = []
        for start in triple_states:
            state = start
            seen = set()
            moves = 0
            while True:
                require(state not in seen, ("lexicographic section cycle", block, start))
                seen.add(state)
                _row, target = lexicographic[state]
                moves += 1
                if target not in triple_set:
                    break
                state = target
            lex_lengths.append(moves)

        record = (
            label,
            block,
            decorated_edge_count,
            sum(len(targets) for targets in adjacency.values()),
            sum(len(choices) for choices in exit_choices.values()),
            sum(not adjacency[state] for state in triple_states),
            longest_all_choice,
            max(lex_lengths),
            stable_digest(tuple(
                (state, lexicographic[state]) for state in triple_states
            )),
        )
        block_records.append(record)
    return tuple(block_records)


def main():
    face_profile = tuple(
        (
            face,
            tuple(sorted(Counter(FACES[face]["poles"]).items())),
            FACES[face]["reset"],
            len(FACES[face]["states"]),
            len(FACES[face]["decisions"]),
        )
        for face in FACE_ORDER
    )
    require(
        tuple((row[0], row[3], row[4]) for row in face_profile)
        == (("F12", 239, 238), ("F13", 4319, 4318), ("F23", 3839, 3838)),
        "inherited face universe drift",
    )
    triple_states = tuple(sorted(
        set.intersection(*(set(FACES[face]["decisions"]) for face in FACE_ORDER)),
        key=lambda state: (len(state), state),
    ))
    require(len(triple_states) == 94, "inherited triple universe is not 94")

    successors = {face: {} for face in FACE_ORDER}
    needed_states = {face: set(triple_states) for face in FACE_ORDER}
    for face in FACE_ORDER:
        for state in triple_states:
            targets = toward_reset_neighbours(
                state, FACES[face]["poles"], FACES[face]["reset"]
            )
            successors[face][state] = targets
            needed_states[face].update(targets)
    tasks = tuple(
        (face, FACES[face]["a"], FACES[face]["b"], state)
        for face in FACE_ORDER
        for state in sorted(needed_states[face], key=lambda item: (len(item), item))
    )
    with Pool(processes=4) as pool:
        responses = dict(pool.imap_unordered(response_task, tasks, chunksize=4))
    require(len(responses) == len(tasks), "a response task was lost")

    decorated = {face: {} for face in FACE_ORDER}
    projected = {face: {} for face in FACE_ORDER}
    for face in FACE_ORDER:
        for state in triple_states:
            source = responses[face, state]
            pairs = tuple(sorted(
                (row, target)
                for target in successors[face][state]
                for row in range(1, 23)
                if responses[face, target][row - 1] > source[row - 1]
            ))
            require(pairs, ("move-decorated face fibre empty", face, state))
            decorated[face][state] = pairs
            projected[face][state] = tuple(sorted({row for row, _target in pairs}))

    projected_triple_records = tuple(
        (
            state,
            projected["F12"][state],
            projected["F13"][state],
            projected["F23"][state],
            common_rows(projected, ("F12", "F13"), state),
            common_rows(projected, ("F12", "F23"), state),
            common_rows(projected, ("F13", "F23"), state),
            common_rows(projected, FACE_ORDER, state),
        )
        for state in triple_states
    )
    require(
        stable_digest(projected_triple_records)
        == "d0ed0a5ace30b771637daaf388bdd794d787de5d2df566210fc2a413bc78ce2c",
        "projection failed THM-3286 frozen triple atlas",
    )
    require(
        all(
            tuple(projected[face][state] for face in FACE_ORDER) == expected
            for state, expected in EXPECTED_HOSTILE_ROWS.items()
        ),
        "projected hostile control rows drift",
    )

    pair_records = {}
    pair_empty = {}
    pair_successor_empty = {}
    pair_projection_empty = {}
    for pair in PAIR_ORDER:
        records = tuple(
            (
                state,
                common_rows(projected, pair, state),
                common_successors(successors, pair, state),
                common_decorated(decorated, pair, state),
            )
            for state in triple_states
        )
        pair_records[pair] = records
        pair_empty[pair] = tuple(record[0] for record in records if not record[3])
        pair_projection_empty[pair] = tuple(
            record[0] for record in records if not record[1]
        )
        pair_successor_empty[pair] = tuple(
            record[0] for record in records if not record[2]
        )

    triple_records = tuple(
        (
            state,
            common_rows(projected, FACE_ORDER, state),
            common_successors(successors, FACE_ORDER, state),
            common_decorated(decorated, FACE_ORDER, state),
            tuple(common_decorated(decorated, pair, state) for pair in PAIR_ORDER),
        )
        for state in triple_states
    )
    triple_empty = tuple(record[0] for record in triple_records if not record[3])
    decorated_proper_helly = tuple(
        record[0] for record in triple_records
        if not record[3] and all(record[4])
    )

    local_widths = tuple(
        (state, minimum_local_labels(decorated, state)) for state in triple_states
    )
    local_histogram = tuple(sorted(Counter(width for _state, width in local_widths).items()))
    local_width_two_locus = tuple(
        state for state, width in local_widths if width == 2
    )

    valid_colourings = {}
    for number in range(1, len(FACE_ORDER) + 1):
        valid_colourings[number] = tuple(
            colouring for colouring in surjective_colourings(number)
            if colouring_valid(decorated, triple_states, colouring)
        )
    minimum_fixed = min(
        number for number, colourings in valid_colourings.items() if colourings
    )
    optimal_colourings = valid_colourings[minimum_fixed]
    canonical_optimal = tuple(sorted({
        canonical_colouring(colouring) for colouring in optimal_colourings
    }))
    uncoloured_section = not triple_empty
    coloured_section = bool(optimal_colourings)
    require(coloured_section, "singleton face colours failed to give a section")

    fixed_pair_extensions = tuple(
        (number, tuple(
            colouring
            for colouring in surjective_colourings(number)
            if colouring_valid(
                decorated, triple_states, colouring, allowed_rows=(16, 17)
            )
        ))
        for number in range(1, len(FACE_ORDER) + 1)
    )

    iteration_records = tuple(
        (colouring, section_iteration_record(decorated, triple_states, colouring))
        for colouring in canonical_optimal
    )

    failure_candidates = []
    incompatible_candidates = []
    for pair_index, pair in enumerate(PAIR_ORDER):
        for record in pair_records[pair]:
            state, rows, targets, pairs = record
            if rows and not pairs:
                mechanism = (
                    "incompatible_reset_successors" if not targets
                    else "row_successor_coupling"
                )
                failure_candidates.append(
                    ((len(state), state, pair_index), pair, state, mechanism, rows, targets)
                )
            if rows and not targets:
                incompatible_candidates.append(
                    ((len(state), state, pair_index), pair, state, rows)
                )
    require(failure_candidates, "decorating successors created no projection failure")
    require(incompatible_candidates, "no incompatible-reset hostile was found")
    first_failure = min(failure_candidates)[1:]
    first_pair, first_state, first_mechanism, first_rows, first_targets = first_failure
    first_record = next(
        record for record in pair_records[first_pair] if record[0] == first_state
    )
    first_row_realisers = tuple(
        (
            face,
            tuple(pair for pair in decorated[face][first_state] if pair[0] in first_rows),
        )
        for face in first_pair
    )
    first_incompatible = min(incompatible_candidates)[1:]

    state5_record = next(record for record in triple_records if record[0] == (5,))
    state5_successors = tuple(
        (face, successors[face][(5,)]) for face in FACE_ORDER
    )
    state5_decorated_counts = tuple(
        (face, len(decorated[face][(5,)])) for face in FACE_ORDER
    )
    state5_pair_summary = tuple(
        (
            pair,
            len(common_decorated(decorated, pair, (5,))),
            stable_digest(common_decorated(decorated, pair, (5,))),
            common_decorated(decorated, pair, (5,))
            if len(common_decorated(decorated, pair, (5,))) <= 6 else "DIGEST_ONLY",
        )
        for pair in PAIR_ORDER
    )
    inherited_conflict_controls = tuple(
        (
            state,
            common_rows(projected, ("F12", "F13"), state),
            common_decorated(decorated, ("F12", "F13"), state),
        )
        for state in ((3, 4, 5), (1, 3, 4, 5))
    )

    pair_summary = tuple(
        (
            pair,
            len(pair_empty[pair]),
            pair_empty[pair],
            pair_projection_empty[pair],
            pair_successor_empty[pair],
            stable_digest(pair_records[pair]),
        )
        for pair in PAIR_ORDER
    )
    decorated_face_digests = tuple(
        (
            face,
            stable_digest(tuple(
                (state, decorated[face][state]) for state in triple_states
            )),
        )
        for face in FACE_ORDER
    )
    triple_digest = stable_digest(triple_records)

    expected_pair_summary = (
        (
            ("F12", "F13"),
            13,
            (
                (4, 5), (1, 3, 4), (1, 4, 5), (2, 4, 5), (3, 4, 5),
                (1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5),
                (2, 2, 4, 5), (2, 3, 4, 5), (1, 2, 2, 4, 5),
                (1, 2, 3, 4, 5), (2, 2, 3, 4, 5),
            ),
            ((3, 4, 5), (1, 3, 4, 5)),
            ((1, 3, 4, 5), (1, 2, 3, 4, 5)),
            "80901a99a0958325589bb2a9fdd07a704645f75fe98fa7d66cc02d01efa78a89",
        ),
        (
            ("F12", "F23"),
            5,
            (
                (2, 3, 4), (1, 2, 3, 4), (2, 3, 4, 5),
                (1, 2, 3, 4, 5), (2, 2, 3, 4, 5),
            ),
            (),
            ((2, 3, 4, 5), (1, 2, 3, 4, 5), (2, 2, 3, 4, 5)),
            "221f1805df7c37f948f70cf4e926a9174eac66facc372d38d1b8b98117879ecf",
        ),
        (
            ("F13", "F23"),
            0,
            (),
            (),
            (),
            "684e83ead4c16ae8fa9d955c529083ec63dd3f9eac67e72b291bfb21a9d222a0",
        ),
    )
    expected_triple_empty = (
        (5,), (3, 4), (4, 5), (1, 3, 4), (1, 4, 5), (2, 3, 4),
        (2, 4, 5), (3, 4, 5), (1, 2, 3, 4), (1, 2, 4, 5),
        (1, 3, 4, 5), (2, 2, 4, 5), (2, 3, 4, 5),
        (1, 2, 2, 4, 5), (1, 2, 3, 4, 5), (2, 2, 3, 4, 5),
    )
    expected_iteration = (
        (
            (0, 1, 1),
            (
                (
                    0, ("F12",), 2859, 268, 69, 7, 5, 6,
                    "7892a85be1005e3070204add40502252fd2142ad3f91ea5610a153a8cd3a0845",
                ),
                (
                    1, ("F13", "F23"), 4595, 215, 2651, 5, 6, 6,
                    "acec2a32647c705bd53839a1bb12fcb32a49a35dcc70295e36c6506e7930425a",
                ),
            ),
        ),
    )
    require(
        stable_digest(triple_states)
        == "5835c0ad5c8d0747bf366a11a14c6b5f80f2326317e4539c13bc16b713711180"
        and len(tasks) == 943,
        "restricted universe or response-task closure drift",
    )
    require(
        decorated_face_digests == (
            ("F12", "c9def3e31568c4abe88deefc2b9835e31939dda60b9c62b7ab92c258407db2bd"),
            ("F13", "f399eb76e5bc9851941d4483295f8a060a921f959a92365fbc97d71617c66df4"),
            ("F23", "44020e284ad310da061b18b9b0bbb4d372bdbf458aa4aafdb3e7f20a6a113618"),
        ),
        "move-decorated face digest drift",
    )
    require(pair_summary == expected_pair_summary, "decorated pair atlas drift")
    require(
        triple_empty == expected_triple_empty
        and decorated_proper_helly == ((5,), (3, 4))
        and triple_digest
        == "4b9e8ecf8548357b19836cd8575a0ecadff3798001a078d3ce52415b6fb8d793",
        "decorated triple atlas drift",
    )
    require(
        local_histogram == ((1, 78), (2, 16))
        and tuple(state for state, width in local_widths if width == 2)
        == expected_triple_empty,
        "decorated local supplied-origin width drift",
    )
    require(
        minimum_fixed == 2
        and optimal_colourings == ((0, 1, 1), (1, 0, 0))
        and canonical_optimal == ((0, 1, 1),)
        and not uncoloured_section
        and coloured_section,
        "decorated fixed supplied-origin classification drift",
    )
    require(iteration_records == expected_iteration, "finite section iteration drift")
    require(
        fixed_pair_extensions == (
            (1, ()),
            (2, ((0, 1, 1), (1, 0, 0))),
            (3, (
                (0, 1, 2), (0, 2, 1), (1, 0, 2),
                (1, 2, 0), (2, 0, 1), (2, 1, 0),
            )),
        ),
        "restricted 94-state (16,17) rescue drift",
    )
    require(
        first_failure
        == (("F12", "F13"), (4, 5), "row_successor_coupling", (5,),
            ((1, 4, 5), (3, 4, 5)))
        and first_incompatible
        == (("F12", "F23"), (2, 3, 4, 5), (5, 12)),
        "minimal hostile witness drift",
    )

    print("THM-3286 MOVE-DECORATED THREE-FACE SECTION PARTIAL SCOUT")
    print("status=FINITE-EXACT PARTIAL;no_theorem_ID;no_promotion")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path.relative_to(ROOT)}")
    print("dependency_hashes=END")
    print("face_profiles=" + repr(face_profile))
    print(
        f"triple_universe={len(triple_states)};"
        f"triple_state_digest={stable_digest(triple_states)};"
        f"response_tasks={len(tasks)}"
    )
    print(
        "typed_map=B_f(n)={(row,successor):strict_ascent_and_reset_directed};"
        "projection_(row,successor)->row_recovers_THM3286_triple_atlas=PASS"
    )
    print("decorated_face_digests=" + repr(decorated_face_digests))
    for summary in pair_summary:
        print("decorated_pair_summary=" + repr(summary))
    print(
        f"decorated_triple_empty_count={len(triple_empty)};"
        f"decorated_triple_empty_locus={triple_empty};"
        f"proper_decorated_Helly_locus={decorated_proper_helly};"
        f"triple_digest={triple_digest}"
    )
    print(
        f"local_minimum_supplied_origin_histogram={local_histogram};"
        f"local_width_two_locus={local_width_two_locus};local_width_three_locus=()"
    )
    print(
        f"minimum_fixed_supplied_origin_alphabet={minimum_fixed};"
        f"optimal_colourings={optimal_colourings};"
        f"canonical_optimal_colourings={canonical_optimal}"
    )
    print(
        f"uncoloured_global_static_section={uncoloured_section};"
        f"minimum_coloured_global_static_section={coloured_section};"
        "iteration_record_schema=(label,face_block,decorated_choices,"
        "distinct_internal_targets,exit_choices,forced_exit_states,"
        "max_arbitrary_internal_edges,max_lex_moves_to_exit,lex_digest);"
        f"iteration_records={iteration_records}"
    )
    print(
        "control_state_(5)="
        + repr((
            state5_successors,
            state5_decorated_counts,
            state5_record[1],
            state5_record[2],
            state5_record[3],
            state5_pair_summary,
        ))
    )
    print("inherited_conflict_controls=" + repr(inherited_conflict_controls))
    print(
        "canonical_row_alphabet_(16,17)_restricted_94_extensions="
        + repr(fixed_pair_extensions)
    )
    print(
        "restriction_boundary=PINNED_THM3286_proves_(16,17)_has_no_extension_"
        "on_all_8394_face_state_decisions;restricted_94_rescue_is_not_a_"
        "full_domain_controller"
    )
    print(
        "minimal_projection_failure="
        + repr((
            first_pair,
            first_state,
            first_mechanism,
            first_rows,
            first_targets,
            first_record[3],
            tuple((face, successors[face][first_state]) for face in first_pair),
            first_row_realisers,
        ))
    )
    print("minimal_incompatible_reset_hostile=" + repr(first_incompatible))
    print(
        "first_failed_THM3286_implication=common_projected_row_does_not_imply_"
        "a_common_(row,physical_successor);incompatible_resets_or_row_move_"
        "coupling_are_the_missing_coordinates"
    )
    print(
        f"source_ast=(assert_nodes={ASSERT_NODES},float_literals={FLOAT_LITERALS})"
    )
    print(
        "scope=exactly_three_named_bank_I2_faces_and_inherited_94_triple_states;"
        "iteration_stops_on_exit_from_inherited_graph;"
        "no_FC3_SFC3_GMC_positivity_or_unbounded_history_claim"
    )
    print("ALL FINITE-EXACT PARTIAL CHECKS PASSED")


if __name__ == "__main__":
    freeze_support()
    main()
