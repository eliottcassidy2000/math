#!/usr/bin/env python3
"""Finite-exact all-choice closure scout for the move-decorated three-face lane.

Starting from the 94 common decision states of the three THM-3286 faces, use
the surviving supplied-origin colouring

    {F12} | {F13,F23}

and close each colour lane under *every* common decorated pair
``(row, reset-directed physical successor)``.  A face leaves its lane only
when that successor is its pinned reset; faces are never silently acquired,
and every retained face is checked against its own physical decision
universe.  This is the smallest lane-wise all-choice closure of the 94-state
carrier, not a controller on the full 8,394 face/state decisions.

The computation also repeats the closure with rows restricted to ``{16,17}``
and records the first reachable obstruction.  It is FINITE-EXACT PARTIAL
evidence only: no FC(3), SFC(3), GMC, positivity, or unbounded-history claim
is made.
"""

from __future__ import annotations

import ast
from collections import Counter, deque
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from multiprocessing import Pool, freeze_support
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
PRIOR_SCRIPT = (
    COMP / "fc3_three_support_face_move_decorated_section_partial_scout_20260803.py"
)
PRIOR_OUTPUT = (
    RESULTS / "fc3_three_support_face_move_decorated_section_partial_scout_20260803.out"
)
THM3286_THEOREM = (
    ROOT
    / "01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md"
)
BASE_SCRIPT = COMP / "gmc_complete_physical_bank_unique_reset_thm3238.py"
BASE_OUTPUT = RESULTS / "gmc_complete_physical_bank_unique_reset_thm3238.out"
PINNED = {
    PRIOR_SCRIPT:
        "227c1558da114429a92f08a41e87f919b667ab42f346b4d4e80c5bdeabfe3183",
    PRIOR_OUTPUT:
        "96d8406a921e255647f5a8b3cd66d34ee8777705f117c74d0a15673cd8f7d69f",
    THM3286_THEOREM:
        "fd48e7e142d44b5c60baec7e7bc1b59e2640f672f8de019d3ed2a2e7d969163c",
    BASE_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    BASE_OUTPUT:
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
}

COLOUR_BLOCKS = (("F12",), ("F13", "F23"))
CANONICAL_ROWS = frozenset((16, 17))
EXPECTED_CLOSURE_BLOCKERS = (
    (3, ("F13", "F23"), (1, 2, 3, 4, 5, 6, 7, 8)),
    (3, ("F13", "F23"), (1, 3, 4, 5, 6, 7, 8)),
    (3, ("F13", "F23"), (2, 3, 4, 5, 6, 7, 8)),
    (3, ("F13", "F23"), (3, 4, 5, 6, 7, 8)),
)
EXPECTED_FULL_UNRESTRICTED_BLOCKERS = (
    (("F13", "F23"), (1, 2, 3, 4, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (1, 2, 3, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (1, 3, 4, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (1, 3, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (2, 3, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (3, 4, 4, 5, 6, 7, 8)),
    (("F13", "F23"), (3, 4, 5, 6, 7, 8)),
)


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
    "pinned THM-3286 full-domain (16,17) control wording changed",
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


PRIOR_SPEC = spec_from_file_location("move_decorated_partial_for_closure", PRIOR_SCRIPT)
prior = module_from_spec(PRIOR_SPEC)
PRIOR_SPEC.loader.exec_module(prior)

FACE_ORDER = prior.FACE_ORDER
FACES = prior.FACES


def response_task(task):
    face, state = task
    data = FACES[face]
    vectors = prior.base.coefficient_vectors(
        1, prior.base.BANK, data["a"], data["b"], state
    )
    return (face, state), prior.row_coordinates(vectors)


def common_decorated(responses, active, state, allowed_rows=None):
    """All literal pairs common to every currently active typed face."""

    fibres = []
    for face in active:
        data = FACES[face]
        require(state in data["decisions"], ("ill-typed source", active, face, state))
        source = responses[face, state]
        pairs = frozenset(
            (row, target)
            for target in prior.toward_reset_neighbours(
                state, data["poles"], data["reset"]
            )
            for row in range(1, 23)
            if responses[face, target][row - 1] > source[row - 1]
            and (allowed_rows is None or row in allowed_rows)
        )
        fibres.append(pairs)
    return tuple(sorted(frozenset.intersection(*fibres)))


def successor_active(active, target):
    """Drop exactly reset faces and check every retained face remains typed."""

    for face in active:
        require(
            target in FACES[face]["states"],
            ("common physical edge lost its face typing", active, face, target),
        )
    retained = tuple(face for face in active if target != FACES[face]["reset"])
    require(
        all(target in FACES[face]["decisions"] for face in retained),
        ("retained successor is not a decision", active, retained, target),
    )
    return retained


def closure(responses, starts, allowed_rows=None):
    """Breadth-first all-choice closure with typed reset-face deletion."""

    distance = {node: 0 for node in starts}
    parents = {node: None for node in starts}
    queue = deque(starts)
    edges = []
    terminal_edges = []
    blockers = []
    potential_drops = []
    while queue:
        active, state = queue.popleft()
        choices = common_decorated(
            responses, active, state, allowed_rows=allowed_rows
        )
        if not choices:
            blockers.append((distance[active, state], active, state))
            continue
        source_potential = sum(
            prior.reset_distance(state, FACES[face]["reset"]) for face in active
        )
        for row, target in choices:
            retained = successor_active(active, target)
            target_potential = sum(
                prior.reset_distance(target, FACES[face]["reset"])
                for face in retained
            )
            require(
                target_potential == source_potential - len(active),
                ("reset potential law failed", active, state, row, target),
            )
            potential_drops.append(source_potential - target_potential)
            edge = ((active, state), row, (retained, target))
            edges.append(edge)
            if not retained:
                terminal_edges.append(edge)
                continue
            node = (retained, target)
            if node not in distance:
                distance[node] = distance[active, state] + 1
                parents[node] = ((active, state), row)
                queue.append(node)

    def path_to(node):
        path = []
        cursor = node
        while parents[cursor] is not None:
            source, row = parents[cursor]
            path.append((source, row, cursor))
            cursor = source
        return tuple(reversed(path))

    blocker_paths = tuple(
        (record, path_to((record[1], record[2]))) for record in sorted(blockers)
    )
    node_records = tuple(sorted(
        ((depth, active, state) for (active, state), depth in distance.items()),
        key=lambda record: (record[0], record[1], len(record[2]), record[2]),
    ))
    edge_records = tuple(sorted(
        edges,
        key=lambda edge: (
            distance[edge[0]], edge[0][0], len(edge[0][1]), edge[0][1],
            edge[1], edge[2][0], len(edge[2][1]), edge[2][1],
        ),
    ))
    return {
        "distance": distance,
        "nodes": node_records,
        "edges": edge_records,
        "terminal_edges": tuple(sorted(terminal_edges, key=repr)),
        "blockers": tuple(sorted(blockers)),
        "blocker_paths": blocker_paths,
        "potential_drops": tuple(potential_drops),
    }


def first_restricted_failure(unrestricted, restricted):
    """First unrestricted-reachable node at which the fixed row pair fails."""

    restricted_nodes = frozenset(
        (active, state) for _depth, active, state in restricted["nodes"]
    )
    candidates = tuple(
        record for record in unrestricted["nodes"]
        if (record[1], record[2]) in restricted_nodes
        and not common_decorated(
            RESPONSES, record[1], record[2], allowed_rows=CANONICAL_ROWS
        )
    )
    require(candidates, "the canonical row pair never fails on its reachable closure")
    first = min(candidates, key=lambda record: (
        record[0], record[1], len(record[2]), record[2]
    ))
    node = (first[1], first[2])
    path = next(
        path for record, path in restricted["blocker_paths"]
        if (record[1], record[2]) == node
    )
    return first, path, common_decorated(RESPONSES, first[1], first[2])


def face_decorated_fibre(responses, face, state, allowed_rows=None):
    """One face's typed decorated fibre, exposed for failure anatomy."""

    data = FACES[face]
    source = responses[face, state]
    return tuple(sorted(
        (row, target)
        for target in prior.toward_reset_neighbours(
            state, data["poles"], data["reset"]
        )
        for row in range(1, 23)
        if responses[face, target][row - 1] > source[row - 1]
        and (allowed_rows is None or row in allowed_rows)
    ))


RESPONSES = {}


def main():
    global RESPONSES

    triple_states = tuple(sorted(
        set.intersection(*(set(FACES[face]["decisions"]) for face in FACE_ORDER)),
        key=lambda state: (len(state), state),
    ))
    require(len(triple_states) == 94, "inherited triple universe drift")

    # Full inherited face universes are the cheapest exact response cache and
    # guarantee that every all-choice frontier target is available without an
    # adaptive or policy-dependent truncation.
    tasks = tuple(
        (face, state)
        for face in FACE_ORDER
        for state in FACES[face]["states"]
    )
    with Pool(processes=4) as pool:
        RESPONSES = dict(pool.imap_unordered(response_task, tasks, chunksize=8))
    require(len(RESPONSES) == sum(len(FACES[f]["states"]) for f in FACE_ORDER),
            "a full-face response task was lost")

    starts = tuple(
        (block, state) for block in COLOUR_BLOCKS for state in triple_states
    )
    unrestricted = closure(RESPONSES, starts)
    restricted = closure(RESPONSES, starts, allowed_rows=CANONICAL_ROWS)
    first_1617 = first_restricted_failure(unrestricted, restricted)

    full_raw_states = set.union(*(set(FACES[face]["decisions"]) for face in FACE_ORDER))
    full_coloured_nodes = tuple(sorted(
        (
            tuple(face for face in block if state in FACES[face]["decisions"]),
            state,
        )
        for state in full_raw_states
        for block in COLOUR_BLOCKS
        if any(state in FACES[face]["decisions"] for face in block)
    ))
    full_unrestricted_blockers = tuple(sorted(
        (active, state) for active, state in full_coloured_nodes
        if not common_decorated(RESPONSES, active, state)
    ))
    full_restricted_blockers = tuple(sorted(
        (active, state) for active, state in full_coloured_nodes
        if not common_decorated(
            RESPONSES, active, state, allowed_rows=CANONICAL_ROWS
        )
    ))

    seam = (3, 4, 5, 6, 7, 8)
    seam_macro_probe = tuple(
        (
            (left_root, right_root),
            target,
            typed_faces,
            tuple(
                row for row in range(1, 23)
                if len(typed_faces) == 2
                and all(
                    RESPONSES[face, target][row - 1]
                    > RESPONSES[face, seam][row - 1]
                    for face in typed_faces
                )
            ),
            tuple(
                prior.reset_distance(target, FACES[face]["reset"])
                for face in typed_faces
            ),
            (("F13", "F23"), target) in full_unrestricted_blockers,
        )
        for left_root in (1, 3)
        for right_root in (2, 4)
        for target in (tuple(sorted(seam + (left_root, right_root))),)
        for typed_faces in (tuple(
            face for face in COLOUR_BLOCKS[1]
            if Counter(target) <= Counter(FACES[face]["poles"])
        ),)
    )

    pair_blockers = tuple(record[2] for record in unrestricted["blockers"])
    blocker_anatomy = tuple(
        (
            state,
            tuple(
                prior.reset_distance(state, FACES[face]["reset"])
                for face in COLOUR_BLOCKS[1]
            ),
            tuple(
                (
                    face,
                    prior.toward_reset_neighbours(
                        state, FACES[face]["poles"], FACES[face]["reset"]
                    ),
                    tuple(sorted({
                        row for row, _target in face_decorated_fibre(
                            RESPONSES, face, state
                        )
                    })),
                )
                for face in COLOUR_BLOCKS[1]
            ),
            tuple(sorted(set.intersection(*(
                {
                    row for row, _target in face_decorated_fibre(
                        RESPONSES, face, state
                    )
                }
                for face in COLOUR_BLOCKS[1]
            )))),
            tuple(sorted(set.intersection(*(
                set(prior.toward_reset_neighbours(
                    state, FACES[face]["poles"], FACES[face]["reset"]
                ))
                for face in COLOUR_BLOCKS[1]
            )))),
        )
        for state in pair_blockers
    )
    pair_start_lengths = tuple(
        (
            state,
            (
                sum(
                    prior.reset_distance(state, FACES[face]["reset"])
                    for face in COLOUR_BLOCKS[1]
                ) - 4
            ) // 2,
        )
        for state in triple_states
    )
    singleton_start_lengths = tuple(
        (state, prior.reset_distance(state, FACES["F12"]["reset"]))
        for state in triple_states
    )
    singleton_length_by_state = dict(singleton_start_lengths)
    pair_vs_singleton_slack = tuple(
        (state, pair_length - singleton_length_by_state[state])
        for state, pair_length in pair_start_lengths
    )

    require(
        set(unrestricted["potential_drops"]) == {1, 2},
        "an edge has an impossible active-face potential drop",
    )

    def summary(data):
        node_type_hist = tuple(sorted(Counter(
            active for _depth, active, _state in data["nodes"]
        ).items()))
        depth_hist = tuple(sorted(Counter(
            depth for depth, _active, _state in data["nodes"]
        ).items()))
        return (
            len(data["nodes"]), len(data["edges"]), len(data["terminal_edges"]),
            len(data["blockers"]), node_type_hist, depth_hist,
            stable_digest(data["nodes"]), stable_digest(data["edges"]),
        )

    unrestricted_summary = summary(unrestricted)
    restricted_summary = summary(restricted)
    active_changes = tuple(
        edge for edge in unrestricted["edges"]
        if edge[0][0] != edge[2][0]
    )
    active_change_histogram = tuple(sorted(Counter(
        (edge[0][0], edge[2][0]) for edge in active_changes
    ).items()))
    singleton_length_histogram = tuple(sorted(Counter(
        length for _state, length in singleton_start_lengths
    ).items()))
    pair_length_histogram = tuple(sorted(Counter(
        length for _state, length in pair_start_lengths
    ).items()))
    schedule_slack_histogram = tuple(sorted(Counter(
        slack for _state, slack in pair_vs_singleton_slack
    ).items()))
    zero_slack_locus = tuple(
        state for state, slack in pair_vs_singleton_slack if slack == 0
    )
    unique_physical_states = frozenset(
        state for _depth, _active, state in unrestricted["nodes"]
    )

    require(
        unrestricted_summary == (
            854, 30_853, 69, 4,
            ((("F12",), 94), (("F13", "F23"), 760)),
            ((0, 188), (1, 283), (2, 285), (3, 97), (4, 1)),
            "052b60f08574faf4094f21583caa9efc7514e7e4d8b25f0be7ea83a8f5b378bc",
            "a298b6048747e77ea59dd09ad7b451f7f08d3c80792e3469c1e59c57acf556a3",
        )
        and len(unique_physical_states) == 760
        and unrestricted["blockers"] == EXPECTED_CLOSURE_BLOCKERS,
        "unrestricted typed closure drift",
    )
    require(
        restricted_summary == (
            854, 3_126, 7, 4,
            ((("F12",), 94), (("F13", "F23"), 760)),
            ((0, 188), (1, 283), (2, 285), (3, 97), (4, 1)),
            "052b60f08574faf4094f21583caa9efc7514e7e4d8b25f0be7ea83a8f5b378bc",
            "b75e8106fdb2053e772fc6b168b8ade76e5e4e814fe0243a671a43c3bf746aec",
        )
        and restricted["nodes"] == unrestricted["nodes"]
        and restricted["blockers"] == unrestricted["blockers"],
        "canonical row-pair closure stopped before the universal seam",
    )
    require(
        tuple(record[0] for record in blocker_anatomy)
        == tuple(record[2] for record in EXPECTED_CLOSURE_BLOCKERS)
        and all(sum(record[1]) == 4 for record in blocker_anatomy)
        and all(16 in record[3] and not record[4] for record in blocker_anatomy),
        "shared-pair reset-direction fork anatomy drift",
    )
    require(
        singleton_length_histogram
        == ((1, 7), (2, 20), (3, 30), (4, 25), (5, 11), (6, 1))
        and pair_length_histogram
        == ((3, 4), (4, 15), (5, 27), (6, 25), (7, 16), (8, 6), (9, 1)),
        "maximal path-length certificate drift",
    )
    require(
        schedule_slack_histogram
        == ((0, 7), (1, 24), (2, 16), (3, 16), (4, 31))
        and zero_slack_locus
        == ((3,), (4,), (5,), (3, 4), (3, 5), (4, 5), (3, 4, 5)),
        "pointwise singleton-before-fork timing law drift",
    )
    require(
        len(full_coloured_nodes) == 6_668
        and full_unrestricted_blockers == EXPECTED_FULL_UNRESTRICTED_BLOCKERS
        and len(full_restricted_blockers) == 38
        and stable_digest(full_restricted_blockers)
        == "40f63cff4f3a116364e280b6e34fc22df1023ab971b6c1ffd491f031145d8516"
        and set(full_unrestricted_blockers) <= set(full_restricted_blockers),
        "full active-context pinned control drift",
    )
    require(
        seam_macro_probe == (
            (
                (1, 2), (1, 2, 3, 4, 5, 6, 7, 8), ("F13", "F23"),
                (2, 5, 8, 9, 11, 12, 14, 16, 18, 22), (2, 2), True,
            ),
            (
                (1, 4), (1, 3, 4, 4, 5, 6, 7, 8), ("F13", "F23"),
                (2, 5, 8, 9, 11, 12, 14, 16, 18, 22), (2, 2), True,
            ),
            ((3, 2), (2, 3, 3, 4, 5, 6, 7, 8), ("F13",), (), (2,), False),
            ((3, 4), (3, 3, 4, 4, 5, 6, 7, 8), ("F13",), (), (2,), False),
        ),
        "two-pole seam macro probe drift",
    )
    require(
        len(active_changes) == 69
        and active_change_histogram == (((('F12',), ()), 69),)
        and stable_digest(active_changes)
        == "ed47f8e6ca02c4aaf3a2197bd16322c8501a09987bc9cab43289010dbb9635b8"
        and all(edge[0][0] == ("F12",) and not edge[2][0]
                for edge in unrestricted["terminal_edges"]),
        "active-face change or terminal typing drift",
    )
    require(
        first_1617 == (
            (3, ("F13", "F23"), (3, 4, 5, 6, 7, 8)),
            (
                ((
                    ("F13", "F23"), (3, 4, 5)
                ), 17, (("F13", "F23"), (3, 4, 5, 6))),
                ((
                    ("F13", "F23"), (3, 4, 5, 6)
                ), 16, (("F13", "F23"), (3, 4, 5, 6, 7))),
                ((
                    ("F13", "F23"), (3, 4, 5, 6, 7)
                ), 16, (("F13", "F23"), (3, 4, 5, 6, 7, 8))),
            ),
            (),
        ),
        "minimal canonical-row chronology obstruction drift",
    )

    print("THM-3286 MOVE-DECORATED ALL-CHOICE CLOSURE PARTIAL SCOUT")
    print("status=FINITE-EXACT PARTIAL;no_theorem_ID;no_promotion")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path.relative_to(ROOT)}")
    print("dependency_hashes=END")
    print(
        "closure_type=lane_wise_all_choice_from_94_common_states;"
        "colour_blocks=((F12),(F13,F23));faces_drop_only_at_their_exact_reset;"
        "faces_are_never_silently_acquired"
    )
    print(
        f"full_face_response_tasks={len(tasks)};"
        f"face_state_counts={tuple((f, len(FACES[f]['states'])) for f in FACE_ORDER)};"
        f"start_nodes={len(starts)};typed_closure_nodes={len(unrestricted['nodes'])};"
        f"underlying_untyped_physical_states={len(unique_physical_states)}"
    )
    print("unrestricted_closure_summary=" + repr(unrestricted_summary))
    print("unrestricted_blockers=" + repr(unrestricted["blockers"]))
    print("unrestricted_blocker_paths=" + repr(unrestricted["blocker_paths"]))
    print("shared_pair_blocker_anatomy=" + repr(blocker_anatomy))
    print(
        "two_pole_seam_macro_probe_at_C=" + repr(seam_macro_probe)
        + ";boundary=macro_has_no_common_first_one_pole_step_and_both_typed_"
        "targets_are_potential_neutral_unrestricted_blockers"
    )
    print(
        f"full_active_coloured_context_universe={len(full_coloured_nodes)};"
        f"full_unrestricted_blocker_count={len(full_unrestricted_blockers)};"
        f"full_unrestricted_blocker_digest={stable_digest(full_unrestricted_blockers)};"
        f"full_unrestricted_first_blockers={full_unrestricted_blockers[:12]}"
    )
    print(
        f"full_canonical_rows_(16,17)_blocker_count={len(full_restricted_blockers)};"
        f"full_canonical_rows_(16,17)_blocker_digest="
        f"{stable_digest(full_restricted_blockers)};"
        f"full_canonical_rows_(16,17)_first_blockers={full_restricted_blockers[:12]}"
    )
    print(
        "maximal_path_length_histograms=(F12_to_reset,"
        + repr(singleton_length_histogram)
        + ";F13_F23_to_blocker,"
        + repr(pair_length_histogram)
        + ")"
    )
    print(
        f"pointwise_pair_fork_time_minus_F12_reset_time={schedule_slack_histogram};"
        f"zero_slack_locus={zero_slack_locus};"
        "interpretation=F12_finishes_no_later_than_the_pair_fork_but_reusing_"
        "its_origin_would_require_an_unproved_rebinding_or_copy_sidecar"
    )
    print(
        f"active_face_change_edge_count={len(active_changes)};"
        f"active_face_change_histogram={active_change_histogram};"
        f"active_face_change_digest={stable_digest(active_changes)};"
        f"active_face_change_sample={active_changes[:12]}"
    )
    print("canonical_rows_(16,17)_closure_summary=" + repr(restricted_summary))
    print("canonical_rows_(16,17)_blockers=" + repr(restricted["blockers"]))
    print("first_canonical_rows_(16,17)_failure=" + repr(first_1617))
    print(
        "acyclicity_certificate=every_edge_drops_sum_of_active_reset_distances_"
        "by_exactly_the_number_of_active_faces;"
        f"observed_drops={tuple(sorted(set(unrestricted['potential_drops'])))}"
    )
    print(
        f"source_ast=(assert_nodes={ASSERT_NODES},float_literals={FLOAT_LITERALS})"
    )
    print(
        "scope=smallest_lane_wise_all_choice_closure_from_the_inherited_94_"
        "states_inside_the_three_pinned_face_universes;not_a_shared_state_"
        "interleaving_or_full_8394_decision_controller;"
        "no_FC3_SFC3_GMC_positivity_or_unbounded_history_claim"
    )
    print("ALL FINITE-EXACT PARTIAL CHECKS PASSED")


if __name__ == "__main__":
    freeze_support()
    main()
