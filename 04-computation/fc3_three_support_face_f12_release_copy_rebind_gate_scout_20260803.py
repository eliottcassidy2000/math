#!/usr/bin/env python3
"""Finite-exact origin-release and explicit copy/rebind gate for three faces.

The inherited all-choice carrier has two separate lanes,

    {F12} | {F13,F23}.

This scout first audits what an F12 reset actually releases in that carrier.
It then tests the four reachable shared-lane blockers without adding an
operation.  Finally, and separately, it adds one explicit sidecar operation:
after a certified F12 reset, copy the *current shared-lane state* into the
freed origin and split {F13,F23} into two singleton obligations.  The copy is
not a response edge, has no row label, and preserves the total reset
potential; every subsequent edge is again a literal row-labelled one-pole
physical successor.

This is a finite-exact partial experiment.  The copy rule is an added carrier,
not inherited mathematics, and no FC(3), SFC(3), GMC, positivity, or shared
physical interleaving conclusion is claimed.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict, deque
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from multiprocessing import Pool, freeze_support
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
PRIOR_SCRIPT = (
    COMP / "fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.py"
)
PRIOR_OUTPUT = (
    RESULTS / "fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.out"
)
THM3286 = (
    ROOT
    / "01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md"
)
THM3278 = (
    ROOT
    / "01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md"
)
PINNED = {
    PRIOR_SCRIPT:
        "74cd5e850d465db028bb810a246cc6fddec03cad64a20760145cedd29e9cb9a6",
    PRIOR_OUTPUT:
        "616758c2e174d2ffa0cedd37d92a78e6e61ac48e876853973901d2e81080893f",
    THM3286:
        "fd48e7e142d44b5c60baec7e7bc1b59e2640f672f8de019d3ed2a2e7d969163c",
    THM3278:
        "d6c0fe74d61b11951c9efa734fadbba0bec7f86a76621812b584ec1590eba89e",
}

PAIR = ("F13", "F23")
F12 = ("F12",)
ROWS_1617 = frozenset((16, 17))
EXPECTED_REACHABLE_BLOCKERS = (
    (1, 2, 3, 4, 5, 6, 7, 8),
    (1, 3, 4, 5, 6, 7, 8),
    (2, 3, 4, 5, 6, 7, 8),
    (3, 4, 5, 6, 7, 8),
)
EXPECTED_FULL_UNIVERSAL_BLOCKERS = (
    (1, 2, 3, 4, 4, 5, 6, 7, 8),
    (1, 2, 3, 4, 5, 6, 7, 8),
    (1, 3, 4, 4, 5, 6, 7, 8),
    (1, 3, 4, 5, 6, 7, 8),
    (2, 3, 4, 5, 6, 7, 8),
    (3, 4, 4, 5, 6, 7, 8),
    (3, 4, 5, 6, 7, 8),
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
    "The origin bit is therefore an exact cut coordinate. It is not a"
    in " ".join(THM3278.read_text(encoding="utf-8").split()),
    "THM-3278 origin-bit typing changed",
)
require(
    "faces are never silently acquired"
    in PRIOR_SCRIPT.read_text(encoding="utf-8"),
    "inherited active-face law changed",
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

PRIOR_SPEC = spec_from_file_location("fc3_all_choice_for_rebind_gate", PRIOR_SCRIPT)
prior = module_from_spec(PRIOR_SPEC)
PRIOR_SPEC.loader.exec_module(prior)

FACES = prior.FACES
FACE_ORDER = prior.FACE_ORDER


def response_task(task):
    face, state = task
    data = FACES[face]
    vectors = prior.prior.base.coefficient_vectors(
        1, prior.prior.base.BANK, data["a"], data["b"], state
    )
    return (face, state), prior.prior.row_coordinates(vectors)


def face_choices(responses, face, state, allowed_rows=None):
    if state == FACES[face]["reset"]:
        return ()
    require(state in FACES[face]["decisions"], ("ill-typed singleton", face, state))
    return prior.face_decorated_fibre(
        responses, face, state, allowed_rows=allowed_rows
    )


def closure_adjacency(data, active):
    adjacency = defaultdict(list)
    for source, row, target in data["edges"]:
        if source[0] == active:
            adjacency[source[1]].append((row, target[0], target[1]))
    return {
        state: tuple(sorted(edges)) for state, edges in adjacency.items()
    }


def pair_endpoint_census(data, starts):
    """All row-labelled blocker endpoints from each start, by DAG recursion."""

    adjacency = closure_adjacency(data, PAIR)
    blockers = frozenset(state for _depth, active, state in data["blockers"]
                         if active == PAIR)

    @lru_cache(maxsize=None)
    def endpoint_counts(state):
        edges = adjacency.get(state, ())
        if not edges:
            require(state in blockers, ("unexpected pair sink", state))
            return ((state, 1),)
        total = Counter()
        for _row, active, target in edges:
            require(active == PAIR, ("pair active set changed", state, target))
            for endpoint, count in endpoint_counts(target):
                total[endpoint] += count
        return tuple(sorted(total.items()))

    records = tuple((state, endpoint_counts(state)) for state in starts)
    return records


def f12_release_census(data, starts):
    """All row-labelled F12 paths end at the unique exact reset."""

    adjacency = closure_adjacency(data, F12)
    reset = FACES["F12"]["reset"]

    @lru_cache(maxsize=None)
    def count_paths(state):
        edges = adjacency.get(state, ())
        require(edges, ("F12 path blocked before reset", state))
        total = 0
        for _row, active, target in edges:
            if active:
                require(active == F12, ("unexpected F12 active change", state, target))
                total += count_paths(target)
            else:
                require(target == reset, ("wrong F12 terminal", state, target))
                total += 1
        return total

    return tuple(
        (
            state,
            prior.prior.reset_distance(state, reset),
            count_paths(state),
        )
        for state in starts
    )


def singleton_completion(responses, face, start, allowed_rows=None):
    """Close every literal singleton choice and count row-labelled paths."""

    reset = FACES[face]["reset"]
    queue = deque((start,))
    seen = {start}
    edges = []
    blockers = []
    while queue:
        state = queue.popleft()
        if state == reset:
            continue
        choices = face_choices(responses, face, state, allowed_rows=allowed_rows)
        if not choices:
            blockers.append(state)
            continue
        for row, target in choices:
            require(
                prior.prior.reset_distance(target, reset) + 1
                == prior.prior.reset_distance(state, reset),
                ("singleton edge is not reset-directed", face, state, row, target),
            )
            edges.append((state, row, target))
            if target not in seen:
                seen.add(target)
                queue.append(target)

    by_source = defaultdict(list)
    for source, row, target in edges:
        by_source[source].append((row, target))

    @lru_cache(maxsize=None)
    def count_paths(state):
        if state == reset:
            return 1
        if state in blockers:
            return 0
        return sum(count_paths(target) for _row, target in by_source[state])

    lex_path = []
    cursor = start
    while cursor != reset and cursor not in blockers:
        row, target = min(by_source[cursor])
        lex_path.append((cursor, row, target))
        cursor = target

    return (
        face,
        start,
        prior.prior.reset_distance(start, reset),
        len(seen),
        len(edges),
        tuple(sorted(blockers, key=lambda state: (len(state), state))),
        count_paths(start),
        tuple(lex_path),
        stable_digest(tuple(sorted(edges))),
    )


def macro_probe(responses):
    seam = (3, 4, 5, 6, 7, 8)
    return tuple(
        (
            (left_root, right_root),
            target,
            typed_faces,
            tuple(
                row for row in range(1, 23)
                if len(typed_faces) == 2
                and all(
                    responses[face, target][row - 1]
                    > responses[face, seam][row - 1]
                    for face in typed_faces
                )
            ),
            tuple(
                prior.prior.reset_distance(target, FACES[face]["reset"])
                for face in typed_faces
            ),
            bool(prior.common_decorated(responses, PAIR, seam)),
            all(
                target not in prior.prior.toward_reset_neighbours(
                    seam, FACES[face]["poles"], FACES[face]["reset"]
                )
                for face in typed_faces
            ),
        )
        for left_root in (1, 3)
        for right_root in (2, 4)
        for target in (tuple(sorted(seam + (left_root, right_root))),)
        for typed_faces in (tuple(
            face for face in PAIR
            if Counter(target) <= Counter(FACES[face]["poles"])
        ),)
    )


def main():
    tasks = tuple(
        (face, state)
        for face in FACE_ORDER
        for state in FACES[face]["states"]
    )
    with Pool(processes=4) as pool:
        responses = dict(pool.imap_unordered(response_task, tasks, chunksize=8))
    require(len(responses) == 8_397, "full response cache drift")

    triple_states = tuple(sorted(
        set.intersection(*(set(FACES[face]["decisions"]) for face in FACE_ORDER)),
        key=lambda state: (len(state), state),
    ))
    require(len(triple_states) == 94, "triple-state universe drift")
    starts = tuple(
        (block, state) for block in prior.COLOUR_BLOCKS for state in triple_states
    )
    unrestricted = prior.closure(responses, starts)
    restricted = prior.closure(responses, starts, allowed_rows=ROWS_1617)
    reachable_blockers = tuple(record[2] for record in unrestricted["blockers"])
    require(
        reachable_blockers == EXPECTED_REACHABLE_BLOCKERS
        and restricted["blockers"] == unrestricted["blockers"],
        "inherited four-blocker closure drift",
    )

    # What the current carrier actually releases.
    release_edges = tuple(edge for edge in unrestricted["terminal_edges"]
                          if edge[0][0] == F12)
    restricted_release_edges = tuple(edge for edge in restricted["terminal_edges"]
                                     if edge[0][0] == F12)
    active_changes = tuple(
        edge for edge in unrestricted["edges"] if edge[0][0] != edge[2][0]
    )
    acquisitions = tuple(
        edge for edge in unrestricted["edges"]
        if not set(edge[2][0]) <= set(edge[0][0])
    )
    pair_active_changes = tuple(edge for edge in active_changes if edge[0][0] == PAIR)
    reset12 = FACES["F12"]["reset"]
    current_carrier_records = tuple(
        (
            blocker,
            prior.common_decorated(responses, PAIR, blocker),
            tuple(
                edge for edge in release_edges if edge[2][1] == blocker
            ),
            reset12 == blocker,
        )
        for blocker in reachable_blockers
    )
    require(
        len(release_edges) == 69
        and len(restricted_release_edges) == 7
        and all(edge[2] == ((), reset12) for edge in release_edges)
        and len(active_changes) == 69
        and not acquisitions
        and not pair_active_changes
        and all(not record[1] and not record[2] and not record[3]
                for record in current_carrier_records),
        "current carrier unexpectedly contains a rebind",
    )

    pair_endpoints = pair_endpoint_census(unrestricted, triple_states)
    restricted_pair_endpoints = pair_endpoint_census(restricted, triple_states)
    f12_releases = f12_release_census(unrestricted, triple_states)
    restricted_f12_releases = f12_release_census(restricted, triple_states)

    def endpoint_summary(records):
        path_counts = tuple(
            count for _state, endpoints in records for _endpoint, count in endpoints
        )
        return (
            tuple(sorted(Counter(len(endpoints) for _state, endpoints in records).items())),
            tuple(
                (blocker, sum(
                    blocker in dict(endpoints) for _state, endpoints in records
                ))
                for blocker in reachable_blockers
            ),
            sum(len(endpoints) for _state, endpoints in records),
            (min(path_counts), max(path_counts)),
            stable_digest(records),
        )

    pair_endpoint_summary = endpoint_summary(pair_endpoints)
    restricted_pair_endpoint_summary = endpoint_summary(restricted_pair_endpoints)
    low_support_endpoint_law = tuple(
        (
            state,
            tuple(endpoint for endpoint, _count in endpoints),
            tuple(sorted(
                (3, 4, 5, 6, 7, 8)
                + tuple(root for root in (1, 2) if root in state)
            )),
        )
        for state, endpoints in pair_endpoints
    )
    restricted_endpoint_assignment = tuple(
        (state, tuple(endpoint for endpoint, _count in endpoints))
        for state, endpoints in restricted_pair_endpoints
    )
    require(
        all(len(actual) == 1 and actual[0] == predicted
            for _state, actual, predicted in low_support_endpoint_law)
        and restricted_endpoint_assignment
        == tuple((state, actual) for state, actual, _predicted
                 in low_support_endpoint_law),
        "row choices changed the Boolean low-support blocker law",
    )
    require(
        pair_endpoint_summary == (
            ((1, 94),),
            (
                ((1, 2, 3, 4, 5, 6, 7, 8), 47),
                ((1, 3, 4, 5, 6, 7, 8), 16),
                ((2, 3, 4, 5, 6, 7, 8), 24),
                ((3, 4, 5, 6, 7, 8), 7),
            ),
            94,
            (4_698, 159_909_458_497_526),
            "0889754e8a4697fd9addf830addb439805f9383839a29e84d35f14541e4e0dcb",
        )
        and restricted_pair_endpoint_summary == (
            ((1, 94),),
            (
                ((1, 2, 3, 4, 5, 6, 7, 8), 47),
                ((1, 3, 4, 5, 6, 7, 8), 16),
                ((2, 3, 4, 5, 6, 7, 8), 24),
                ((3, 4, 5, 6, 7, 8), 7),
            ),
            94,
            (6, 95_694),
            "bd779ecf5fb364237921a38e510af6447a548bad704cb95ffe210571134f9988",
        )
        and stable_digest(tuple(
            (state, actual) for state, actual, _predicted in low_support_endpoint_law
        )) == "2ac166ad386167b3787ed4e6dc792034416f0b8350034b8ee2965e6b428a1725",
        "all-start endpoint census drift",
    )

    release_time = {state: time for state, time, _count in f12_releases}
    schedule_records = []
    for state, endpoints in pair_endpoints:
        phi = sum(
            prior.prior.reset_distance(state, FACES[face]["reset"])
            for face in PAIR
        )
        require((phi - 4) % 2 == 0, ("bad pair parity", state, phi))
        pair_time = (phi - 4) // 2
        schedule_records.append(
            (state, release_time[state], pair_time, pair_time - release_time[state],
             tuple(endpoint for endpoint, _count in endpoints))
        )
    schedule_records = tuple(schedule_records)
    require(
        all(record[3] >= 0 for record in schedule_records)
        and stable_digest(f12_releases)
        == "9ae95fba446c10f1e2652281aef2d8a0e0b8b0ac1a1356790eec2a3bf571adc4"
        and stable_digest(restricted_f12_releases)
        == "828151f4bc0eaf57fde15199c0c87cf0d0f186c9df86c3e63e36bb4b8f3d5740"
        and stable_digest(schedule_records)
        == "1fbf8fadfe5e9878899697ea8ef762b80b0c8b2d79b749a22f4dc2a6c58e044f"
        and tuple(sorted(Counter(record[3] for record in schedule_records).items()))
        == ((0, 7), (1, 24), (2, 16), (3, 16), (4, 31)),
        "F12 release or schedule census drift",
    )

    # Smallest tested extension: a zero-potential, non-response copy at the
    # current pair blocker, followed only by inherited singleton edges.
    copy_records = []
    completion_records = []
    completion_1617_records = []
    for blocker in reachable_blockers:
        per_face = {}
        per_face_1617 = {}
        for face in PAIR:
            record = singleton_completion(responses, face, blocker)
            record_1617 = singleton_completion(
                responses, face, blocker, allowed_rows=ROWS_1617
            )
            completion_records.append(record)
            completion_1617_records.append(record_1617)
            per_face[face] = record
            per_face_1617[face] = record_1617
        phi_before = sum(
            prior.prior.reset_distance(blocker, FACES[face]["reset"])
            for face in PAIR
        )
        for released_origin_face in PAIR:
            old_origin_face = next(face for face in PAIR
                                   if face != released_origin_face)
            copy_records.append((
                blocker,
                ("released_F12_origin", released_origin_face, blocker,
                 "COPY_FROM_SHARED_ANCESTRY"),
                ("old_shared_origin", old_origin_face, blocker,
                 "RESTRICT_SHARED_ANCESTRY"),
                None,
                (phi_before, phi_before),
                per_face[released_origin_face][6],
                per_face[old_origin_face][6],
                not per_face[released_origin_face][5]
                and not per_face[old_origin_face][5],
                not per_face_1617[released_origin_face][5]
                and not per_face_1617[old_origin_face][5],
            ))
    copy_records = tuple(copy_records)
    completion_records = tuple(completion_records)
    completion_1617_records = tuple(completion_1617_records)
    require(
        len(copy_records) == 8
        and all(record[4] == (4, 4) and record[7] and record[8]
                for record in copy_records)
        and all(not record[5] for record in completion_records)
        and all(not record[5] for record in completion_1617_records)
        and all(record[2] >= 1 for record in completion_records)
        and all(record[2] >= 1 for record in completion_1617_records),
        "explicit unrestricted copy extension failed",
    )

    # Full-domain hostile: four-state copy permission is not a global repair.
    full_raw_states = set.union(*(set(FACES[face]["decisions"])
                                 for face in FACE_ORDER))
    full_coloured_nodes = tuple(sorted(
        (
            tuple(face for face in block if state in FACES[face]["decisions"]),
            state,
        )
        for state in full_raw_states
        for block in prior.COLOUR_BLOCKS
        if any(state in FACES[face]["decisions"] for face in block)
    ))
    full_unrestricted_blockers = tuple(sorted(
        (active, state) for active, state in full_coloured_nodes
        if not prior.common_decorated(responses, active, state)
    ))
    full_restricted_blockers = tuple(sorted(
        (active, state) for active, state in full_coloured_nodes
        if not prior.common_decorated(
            responses, active, state, allowed_rows=ROWS_1617
        )
    ))
    require(
        tuple(state for active, state in full_unrestricted_blockers if active == PAIR)
        == EXPECTED_FULL_UNIVERSAL_BLOCKERS
        and len(full_restricted_blockers) == 38,
        "full-domain blocker control drift",
    )
    broad_copy_1617 = tuple(
        (
            state,
            tuple(
                singleton_completion(
                    responses, face, state, allowed_rows=ROWS_1617
                )
                for face in PAIR
            ),
        )
        for state in EXPECTED_FULL_UNIVERSAL_BLOCKERS
    )
    broad_copy_1617_success = tuple(
        state for state, records in broad_copy_1617
        if all(not record[5] for record in records)
    )
    four_gate_1617_success = tuple(
        blocker for blocker in reachable_blockers
        if all(
            not singleton_completion(
                responses, face, blocker, allowed_rows=ROWS_1617
            )[5]
            for face in PAIR
        )
    )

    seam_macro_probe = macro_probe(responses)
    q12_face_typing = tuple(
        (
            face,
            reset12 in FACES[face]["decisions"],
            prior.prior.reset_distance(reset12, FACES[face]["reset"]),
        )
        for face in PAIR
    )
    require(
        four_gate_1617_success == EXPECTED_REACHABLE_BLOCKERS
        and broad_copy_1617_success == EXPECTED_FULL_UNIVERSAL_BLOCKERS
        and len(set(full_restricted_blockers) - set(full_unrestricted_blockers))
        == 31
        and q12_face_typing == (("F13", True, 6), ("F23", True, 6)),
        "full-domain or terminal-state hostile control drift",
    )
    require(
        seam_macro_probe == (
            (
                (1, 2), (1, 2, 3, 4, 5, 6, 7, 8), PAIR,
                (2, 5, 8, 9, 11, 12, 14, 16, 18, 22), (2, 2), False, True,
            ),
            (
                (1, 4), (1, 3, 4, 4, 5, 6, 7, 8), PAIR,
                (2, 5, 8, 9, 11, 12, 14, 16, 18, 22), (2, 2), False, True,
            ),
            ((3, 2), (2, 3, 3, 4, 5, 6, 7, 8), ("F13",), (), (2,), False, True),
            ((3, 4), (3, 3, 4, 4, 5, 6, 7, 8), ("F13",), (), (2,), False, True),
        ),
        "two-pole macro hostile drift",
    )

    print("FC3 THREE-FACE F12 RELEASE / COPY-REBIND GATE PARTIAL SCOUT")
    print("status=FINITE-EXACT PARTIAL;no_theorem_ID;no_promotion")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path.relative_to(ROOT)}")
    print("dependency_hashes=END")
    print(
        "inherited_carrier=node_(active_faces,state);edge_requires_one_common_"
        "literal_(row,one_pole_successor);target_active_drops_exact_reset_faces;"
        "faces_never_acquired;THM3278_origin_bit_is_a_static_cut_coordinate_"
        "not_an_occupancy_token_or_state_register"
    )
    print(
        f"F12_release_edges={len(release_edges)};"
        f"F12_release_edges_rows_16_17={len(restricted_release_edges)};"
        f"released_terminal=(empty,{reset12});active_change_edges={len(active_changes)};"
        f"pair_active_change_edges={len(pair_active_changes)};"
        f"face_acquisition_edges={len(acquisitions)};"
        f"Q12_is_separately_typed_on_pair={q12_face_typing}"
    )
    print(
        f"current_carrier_release_x_four_blocker_combinations="
        f"{len(release_edges) * len(reachable_blockers)};"
        f"lawful_current_carrier_rebinds=0;records={current_carrier_records};"
        "mechanism=F12_reset_deletes_only_its_own_obligation_and_ends_at_Q12_"
        "in_a_separate_lane;it_creates_neither_a_pair_state_copy_nor_a_face_"
        "acquisition_edge;although_Q12_is_facewise_typed_retagging_it_would_"
        "replace_the_current_pair_ancestry_and_is_not_an_inherited_edge"
    )
    print(
        f"all_94_F12_release_records_digest={stable_digest(f12_releases)};"
        f"all_94_F12_row_labelled_completion_count_range="
        f"{(min(record[2] for record in f12_releases), max(record[2] for record in f12_releases))};"
        f"rows_16_17_release_records_digest={stable_digest(restricted_f12_releases)};"
        f"rows_16_17_completion_count_range="
        f"{(min(record[2] for record in restricted_f12_releases), max(record[2] for record in restricted_f12_releases))}"
    )
    print("all_94_pair_endpoint_summary=" + repr(pair_endpoint_summary))
    print(
        "all_94_unique_blocker_law=blocker_C_union_support(start)_intersect_{1,2};"
        f"assignment_digest={stable_digest(tuple((state, actual) for state, actual, _predicted in low_support_endpoint_law))};"
        "mechanism=common_reset_directions_reduce_low_root_multiplicity_to_"
        "its_presence_bit_then_freeze_at_the_opposed_0/1_targets_while_roots_"
        "3_through_8_normalize_to_one"
    )
    print("rows_16_17_pair_endpoint_summary=" + repr(restricted_pair_endpoint_summary))
    print(
        f"schedule_records_digest={stable_digest(schedule_records)};"
        f"slack_histogram={tuple(sorted(Counter(record[3] for record in schedule_records).items()))};"
        f"zero_slack_locus={tuple(record[0] for record in schedule_records if record[3] == 0)};"
        "typing=no_shared_interleaving;comparison_is_only_between_two_separate_"
        "reset_directed_DAG_lengths_from_the_same_start"
    )
    print(
        "explicit_extension=COPY_AT_REACHABLE_BLOCKER_after_certified_F12_release;"
        "copy_source_is_current_shared_state_not_Q12;copy_row=None;copy_changes_"
        "one_shared_state_register_into_two_ancestry_typed_singleton_registers;"
        "copy_potential_drop=0;all_following_moves_are_inherited_row_labelled_"
        "one_pole_singleton_edges"
    )
    print("unrestricted_singleton_completion_records=" + repr(completion_records))
    print("rows_16_17_singleton_completion_records=" + repr(completion_1617_records))
    print("copy_assignment_records=" + repr(copy_records))
    print(
        f"copy_extension_all_94_start_blocker_pairs="
        f"{sum(len(endpoints) for _state, endpoints in pair_endpoints)};"
        "unrestricted_success=ALL;rows_16_17_success="
        + ("ALL" if all(record[8] for record in copy_records) else "NOT_ALL")
        + ";post_copy_total_singleton_steps=4"
    )
    print(
        f"full_domain_coloured_contexts={len(full_coloured_nodes)};"
        f"full_unrestricted_pair_blockers={len(full_unrestricted_blockers)};"
        f"full_rows_16_17_blockers={len(full_restricted_blockers)};"
        f"row_specific_rows_16_17_blockers="
        f"{len(set(full_restricted_blockers) - set(full_unrestricted_blockers))};"
        f"four_gate_rows_16_17_locally_copyable={four_gate_1617_success};"
        f"remaining_after_only_four_gate="
        f"{len(full_restricted_blockers) - len(four_gate_1617_success)};"
        f"all_seven_universal_rows_16_17_locally_copyable={broad_copy_1617_success};"
        f"remaining_even_after_broad_seven_gate="
        f"{len(full_restricted_blockers) - len(broad_copy_1617_success)};"
        "boundary=local_copy_continuation_does_not_make_(16,17)_a_global_"
        "three_face_selector"
    )
    print("minimal_blocker_C_copy_records=" + repr(tuple(
        record for record in copy_records if record[0] == (3, 4, 5, 6, 7, 8)
    )))
    print(
        "two_pole_macro_negative_control=" + repr(seam_macro_probe)
        + ";macro_is_not_an_extension_edge;it_has_no_common_first_one_pole_"
        "step_and_the_jointly_typed_targets_preserve_pair_potential"
    )
    print(
        f"source_ast=(assert_nodes={ASSERT_NODES},float_literals={FLOAT_LITERALS})"
    )
    print(
        "scope=exact_three_pinned_faces_and_94_common_starts;current_carrier_"
        "obstruction_plus_one_explicit_copy_sidecar_at_four_reachable_blockers;"
        "copy_is_not_inherited_and_no_shared_state_interleaving_FC3_SFC3_GMC_"
        "positivity_or_unbounded_controller_claim"
    )
    print("ALL FINITE-EXACT PARTIAL CHECKS PASSED")


if __name__ == "__main__":
    freeze_support()
    main()
