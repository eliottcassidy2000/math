#!/usr/bin/env python3
"""Cheap profile-level probe for the two-current-digit period-three ranks.

This script deliberately does not repeat the endpoint event sweep, character
bank, or relation inversion.  It compares the candidate source profiles with
the independently reconstructed two-digit profiles entrywise, applies the
canonical owner-visible chamber/state mask independently, and asks for the row
ranks of the 13 older-digit sections on the six pointed owner cells.
"""

from __future__ import annotations

from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

TWO = importlib.import_module(
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_probe_20260816"
)
AUDIT = importlib.import_module(
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_independent_audit_20260816"
)

P = 13
PRIME = TWO.PRIME
POINTED = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
EXPECTED_PERIOD3 = (4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4)
EXPECTED_PROFILE_RANKS = (6, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 6)
EXPECTED_SUPPORT_RANKS = (3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3)
EXPECTED_CLASS_UNION_RANKS = (3, 3, 6, 10)
EXPECTED_CLASS_PAIR_UNION_RANKS = (
    (3, 3, 6, 10),
    (3, 3, 6, 10),
    (6, 6, 6, 12),
    (10, 10, 12, 10),
)
EXPECTED_REFLECTION_RECORDS = (
    (3, 3, 2, 1),
    (6, 6, 3, 3),
    (10, 10, 5, 5),
    (12, 12, 6, 6),
)
EXPECTED_SLICE_REFLECTION_RECORDS = (
    (4, 4, 2, 2),
    (5, 5, 3, 2),
)
EXPECTED_FLAG_LEDGER = (
    (3, 6, 10, 12),
    (4, 4, 4),
    (5, 5, 5),
    (7, 7),
    (10, 12, 6),
    (9, 9),
    (11, 10, 11),
)
EXPECTED_STATE_PROJECTION_RANKS = (
    (1, 2, 1, 2),
    (1, 3, 1, 3),
    (3, 3, 3, 3),
    (3, 4, 3, 4),
)
TRANSLATION3_ORBIT = (0, 3, 6, 9, 12, 2, 5, 8, 11, 1, 4, 7, 10)
EXPECTED_PROFILE_DIGESTS = (
    "f7d1ac1bba79a25b232dd4f0539b9a236f047b5224e49050a44c70f4d2544b68",
    "6fd60baa0f82c5f234164b2869cb369468fcbcf128e995603c49e5c07ed7ea68",
)
EXPECTED_AUDITED_SOURCE_SHA256 = (
    "dd2f48375fc38b419babdd0ed13365fb56918829a663270c4b23d56635d6e097"
)
EXPECTED_SEMANTIC_SHA256 = (
    "0d5798935986c13b1aa70fc1db7abac162994dabaf7dd3ab5f323ffeb1d1d63e"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank(rows) -> int:
    return TWO.C.rank_mod(tuple(tuple(value % PRIME for value in row) for row in rows))


def profile_rows(cells, r0: int, support_only: bool):
    rows = []
    for r1 in range(P):
        row = []
        for state, _u_values, _v_values, _mask, address_values in cells:
            for point_state, root in POINTED:
                value = address_values[root][r0][r1] if state == point_state else 0
                row.append(int(value != 0) if support_only else value)
        rows.append(tuple(row))
    return tuple(rows)


def reflect_profile_row(row, cell_count: int, point_reflection):
    return tuple(
        row[(cell_count - 1 - cell) * len(POINTED) + point_reflection[point]]
        for cell in range(cell_count)
        for point in range(len(POINTED))
    )


def reflection_record(rows, cell_count: int, point_reflection):
    reflected = tuple(
        reflect_profile_row(row, cell_count, point_reflection) for row in rows
    )
    plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    return rank(rows), rank(tuple(rows) + reflected), rank(plus), rank(minus)


def state_projection_ranks(rows, cell_count: int):
    answer = []
    for state in range(4):
        columns = tuple(
            cell * len(POINTED) + point
            for cell in range(cell_count)
            for point, (point_state, _root) in enumerate(POINTED)
            if point_state == state
        )
        answer.append(rank(tuple(tuple(row[column] for column in columns)
                                 for row in rows)))
    return tuple(answer)


def independently_masked_cells():
    values = AUDIT.source_context()
    profiles, boundaries, cells = values[:3]
    source_sha = values[-2]
    require(source_sha == EXPECTED_AUDITED_SOURCE_SHA256,
            ("independent source digest", source_sha))
    masked = []
    typed_count = 0
    for (left, right), cell in zip(zip(boundaries, boundaries[1:]), cells):
        state = cell[0]
        chamber = AUDIT.ONE.R.chamber(left, right)
        owner_visible = (
            (chamber == "left" and state in (0, 1))
            or (chamber == "right" and state in (2, 3))
        )
        state = state if owner_visible else None
        typed_count += int(state is not None)
        masked.append((state,) + cell[1:])
    require((len(masked), typed_count) == (33, 16),
            ("independent visibility mask", len(masked), typed_count))
    return profiles, boundaries, tuple(masked), source_sha, typed_count


def main() -> None:
    profiles, boundaries, candidate_cells, source_record = TWO.two_digit_context()
    (audit_profiles, audit_boundaries, cells,
     audit_source_sha, typed_count) = independently_masked_cells()
    require(profiles == audit_profiles, "candidate/independent profile mismatch")
    require(boundaries == audit_boundaries,
            "candidate/independent boundary mismatch")
    require(candidate_cells == cells,
            "candidate/independent endpoint-masked cell mismatch")
    require(len(cells) == len(boundaries) - 1, "cell/boundary mismatch")

    weighted_rows = tuple(profile_rows(cells, r0, False) for r0 in range(P))
    support_rows = tuple(profile_rows(cells, r0, True) for r0 in range(P))
    weighted_ranks = tuple(rank(rows) for rows in weighted_rows)
    support_ranks = tuple(rank(rows) for rows in support_rows)
    global_weighted_rank = rank(tuple(row for block in weighted_rows for row in block))
    global_support_rank = rank(tuple(row for block in support_rows for row in block))

    point_lookup = {pair: index for index, pair in enumerate(POINTED)}
    point_reflection = tuple(
        point_lookup[(state ^ 2, P - 1 - root)] for state, root in POINTED
    )
    for r0 in range(P):
        for r1 in range(P):
            row = weighted_rows[r0][r1]
            mate = weighted_rows[P - 1 - r0][P - 1 - r1]
            reflected = reflect_profile_row(mate, len(cells), point_reflection)
            require(row == reflected, ("typed profile reflection", r0, r1))

    digit_classes = ((1, 4, 7, 10), (2, 5, 8, 11), (3, 6, 9), (0, 12))
    class_union_ranks = tuple(
        rank(tuple(row for r0 in digit_class for row in weighted_rows[r0]))
        for digit_class in digit_classes
    )
    class_pair_union_ranks = tuple(
        tuple(
            rank(tuple(
                row
                for class_index in sorted({left_class, right_class})
                for r0 in digit_classes[class_index]
                for row in weighted_rows[r0]
            ))
            for right_class in range(len(digit_classes))
        )
        for left_class in range(len(digit_classes))
    )
    class_rows = tuple(
        tuple(row for r0 in digit_class for row in weighted_rows[r0])
        for digit_class in digit_classes
    )
    interior_nonmultiple_rows = class_rows[0] + class_rows[1]
    interior_multiple_rows = class_rows[2]
    boundary_rows = class_rows[3]
    global_rows = tuple(row for block in weighted_rows for row in block)
    middle_slices = tuple(weighted_rows[r0] for r0 in (3, 6, 9))
    boundary_slices = tuple(weighted_rows[r0] for r0 in (0, 12))
    reflection_records = (
        reflection_record(interior_nonmultiple_rows, len(cells), point_reflection),
        reflection_record(interior_multiple_rows, len(cells), point_reflection),
        reflection_record(boundary_rows, len(cells), point_reflection),
        reflection_record(global_rows, len(cells), point_reflection),
    )
    slice_reflection_records = (
        reflection_record(weighted_rows[6], len(cells), point_reflection),
        reflection_record(weighted_rows[3] + weighted_rows[9],
                          len(cells), point_reflection),
    )
    flag_ledger = (
        (rank(interior_nonmultiple_rows), rank(interior_multiple_rows),
         rank(boundary_rows), rank(global_rows)),
        tuple(rank(interior_nonmultiple_rows + rows) for rows in middle_slices),
        tuple(rank(left + right) for left, right in
              ((middle_slices[0], middle_slices[1]),
               (middle_slices[0], middle_slices[2]),
               (middle_slices[1], middle_slices[2]))),
        tuple(rank(interior_nonmultiple_rows + rows) for rows in boundary_slices),
        (rank(interior_nonmultiple_rows + boundary_rows),
         rank(interior_multiple_rows + boundary_rows),
         rank(interior_nonmultiple_rows + interior_multiple_rows)),
        tuple(rank(interior_multiple_rows + rows) for rows in boundary_slices),
        tuple(rank(boundary_rows + rows) for rows in middle_slices),
    )
    state_projection_ledger = tuple(
        state_projection_ranks(rows, len(cells))
        for rows in (interior_nonmultiple_rows, interior_multiple_rows,
                     boundary_rows, global_rows)
    )
    profile_translation3_ranks = tuple(weighted_ranks[r0]
                                       for r0 in TRANSLATION3_ORBIT)
    endpoint_translation3_ranks = tuple(EXPECTED_PERIOD3[r0]
                                        for r0 in TRANSLATION3_ORBIT)

    require(weighted_ranks == EXPECTED_PROFILE_RANKS,
            ("weighted profile ranks", weighted_ranks))
    require(support_ranks == EXPECTED_SUPPORT_RANKS,
            ("support profile ranks", support_ranks))
    require((global_weighted_rank, global_support_rank) == (12, 3),
            ("global profile ranks", global_weighted_rank, global_support_rank))
    require(class_union_ranks == EXPECTED_CLASS_UNION_RANKS,
            ("class union ranks", class_union_ranks))
    require(class_pair_union_ranks == EXPECTED_CLASS_PAIR_UNION_RANKS,
            ("class pair union ranks", class_pair_union_ranks))
    require(reflection_records == EXPECTED_REFLECTION_RECORDS,
            ("reflection representations", reflection_records))
    require(slice_reflection_records == EXPECTED_SLICE_REFLECTION_RECORDS,
            ("slice reflection representations", slice_reflection_records))
    require(flag_ledger == EXPECTED_FLAG_LEDGER,
            ("profile flag ledger", flag_ledger))
    require(state_projection_ledger == EXPECTED_STATE_PROJECTION_RANKS,
            ("state projection ranks", state_projection_ledger))
    require(len(set(TRANSLATION3_ORBIT)) == P
            and profile_translation3_ranks[0] != profile_translation3_ranks[-1]
            and endpoint_translation3_ranks[0] != endpoint_translation3_ranks[-1],
            ("translation-by-three hostile", profile_translation3_ranks,
             endpoint_translation3_ranks))
    require((digest_json(weighted_rows), digest_json(support_rows))
            == EXPECTED_PROFILE_DIGESTS, "profile digest drift")

    records = (
        (audit_source_sha, len(cells), typed_count, digest_json(cells)),
        weighted_ranks,
        support_ranks,
        global_weighted_rank,
        global_support_rank,
        digit_classes,
        class_union_ranks,
        class_pair_union_ranks,
        reflection_records,
        slice_reflection_records,
        flag_ledger,
        state_projection_ledger,
        TRANSLATION3_ORBIT,
        profile_translation3_ranks,
        endpoint_translation3_ranks,
        point_reflection,
        digest_json(weighted_rows),
        digest_json(support_rows),
        source_record[-2:],
    )
    semantic = digest_json(records)
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 two-digit period-three pointed-profile probe ==")
    print(f"independent_source_gate=(sha256={audit_source_sha},cells={len(cells)},owner_visible={typed_count},masked_cell_sha256={records[0][3]}): PASS")
    print(f"pointed_states={POINTED}")
    print(f"profile_r1_conditional_ranks_by_r0={weighted_ranks}")
    print(f"support_conditional_ranks_by_r0={support_ranks}")
    print(f"global_address_profile_ranks_(weighted,support)=({global_weighted_rank},{global_support_rank})")
    print(f"digit_classes={digit_classes};class_union_ranks={class_union_ranks}")
    print(f"class_pair_union_ranks={class_pair_union_ranks}")
    print("reflection_record_order=(rank,rank_with_reflection,plus_rank,minus_rank)")
    print(f"reflection_records_(interior_nonmultiple,interior_multiple,boundary,global)={reflection_records}")
    print(f"slice_reflection_records_(r0_6,r0_3_union_r0_9)={slice_reflection_records}")
    print(f"flag_ledger={flag_ledger}")
    print(f"state_projection_ranks_(L,M,B,G)={state_projection_ledger}")
    print(f"translation_by_3_orbit={TRANSLATION3_ORBIT}")
    print(f"translation_by_3_ranks_(profile,endpoint)=({profile_translation3_ranks},{endpoint_translation3_ranks}): NONCONSTANT")
    print(f"typed_pointed_reflection={point_reflection}: PASS")
    print(f"digests_(weighted,support)={(digest_json(weighted_rows), digest_json(support_rows))}")
    print(f"parent_profile_digests={source_record[-2:]}")
    print(f"semantic_sha256={semantic}")
    print(f"endpoint_conditional_ranks={EXPECTED_PERIOD3};profile_to_endpoint_rank_drop={tuple(left-right for left,right in zip(weighted_ranks,EXPECTED_PERIOD3))}")
    print("endpoint_address_rank=4;profile_address_rank=12;endpoint_projection_drop=8")
    print("scope=owner-visible endpoint-admissible nested-window profiles only;no endpoint sweep,character bank,relation inversion,or physical current")


if __name__ == "__main__":
    main()
