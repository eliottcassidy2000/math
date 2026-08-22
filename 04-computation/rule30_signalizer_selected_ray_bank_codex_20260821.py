#!/usr/bin/env python3
"""Exact selected-ray observer bank for the first Rule 30 signalizers.

This is a finite scout downstream of proved THM-3511.  It asks how many
point evaluations of the depth-four portrait separate the 23 physical
signalizers frozen there, and then checks whether the smallest bank is a
congruence for the native first-return renormalization.  It is not a finite
signalizer graph and has no Rule 30 prize consequence.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
from pathlib import Path


DEPENDENCY = Path(__file__).with_name("rule30_orbit_signalizer_gap_thm3511.py")
DEPENDENCY_SHA256 = "2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def load_dependency():
    payload = DEPENDENCY.read_bytes()
    require(hashlib.sha256(payload).hexdigest() == DEPENDENCY_SHA256,
            "THM-3511 dependency hash")
    spec = importlib.util.spec_from_file_location("rule30_thm3511", DEPENDENCY)
    require(spec is not None and spec.loader is not None, "dependency loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def word(module, text: str) -> bytearray:
    table = {"A": module.A, "B": module.B, "C": module.C}
    return bytearray(table[letter] for letter in text)


def portrait(module, state: bytearray, depth: int) -> tuple[int, ...]:
    actions = module.generator_actions(depth)
    return module.permutation(state, actions, {})


def ray_key(permutation: tuple[int, ...], rays: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(permutation[ray] for ray in rays)


def renormalized_record(module, text: str, depth: int,
                        rays: tuple[int, ...]) -> tuple[object, ...]:
    state = word(module, text)
    require(module.activity(state) == 1, (text, "active"))
    current = portrait(module, state, depth)
    successor, gap = module.renormalize(state)
    following = portrait(module, successor, depth)
    return (
        text,
        ray_key(current, rays),
        gap,
        ray_key(following, rays),
        current,
        following,
    )


def recover_bounded_first_return(module, state: bytearray, target_depth: int,
                                 gap_budget: int) -> tuple[object, ...]:
    """Recover a first return from one sufficiently deep source portrait.

    The depth-(target_depth + gap_budget) action determines the square action
    on that finite tree.  If the first moved zero-ray bit has index at most
    gap_budget, restriction below the corresponding fixed zero prefix gives
    the complete target-depth portrait of the first-return section.  If no
    such bit is visible, the only sound conclusion is overflow.
    """
    require(target_depth >= 1 and gap_budget >= 1, "positive portrait budget")
    source_depth = target_depth + gap_budget
    source = portrait(module, state, source_depth)
    square = tuple(source[source[value]] for value in range(1 << source_depth))

    visible_zero_image = square[0] & ((1 << (gap_budget + 1)) - 1)
    if visible_zero_image == 0:
        return ("overflow", gap_budget, source_depth)

    gap = module.nu2(visible_zero_image)
    require(1 <= gap <= gap_budget, ("visible gap", gap, gap_budget))
    target_mask = (1 << target_depth) - 1
    recovered = []
    for value in range(1 << target_depth):
        lifted_input = value << gap
        lifted_output = square[lifted_input]
        require(lifted_output & ((1 << gap) - 1) == 0,
                ("fixed zero prefix", gap, value))
        recovered.append((lifted_output >> gap) & target_mask)
    return ("exact", gap, tuple(recovered), source_depth)


def main() -> None:
    module = load_dependency()
    physical_words, gaps, transcript_digest, direct_digest = module.signalizer_words()
    require(len(physical_words) == len(gaps) == 23, "physical universe")

    depth = 4
    actions = module.generator_actions(depth)
    block_cache: dict[bytes, tuple[int, ...]] = {}
    physical_portraits = tuple(
        module.permutation(state, actions, block_cache) for state in physical_words
    )
    require(len(set(physical_portraits)) == 23, "THM-3511 portrait control")

    separating: dict[int, tuple[tuple[int, ...], ...]] = {}
    for size in (1, 2, 3):
        banks = []
        for rays in itertools.combinations(range(1 << depth), size):
            keys = tuple(ray_key(item, rays) for item in physical_portraits)
            if len(set(keys)) == len(keys):
                banks.append(rays)
        separating[size] = tuple(banks)

    require(separating[1] == (), "no one-ray separator")
    require(separating[2] == (), "no two-ray separator")
    require(len(separating[3]) == 104, "three-ray separator count")
    require(separating[3][0] == (0, 1, 2), "lexicographic bank")
    rays = separating[3][0]

    physical_records = tuple(
        (index, gaps[index], ray_key(physical_portraits[index], rays))
        for index in range(23)
    )
    require(len({item[2] for item in physical_records}) == 23,
            "selected bank injective")

    # The selected bank separates the frozen physical sample but is not a
    # congruence for the first-return operation on ambient active words.
    ca = renormalized_record(module, "CA", depth, rays)
    ccac = renormalized_record(module, "CCAC", depth, rays)
    require(ca[1] == ccac[1] == (7, 10, 1), "selected-ray collision")
    require((ca[2], ccac[2]) == (3, 10), "selected-ray gap split")
    require(ca[3] != ccac[3], "selected-ray successor split")

    # Even the complete fixed-depth portrait is not an operation-congruence.
    full_hostiles = []
    for portrait_depth, left_text, right_text, expected_gaps in (
        (4, "CBC", "AAAACBC", (6, 8)),
        (5, "AACAAACC", "CAAACAAC", (6, 5)),
    ):
        left = renormalized_record(module, left_text, portrait_depth,
                                   tuple(range(1 << portrait_depth)))
        right = renormalized_record(module, right_text, portrait_depth,
                                    tuple(range(1 << portrait_depth)))
        require(left[4] == right[4], (portrait_depth, "full portrait collision"))
        require((left[2], right[2]) == expected_gaps,
                (portrait_depth, "full portrait gap split"))
        full_hostiles.append(
            (portrait_depth, left_text, right_text, expected_gaps,
             hashlib.sha256(bytes(left[4])).hexdigest())
        )

    # A fixed portrait is not a first-return congruence, but there is an exact
    # bounded-gap repair.  A depth-(D+B) source portrait recovers the gap and
    # depth-D successor portrait when d<=B, and otherwise reports overflow.
    # The smaller budgets below separate each collision by exact/overflow;
    # the larger budgets recover both successors exactly.
    adaptive_recoveries = []
    for portrait_depth, left_text, right_text, split_budget, full_budget in (
        (4, "CA", "CCAC", 3, 10),
        (4, "CBC", "AAAACBC", 6, 8),
        (5, "AACAAACC", "CAAACAAC", 5, 6),
    ):
        all_rays = tuple(range(1 << portrait_depth))
        left_state = word(module, left_text)
        right_state = word(module, right_text)
        left_direct = renormalized_record(
            module, left_text, portrait_depth, all_rays,
        )
        right_direct = renormalized_record(
            module, right_text, portrait_depth, all_rays,
        )

        left_split = recover_bounded_first_return(
            module, left_state, portrait_depth, split_budget,
        )
        right_split = recover_bounded_first_return(
            module, right_state, portrait_depth, split_budget,
        )
        split_statuses = (left_split[0], right_split[0])
        require(set(split_statuses) == {"exact", "overflow"},
                (left_text, right_text, "split statuses"))
        for recovered, direct in ((left_split, left_direct),
                                  (right_split, right_direct)):
            if recovered[0] == "exact":
                require((recovered[1], recovered[2]) == (direct[2], direct[5]),
                        (left_text, right_text, "split recovery"))
            else:
                require(direct[2] > split_budget,
                        (left_text, right_text, "sound overflow"))

        left_full = recover_bounded_first_return(
            module, left_state, portrait_depth, full_budget,
        )
        right_full = recover_bounded_first_return(
            module, right_state, portrait_depth, full_budget,
        )
        require(left_full[0] == right_full[0] == "exact",
                (left_text, right_text, "full recovery statuses"))
        require((left_full[1], left_full[2]) == (left_direct[2], left_direct[5]),
                (left_text, "full recovery"))
        require((right_full[1], right_full[2]) == (right_direct[2], right_direct[5]),
                (right_text, "full recovery"))
        adaptive_recoveries.append(
            (
                portrait_depth,
                left_text,
                right_text,
                split_budget,
                split_statuses,
                full_budget,
                (left_full[1], left_full[3]),
                (right_full[1], right_full[3]),
            )
        )

    semantic = (
        DEPENDENCY_SHA256,
        transcript_digest,
        direct_digest,
        tuple((size, separating[size]) for size in (1, 2, 3)),
        physical_records,
        ca[:4],
        ccac[:4],
        tuple(full_hostiles),
        tuple(adaptive_recoveries),
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_SIGNALIZER_SELECTED_RAY_BANK_20260821")
    print("status=PROVED_BOUNDED_GAP_RECOVERY+FINITE-EXACT_BANK;no_finite_graph;no_rule30_prize")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print(f"physical_universe=s_0..s_22;count={len(physical_records)};depth={depth}")
    print("separator_counts=" + repr(tuple((k, len(separating[k])) for k in (1, 2, 3))))
    print(f"minimum_bank_size=3;lexicographic_bank={rays};three_ray_banks=104")
    print("physical_records=(m,gap,ray_outputs)=" + repr(physical_records).replace(" ", ""))
    print("selected_bank_hostile=" + repr((ca[:4], ccac[:4])).replace(" ", ""))
    print("full_portrait_hostiles=" + repr(tuple(full_hostiles)).replace(" ", ""))
    print("adaptive_recoveries=" + repr(tuple(adaptive_recoveries)).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print("interpretation=frozen_separator_not_operation_closed;"
          "depth_D_plus_gap_budget_recovers_exact_transition_or_overflow")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
