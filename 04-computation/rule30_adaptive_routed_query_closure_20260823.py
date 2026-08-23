#!/usr/bin/env python3
"""Exact hostile probe for a bounded Rule 30 adaptive-query compiler.

This is deliberately a scratch-only finite experiment.  It studies the
literal point-evaluation compiler for the first-return operation downstream
of THM-3511 and the 2026-08-21 adaptive-observer signal.  It does not claim a
finite signalizer graph, a uniform gap bound, or any Rule 30 prize.
"""

from __future__ import annotations

import collections
import hashlib
import importlib.util
import itertools
import json
import sys
from pathlib import Path


sys.stdout.reconfigure(newline="\n")


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / "04-computation" / "rule30_orbit_signalizer_gap_thm3511.py"
DEPENDENCY_SHA256 = "2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445"

TARGET_DEPTH = 4
GAP_BUDGET = 4
SOURCE_DEPTH = TARGET_DEPTH + GAP_BUDGET
TARGET_RAYS = (0, 1, 2)
MAX_WORD_LENGTH = 9
EXPECTED_SEMANTIC_SHA256 = "0608bc5298aa7855dabb4784dd56df3c2df6cba06a7bfe0e2a60b754f1a3ca63"


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


def text(module, word: bytes | bytearray) -> str:
    return "".join(module.LETTERS[letter] for letter in word)


def active_orbit(module, start: int, depth: int) -> tuple[set[int], int, dict[int, str]]:
    """Orbit of one ray under words with odd root activity.

    The BFS state retains root activity because the desired operation is
    defined on active signalizers.  Representatives are shortest, then
    lexicographically first in the generator order A,B,C.
    """
    actions = module.generator_actions(depth)
    queue = collections.deque(((start, 0),))
    distance = {(start, 0): 0}
    representative = {(start, 0): ""}
    while queue:
        value, activity = queue.popleft()
        prefix = representative[(value, activity)]
        for letter in (module.A, module.B, module.C):
            state = (actions[letter][value], activity ^ (letter != module.A))
            if state in distance:
                continue
            distance[state] = distance[(value, activity)] + 1
            representative[state] = prefix + module.LETTERS[letter]
            queue.append(state)
    active_values = {value for value, activity in distance if activity == 1}
    active_representatives = {
        value: representative[(value, 1)] for value in active_values
    }
    return active_values, max(distance.values()), active_representatives


def first_stage_inputs(budget: int) -> tuple[int, ...]:
    """Inputs needed before the routed second query is known.

    Ray zero is needed for the gap.  If the gap is d, target rays one and two
    lift below 0^d to 2^d and 2^(d+1).  Taking all 1<=d<=B leaves zero and the
    powers 2^1,...,2^(B+1).
    """
    return (0,) + tuple(1 << exponent for exponent in range(1, budget + 2))


def operation_record(module, permutation: tuple[int, ...]) -> tuple[object, ...]:
    square = tuple(permutation[permutation[value]] for value in range(len(permutation)))
    target_mask = (1 << TARGET_DEPTH) - 1
    current_bank = tuple(permutation[ray] & target_mask for ray in TARGET_RAYS)
    visible = square[0] & ((1 << (GAP_BUDGET + 1)) - 1)
    zero_chain = (permutation[0], square[0])
    if visible == 0:
        return ("overflow", current_bank, zero_chain)
    gap = module.nu2(visible)
    require(1 <= gap <= GAP_BUDGET, ("visible gap", gap))
    chains = []
    outputs = []
    for ray in TARGET_RAYS:
        lifted = ray << gap
        routed = permutation[lifted]
        twice = permutation[routed]
        require(twice & ((1 << gap) - 1) == 0,
                ("fixed prefix", gap, ray, twice))
        chains.append((lifted, routed, twice))
        outputs.append((twice >> gap) & target_mask)
    require(chains[0][1:] == zero_chain, "zero chain agreement")
    return ("exact", gap, current_bank, tuple(chains), tuple(outputs))


def direct_target(module, word: bytearray, gap: int) -> tuple[int, ...]:
    successor, direct_gap = module.renormalize(word)
    require(direct_gap == gap, ("direct gap", text(module, word), gap, direct_gap))
    actions = module.generator_actions(TARGET_DEPTH)
    portrait = module.permutation(successor, actions, {})
    return tuple(portrait[ray] for ray in TARGET_RAYS)


def find_hostile(
    records: tuple[tuple[str, tuple[object, ...]], ...],
    kept_off_rays: tuple[int, ...],
    changed_ray: int,
) -> tuple[object, ...]:
    """Find a collision after retaining only the named off-ray chains."""
    buckets: dict[tuple[object, ...], tuple[str, tuple[object, ...]]] = {}
    for word_text, record in records:
        if record[0] != "exact":
            continue
        _, gap, current_bank, chains, outputs = record
        # chains[0] is the zero-ray chain.  The input coordinate is determined
        # by gap, but retaining it makes the observation contract explicit.
        key = (current_bank, gap, chains[0]) + tuple(
            chains[ray] for ray in kept_off_rays
        )
        previous = buckets.get(key)
        if previous is None:
            buckets[key] = (word_text, record)
            continue
        old_text, old_record = previous
        old_outputs = old_record[4]
        if old_outputs[changed_ray] != outputs[changed_ray]:
            return (
                old_text,
                word_text,
                gap,
                key,
                old_outputs,
                outputs,
            )
    raise RuntimeError(
        f"no hostile through length {MAX_WORD_LENGTH}: "
        f"kept={kept_off_rays}, changed={changed_ray}"
    )


def enumerate_ambient_words(module) -> tuple[
    tuple[tuple[str, tuple[object, ...]], ...], int, int
]:
    actions = module.generator_actions(SOURCE_DEPTH)
    block_cache: dict[bytes, tuple[int, ...]] = {}
    records = []
    active_word_count = 0
    distinct_permutations: dict[tuple[int, ...], str] = {}

    for length in range(1, MAX_WORD_LENGTH + 1):
        for letters in itertools.product((module.A, module.B, module.C), repeat=length):
            word = bytearray(letters)
            if module.activity(word) != 1:
                continue
            active_word_count += 1
            permutation = module.permutation(word, actions, block_cache)
            word_text = text(module, word)
            if permutation in distinct_permutations:
                continue
            distinct_permutations[permutation] = word_text
            record = operation_record(module, permutation)
            if record[0] == "exact":
                require(direct_target(module, word, record[1]) == record[4],
                        ("adaptive/direct target", word_text))
            records.append((word_text, record))
    return tuple(records), active_word_count, len(distinct_permutations)


def selected_bank_hostile(module) -> tuple[object, ...]:
    actions = module.generator_actions(TARGET_DEPTH)
    result = []
    for word_text in ("CA", "CCAC"):
        word = bytearray(module.LETTERS.index(letter) for letter in word_text)
        portrait = module.permutation(word, actions, {})
        successor, gap = module.renormalize(word)
        successor_portrait = module.permutation(successor, actions, {})
        result.append(
            (
                word_text,
                tuple(portrait[ray] for ray in TARGET_RAYS),
                gap,
                tuple(successor_portrait[ray] for ray in TARGET_RAYS),
            )
        )
    require(result[0][1] == result[1][1] == (7, 10, 1),
            "inherited selected-bank collision")
    require((result[0][2], result[1][2]) == (3, 10),
            "inherited gap split")
    return tuple(result)


def main() -> None:
    module = load_dependency()

    # Static two-stage point-query closure.  At each tested budget, every
    # required even first-stage input can be routed to every odd ray by an
    # active word.  Thus a nonadaptive literal composer must retain the whole
    # odd half of the source portrait, in addition to its first-stage inputs.
    closure_rows = []
    for budget in range(1, 7):
        depth = TARGET_DEPTH + budget
        inputs = first_stage_inputs(budget)
        odd_values = set(range(1, 1 << depth, 2))
        orbit_diameters = []
        witnesses = []
        for start in inputs:
            values, diameter, representatives = active_orbit(module, start, depth)
            require(values == odd_values,
                    ("active orbit is odd half", budget, start, len(values)))
            orbit_diameters.append(diameter)
            witnesses.append((start, representatives[max(odd_values)]))
        closure = set(inputs) | odd_values
        require(len(closure) == (1 << (depth - 1)) + budget + 2,
                ("closure size", budget))
        closure_rows.append(
            (
                budget,
                depth,
                inputs,
                len(odd_values),
                len(closure),
                max(orbit_diameters),
                tuple(witnesses),
            )
        )

    records, active_word_count, distinct_permutations = enumerate_ambient_words(module)
    exact_count = sum(record[0] == "exact" for _, record in records)
    overflow_count = len(records) - exact_count
    require(exact_count > 0 and overflow_count > 0, "positive and overflow controls")

    zero_only_hostile = find_hostile(records, (), 1)
    omit_ray_one_hostile = find_hostile(records, (2,), 1)
    omit_ray_two_hostile = find_hostile(records, (1,), 2)
    inherited_hostile = selected_bank_hostile(module)

    semantic = (
        DEPENDENCY_SHA256,
        TARGET_DEPTH,
        GAP_BUDGET,
        SOURCE_DEPTH,
        TARGET_RAYS,
        MAX_WORD_LENGTH,
        tuple(closure_rows),
        active_word_count,
        distinct_permutations,
        exact_count,
        overflow_count,
        zero_only_hostile,
        omit_ray_one_hostile,
        omit_ray_two_hostile,
        inherited_hostile,
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "frozen semantic transcript")

    print("RULE30_ADAPTIVE_QUERY_CLOSURE_20260823")
    print("status=FINITE-EXACT;literal_point-query_model;no_rule30_prize")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print(
        f"target_depth={TARGET_DEPTH};gap_budget={GAP_BUDGET};"
        f"source_depth={SOURCE_DEPTH};target_rays={TARGET_RAYS}"
    )
    print("static_closure_rows=" + repr(tuple(closure_rows)).replace(" ", ""))
    print(
        f"ambient_words=active_length_1_to_{MAX_WORD_LENGTH}:{active_word_count};"
        f"distinct_depth_{SOURCE_DEPTH}_permutations={distinct_permutations};"
        f"exact={exact_count};overflow={overflow_count}"
    )
    print("zero_only_hostile=" + repr(zero_only_hostile).replace(" ", ""))
    print("omit_ray_one_hostile=" + repr(omit_ray_one_hostile).replace(" ", ""))
    print("omit_ray_two_hostile=" + repr(omit_ray_two_hostile).replace(" ", ""))
    print("inherited_bank_hostile=" + repr(inherited_hostile).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "interpretation=static_literal_composition_saturates_odd_half;"
        "adaptive_zero_chain_plus_two_off_ray_routed_chains_is_exact;"
        "each_off_ray_chain_has_a_finite_hostile_when_omitted"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
