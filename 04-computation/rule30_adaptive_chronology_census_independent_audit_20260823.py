#!/usr/bin/env python3
"""Clean-room finite audit of the Rule 30 conditional-query selector.

This companion does not import the primary chronology census or THM-3511's
executable.  It rebuilds the three Mealy generators directly from their
wreath recursion and enumerates the bounded word universe independently.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib
import itertools
import json
import sys


sys.stdout.reconfigure(newline="\n")

A, B, C = 0, 1, 2
LETTERS = "ABC"
SOURCE_DEPTH = 8
TARGET_DEPTH = 4
GAP_BUDGET = 4
MAX_WORD_LENGTH = 9
EXPECTED_SEMANTIC_SHA256 = "dfaa16be3db0361f192588a50eebc53fded33e3c8df5c9c708198fd0a8cd191e"


def gate(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"gate failed: {label}")


def generator_image(letter: int, value: int, depth: int) -> int:
    """Evaluate A=(A,B), B=(C,B)sigma, C=(A,B)sigma, low bit first."""
    state = letter
    output = 0
    for level in range(depth):
        bit = (value >> level) & 1
        if state == A:
            out_bit = bit
            state = A if bit == 0 else B
        elif state == B:
            out_bit = bit ^ 1
            state = C if bit == 0 else B
        else:
            out_bit = bit ^ 1
            state = A if bit == 0 else B
        output |= out_bit << level
    return output


def compose_word(word: tuple[int, ...], actions) -> tuple[int, ...]:
    permutation = tuple(range(1 << SOURCE_DEPTH))
    for letter in word:
        permutation = tuple(actions[letter][value] for value in permutation)
    return permutation


def valuation_two(value: int) -> int:
    gate(value != 0, "nonzero valuation input")
    return (value & -value).bit_length() - 1


def exact_record(permutation: tuple[int, ...]):
    square = tuple(permutation[permutation[value]] for value in range(1 << SOURCE_DEPTH))
    visible = square[0] & ((1 << (GAP_BUDGET + 1)) - 1)
    if visible == 0:
        return None
    gap = valuation_two(visible)
    gate(1 <= gap <= GAP_BUDGET, ("gap", gap))
    current = tuple(permutation[ray] & ((1 << TARGET_DEPTH) - 1) for ray in (0, 1, 2))
    chains = []
    outputs = []
    for ray in (0, 1, 2):
        lifted = ray << gap
        routed = permutation[lifted]
        twice = permutation[routed]
        gate(twice & ((1 << gap) - 1) == 0,
             ("fixed prefix", ray, gap, twice))
        chains.append((lifted, routed, twice))
        outputs.append((twice >> gap) & ((1 << TARGET_DEPTH) - 1))
    gate(chains[0][1:] == (permutation[0], square[0]), "zero chain")
    return current, gap, tuple(chains), tuple(outputs)


def coordinate_suffices(rows, coordinate: int) -> bool:
    fibres = defaultdict(set)
    for _, observations, output in rows:
        fibres[observations[coordinate]].add(output)
    return all(len(values) == 1 for values in fibres.values())


def first_collision(rows, coordinate: int):
    seen = {}
    for word, observations, output in rows:
        observation = observations[coordinate]
        old = seen.get(observation)
        if old is not None and old[1] != output:
            return old[0], word, observation, old[1], output
        seen.setdefault(observation, (word, output))
    raise RuntimeError(("collision absent", coordinate))


def main() -> None:
    actions = tuple(
        tuple(generator_image(letter, value, SOURCE_DEPTH)
              for value in range(1 << SOURCE_DEPTH))
        for letter in (A, B, C)
    )
    for action in actions:
        gate(len(set(action)) == 1 << SOURCE_DEPTH, "generator permutation")

    representatives = {}
    active_word_count = 0
    for length in range(1, MAX_WORD_LENGTH + 1):
        for word in itertools.product((A, B, C), repeat=length):
            if sum(letter != A for letter in word) % 2 == 0:
                continue
            active_word_count += 1
            permutation = compose_word(word, actions)
            representatives.setdefault(permutation, "".join(LETTERS[x] for x in word))

    classes = defaultdict(list)
    exact_count = 0
    overflow_count = 0
    for permutation, word in representatives.items():
        record = exact_record(permutation)
        if record is None:
            overflow_count += 1
            continue
        exact_count += 1
        current, gap, chains, output = record
        classes[(current, gap, chains[0])].append(
            (word, (chains[1], chains[2]), output)
        )

    labels = (
        "zero_queries",
        "either_one_query",
        "ray_one_query",
        "ray_two_query",
        "both_queries",
    )
    class_census = defaultdict(int)
    row_census = defaultdict(int)
    ray_one_controls = []
    ray_two_controls = []
    hard_controls = []
    maximum_class_size = 0
    maximum_output_count = 0
    pair_ambiguities = 0

    for zero_key, rows in classes.items():
        maximum_class_size = max(maximum_class_size, len(rows))
        outputs = {row[2] for row in rows}
        maximum_output_count = max(maximum_output_count, len(outputs))
        if len(outputs) == 1:
            label = "zero_queries"
        else:
            ray_one = coordinate_suffices(rows, 0)
            ray_two = coordinate_suffices(rows, 1)
            if ray_one and ray_two:
                label = "either_one_query"
            elif ray_one:
                label = "ray_one_query"
                ray_one_controls.append(
                    (zero_key, len(rows), len(outputs), first_collision(rows, 1))
                )
            elif ray_two:
                label = "ray_two_query"
                ray_two_controls.append(
                    (zero_key, len(rows), len(outputs), first_collision(rows, 0))
                )
            else:
                label = "both_queries"
                hard_controls.append(zero_key)

        pair_fibres = defaultdict(set)
        for _, observation_pair, output in rows:
            pair_fibres[observation_pair].add(output)
        pair_ambiguities += sum(len(values) > 1 for values in pair_fibres.values())
        class_census[label] += 1
        row_census[label] += len(rows)

    gate(active_word_count == 14762, "active word census")
    gate(len(representatives) == 14762, "distinct permutation census")
    gate((exact_count, overflow_count) == (13853, 909), "exact/overflow census")
    gate(len(classes) == 11927, "zero class census")
    gate(pair_ambiguities == 0, "two-chain determination")
    gate(not hard_controls, "one conditional chain suffices")
    gate(ray_one_controls and ray_two_controls, "fixed-chain hostiles")

    def shortest(controls):
        return min(
            controls,
            key=lambda item: (
                max(len(item[3][0]), len(item[3][1])),
                len(item[3][0]) + len(item[3][1]),
                repr(item),
            ),
        )

    ray_one_control = shortest(ray_one_controls)
    ray_two_control = shortest(ray_two_controls)
    census = tuple((label, class_census[label], row_census[label]) for label in labels)
    gate(
        census
        == (
            ("zero_queries", 10502, 10705),
            ("either_one_query", 1318, 2889),
            ("ray_one_query", 85, 211),
            ("ray_two_query", 22, 48),
            ("both_queries", 0, 0),
        ),
        "conditional query census",
    )
    gate((maximum_class_size, maximum_output_count) == (7, 7), "class maxima")

    semantic = (
        active_word_count,
        len(representatives),
        exact_count,
        overflow_count,
        len(classes),
        census,
        maximum_class_size,
        maximum_output_count,
        ray_one_control,
        ray_two_control,
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    gate(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic transcript")

    print("RULE30_ADAPTIVE_CHRONOLOGY_INDEPENDENT_AUDIT_20260823")
    print("status=FINITE-EXACT;clean_room_mealy_enumeration;no_uniform_claim")
    print(
        f"active_words={active_word_count};distinct_permutations={len(representatives)};"
        f"exact={exact_count};overflow={overflow_count};classes={len(classes)}"
    )
    print("query_census=" + repr(census).replace(" ", ""))
    print(
        f"maximum_class_size={maximum_class_size};"
        f"maximum_output_count={maximum_output_count};pair_ambiguities={pair_ambiguities}"
    )
    print("ray_one_needed_control=" + repr(ray_one_control).replace(" ", ""))
    print("ray_two_needed_control=" + repr(ray_two_control).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
