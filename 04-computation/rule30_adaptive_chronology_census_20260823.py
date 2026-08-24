#!/usr/bin/env python3
"""Finite exact census of conditional routed queries in the Rule 30 audit.

This scratch probe reuses the frozen ambient enumeration but computes a new
quotient: after the current portrait and the zero-chain/gap are known, how
many of the two off-ray chains are genuinely needed?  It is deliberately
scoped to active words of length at most nine at source depth eight.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import runpy
import sys


sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation" / "rule30_adaptive_routed_query_closure_20260823.py"
SOURCE_SHA256 = "ba00dbf59505ccb1373de0fde4a456e7d839ace78def29ffb27642a1f6b8b849"
EXPECTED_SEMANTIC_SHA256 = "41dfbb528276dbb10077c3d30c2a88fb668585a6c248b1d3f6908347d0d377e0"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def same_observation_implies_output(rows, coordinate: int) -> bool:
    fibres = defaultdict(set)
    for _, q1, q2, output in rows:
        fibres[(q1, q2)[coordinate]].add(output)
    return all(len(values) == 1 for values in fibres.values())


def collision(rows, coordinate: int):
    seen = {}
    for word, q1, q2, output in rows:
        observation = (q1, q2)[coordinate]
        if observation in seen and seen[observation][1] != output:
            old_word, old_output = seen[observation]
            return old_word, word, observation, old_output, output
        seen.setdefault(observation, (word, output))
    raise RuntimeError(("missing collision", coordinate))


def main() -> None:
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256,
            "source hash")
    namespace = runpy.run_path(str(SOURCE))
    module = namespace["load_dependency"]()
    records, active_word_count, distinct_permutations = namespace[
        "enumerate_ambient_words"
    ](module)

    classes = defaultdict(list)
    exact_count = 0
    for word, record in records:
        if record[0] != "exact":
            continue
        exact_count += 1
        _, gap, current, chains, output = record
        zero_key = (current, gap, chains[0])
        classes[zero_key].append((word, chains[1], chains[2], output))

    census = defaultdict(int)
    row_census = defaultdict(int)
    ray_one_controls = []
    ray_two_controls = []
    hard_controls = []
    maximum_class_size = 0
    maximum_output_count = 0
    full_pair_collisions = 0
    for zero_key, rows in classes.items():
        maximum_class_size = max(maximum_class_size, len(rows))
        outputs = {row[3] for row in rows}
        maximum_output_count = max(maximum_output_count, len(outputs))
        if len(outputs) == 1:
            label = "zero_queries"
        else:
            q1 = same_observation_implies_output(rows, 0)
            q2 = same_observation_implies_output(rows, 1)
            if q1 and q2:
                label = "either_one_query"
            elif q1:
                label = "ray_one_query"
                ray_one_controls.append(
                    (zero_key, len(rows), len(outputs), collision(rows, 1))
                )
            elif q2:
                label = "ray_two_query"
                ray_two_controls.append(
                    (zero_key, len(rows), len(outputs), collision(rows, 0))
                )
            else:
                label = "both_queries"
                hard_controls.append(
                    (
                        zero_key,
                        len(rows),
                        len(outputs),
                        collision(rows, 0),
                        collision(rows, 1),
                    )
                )

        pair_fibres = defaultdict(set)
        for _, q1_value, q2_value, output in rows:
            pair_fibres[(q1_value, q2_value)].add(output)
        full_pair_collisions += sum(len(values) > 1 for values in pair_fibres.values())
        census[label] += 1
        row_census[label] += len(rows)

    require(full_pair_collisions == 0, "two routed chains determine output")
    require(not hard_controls, "one conditional off-ray query always suffices")
    require(ray_one_controls and ray_two_controls,
            "both selector directions have positive controls")

    def shortest_control(controls):
        return min(
            controls,
            key=lambda item: (
                max(len(item[3][0]), len(item[3][1])),
                len(item[3][0]) + len(item[3][1]),
                repr(item),
            ),
        )

    shortest_ray_one = shortest_control(ray_one_controls)
    shortest_ray_two = shortest_control(ray_two_controls)

    ordered_labels = (
        "zero_queries",
        "either_one_query",
        "ray_one_query",
        "ray_two_query",
        "both_queries",
    )
    semantic = (
        SOURCE_SHA256,
        active_word_count,
        distinct_permutations,
        exact_count,
        len(classes),
        tuple((label, census[label], row_census[label]) for label in ordered_labels),
        maximum_class_size,
        maximum_output_count,
        len(hard_controls),
        shortest_ray_one,
        shortest_ray_two,
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "frozen semantic transcript")

    print("RULE30_ADAPTIVE_CHRONOLOGY_CENSUS_20260823")
    print("status=FINITE-EXACT;active_words_length_at_most_9;no_uniform_claim")
    print(f"source_sha256={SOURCE_SHA256}")
    print(
        f"active_words={active_word_count};distinct_permutations={distinct_permutations};"
        f"exact={exact_count};zero_observation_classes={len(classes)}"
    )
    print(
        "query_census="
        + repr(tuple((label, census[label], row_census[label]) for label in ordered_labels)).replace(" ", "")
    )
    print(
        f"maximum_class_size={maximum_class_size};"
        f"maximum_output_count={maximum_output_count};"
        f"hard_two_query_classes={len(hard_controls)}"
    )
    print("ray_one_needed_control=" + repr(shortest_ray_one).replace(" ", ""))
    print("ray_two_needed_control=" + repr(shortest_ray_two).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "interpretation=the_zero_observation_selects_at_most_one_off_ray_chain;"
        "neither_fixed_off_ray_chain_is_universal;two_chains_are_never_needed_"
        "in_this_finite_ambient"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
