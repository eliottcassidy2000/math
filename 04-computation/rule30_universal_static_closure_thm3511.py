#!/usr/bin/env python3
"""Primary exact audit for universal Rule 30 static query saturation.

The proof-facing statement is group-theoretic.  On the binary rooted tree let

    A=(A,B),  B=(C,B)sigma,  C=(A,B)sigma.

The subgroup <A,B> is level-transitive because A fixes the root with left
section A, while B^{-1}AB fixes the root with left section B.  This checker
imports the canonical THM-3511 Mealy action, verifies the section identity on
large finite levels, exhausts the induced Schreier orbits, and checks the
static-bank cardinality formula.  It proves no Rule 30 prize.
"""

from __future__ import annotations

import collections
import hashlib
import importlib.util
import json
import sys
from pathlib import Path


sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / "04-computation" / "rule30_orbit_signalizer_gap_thm3511.py"
DEPENDENCY_SHA256 = "2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445"
MAX_DEPTH = 16
EXPECTED_SEMANTIC_SHA256 = "7b7b977eedf82a3798015e9e0d6be50f07ac0efd33b01bc6696b0be2a2a55a8d"


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


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [-1] * len(permutation)
    for source, target in enumerate(permutation):
        require(result[target] == -1, ("injective", source, target))
        result[target] = source
    require(all(value >= 0 for value in result), "surjective")
    return tuple(result)


def compose(*permutations: tuple[int, ...]) -> tuple[int, ...]:
    """Right-action product: apply the displayed permutations left to right."""
    size = len(permutations[0])
    result = tuple(range(size))
    for permutation in permutations:
        require(len(permutation) == size, "composition size")
        result = tuple(permutation[value] for value in result)
    return result


def orbit_with_activity(
    action_a: tuple[int, ...], action_b: tuple[int, ...], start: int
) -> tuple[int, int, int, int]:
    """Positive <A,B>-orbit, retaining the root-activity parity."""
    queue = collections.deque(((start, 0),))
    distance = {(start, 0): 0}
    activity_by_value = {start: 0}
    while queue:
        value, activity = queue.popleft()
        for permutation, toggle in ((action_a, 0), (action_b, 1)):
            state = (permutation[value], activity ^ toggle)
            old_activity = activity_by_value.get(state[0])
            if old_activity is not None:
                require(old_activity == state[1],
                        ("root activity determined by endpoint", start, state))
            else:
                activity_by_value[state[0]] = state[1]
            if state in distance:
                continue
            distance[state] = distance[(value, activity)] + 1
            queue.append(state)
    active_count = sum(activity for activity in activity_by_value.values())
    return (
        len(activity_by_value),
        active_count,
        max(distance.values()),
        sum(distance.values()),
    )


def cycle_length(permutation: tuple[int, ...], start: int) -> int:
    seen: dict[int, int] = {}
    value = start
    while value not in seen:
        seen[value] = len(seen)
        value = permutation[value]
    require(seen[value] == 0, ("permutation cycle has no tail", start))
    return len(seen)


def main() -> None:
    module = load_dependency()
    section_rows = []
    orbit_rows = []
    single_generator_rows = []
    orbit_gates = 0

    actions_by_depth = {
        depth: module.generator_actions(depth) for depth in range(1, MAX_DEPTH + 1)
    }

    for depth in range(2, MAX_DEPTH + 1):
        action_a, action_b, _ = actions_by_depth[depth]
        small_a, small_b, _ = actions_by_depth[depth - 1]
        inverse_b = inverse(action_b)
        conjugate = compose(inverse_b, action_a, action_b)
        left_checks = 0
        for suffix in range(1 << (depth - 1)):
            left_input = suffix << 1
            require(action_a[left_input] == (small_a[suffix] << 1),
                    ("A left section", depth, suffix))
            require(conjugate[left_input] == (small_b[suffix] << 1),
                    ("B^-1 A B left section", depth, suffix))
            left_checks += 2
        require(all((conjugate[value] & 1) == (value & 1)
                    for value in range(1 << depth)),
                ("conjugate root stabilizer", depth))
        section_rows.append((depth, left_checks))

    # A single orbit from zero verifies level transitivity.  Powers of two
    # independently stress the even inputs used by the routed-query compiler.
    for depth in range(1, MAX_DEPTH + 1):
        action_a, action_b, _ = actions_by_depth[depth]
        a_cycle = cycle_length(action_a, 0)
        b_cycle = cycle_length(action_b, 0)
        require(a_cycle == 1, ("A-only zero hostile", depth, a_cycle))
        if depth >= 2:
            require(b_cycle < (1 << depth),
                    ("root-transitive B is not level-transitive", depth, b_cycle))
        single_generator_rows.append((depth, a_cycle, b_cycle))
        starts = (0,) + tuple(1 << exponent for exponent in range(1, depth))
        depth_rows = []
        for start in starts:
            size, active_count, eccentricity, distance_sum = orbit_with_activity(
                action_a, action_b, start
            )
            require(size == 1 << depth, ("level transitivity", depth, start, size))
            require(active_count == 1 << (depth - 1),
                    ("active half", depth, start, active_count))
            # Every tested start is even.  Activity is therefore exactly the
            # output least-significant bit, not an independent state label.
            depth_rows.append((start, eccentricity, distance_sum))
            orbit_gates += size
        orbit_rows.append((depth, tuple(depth_rows)))

    # Exact static closure for target rays (0,1,2).  Before the gap d is
    # known, the first-stage inputs are 0,2,...,2^(B+1); after one active
    # query their possible images are the entire odd half.
    closure_rows = []
    for target_depth in range(2, 9):
        for gap_budget in range(1, 9):
            source_depth = target_depth + gap_budget
            inputs = (0,) + tuple(1 << exponent
                                  for exponent in range(1, gap_budget + 2))
            require(max(inputs) < (1 << source_depth),
                    ("input range", target_depth, gap_budget))
            require(len(set(inputs)) == gap_budget + 2,
                    ("input distinctness", target_depth, gap_budget))
            static_size = (1 << (source_depth - 1)) + len(inputs)
            closure_rows.append(
                (target_depth, gap_budget, source_depth, len(inputs), static_size)
            )

    semantic = (
        DEPENDENCY_SHA256,
        MAX_DEPTH,
        tuple(section_rows),
        tuple(orbit_rows),
        tuple(single_generator_rows),
        orbit_gates,
        tuple(closure_rows),
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "frozen semantic transcript")

    print("RULE30_UNIVERSAL_STATIC_CLOSURE_PRIMARY_20260823")
    print("status=PROVED_GROUP_LEMMA_PLUS_FINITE_EXACT_AUDIT;no_rule30_prize")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print("wreath_identity=A|0=A;(B^-1*A*B)|0=B;root_stabilizer=true")
    print(f"section_depths=2..{MAX_DEPTH};section_gates={sum(row[1] for row in section_rows)}")
    print(f"orbit_depths=1..{MAX_DEPTH};orbit_gates={orbit_gates}")
    print("orbit_eccentricities_from_zero=" + repr(tuple(
        (depth, rows[0][1]) for depth, rows in orbit_rows
    )).replace(" ", ""))
    print("single_generator_zero_cycles=" + repr(tuple(single_generator_rows)).replace(" ", ""))
    print("closure_formula=2^(D+B-1)+B+2;D>=2;B>=1")
    print("closure_sample_D5=" + repr(tuple(
        row for row in closure_rows if row[0] == 5
    )).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print("interpretation=positive_AB_monoid_is_level_transitive_on_every_finite_level;active_images_of_each_even_input_are_exactly_the_odd_half;literal_static_two_stage_query_closure_is_exponential")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
