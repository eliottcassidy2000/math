#!/usr/bin/env python3
"""Validate exact minimum-seed search on scalar-closed states.

For a fixed apex gate A, let cl(P) be the least scalar-activation fixed
point above P.  This is a finite closure operator.  The directed graph

    X  -->  cl(X union {a}),       X=cl(X), a not in X,

with unit edge costs has distance from cl(empty) to A equal to the exact
minimum seed size.  It also returns an interleaved order: after paying for
one apex on the actual current prefix X, append the scalar cascade to the
next closed state.  This avoids treating future seed labels as already
excluded.

The probe reproduces the four committed exact minima and reports how many
closed states/edges replace raw subset enumeration.
"""

from __future__ import annotations

import importlib.util
from collections import deque
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "04-computation/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.py"
)


def load(path: Path):
    spec = importlib.util.spec_from_file_location("closed_state_minseed", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load(SOURCE)


def scalar_closure(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    start: int,
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    """Return the scalar closure mask and deterministic simultaneous rounds."""

    index = {speed: bit for bit, speed in enumerate(gate)}
    active = start
    rounds: list[tuple[int, ...]] = []
    while True:
        additions: list[int] = []
        for bit, apex in enumerate(gate):
            if active & (1 << bit):
                continue
            values: list[F] = []
            for value, speed in profiles[apex]["ranked_candidates"]:
                speed_bit = index.get(speed)
                if speed_bit is not None and active & (1 << speed_bit):
                    continue
                values.append(value)
                if len(values) == 5:
                    break
            if len(values) != 5:
                raise RuntimeError("activation candidate set too small")
            if profiles[apex]["m"] - sum(values, F(0)) > 0:
                additions.append(bit)
        if not additions:
            return active, tuple(rounds)
        rounds.append(tuple(gate[bit] for bit in additions))
        for bit in additions:
            active |= 1 << bit


def minimum_closed_state_path(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
) -> dict[str, object]:
    full = (1 << len(gate)) - 1
    start, start_rounds = scalar_closure(gate, profiles, 0)
    queue = deque([start])
    distance = {start: 0}
    previous: dict[int, tuple[int, int, tuple[tuple[int, ...], ...]]] = {}
    edges = 0
    cache: dict[int, tuple[int, tuple[tuple[int, ...], ...]]] = {
        start: (start, start_rounds)
    }
    found = start == full
    while queue and not found:
        state = queue.popleft()
        for bit, apex in enumerate(gate):
            if state & (1 << bit):
                continue
            raw = state | (1 << bit)
            if raw not in cache:
                cache[raw] = scalar_closure(gate, profiles, raw)
            target, rounds = cache[raw]
            edges += 1
            if target in distance:
                continue
            distance[target] = distance[state] + 1
            previous[target] = (state, apex, rounds)
            if target == full:
                found = True
                break
            queue.append(target)
    if full not in distance:
        raise RuntimeError("full gate is unreachable")

    steps: list[tuple[int, int, tuple[tuple[int, ...], ...]]] = []
    state = full
    while state != start:
        old, apex, rounds = previous[state]
        steps.append((old, apex, rounds))
        state = old
    steps.reverse()
    return {
        "minimum": distance[full],
        "start": start,
        "steps": tuple(steps),
        "states": len(distance),
        "edges": edges,
        "closure_calls": len(cache),
    }


def labels_from_mask(gate: tuple[int, ...], mask: int) -> tuple[int, ...]:
    return tuple(
        speed for bit, speed in enumerate(gate) if mask & (1 << bit)
    )


def main() -> None:
    expected = {
        row[0]: row[4] for row in M.EXPECTED_ROOT_RESULTS
    }
    totals = [0, 0, 0]
    for body in M.M.BATTERY_ROOTS:
        root = M.M.global_root(body)
        gate = tuple(speed for _, speed in root["top"][: root["K"]])
        profiles = {
            apex: M.profile_apex(root, gate, apex) for apex in gate
        }
        result = minimum_closed_state_path(gate, profiles)
        if result["minimum"] != expected[body]:
            raise RuntimeError(
                f"closed-state minimum changed at {body}: {result['minimum']}"
            )
        proof_steps = tuple(
            (
                labels_from_mask(gate, state),
                apex,
                rounds,
            )
            for state, apex, rounds in result["steps"]
        )
        print(
            f"E={body};K={len(gate)};minimum={result['minimum']};"
            f"closed_states={result['states']};edges={result['edges']};"
            f"closure_calls={result['closure_calls']};"
            f"interleaved_steps={proof_steps}"
        )
        totals[0] += result["states"]
        totals[1] += result["edges"]
        totals[2] += result["closure_calls"]
    print(
        f"TOTAL closed_states={totals[0]};edges={totals[1]};"
        f"closure_calls={totals[2]}"
    )
    print("all_exact_minima_reproduced=PASS")


if __name__ == "__main__":
    main()
