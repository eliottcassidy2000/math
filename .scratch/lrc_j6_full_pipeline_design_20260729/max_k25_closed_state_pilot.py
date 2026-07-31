#!/usr/bin/env python3
"""Run the exact scalar-closed-state seed search on the unique K=25 root.

This is the smallest hostile all-root successor experiment after the
top-forty gate atlas: construct all 25 reusable apex profiles, find an exact
minimum scalar-bootstrap path, and audit whether every paid edge satisfies
the first singleton-complement condition B1 < 3h/7.  For each such edge it
also enumerates the globally sealed H4 core and reports the resulting pair
workload, without running the deeper pair-residual certificates.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATLAS_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.py"
)
MINSEED_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.py"
)
BFS_PATH = (
    ROOT
    / ".scratch/lrc_j6_full_pipeline_design_20260729/closed_state_bfs_probe.py"
)
BODY = (1, 8, 10, 11, 12, 13, 14)


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load(ATLAS_PATH, "maxk_atlas")
M = load(MINSEED_PATH, "maxk_minseed")
B = load(BFS_PATH, "maxk_bfs")


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def labels(gate: tuple[int, ...], mask: int) -> tuple[int, ...]:
    return tuple(
        speed for bit, speed in enumerate(gate) if mask & (1 << bit)
    )


def q1_after(
    profile: dict[str, object],
    excluded: set[int],
) -> tuple[F, int]:
    for value, speed in profile["ranked_candidates"]:
        if speed not in excluded:
            return value, speed
    raise RuntimeError("empty apex candidate ledger")


def main() -> None:
    base = A.profile_body(BODY)
    if base["adaptive_k"] != 25:
        raise RuntimeError("hostile body is no longer the unique K=25 root")
    good, components, mass = A.T.CORE.good_norm(BODY)
    root = {
        **base,
        "good": good,
        "K": base["adaptive_k"],
    }
    if components != root["r"] or mass != root["m"]:
        raise RuntimeError("root reconstruction changed")
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    profiles = {
        apex: M.profile_apex(root, gate, apex) for apex in gate
    }
    search = B.minimum_closed_state_path(gate, profiles)
    print(
        f"E={BODY};K={len(gate)};minimum={search['minimum']};"
        f"closed_states={search['states']};edges={search['edges']};"
        f"closure_calls={search['closure_calls']}"
    )

    total_pairs = 0
    failures = 0
    for index, (state, apex, rounds) in enumerate(search["steps"], start=1):
        excluded = set(labels(gate, state))
        excluded.add(apex)
        profile = profiles[apex]
        q1, q1_speed = q1_after(profile, excluded)
        singleton_margin = F(3, 7) * profile["m"] - q1
        if singleton_margin <= 0:
            failures += 1
            print(
                f"paid={index};P={state.bit_count()};a={apex};"
                f"cascade={rounds};B1_condition=FAIL;"
                f"3h7-q1={M.ftext(singleton_margin)};q1speed={q1_speed}"
            )
            continue

        level = (profile["m"] - q1) / 4
        if level <= profile["m"] / 7:
            raise RuntimeError("positive singleton margin did not finite H4")
        threshold = (
            M.S2
            * profile["r"]
            / (7 * (level - profile["m"] / 7))
        )
        tail_first = max(M.FIRST_EXTERNAL, ceiling(threshold))
        literal, direct_r, direct_m = M.M.T.CORE.good_norm(
            tuple(sorted((*BODY, apex)))
        )
        if (
            direct_r != profile["r"]
            or direct_m != profile["m"]
        ):
            raise RuntimeError("paid carrier reconstruction changed")
        rows = M.M.T.coverages_many(
            literal,
            [
                speed
                for speed in range(M.FIRST_EXTERNAL, tail_first)
                if speed not in excluded
            ],
        )
        if (
            profile["m"] / 7
            + M.S2 * profile["r"] / (7 * tail_first)
            > level
        ):
            raise RuntimeError("H4 tail failed to seal")
        core = tuple(
            speed for value, speed in rows if value >= level
        )
        pair_count = len(tuple(combinations(core, 2)))
        total_pairs += pair_count
        print(
            f"paid={index};P={state.bit_count()};a={apex};"
            f"cascade={rounds};B1_condition=PASS;"
            f"3h7-q1={M.ftext(singleton_margin)};q1speed={q1_speed};"
            f"Htail={tail_first};H={len(core)};pairs={pair_count}"
        )

    print(
        f"paid_edges={search['minimum']};B1_failures={failures};"
        f"predicted_H4_pairs={total_pairs}"
    )
    print("scope=unique K25 root;H4 structural layer only")


if __name__ == "__main__":
    main()
