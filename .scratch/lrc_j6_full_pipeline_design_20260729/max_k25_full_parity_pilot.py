#!/usr/bin/env python3
"""Attempt the full singleton-complement certificate on the K=25 BFS path.

This composes the unique hostile root's exact minimum closed-state path with
the existing suffix-H4 pair and recursive H2/singleton routines.  It is a
scratch discovery computation; it does not promote the reserved parity
theorem or claim a uniform j6 result.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = (
    ROOT
    / ".scratch/lrc_j6_full_pipeline_design_20260729/max_k25_closed_state_pilot.py"
)
PARITY_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py"
)


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


X = load(PILOT_PATH, "maxk_full_source")
Q = load(PARITY_PATH, "maxk_full_parity")


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def main() -> None:
    base = X.A.profile_body(X.BODY)
    good, components, root_mass = X.A.T.CORE.good_norm(X.BODY)
    root = {**base, "good": good, "K": base["adaptive_k"]}
    if components != root["r"] or root_mass != root["m"]:
        raise RuntimeError("root reconstruction changed")
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    profiles = {
        apex: X.M.profile_apex(root, gate, apex) for apex in gate
    }
    search = X.B.minimum_closed_state_path(gate, profiles)

    all_pairs: list[tuple[dict[str, object], dict[str, object]]] = []
    branch_rows: list[dict[str, object]] = []
    for index, (state, apex, rounds) in enumerate(search["steps"], start=1):
        prior = X.labels(gate, state)
        excluded = set(prior) | {apex}
        profile = profiles[apex]
        q1, _ = X.q1_after(profile, excluded)
        singleton_margin = F(3, 7) * profile["m"] - q1
        if singleton_margin <= 0:
            raise RuntimeError(f"paid branch misses H4 condition: {apex}")
        level = (profile["m"] - q1) / 4
        threshold = (
            X.M.S2
            * profile["r"]
            / (7 * (level - profile["m"] / 7))
        )
        tail_first = max(X.M.FIRST_EXTERNAL, X.ceiling(threshold))
        carrier, direct_r, direct_m = X.M.M.T.CORE.good_norm(
            tuple(sorted((*X.BODY, apex)))
        )
        if (
            len(carrier) != direct_r
            or direct_r != profile["r"]
            or direct_m != profile["m"]
            or mass(carrier) != direct_m
        ):
            raise RuntimeError("paid carrier reconstruction changed")
        rows = X.M.M.T.coverages_many(
            carrier,
            [
                speed
                for speed in range(X.M.FIRST_EXTERNAL, tail_first)
                if speed not in excluded
            ],
        )
        if (
            profile["m"] / 7
            + X.M.S2 * profile["r"] / (7 * tail_first)
            > level
        ):
            raise RuntimeError("H4 tail failed to seal")
        core = tuple(
            speed for value, speed in rows if value >= level
        )
        branch = {
            "body": X.BODY,
            "rank": index,
            "root_rank": gate.index(apex) + 1,
            "apex": apex,
            "excluded_prefix": (*prior, apex),
            "m": profile["m"],
            "r": profile["r"],
        }
        row = {
            "root": root,
            "branch": branch,
            "carrier": carrier,
            "excluded": excluded,
            "q1": q1,
            "singleton_margin": singleton_margin,
            "level": level,
            "Htail": tail_first,
            "H": core,
            "cascade": rounds,
        }
        branch_rows.append(row)
        local = [
            Q.pair_residual(row, hpair)
            for hpair in combinations(core, 2)
        ]
        all_pairs.extend((row, pair) for pair in local)
        hard_count = sum(not pair["adaptive_closed"] for pair in local)
        print(
            f"paid={index};P={len(prior)};a={apex};H={len(core)};"
            f"pairs={len(local)};"
            f"top3={sum(pair['direct_margin'] > 0 for pair in local)};"
            f"B2B1={sum(pair['pair_margin'] > 0 for pair in local)};"
            f"union={len(local)-hard_count};hard={hard_count};"
            f"paid_pair_unions={sum(pair['cap']['paid'] for pair in local)}"
        )

    hard = [
        (row, pair)
        for row, pair in all_pairs
        if not pair["adaptive_closed"]
    ]
    recursive_rows = [
        (row, pair, Q.recursive_k3_close(row, pair))
        for row, pair in hard
    ]
    for row, pair, recursive in recursive_rows:
        print(
            f"recursive=a:{row['branch']['apex']};"
            f"H4pair={pair['hpair']};H2={recursive['H']};"
            f"heavy={recursive['heavy']};W={recursive['horizons']};"
            f"checks={recursive['checks']};covers={recursive['covers']}"
        )
    print(
        f"TOTAL paid_edges={len(branch_rows)};pairs={len(all_pairs)};"
        f"top3={sum(pair['direct_margin'] > 0 for _, pair in all_pairs)};"
        f"B2B1={sum(pair['pair_margin'] > 0 for _, pair in all_pairs)};"
        f"union={len(all_pairs)-len(hard)};hard={len(hard)};"
        f"recursive_closed="
        f"{sum(recursive['closed'] for _, _, recursive in recursive_rows)};"
        f"recursive_open="
        f"{sum(not recursive['closed'] for _, _, recursive in recursive_rows)};"
        f"singleton_checks="
        f"{sum(recursive['checks'] for _, _, recursive in recursive_rows)}"
    )
    if not all(
        recursive["closed"] for _, _, recursive in recursive_rows
    ):
        raise RuntimeError("max-K path has an open recursive row")
    print(
        "candidate_consequence=all paid edges on the exact-minimum path close;"
        "scalar cascades close the K25 gate"
    )
    print("scope=one root;reserved parity mechanism;not uniform;not LRC14")


if __name__ == "__main__":
    main()
