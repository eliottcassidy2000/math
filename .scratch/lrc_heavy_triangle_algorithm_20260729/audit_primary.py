#!/usr/bin/env python3
"""Independent bounded replay of the concurrent heavy-triangle primary."""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = Path(__file__).resolve().with_name("pilot.py")
PRIMARY_PATH = Path(
    "/private/tmp/math-wt-lrc-j5-recursion-20260729/"
    "04-computation/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.py"
)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


I = load("heavy_triangle_independent_pilot", PILOT_PATH)
P = load("heavy_triangle_concurrent_primary", PRIMARY_PATH)


def main() -> None:
    sample = I.S.collect_sample()
    independent_profiles = [I.S.coarse_profile(row) for row in sample]
    primary_profiles = []
    for pair_row in sample:
        primary_profiles.append(
            P.profile_coarse(
                (
                    pair_row["body"],
                    pair_row["apex_rank"],
                    pair_row["apex"],
                    pair_row["m"],
                    pair_row["r"],
                    pair_row["pair_cap"],
                    None,
                    None,
                )
            )
        )
    for pair_row, independent, primary in zip(
        sample,
        independent_profiles,
        primary_profiles,
    ):
        independent_edges = tuple(
            (row[4], row[5]) for row in independent["coarse_rows"]
        )
        independent_triangles = I.triangles_from_edges(
            independent["head"],
            independent_edges,
        )
        if (
            primary["head_size"] != len(independent["head"])
            or primary["edge_count"] != len(independent_edges)
            or primary["triangles"] != independent_triangles
        ):
            raise RuntimeError(
                "coarse primary/independent mismatch: "
                f"{pair_row['body']}, {pair_row['apex']}"
            )

    triangle_tasks = [
        (
            row["body"],
            row["apex_rank"],
            row["apex"],
            triangle,
            residual_mass,
        )
        for row in primary_profiles
        for triangle, residual_mass in zip(
            row["triangles"],
            row["triangle_masses"],
        )
    ]
    extensions = [P.seal_triangle_extension(task) for task in triangle_tasks]
    if any(not row["closed"] or row["empty"] or row["extensions"] for row in extensions):
        raise RuntimeError("primary found an unexpected sample extension")

    # Replay the independent longest-component predicate on every same literal
    # residual; it must agree that there is no fourth-speed cover.
    geometric_scans = 0
    max_geometric_tail_first = 0
    max_primary_horizon = 0
    for row in extensions:
        root_good, _, _ = P.T.CORE.good_norm(row["body"])
        residual = P.T.THM2883.subtract_local_multi(
            root_good,
            (row["apex"], *row["triangle"]),
        )
        singleton = I.singleton_profile(
            row["body"],
            (row["apex"], *row["triangle"]),
            residual,
        )
        if singleton["covering"]:
            raise RuntimeError("geometric predicate found an unexpected sample cover")
        geometric_scans += singleton["scanned"]
        max_geometric_tail_first = max(
            max_geometric_tail_first,
            singleton["geometric_tail_first"],
        )
        max_primary_horizon = max(max_primary_horizon, row["horizon"])

    print("PRIMARY HEAVY-TRIANGLE 50-CARRIER INDEPENDENT AUDIT")
    print(
        f"carriers={len(sample)};"
        f"vertices={sum(row['head_size'] for row in primary_profiles)};"
        f"edges={sum(row['edge_count'] for row in primary_profiles)};"
        f"triangles={len(triangle_tasks)};"
        f"extensions_closed={sum(row['closed'] for row in extensions)}"
    )
    print(
        f"geometric_scans={geometric_scans};"
        f"maximum_geometric_tail_first={max_geometric_tail_first};"
        f"maximum_primary_discrepancy_horizon={max_primary_horizon}"
    )
    print(
        "controls=independent-H/edge/triangle-generation;"
        "primary-discrepancy-vs-geometric-singleton-seal"
    )
    print("verdict=PASS")


if __name__ == "__main__":
    main()
