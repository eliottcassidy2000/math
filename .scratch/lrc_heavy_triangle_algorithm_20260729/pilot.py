#!/usr/bin/env python3
"""Pilot the THM-2888/2893 heavy-triangle -> singleton residual route.

This is deliberately a scratch-only design probe.  Its universe is the
deterministic 50-carrier sample already frozen by
``lrc14_j5_refined_residual_graph_sample_codex_20260729.py``; it is not the
10,202-carrier theorem universe.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SAMPLE_PATH = (
    ROOT
    / "04-computation/lrc14_j5_refined_residual_graph_sample_codex_20260729.py"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_sample():
    spec = importlib.util.spec_from_file_location("heavy_triangle_sample", SAMPLE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load sample")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_sample()
A = S.A


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def interval_is_tooth_covered(left: F, right: F, speed: int) -> bool:
    """Whether one interval lies in one radius-1/(14w) danger tooth."""

    lower_center = ceil_fraction(speed * right - F(1, 14))
    upper_center = floor_fraction(speed * left + F(1, 14))
    return lower_center <= upper_center


def singleton_covers(carrier: list[tuple[F, F]], speed: int) -> bool:
    """Exact early-exit predicate for ``carrier subset D_speed`` a.e."""

    return all(
        interval_is_tooth_covered(left, right, speed)
        for left, right in sorted(
            carrier,
            key=lambda interval: interval[1] - interval[0],
            reverse=True,
        )
    )


def singleton_profile(
    body: tuple[int, ...],
    prefix: tuple[int, ...],
    carrier: list[tuple[F, F]],
) -> dict[str, object]:
    """Seal every singleton by the longest-component tooth-width bound."""

    residual_mass = mass(carrier)
    require(residual_mass > 0 and carrier, "empty triangle residual")
    longest = max(right - left for left, right in carrier)
    geometric_last = floor_fraction(F(1, 7) / longest)
    geometric_tail_first = geometric_last + 1

    components = len(carrier)
    discrepancy_ratio = A.S2 * components / (6 * residual_mass)
    discrepancy_tail_first = ceil_fraction(discrepancy_ratio)
    require(
        residual_mass / 7
        + A.S2 * components / (7 * max(1, discrepancy_tail_first))
        <= residual_mass,
        "discrepancy singleton tail arithmetic failed",
    )
    require(
        geometric_tail_first <= discrepancy_tail_first,
        "longest-component seal should dominate total-mass discrepancy",
    )

    excluded = set(prefix)
    scanned = [
        speed
        for speed in range(A.FIRST_EXTERNAL, geometric_tail_first)
        if speed not in excluded
    ]
    covering = tuple(
        speed for speed in scanned if singleton_covers(carrier, speed)
    )

    # Compare the early-exit containment test to literal subtraction at both
    # scan boundaries and at every positive candidate.
    controls = []
    if scanned:
        controls.extend((scanned[0], scanned[-1]))
    controls.extend(covering)
    for speed in sorted(set(controls)):
        literal_empty = not A.P1.THM2883.subtract_local(carrier, speed)
        require(
            singleton_covers(carrier, speed) == literal_empty,
            f"containment/literal mismatch: {body}, {prefix}, {speed}",
        )
        scalar_full = A.P1.THM2885.coverage(carrier, speed) == residual_mass
        require(
            literal_empty == scalar_full,
            f"literal/scalar mismatch: {body}, {prefix}, {speed}",
        )

    return {
        "mass": residual_mass,
        "components": components,
        "longest": longest,
        "geometric_tail_first": geometric_tail_first,
        "discrepancy_tail_first": discrepancy_tail_first,
        "scanned": len(scanned),
        "covering": covering,
    }


def triangles_from_edges(
    vertices: tuple[int, ...],
    edges: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int, int], ...]:
    adjacency = {vertex: set() for vertex in vertices}
    for x, y in edges:
        adjacency[x].add(y)
        adjacency[y].add(x)
    triangles = []
    for x in vertices:
        for y in sorted(vertex for vertex in adjacency[x] if vertex > x):
            for z in sorted(
                vertex
                for vertex in adjacency[x].intersection(adjacency[y])
                if vertex > y
            ):
                triangles.append((x, y, z))
    brute = tuple(
        triple
        for triple in combinations(vertices, 3)
        if all(tuple(sorted(pair)) in set(edges) for pair in combinations(triple, 2))
    )
    require(tuple(triangles) == brute, "adjacency/brute triangle census differs")
    return tuple(triangles)


def main() -> None:
    sample = S.collect_sample()
    profiles = [S.coarse_profile(row) for row in sample]
    summary_rows = []
    triangle_rows = []
    total_triangles = 0
    total_singleton_scans = 0
    covering_residuals = 0
    empty_triangle_residuals = 0
    maximum_geom = None
    maximum_disc = None
    maximum_triangles = None

    for pair_row, profile in zip(sample, profiles):
        body = profile["body"]
        apex = profile["apex"]
        carrier = A.reconstruct_carrier(body, (apex,))
        edges = tuple((row[4], row[5]) for row in profile["coarse_rows"])
        triangles = triangles_from_edges(profile["head"], edges)
        total_triangles += len(triangles)
        maximum_triangles = max(
            maximum_triangles or (0, body, apex),
            (len(triangles), body, apex),
        )
        local_covering = 0
        for triangle in triangles:
            residual = A.P1.THM2883.subtract_local_multi(carrier, triangle)
            sequential = carrier
            for speed in triangle:
                sequential = A.P1.THM2883.subtract_local(sequential, speed)
            require(
                residual == sequential,
                f"triple simultaneous/sequential mismatch: {body}, {apex}, {triangle}",
            )
            if not residual:
                empty_triangle_residuals += 1
                triangle_rows.append(
                    f"{body};a={apex};T={triangle};empty=1\n"
                )
                continue
            singleton = singleton_profile(
                body,
                (apex, *triangle),
                residual,
            )
            total_singleton_scans += singleton["scanned"]
            if singleton["covering"]:
                covering_residuals += 1
                local_covering += 1
            maximum_geom = max(
                maximum_geom or (0, body, apex, triangle),
                (
                    singleton["geometric_tail_first"],
                    body,
                    apex,
                    triangle,
                ),
            )
            maximum_disc = max(
                maximum_disc or (0, body, apex, triangle),
                (
                    singleton["discrepancy_tail_first"],
                    body,
                    apex,
                    triangle,
                ),
            )
            triangle_rows.append(
                f"{body};a={apex};T={triangle};"
                f"m={singleton['mass']};r={singleton['components']};"
                f"L={singleton['longest']};"
                f"G={singleton['geometric_tail_first']};"
                f"D={singleton['discrepancy_tail_first']};"
                f"scan={singleton['scanned']};cover={singleton['covering']}\n"
            )
        summary_rows.append(
            f"{body};rank={pair_row['apex_rank']};a={apex};"
            f"H={len(profile['head'])};E={len(edges)};"
            f"T={len(triangles)};cover={local_covering}\n"
        )

    digest = hashlib.sha256(
        (
            "LRC14/heavy-triangle-singleton-pilot/v1\n"
            + "".join(summary_rows)
            + "".join(triangle_rows)
        ).encode()
    ).hexdigest()
    print("LRC14 HEAVY-TRIANGLE SINGLETON PILOT")
    print("scope=existing deterministic 50-carrier sample;not theorem universe")
    print(
        f"carriers={len(sample)};"
        f"high_vertices={sum(len(profile['head']) for profile in profiles)};"
        f"heavy_edges={sum(len(profile['coarse_rows']) for profile in profiles)};"
        f"heavy_triangles={total_triangles};"
        f"max_triangles={maximum_triangles};"
        f"empty_triangle_residuals={empty_triangle_residuals}"
    )
    print(
        f"singleton_scans={total_singleton_scans};"
        f"covering_residuals={covering_residuals};"
        f"maximum_geometric_tail_first={maximum_geom};"
        f"maximum_discrepancy_tail_first={maximum_disc};"
        f"digest={digest}"
    )
    print("controls=adjacency-vs-brute;triple-path;containment-vs-literal-vs-scalar")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
