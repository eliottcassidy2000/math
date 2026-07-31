#!/usr/bin/env python3
"""Independent exact checks of the discovery-run sharp margins.

The primary heavy-triangle module is intentionally not imported.  This loads
only the already-canonical THM-2888 arithmetic dependencies, reconstructs the
literal residual from ``CORE.good_norm``, and uses the scalar tooth primitive
at every finite speed.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATLAS_PATH = (
    ROOT
    / "04-computation/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.py"
)
ATLAS_SHA256 = "cba433bce508ca8cc1e90c813e1988bb73c765ffafe350b74c3bad240eeca10f"
FIRST_EXTERNAL = 15
S2 = F(99, 70)


def load_atlas():
    if hashlib.sha256(ATLAS_PATH.read_bytes()).hexdigest() != ATLAS_SHA256:
        raise RuntimeError("canonical THM-2888 arithmetic dependency changed")
    spec = importlib.util.spec_from_file_location(
        "heavy_triangle_extrema_atlas",
        ATLAS_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load canonical atlas")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_atlas()


def mass(carrier):
    return sum((right - left for left, right in carrier), F(0))


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def main() -> None:
    body = (1, 5, 6, 7, 8, 9, 11, 13)
    rank = 2
    apex = 19
    triangle = (37, 121, 130)
    root_good, _, _ = T.CORE.good_norm(body)
    residual = T.THM2883.subtract_local_multi(
        root_good,
        (apex, *triangle),
    )
    residual_mass = mass(residual)
    full_family = tuple(sorted((*body, apex, *triangle)))
    direct_good, direct_r, direct_mass = T.CORE.good_norm(full_family)
    if (
        residual != direct_good
        or residual_mass != direct_mass
        or len(residual) != direct_r
    ):
        raise RuntimeError("literal/direct sharp residual reconstruction failed")
    if (
        residual_mass != F(1_099_999, 24_444_420)
        or len(residual) != 42
    ):
        raise RuntimeError(
            f"sharp residual mass/components changed: {residual_mass}, {len(residual)}"
        )

    # Independently rebuild the heavy core and graph for the sharp carrier
    # using scalar coverages at every candidate speed.
    apex_rows = T.profile_body(body)
    apex_row = next(row for row in apex_rows if row["apex"] == apex)
    hostile_keys = {
        (hostile_body, hostile_apex)
        for hostile_body, _, hostile_apex, _, _, _ in T.EXPECTED_HOSTILES
    }
    if (
        apex_row["apex_rank"] != rank
        or apex_row["scalar_class"] == "direct"
        or apex_row["pairpair_direct"]
        or (body, apex) in hostile_keys
        or apex_row["margin"] <= 0
    ):
        raise RuntimeError("sharp carrier is not an honest nonterminal target")
    sharp_terminals = {
        row["apex"]
        for row in apex_rows
        if row["scalar_class"] == "direct"
        or row["pairpair_direct"]
        or (body, row["apex"]) in hostile_keys
    }
    sharp_root = T.THM2885.profile_body(body)
    sharp_allowed = tuple(
        (value, speed)
        for rank_index, (value, speed) in enumerate(
            sharp_root["top15"],
            start=1,
        )
        if rank_index >= 11 or speed not in sharp_terminals
    )
    sharp_root_margin = sharp_root["m"] - sum(
        (value for value, _ in sharp_allowed[:5]),
        F(0),
    )
    if sharp_root_margin > 0:
        raise RuntimeError("sharp carrier belongs to an initially closed root")
    carrier = T.THM2883.subtract_local(root_good, apex)
    carrier_mass = mass(carrier)
    theta = carrier_mass - apex_row["global_cap"]
    level = theta / 2
    core_tail_first = ceil_fraction(
        S2 * len(carrier) / (7 * (level - carrier_mass / 7))
    )
    scalar_core_rows = [
        (T.THM2885.coverage(carrier, speed), speed)
        for speed in range(FIRST_EXTERNAL, core_tail_first)
        if speed != apex
    ]
    high = tuple(
        sorted(speed for value, speed in scalar_core_rows if value >= level)
    )
    high_coverage = {speed: value for value, speed in scalar_core_rows if speed in high}
    edges = set()
    for x, y in combinations(high, 2):
        if high_coverage[x] + high_coverage[y] < theta:
            continue
        pair_residual = T.THM2883.subtract_local_multi(carrier, (x, y))
        if carrier_mass - mass(pair_residual) >= theta:
            edges.add((x, y))
    triangles = tuple(
        triple
        for triple in combinations(high, 3)
        if all(tuple(sorted(pair)) in edges for pair in combinations(triple, 2))
    )
    if triangle not in triangles:
        raise RuntimeError("sharp triangle missing from independent heavy graph")

    tail_threshold = S2 * len(residual) / (6 * residual_mass)
    horizon = max(
        FIRST_EXTERNAL - 1,
        tail_threshold.numerator // tail_threshold.denominator,
    )
    excluded = {apex, *triangle}
    scalar_rows = [
        (T.THM2885.coverage(residual, speed), speed)
        for speed in range(FIRST_EXTERNAL, horizon + 1)
        if speed not in excluded
    ]
    finite_max, maximizer = max(
        scalar_rows,
        key=lambda item: (item[0], -item[1]),
    )
    tail_cap = (
        residual_mass / 7
        + S2 * len(residual) / (7 * (horizon + 1))
    )
    cap = max(finite_max, tail_cap)
    margin = residual_mass - cap
    if (
        horizon != 219
        or finite_max != F(7_111, 459_800)
        or maximizer != 25
        or tail_cap != F(19_249_981, 427_777_350)
        or cap != tail_cap
        or margin != F(1, 285_184_900)
        or any(value == residual_mass for value, _ in scalar_rows)
        or not tail_cap < residual_mass
        or not T.THM2883.subtract_local(residual, maximizer)
    ):
        raise RuntimeError(
            "independent sharp triangle scalar replay failed: "
            f"horizon={horizon}, finite={finite_max}, cap={cap}, "
            f"maximizer={maximizer}, "
            f"finite_margin={residual_mass-finite_max}, "
            f"margin={margin}, tail={tail_cap}, "
            f"literal_nonempty={bool(T.THM2883.subtract_local(residual, maximizer))}"
        )

    root_body = (2, 5, 7, 8, 10, 12, 13, 14)
    root = T.THM2885.profile_body(root_body)
    body_rows = T.profile_body(root_body)
    root_hostile_keys = {
        (hostile_body, hostile_apex)
        for hostile_body, _, hostile_apex, _, _, _ in T.EXPECTED_HOSTILES
    }
    terminals = {
        row["apex"]
        for row in body_rows
        if row["scalar_class"] == "direct"
        or row["pairpair_direct"]
        or (root_body, row["apex"]) in root_hostile_keys
    }
    allowed = tuple(
        (value, speed)
        for rank_index, (value, speed) in enumerate(root["top15"], start=1)
        if rank_index >= 11 or speed not in terminals
    )
    final_top5 = allowed[:5]
    root_good, _, root_mass = T.CORE.good_norm(root_body)
    scalar_top5 = tuple(
        (T.THM2885.coverage(root_good, speed), speed)
        for _, speed in final_top5
    )
    if scalar_top5 != final_top5 or root_mass != root["m"]:
        raise RuntimeError("sharp final root scalar reconstruction failed")
    final_margin = root["m"] - sum(
        (value for value, _ in scalar_top5),
        F(0),
    )
    if final_margin != F(3_973, 89_664_120):
        raise RuntimeError(f"sharp final root margin changed: {final_margin}")

    print("HEAVY-TRIANGLE EXTREMA INDEPENDENT AUDIT")
    print(
        f"triangle_body={body};apex={apex};triangle={triangle};"
        f"initial_root_margin={sharp_root_margin};"
        f"core_tail_first={core_tail_first};H={len(high)};"
        f"heavy_edges={len(edges)};heavy_triangles={len(triangles)};"
        f"mass={residual_mass};components={len(residual)};"
        f"tail_threshold={tail_threshold};horizon={horizon};maximizer={maximizer};"
        f"finite_max={finite_max};tail_cap={tail_cap};"
        f"margin={margin};finite_speeds_scanned={len(scalar_rows)}"
    )
    print(
        f"root_body={root_body};"
        f"initial_terminals={tuple(sorted(terminals))};"
        f"allowed_top5={tuple(speed for _, speed in final_top5)};"
        f"margin={final_margin}"
    )
    print("controls=scalar-not-vectorized-sharp-row;literal-nonempty-maximizer;root-recompute")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
