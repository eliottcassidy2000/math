#!/usr/bin/env python3
r"""Exact heavy-triangle census on the honest residual j=5 apex branches.

THM-2888 gives a global pair-union cap ``B<5h/7`` for every nonterminal
first-apex carrier left by its genuine weighted root gate.  Put

    theta = h-B > 2h/7,       H = {x : c(x) >= theta/2}.

If four further dangers cover the carrier, every pair among their speeds has
union at least ``theta``: the complementary pair has union at most ``B``.
Every theta-heavy pair meets the finite set ``H``, so at least three vertices
of the hypothetical K4 lie in H.  Those three form a theta-heavy triangle.

This verifier enumerates every such triangle in H.  For each one it forms the
literal residual after its three dangers and globally excludes a fourth danger
by an exact finite scan plus the strict THM-2885 discrepancy tail

    c_R(w) < |R|/7 + (99/70) components(R)/(7w).

Thus a triangle is nonextendible as soon as the finite maximum and the sealed
tail cap are both strictly below ``|R|``.  This keeps the possible fourth
speed outside H; an H-only K4 census would not be sufficient.

The universe is restricted to the nonterminal apices of the active roots
which the honest THM-2888 root gate does not already close.  Rank-three
finiteness is never treated as literal positivity.  The 495 bodies contained
in {1,...,12} remain the separate THM-885 terminal chamber.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = (
    ROOT
    / "04-computation/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.py"
)
ATLAS_SHA256 = (
    "cba433bce508ca8cc1e90c813e1988bb73c765ffafe350b74c3bad240eeca10f"
)
ATLAS_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j5_first_apex_pair_cap_atlas_codex_20260729.out"
)
ATLAS_OUTPUT_SHA256 = (
    "b6bcbf90257942523e7d26d33a2c06bde67805e1e008bfe190ccca7df0d83669"
)

FIRST_EXTERNAL = 15
S2 = F(99, 70)
EXPECTED_ATLAS_COUNTS = (13_969, 10_939, 5_122, 13_802, 16_122)
EXPECTED_ROOT_COUNTS = (2_508, 1_064, 1_444, 10_202)

EXPECTED_TARGET_DIGEST = (
    "4040369ff36500cc804c6c2ffa54fdb1f3072e2aa775c9db162b459209d10316"
)
EXPECTED_NORMALIZATION_COUNTS = (26, 130_827, 1_497, 282, 321, 2_409)
EXPECTED_NORMALIZATION_DIGEST = (
    "adb9881d84a7d209ec333ef5ff7ac4548d15f11fd6d715b7f4651ba3aacf8c54"
)
EXPECTED_COARSE_COUNTS = (
    55_159,
    65,
    140_846,
    1_021,
    372_209,
    7_814,
    4_578,
    1_007_764,
    36_265,
    20_404,
    2_410,
    4_578,
    26,
    0,
)
EXPECTED_COARSE_DIGEST = (
    "2b00dd59fd7320c95a2e29226e540df5eb890fe1fff58e7b6fbb708bd1cffd92"
)
EXPECTED_TRIANGLE_COUNTS = (372_209, 372_209, 0, 0, 744_418, 476)
EXPECTED_TRIANGLE_DIGEST = (
    "a49a961b66ce2e5fee114313ef58eb2c2808f5c73b70b2e61839b9a2eb3c96bf"
)
EXPECTED_MINIMUM_TRIANGLE_MARGIN = (
    F(1, 285_184_900),
    (1, 5, 6, 7, 8, 9, 11, 13),
    2,
    19,
    (37, 121, 130),
    25,
)
EXPECTED_HYPEREDGE_COUNTS = (0, 0, 0, 0, 0)
EXPECTED_HYPEREDGE_DIGEST = (
    "07de7784ad51a6d10d450dc9bb99edf738c2f6bcdcf9d78b2fd59803dedf8508"
)
EXPECTED_CARRIER_COUNTS = (10_202, 5_624, 4_578, 10_202, 0)
EXPECTED_CARRIER_DIGEST = (
    "1f028c236253428f3459c34b9efe9428410fd4a3858343ef2383f9679f1f764d"
)
EXPECTED_FINAL_ROOT_COUNTS = (2_508, 2_508, 0)
EXPECTED_ROOT_DIGEST = (
    "143f995d7aaf4afedba4208b9e09d7611bd6ae8b86a3e6df465a57c97837005f"
)
EXPECTED_MINIMUM_ROOT_MARGIN = (
    F(3_973, 89_664_120),
    (2, 5, 7, 8, 10, 12, 13, 14),
)
EXPECTED_STRATUM_COUNTS = (
    1_584,
    718,
    866,
    6_032,
    32_715,
    86_017,
    237_601,
    237_601,
    0,
    6_032,
    0,
    1_584,
    0,
    0,
    924,
    346,
    578,
    4_170,
    22_444,
    54_829,
    134_608,
    134_608,
    0,
    4_170,
    0,
    924,
    0,
    0,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_atlas():
    require(file_sha256(ATLAS_PATH) == ATLAS_SHA256, "THM-2888 verifier changed")
    require(
        file_sha256(ATLAS_OUTPUT_PATH) == ATLAS_OUTPUT_SHA256,
        "THM-2888 transcript changed",
    )
    transcript = ATLAS_OUTPUT_PATH.read_text()
    require(
        "active_residual=1444;active_residual_apex_branches=10202"
        in transcript
        and transcript.endswith("all_exact_controls=PASS\n"),
        "THM-2888 root-gate verdict changed",
    )
    spec = importlib.util.spec_from_file_location("j5_postgate_thm2888", ATLAS_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-2888")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_atlas()
HOSTILE_KEYS = {
    (body, apex)
    for body, _, apex, _, _, _ in T.EXPECTED_HOSTILES
}


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def digest_text(header: str, rows: list[str]) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def terminal(row: dict[str, object]) -> bool:
    return bool(
        row["scalar_class"] == "direct"
        or row["pairpair_direct"]
        or (row["body"], row["apex"]) in HOSTILE_KEYS
    )


def profile_atlas_body(body: tuple[int, ...]) -> list[dict[str, object]]:
    """Pickle-stable wrapper around the hash-loaded THM-2888 worker."""

    return T.profile_body(body)


def normalize_atlas_row(row: dict[str, object]) -> dict[str, object]:
    """Pickle-stable wrapper for THM-2888's proved heavy-edge deletion."""

    return T.normalize_heavy_edge(row)


def body_stratum(body: tuple[int, ...]) -> str:
    high_count = int(13 in body) + int(14 in body)
    require(high_count in (1, 2), f"inactive body reached stratum: {body}")
    return "j6_one_high" if high_count == 1 else "j7_both_high"


def complement_margin(
    top15: tuple[tuple[F, int], ...],
    terminals: set[int],
) -> tuple[F, tuple[tuple[F, int], ...]]:
    allowed = tuple(
        (value, speed)
        for rank, (value, speed) in enumerate(top15, start=1)
        if rank >= 11 or speed not in terminals
    )
    require(len(allowed) >= 5, "weighted complement has fewer than five speeds")
    top5 = allowed[:5]
    return (
        -sum((value for value, _ in top5), F(0)),
        top5,
    )


def profile_coarse(
    task: tuple[
        tuple[int, ...],
        int,
        int,
        F,
        int,
        F,
        tuple[int, int] | None,
        F | None,
    ],
) -> dict[str, object]:
    (
        body,
        apex_rank,
        apex,
        expected_mass,
        expected_components,
        pair_cap,
        excluded_edge,
        forced_margin,
    ) = task
    root_good, _, _ = T.CORE.good_norm(body)
    carrier = T.THM2883.subtract_local(root_good, apex)
    mass = interval_mass(carrier)
    components = len(carrier)
    direct_good, direct_r, direct_m = T.CORE.good_norm((*body, apex))
    require(
        carrier == direct_good
        and mass == direct_m == expected_mass
        and components == direct_r == expected_components,
        f"literal first-apex reconstruction changed: {body}, {apex}",
    )
    theta = mass - pair_cap
    if 2 * pair_cap < mass:
        require(
            excluded_edge is not None and forced_margin is not None,
            f"unexpected undeclared pair-direct target: {body}, {apex}",
        )
        local_digest = digest_text(
            "LRC14/j5/postgate-heavy-triangle/coarse-local/v1\n",
            [
                "key="
                + ",".join(map(str, body))
                + f";rank={apex_rank};a={apex};h={ftext(mass)};"
                + f"B={ftext(pair_cap)};theta={ftext(theta)};"
                + f"excluded={excluded_edge};"
                + f"forced_margin={ftext(forced_margin)};"
                + "deleted_pair_direct=1\n"
            ],
        )
        return {
            "body": body,
            "apex_rank": apex_rank,
            "apex": apex,
            "head_size": 0,
            "edge_count": 0,
            "triangle_count": 0,
            "k4_count": 0,
            "triangles": (),
            "triangle_masses": (),
            "tail_first": 0,
            "controls": 0,
            "triple_controls": 0,
            "normalized": True,
            "normalized_pair_direct": True,
            "forced_closed": forced_margin > 0,
            "digest": local_digest,
        }
    require(
        theta > 2 * mass / 7 and theta <= mass / 2,
        f"target is outside the honest pair-cap residual: {body}, {apex}",
    )
    level = theta / 2
    threshold = S2 * components / (7 * (level - mass / 7))
    tail_first = ceiling(threshold)
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= level,
        f"finite H tail did not seal: {body}, {apex}",
    )
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    coverage_rows = T.THM2885.coverages_many(carrier, speeds)
    require(
        all(F(0) <= value <= mass for value, _ in coverage_rows),
        f"H coverage escaped carrier mass: {body}, {apex}",
    )
    by_speed = {speed: value for value, speed in coverage_rows}
    controls = 0
    if speeds:
        for speed in dict.fromkeys((speeds[0], speeds[-1])):
            require(
                by_speed[speed] == T.THM2885.coverage(carrier, speed),
                f"H vector/scalar control failed: {body}, {apex}, {speed}",
            )
            controls += 1
    head = tuple(
        sorted(speed for speed, value in by_speed.items() if value >= level)
    )
    edge_union: dict[tuple[int, int], F] = {}
    for x, y in combinations(head, 2):
        if excluded_edge is not None and (x, y) == excluded_edge:
            continue
        if by_speed[x] + by_speed[y] < theta:
            continue
        residual = T.THM2883.subtract_local_multi(carrier, (x, y))
        union = mass - interval_mass(residual)
        if union >= theta:
            edge_union[(x, y)] = union
    adjacency = {speed: set() for speed in head}
    for x, y in edge_union:
        adjacency[x].add(y)
        adjacency[y].add(x)
    triangles: list[tuple[int, int, int]] = []
    triangle_masses: list[F] = []
    triple_controls = 0
    for x, y, z in combinations(head, 3):
        if y in adjacency[x] and z in adjacency[x] and z in adjacency[y]:
            residual = T.THM2883.subtract_local_multi(carrier, (x, y, z))
            if triple_controls == 0:
                sequential = carrier
                for speed in (x, y, z):
                    sequential = T.THM2883.subtract_local(sequential, speed)
                require(
                    sequential == residual,
                    f"triple subtraction order changed: {body}, {apex}, {(x, y, z)}",
                )
                triple_controls += 1
            triangles.append((x, y, z))
            triangle_masses.append(interval_mass(residual))
    k4_count = 0
    for x, y, z in triangles:
        k4_count += sum(
            w in adjacency[x] and w in adjacency[y] and w in adjacency[z]
            for w in head
            if w > z
        )
    local_rows = [
        "key="
        + ",".join(map(str, body))
        + f";rank={apex_rank};a={apex};h={ftext(mass)};"
        + f"B={ftext(pair_cap)};theta={ftext(theta)};"
        + f"excluded={excluded_edge};forced_margin="
        + ("NA" if forced_margin is None else ftext(forced_margin))
        + ";"
        + f"Htail={tail_first};H="
        + ",".join(f"{speed}:{ftext(by_speed[speed])}" for speed in head)
        + "\n"
    ]
    local_rows.extend(
        f"edge={x},{y};U={ftext(union)}\n"
        for (x, y), union in sorted(edge_union.items())
    )
    local_rows.extend(
        f"triangle={x},{y},{z};L={ftext(residual_mass)}\n"
        for (x, y, z), residual_mass in zip(triangles, triangle_masses)
    )
    local_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/coarse-local/v1\n",
        local_rows,
    )
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "head_size": len(head),
        "edge_count": len(edge_union),
        "triangle_count": len(triangles),
        "k4_count": k4_count,
        "triangles": tuple(triangles),
        "triangle_masses": tuple(triangle_masses),
        "tail_first": tail_first,
        "controls": controls,
        "triple_controls": triple_controls,
        "normalized": excluded_edge is not None,
        "normalized_pair_direct": False,
        "forced_closed": forced_margin is None or forced_margin > 0,
        "digest": local_digest,
    }


def seal_triangle_extension(
    task: tuple[tuple[int, ...], int, int, tuple[int, int, int], F],
) -> dict[str, object]:
    body, apex_rank, apex, triangle, expected_mass = task
    root_good, _, _ = T.CORE.good_norm(body)
    residual = T.THM2883.subtract_local_multi(
        root_good,
        (apex, *triangle),
    )
    mass = interval_mass(residual)
    family = tuple(sorted((*body, apex, *triangle)))
    direct_good, direct_r, direct_m = T.CORE.good_norm(family)
    require(
        direct_good == residual and direct_m == mass == expected_mass,
        f"literal triangle residual changed: {body}, {apex}, {triangle}",
    )
    if mass == 0:
        local_digest = digest_text(
            "LRC14/j5/postgate-heavy-triangle/extension-local/v1\n",
            [
                "key="
                + ",".join(map(str, body))
                + f";rank={apex_rank};a={apex};T={triangle};empty=1\n"
            ],
        )
        return {
            "body": body,
            "apex_rank": apex_rank,
            "apex": apex,
            "triangle": triangle,
            "closed": False,
            "empty": True,
            "mass": mass,
            "finite_max": mass,
            "tail_cap": F(0),
            "cap": mass,
            "margin": F(0),
            "horizon": 0,
            "maximizer": None,
            "extensions": (),
            "controls": 0,
            "digest": local_digest,
        }
    components = direct_r
    # For w >= horizon+1 the strict discrepancy bound is below |R|.
    tail_threshold = S2 * components / (6 * mass)
    horizon = max(
        FIRST_EXTERNAL - 1,
        tail_threshold.numerator // tail_threshold.denominator,
    )
    tail_cap = mass / 7 + S2 * components / (7 * (horizon + 1))
    require(
        tail_cap < mass,
        f"one-comb tail did not seal: {body}, {apex}, {triangle}",
    )
    excluded = {apex, *triangle}
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, horizon + 1)
        if speed not in excluded
    ]
    coverage_rows = T.THM2885.coverages_many(residual, speeds)
    require(
        all(F(0) <= value <= mass for value, _ in coverage_rows),
        f"triangle coverage escaped residual mass: {family}",
    )
    by_speed = {speed: value for value, speed in coverage_rows}
    controls = 0
    if speeds:
        for speed in dict.fromkeys((speeds[0], speeds[-1])):
            require(
                by_speed[speed] == T.THM2885.coverage(residual, speed),
                f"triangle vector/scalar control failed: {family}, {speed}",
            )
            controls += 1
    ranked = sorted(coverage_rows, key=lambda item: (-item[0], item[1]))
    finite_max, maximizer = ranked[0] if ranked else (F(0), None)
    cap = max(finite_max, tail_cap)
    extensions = tuple(
        speed for value, speed in coverage_rows if value == mass
    )
    closed = cap < mass
    require(
        closed == (not extensions),
        f"triangle cap/witness mismatch: {body}, {apex}, {triangle}",
    )
    margin = mass - cap
    local_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/extension-local/v1\n",
        [
            "key="
            + ",".join(map(str, body))
            + f";rank={apex_rank};a={apex};"
            + f"T={triangle[0]},{triangle[1]},{triangle[2]};"
            + f"L={ftext(mass)};r={components};W={horizon};"
            + f"finite={ftext(finite_max)};tail={ftext(tail_cap)};"
            + f"cap={ftext(cap)};margin={ftext(margin)};"
            + f"max={maximizer};extensions={','.join(map(str, extensions))}\n"
        ],
    )
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "triangle": triangle,
        "closed": closed,
        "empty": False,
        "mass": mass,
        "finite_max": finite_max,
        "tail_cap": tail_cap,
        "cap": cap,
        "margin": margin,
        "horizon": horizon,
        "maximizer": maximizer,
        "extensions": extensions,
        "controls": controls,
        "digest": local_digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, os.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    context = mp.get_context("spawn")

    if args.workers == 1:
        nested = list(map(profile_atlas_body, T.BODIES))
    else:
        with context.Pool(args.workers) as pool:
            nested = pool.map(profile_atlas_body, T.BODIES, chunksize=1)
    rows = [row for body_rows in nested for row in body_rows]
    require(
        len(rows) == 30_030
        and tuple(row["body"] for row in rows[::10]) == T.BODIES,
        "THM-2888 atlas order changed",
    )
    scalar_counts = tuple(
        sum(row["scalar_class"] == name for row in rows)
        for name in ("direct", "rank3", "failure")
    )
    terminal_count = sum(terminal(row) for row in rows)
    require(
        (*scalar_counts, sum(row["pairpair_direct"] for row in rows), terminal_count)
        == EXPECTED_ATLAS_COUNTS,
        "THM-2888 atlas census changed",
    )

    root_metadata: dict[tuple[int, ...], dict[str, object]] = {}
    target_rows: list[dict[str, object]] = []
    active_roots = 0
    initial_closed = 0
    initial_residual = 0
    for body_index, body in enumerate(T.BODIES):
        if max(body) < 13:
            continue
        active_roots += 1
        body_rows = rows[10 * body_index : 10 * body_index + 10]
        root = T.THM2885.profile_body(body)
        terminals = {row["apex"] for row in body_rows if terminal(row)}
        neg_sum, top5 = complement_margin(root["top15"], terminals)
        margin = root["m"] + neg_sum
        closed = margin > 0
        initial_closed += closed
        initial_residual += not closed
        root_metadata[body] = {
            "m": root["m"],
            "top15": root["top15"],
            "initial_terminals": terminals,
            "initial_margin": margin,
            "initial_top5": top5,
            "initial_closed": closed,
        }
        if not closed:
            target_rows.extend(row for row in body_rows if not terminal(row))
    require(
        (active_roots, initial_closed, initial_residual, len(target_rows))
        == EXPECTED_ROOT_COUNTS,
        "honest post-gate target universe changed",
    )
    normalization_sources = [
        row
        for row in target_rows
        if row["scalar_class"] == "failure"
        and row["cutoff"] is not None
        and row["cutoff"] > T.PAIR_HORIZON
    ]
    normalization_keys = {
        (row["body"], row["apex"]) for row in normalization_sources
    }
    unnormalized_high_cutoff = [
        row
        for row in target_rows
        if (row["body"], row["apex"]) not in normalization_keys
        and row["cutoff"] is not None
        and row["cutoff"] > T.PAIR_HORIZON
    ]
    require(
        not unnormalized_high_cutoff,
        "high-cutoff target escaped heavy-edge normalization",
    )
    if args.workers == 1:
        normalized = list(map(normalize_atlas_row, normalization_sources))
    else:
        with context.Pool(args.workers) as pool:
            normalized = pool.map(
                normalize_atlas_row,
                normalization_sources,
                chunksize=1,
            )
    normalized_by_key = {
        (row["body"], row["apex"]): row
        for row in normalized
    }
    require(
        len(normalized_by_key) == len(normalization_sources)
        and all(row["forced"]["margin"] > 0 for row in normalized),
        "heavy-edge normalization changed",
    )
    normalization_counts = (
        len(normalized),
        max((row["old_cutoff"] for row in normalized), default=0),
        max((row["deleted_cutoff"] for row in normalized), default=0),
        sum(row["forced"]["paid"] for row in normalized),
        sum(row["deleted_paid"] for row in normalized),
        max(
            (
                row["cutoff"]
                for row in target_rows
                if (row["body"], row["apex"]) not in normalization_keys
            ),
            default=0,
        ),
    )
    normalization_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/normalizations/v1\n",
        [
            ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"edge={row['edge'][0]},{row['edge'][1]};"
            + f"forced_margin={ftext(row['forced']['margin'])};"
            + f"forced_sha={row['forced']['paid_digest']};"
            + f"deleted_B={ftext(row['deleted_cap'])};"
            + f"deleted_W={row['deleted_cutoff']};"
            + f"deleted_sha={row['deleted_digest']}\n"
            for row in normalized
        ],
    )
    if EXPECTED_NORMALIZATION_COUNTS is not None:
        require(
            normalization_counts == EXPECTED_NORMALIZATION_COUNTS
            and normalization_digest == EXPECTED_NORMALIZATION_DIGEST,
            "target heavy-edge normalization changed",
        )
    targets = []
    for row in target_rows:
        normalized_row = normalized_by_key.get((row["body"], row["apex"]))
        targets.append(
            (
                row["body"],
                row["apex_rank"],
                row["apex"],
                row["m"],
                row["r"],
                row["global_cap"]
                if normalized_row is None
                else normalized_row["deleted_cap"],
                None if normalized_row is None else normalized_row["edge"],
                None
                if normalized_row is None
                else normalized_row["forced"]["margin"],
            )
        )
    target_lines = [
        ",".join(map(str, body))
        + f";rank={rank};a={apex};h={ftext(mass)};"
        + f"r={components};B={ftext(pair_cap)};"
        + f"excluded={excluded_edge};forced_margin="
        + ("NA" if forced_margin is None else ftext(forced_margin))
        + "\n"
        for (
            body,
            rank,
            apex,
            mass,
            components,
            pair_cap,
            excluded_edge,
            forced_margin,
        ) in targets
    ]
    target_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/targets/v1\n",
        target_lines,
    )
    if EXPECTED_TARGET_DIGEST is not None:
        require(target_digest == EXPECTED_TARGET_DIGEST, "target digest changed")

    def predicted_tail_first(task):
        mass = task[3]
        components = task[4]
        theta = mass - task[5]
        level = theta / 2
        return ceiling(S2 * components / (7 * (level - mass / 7)))

    work_targets = sorted(
        targets,
        key=lambda task: (
            -predicted_tail_first(task),
            task[0],
            task[1],
            task[2],
        ),
    )
    if args.workers == 1:
        coarse_work = list(map(profile_coarse, work_targets))
    else:
        with context.Pool(args.workers) as pool:
            coarse_work = pool.map(profile_coarse, work_targets, chunksize=1)
    coarse_by_key = {
        (row["body"], row["apex"]): row
        for row in coarse_work
    }
    require(len(coarse_by_key) == len(targets), "duplicate coarse carrier key")
    coarse = [
        coarse_by_key[(body, apex)]
        for body, _, apex, _, _, _, _, _ in targets
    ]
    require(
        tuple((row["body"], row["apex_rank"], row["apex"]) for row in coarse)
        == tuple(
            (body, rank, apex)
            for body, rank, apex, _, _, _, _, _ in targets
        ),
        "coarse carrier order changed",
    )
    coarse_counts = (
        sum(row["head_size"] for row in coarse),
        max(row["head_size"] for row in coarse),
        sum(row["edge_count"] for row in coarse),
        max(row["edge_count"] for row in coarse),
        sum(row["triangle_count"] for row in coarse),
        max(row["triangle_count"] for row in coarse),
        sum(row["triangle_count"] > 0 for row in coarse),
        sum(row["k4_count"] for row in coarse),
        max(row["k4_count"] for row in coarse),
        sum(row["controls"] for row in coarse),
        max(row["tail_first"] for row in coarse),
        sum(row["triple_controls"] for row in coarse),
        sum(row["normalized"] for row in coarse),
        sum(row["normalized_pair_direct"] for row in coarse),
    )
    maximum_head = max(
        (row["head_size"], row["body"], row["apex_rank"], row["apex"])
        for row in coarse
    )
    maximum_edges = max(
        (row["edge_count"], row["body"], row["apex_rank"], row["apex"])
        for row in coarse
    )
    maximum_triangles = max(
        (row["triangle_count"], row["body"], row["apex_rank"], row["apex"])
        for row in coarse
    )
    maximum_k4 = max(
        (row["k4_count"], row["body"], row["apex_rank"], row["apex"])
        for row in coarse
    )
    maximum_h_tail = max(
        (row["tail_first"], row["body"], row["apex_rank"], row["apex"])
        for row in coarse
    )
    coarse_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/coarse-aggregate/v1\n",
        [row["digest"] + "\n" for row in coarse],
    )
    if EXPECTED_COARSE_COUNTS is not None:
        require(
            coarse_counts == EXPECTED_COARSE_COUNTS
            and coarse_digest == EXPECTED_COARSE_DIGEST,
            "coarse heavy-triangle census changed",
        )

    triangle_tasks = [
        (
            row["body"],
            row["apex_rank"],
            row["apex"],
            triangle,
            residual_mass,
        )
        for row in coarse
        for triangle, residual_mass in zip(
            row["triangles"],
            row["triangle_masses"],
        )
    ]
    if args.workers == 1:
        extensions = list(map(seal_triangle_extension, triangle_tasks))
    else:
        with context.Pool(args.workers) as pool:
            extensions = pool.map(
                seal_triangle_extension,
                triangle_tasks,
                chunksize=4,
            )
    require(
        tuple(
            (row["body"], row["apex_rank"], row["apex"], row["triangle"])
            for row in extensions
        )
        == tuple((body, rank, apex, triangle) for body, rank, apex, triangle, _ in triangle_tasks),
        "triangle extension order changed",
    )
    closed_extensions = [row for row in extensions if row["closed"]]
    minimum_margin_row = min(
        closed_extensions,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["triangle"],
            row["maximizer"],
        ),
    ) if closed_extensions else None
    minimum_margin = (
        (
            minimum_margin_row["margin"],
            minimum_margin_row["body"],
            minimum_margin_row["apex_rank"],
            minimum_margin_row["apex"],
            minimum_margin_row["triangle"],
            minimum_margin_row["maximizer"],
        )
        if minimum_margin_row is not None
        else None
    )
    triangle_counts = (
        len(extensions),
        len(closed_extensions),
        sum(not row["closed"] for row in extensions),
        sum(row["empty"] for row in extensions),
        sum(row["controls"] for row in extensions),
        max((row["horizon"] for row in extensions), default=0),
    )
    maximum_scalar_horizon = max(
        (
            row["horizon"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["triangle"],
        )
        for row in extensions
    ) if extensions else None
    triangle_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/extensions-aggregate/v1\n",
        [row["digest"] + "\n" for row in extensions],
    )
    if EXPECTED_TRIANGLE_COUNTS is not None:
        require(
            triangle_counts == EXPECTED_TRIANGLE_COUNTS
            and triangle_digest == EXPECTED_TRIANGLE_DIGEST,
            "triangle extension census changed",
        )
    require(
        minimum_margin == EXPECTED_MINIMUM_TRIANGLE_MARGIN,
        "minimum triangle margin changed",
    )

    hyperedge_sources: dict[
        tuple[tuple[int, ...], tuple[int, ...]],
        list[tuple[int, tuple[int, int, int], int | None]],
    ] = {}
    target_apex_keys = {
        (body, apex)
        for body, _, apex, _, _, _, _, _ in targets
    }
    raw_hyperedge_sources = 0
    for row in extensions:
        if row["closed"]:
            continue
        if row["empty"]:
            core = {row["apex"], *row["triangle"]}
            fourth = FIRST_EXTERNAL
            while fourth in core:
                fourth += 1
            extension_speeds: tuple[int | None, ...] = (fourth,)
        else:
            extension_speeds = row["extensions"]
        for fourth in extension_speeds:
            require(fourth is not None, "missing canonical fourth speed")
            speeds = tuple(sorted((row["apex"], *row["triangle"], fourth)))
            require(
                len(speeds) == 5 and len(set(speeds)) == 5,
                f"bad exact-cover hyperedge: {row['body']}, {speeds}",
            )
            key = (row["body"], speeds)
            hyperedge_sources.setdefault(key, []).append(
                (row["apex"], row["triangle"], fourth)
            )
            raw_hyperedge_sources += 1
    hyperedge_lines: list[str] = []
    for (body, speeds), sources in sorted(hyperedge_sources.items()):
        full_good, _, full_mass = T.CORE.good_norm(
            tuple(sorted((*body, *speeds)))
        )
        require(
            full_mass == 0 and not full_good,
            f"claimed hyperedge has positive literal survivor: {body}, {speeds}",
        )
        top10 = {
            speed
            for _, speed in T.THM2885.profile_body(body)["top15"][:10]
        }
        incidence = tuple(speed for speed in speeds if speed in top10)
        source_apices = tuple(sorted({source[0] for source in sources}))
        require(
            not root_metadata[body]["initial_closed"]
            and incidence
            and set(source_apices) == set(incidence)
            and all((body, apex) in target_apex_keys for apex in incidence),
            f"hyperedge misses its ranked-apex incidence: {body}, {speeds}",
        )
        hyperedge_lines.append(
            ",".join(map(str, body))
            + ";S="
            + ",".join(map(str, speeds))
            + ";top10="
            + ",".join(map(str, incidence))
            + ";source_apices="
            + ",".join(map(str, source_apices))
            + f";multiplicity={len(sources)}\n"
        )
    hyperedge_counts = (
        raw_hyperedge_sources,
        len(hyperedge_sources),
        max((len(sources) for sources in hyperedge_sources.values()), default=0),
        sum(body_stratum(body) == "j6_one_high" for body, _ in hyperedge_sources),
        sum(body_stratum(body) == "j7_both_high" for body, _ in hyperedge_sources),
    )
    hyperedge_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/exact-cover-hyperedges/v1\n",
        hyperedge_lines,
    )
    if EXPECTED_HYPEREDGE_COUNTS is not None:
        require(
            hyperedge_counts == EXPECTED_HYPEREDGE_COUNTS
            and hyperedge_digest == EXPECTED_HYPEREDGE_DIGEST,
            "exact-cover hypergraph changed",
        )

    extension_by_carrier: dict[
        tuple[tuple[int, ...], int],
        list[dict[str, object]],
    ] = {}
    for row in extensions:
        extension_by_carrier.setdefault((row["body"], row["apex"]), []).append(row)
    carrier_lines: list[str] = []
    closed_carriers: set[tuple[tuple[int, ...], int]] = set()
    unresolved_carriers: list[
        tuple[tuple[int, ...], int, int, tuple[int, int, int] | None]
    ] = []
    triangle_free_carriers = 0
    for row in coarse:
        key = (row["body"], row["apex"])
        branch_rows = extension_by_carrier.get(key, [])
        require(
            tuple(extension["triangle"] for extension in branch_rows)
            == row["triangles"],
            f"carrier triangle block changed: {key}",
        )
        closed = row["forced_closed"] and all(
            extension["closed"] for extension in branch_rows
        )
        triangle_free_carriers += not branch_rows
        if closed:
            closed_carriers.add(key)
            witness = None
        else:
            witness = next(
                extension["triangle"]
                for extension in branch_rows
                if not extension["closed"]
            )
            unresolved_carriers.append(
                (row["body"], row["apex_rank"], row["apex"], witness)
            )
        carrier_lines.append(
            ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"triangles={len(branch_rows)};closed={int(closed)};"
            + f"forced_closed={int(row['forced_closed'])};"
            + f"first_survivor={witness}\n"
        )
    carrier_counts = (
        len(coarse),
        triangle_free_carriers,
        len(closed_carriers) - triangle_free_carriers,
        len(closed_carriers),
        len(unresolved_carriers),
    )
    carrier_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/carriers/v1\n",
        carrier_lines,
    )
    if EXPECTED_CARRIER_COUNTS is not None:
        require(
            carrier_counts == EXPECTED_CARRIER_COUNTS
            and carrier_digest == EXPECTED_CARRIER_DIGEST,
            "carrier closure census changed",
        )

    final_closed = 0
    final_residual: list[tuple[int, ...]] = []
    final_lines: list[str] = []
    final_margins: list[tuple[F, tuple[int, ...]]] = []
    final_status: dict[tuple[int, ...], bool] = {}
    closed_by_body: dict[tuple[int, ...], set[int]] = {}
    for body, apex in closed_carriers:
        closed_by_body.setdefault(body, set()).add(apex)
    for body, metadata in root_metadata.items():
        terminals = set(metadata["initial_terminals"])
        terminals.update(closed_by_body.get(body, set()))
        neg_sum, top5 = complement_margin(metadata["top15"], terminals)
        margin = metadata["m"] + neg_sum
        closed = margin > 0
        require(
            not metadata["initial_closed"] or closed,
            f"initially closed root regressed: {body}",
        )
        final_status[body] = closed
        final_closed += closed
        if closed:
            final_margins.append((margin, body))
        else:
            final_residual.append(body)
        final_lines.append(
            ",".join(map(str, body))
            + ";terminals="
            + ",".join(map(str, sorted(terminals)))
            + f";margin={ftext(margin)};top5="
            + ",".join(
                f"{speed}:{ftext(value)}"
                for value, speed in top5
            )
            + f";closed={int(closed)}\n"
        )
    minimum_root_margin = min(final_margins) if final_margins else None
    require(
        minimum_root_margin == EXPECTED_MINIMUM_ROOT_MARGIN,
        "minimum final root margin changed",
    )
    final_root_counts = (
        len(root_metadata),
        final_closed,
        len(final_residual),
    )
    root_digest = digest_text(
        "LRC14/j5/postgate-heavy-triangle/final-roots/v1\n",
        final_lines,
    )
    if EXPECTED_FINAL_ROOT_COUNTS is not None:
        require(
            final_root_counts == EXPECTED_FINAL_ROOT_COUNTS
            and root_digest == EXPECTED_ROOT_DIGEST,
            "final weighted root census changed",
        )

    stratum_counts_list: list[int] = []
    stratum_summaries: list[tuple[str, tuple[int, ...]]] = []
    for label in ("j6_one_high", "j7_both_high"):
        bodies = {
            body for body in root_metadata if body_stratum(body) == label
        }
        coarse_rows = [row for row in coarse if row["body"] in bodies]
        extension_rows = [row for row in extensions if row["body"] in bodies]
        target_keys = {
            (body, apex)
            for body, _, apex, _, _, _, _, _ in targets
            if body in bodies
        }
        counts = (
            len(bodies),
            sum(root_metadata[body]["initial_closed"] for body in bodies),
            sum(not root_metadata[body]["initial_closed"] for body in bodies),
            len(target_keys),
            sum(row["head_size"] for row in coarse_rows),
            sum(row["edge_count"] for row in coarse_rows),
            sum(row["triangle_count"] for row in coarse_rows),
            sum(row["closed"] for row in extension_rows),
            sum(not row["closed"] for row in extension_rows),
            sum(key in closed_carriers for key in target_keys),
            sum(key not in closed_carriers for key in target_keys),
            sum(final_status[body] for body in bodies),
            sum(not final_status[body] for body in bodies),
            sum(
                hyperedge_body in bodies
                for hyperedge_body, _ in hyperedge_sources
            ),
        )
        stratum_summaries.append((label, counts))
        stratum_counts_list.extend(counts)
    stratum_counts = tuple(stratum_counts_list)
    if EXPECTED_STRATUM_COUNTS is not None:
        require(
            stratum_counts == EXPECTED_STRATUM_COUNTS,
            "j6/j7 stratum census changed",
        )

    print("LRC14 J=5 POST-GATE HEAVY-TRIANGLE CENSUS")
    print("status=FINITE-EXACT+GLOBAL-SCALAR-TAIL-SEALED")
    print(
        "universe=active_roots:2508;"
        "honest_initial_closed:1064;residual_roots:1444;"
        f"residual_nonterminal_apex_branches:{len(targets)}"
    )
    print(f"target_digest_sha256={target_digest}")
    print(
        f"heavy_edge_normalizations={normalization_counts[0]};"
        f"maximum_old_cutoff={normalization_counts[1]};"
        f"maximum_deleted_cutoff={normalization_counts[2]};"
        f"forced_pair_checks={normalization_counts[3]};"
        f"deleted_pair_checks={normalization_counts[4]};"
        f"maximum_unnormalized_cutoff={normalization_counts[5]};"
        f"normalization_digest_sha256={normalization_digest}"
    )
    print(
        f"finite_H_vertices_total={coarse_counts[0]};"
        f"max_per_carrier={coarse_counts[1]};"
        f"heavy_edges_total={coarse_counts[2]};"
        f"max_per_carrier={coarse_counts[3]};"
        f"heavy_triangles_total={coarse_counts[4]};"
        f"max_per_carrier={coarse_counts[5]};"
        f"triangle_carriers={coarse_counts[6]};"
        f"H_only_K4_total={coarse_counts[7]};"
        f"max_per_carrier={coarse_counts[8]}"
    )
    print(
        f"maximum_H={maximum_head};maximum_edges={maximum_edges};"
        f"maximum_triangles={maximum_triangles};maximum_H_only_K4={maximum_k4};"
        f"maximum_H_tail={maximum_h_tail}"
    )
    print(
        f"H_vector_scalar_controls={coarse_counts[9]};"
        f"maximum_H_tail_first={coarse_counts[10]};"
        f"triple_order_controls={coarse_counts[11]};"
        f"normalized_carriers={coarse_counts[12]};"
        f"normalized_pair_direct={coarse_counts[13]};"
        f"coarse_digest_sha256={coarse_digest}"
    )
    print(
        f"literal_triangle_residuals={triangle_counts[0]};"
        f"one_comb_closed={triangle_counts[1]};"
        f"extendible={triangle_counts[2]};"
        f"empty_triples={triangle_counts[3]};"
        f"triangle_vector_scalar_controls={triangle_counts[4]};"
        f"maximum_scalar_horizon={maximum_scalar_horizon}"
    )
    print(
        "minimum_triangle_margin="
        + (
            "NA"
            if minimum_margin is None
            else (
                f"{ftext(minimum_margin[0])};body={minimum_margin[1]};"
                f"rank={minimum_margin[2]};apex={minimum_margin[3]};"
                f"triangle={minimum_margin[4]};"
                f"finite_head_maximizer={minimum_margin[5]};"
                f"residual_mass={ftext(minimum_margin_row['mass'])};"
                f"finite_max={ftext(minimum_margin_row['finite_max'])};"
                f"tail_cap={ftext(minimum_margin_row['tail_cap'])};"
                f"cap_source="
                + (
                    "tail"
                    if minimum_margin_row["tail_cap"]
                    > minimum_margin_row["finite_max"]
                    else "finite"
                )
            )
        )
    )
    print(f"triangle_digest_sha256={triangle_digest}")
    print(
        f"exact_cover_hypergraph_raw={hyperedge_counts[0]};"
        f"unique={hyperedge_counts[1]};"
        f"maximum_multiplicity={hyperedge_counts[2]};"
        f"j6={hyperedge_counts[3]};j7={hyperedge_counts[4]};"
        f"hyperedge_digest_sha256={hyperedge_digest}"
    )
    print(
        f"carrier_tier0_triangle_free={carrier_counts[1]};"
        f"tier1_literal_scalar={carrier_counts[2]};"
        f"closed={carrier_counts[3]}/{carrier_counts[0]};"
        f"unresolved={carrier_counts[4]};"
        f"first_unresolved={unresolved_carriers[0] if unresolved_carriers else None}"
    )
    print(f"carrier_digest_sha256={carrier_digest}")
    print(
        f"final_active_weighted_gate={final_root_counts[1]}/"
        f"{final_root_counts[0]};residual={final_root_counts[2]};"
        "minimum_positive_margin="
        + (
            "NA"
            if minimum_root_margin is None
            else f"{ftext(minimum_root_margin[0])};body={minimum_root_margin[1]}"
        )
    )
    print(f"root_digest_sha256={root_digest}")
    for label, counts in stratum_summaries:
        print(
            f"{label}=roots:{counts[0]},"
            f"initial_closed:{counts[1]},initial_residual:{counts[2]},"
            f"target_apices:{counts[3]},H_vertices:{counts[4]},"
            f"heavy_edges:{counts[5]},heavy_triangles:{counts[6]},"
            f"triangle_residuals_closed:{counts[7]},extendible:{counts[8]},"
            f"carriers_closed:{counts[9]},carriers_unresolved:{counts[10]},"
            f"final_closed:{counts[11]},final_residual:{counts[12]},"
            f"exact_hyperedges:{counts[13]}"
        )
    print(
        "scope="
        + (
            "complete_THM2885/2888_eight_body_five_slot_rung"
            if final_root_counts[2] == 0
            else "partial_THM2885/2888_eight_body_five_slot_rung"
        )
        + ";495_THM885_low_bodies_separate;"
        "v8_ge_15_at_most_seven_in_window_unresolved;"
        "not_unrestricted_LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
