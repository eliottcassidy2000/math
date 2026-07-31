#!/usr/bin/env python3
r"""Scoped exact feasibility battery for the seven-body / six-slot LRC rung.

This is a research probe, not a uniform theorem.  It tests two successor
ideas suggested by THM-2888.

First, for a six-speed external set, a root hitting gate retaining ``K``
apices is valid when global coverage ranks ``K+1,...,K+6`` sum to less than
the seven-body good-set mass.  A deterministic 36-root stratified sample
tests the proposed fixed value ``K=12`` and finds the least valid adaptive
``K`` through 24.  Global top thirty ranks are sealed by the strict
THM-735 discrepancy tail.

Second, fixing one root apex leaves five external speeds, not four.  A
two-comb cap ``B2<5h/7`` is therefore not enough.  Let ``q1`` be the largest
single-comb coverage and let ``T3`` be a global triple-union cap.  The cheap
bound is ``T3<=B2+q1``; where it is inconclusive, this verifier computes the
exact finite-head triple maximum and seals the tail by

    U(x,y,z) < B2 + h/7 + (99/70) r/(7*2501)

whenever one endpoint exceeds 2500.

THM-2893 (the complement-cap finite-core flag lemma) supplies the proof
frame.  If ``T3<5h/7``, put ``theta=h-T3>2h/7`` and
``H={x:c(x)>=theta/2}``.  A five-cover makes every pair theta-heavy: its
complementary triple has union at most ``T3``.  Thus it induces a heavy K5.
Since H is a vertex cover of the heavy graph, at least four vertices lie in
H.  In particular it contains a heavy H-triangle.  For three selected
hostile root/apex cases, every such triangle residual is tested against a
global exact two-comb cap.

The same cases also test THM-2893's alternative ``(k,s)=(5,3)`` flag.  If
the stronger pair inequality ``B2<4h/7`` holds, triple edges are declared
heavy at ``h-B2`` and the high core uses threshold ``(h-B2)/3``.  Every
heavy high triple again leaves a residual pair problem.

The 792 seven-body roots contained in {1,...,12} are *not* inherited
terminals: THM-885's j=6 sweep is explicitly unfinished.  Even a successful
uniform seven-body computation would only close the seven-body/six-slot
near-AP rung, not unrestricted LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
ROOT_BASE_HORIZON = 1600
ROOT_TOP_COUNT = 30
PAIR_HORIZON = 2500
SAMPLE_PER_STRATUM = 12
S2 = F(99, 70)

BATTERY_ROOTS = (
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
)

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_SAMPLE_COUNTS: tuple[tuple[str, int, int, int], ...] | None = (
    ("low", 12, 10, 16),
    ("one", 12, 5, 20),
    ("both", 12, 5, 21),
)
EXPECTED_SAMPLE_DIGEST: str | None = (
    "ce3be5476ee287f7d972d1853f4e9fe52175583347832a75f103fd6e9172ba51"
)
EXPECTED_BATTERY_COUNTS: tuple[int, ...] | None = (
    3,
    89,
    54,
    957,
    6_106,
    6_106,
    49_472,
    69,
    33,
    7_700,
    2_207,
    1_072,
    25_143,
)
EXPECTED_BATTERY_DIGEST: str | None = (
    "865478ea6969df9633ae2d08346f90f93dbac10d4d945994a07988202d4c8cc1"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_atlas():
    require(file_sha256(ATLAS_PATH) == ATLAS_SHA256, "THM-2888 script changed")
    require(
        file_sha256(ATLAS_OUTPUT_PATH) == ATLAS_OUTPUT_SHA256,
        "THM-2888 output changed",
    )
    spec = importlib.util.spec_from_file_location("j6_probe_thm2888", ATLAS_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-2888")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_atlas()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def digest_text(header: str, rows: list[str]) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def quantile_sample(
    rows: list[tuple[int, ...]],
    count: int,
) -> tuple[tuple[int, ...], ...]:
    require(count >= 2 and len(rows) >= count, "bad quantile sample")
    return tuple(
        rows[(len(rows) - 1) * index // (count - 1)]
        for index in range(count)
    )


def global_top(
    body: tuple[int, ...],
    count: int = ROOT_TOP_COUNT,
) -> dict[str, object]:
    """Seal exact global root coverages through the requested rank."""

    good, components, mass = T.CORE.good_norm(body)
    rows = T.THM2885.coverages_many(
        good,
        range(FIRST_EXTERNAL, ROOT_BASE_HORIZON + 1),
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q_count = ranked[count - 1][0]
    require(q_count > mass / 7, f"rank {count} misses limit: {body}")
    threshold = S2 * components / (7 * (q_count - mass / 7))
    threshold_first = ceiling(threshold)
    tail_first = max(ROOT_BASE_HORIZON + 1, threshold_first)
    if tail_first > ROOT_BASE_HORIZON + 1:
        rows.extend(
            T.THM2885.coverages_many(
                good,
                range(ROOT_BASE_HORIZON + 1, tail_first),
            )
        )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q_count,
        f"rank {count} tail did not seal: {body}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    by_speed = {speed: value for value, speed in rows}
    controls = 0
    for speed in dict.fromkeys(
        (
            ranked[0][1],
            ranked[count - 1][1],
            FIRST_EXTERNAL,
            tail_first - 1,
        )
    ):
        if speed >= FIRST_EXTERNAL:
            require(
                by_speed[speed] == T.THM2885.coverage(good, speed),
                f"root vector/scalar mismatch: {body}, {speed}",
            )
            controls += 1
    return {
        "body": body,
        "good": good,
        "m": mass,
        "r": components,
        "top": tuple(ranked[:count]),
        "threshold_first": threshold_first,
        "tail_first": tail_first,
        "controls": controls,
    }


def sample_universes() -> tuple[tuple[str, tuple[tuple[int, ...], ...]], ...]:
    low = list(combinations(range(1, 13), 7))
    one = [
        tuple(sorted((*base, special)))
        for base in combinations(range(1, 13), 6)
        for special in (13, 14)
    ]
    both = [
        (*base, 13, 14)
        for base in combinations(range(1, 13), 5)
    ]
    require(
        (len(low), len(one), len(both)) == (792, 1848, 792),
        "seven-body strata changed",
    )
    return (
        ("low", quantile_sample(low, SAMPLE_PER_STRATUM)),
        ("one", quantile_sample(one, SAMPLE_PER_STRATUM)),
        ("both", quantile_sample(both, SAMPLE_PER_STRATUM)),
    )


def root_gate_sample() -> tuple[list[dict[str, object]], str]:
    rows: list[dict[str, object]] = []
    digest_rows: list[str] = []
    for stratum, bodies in sample_universes():
        for body in bodies:
            root = global_top(body)
            top = root["top"]
            fixed_margin = root["m"] - sum(
                (value for value, _ in top[12:18]),
                F(0),
            )
            adaptive_k = None
            adaptive_margin = None
            for k in range(0, ROOT_TOP_COUNT - 5):
                margin = root["m"] - sum(
                    (value for value, _ in top[k : k + 6]),
                    F(0),
                )
                if margin > 0:
                    adaptive_k = k
                    adaptive_margin = margin
                    break
            require(
                adaptive_k is not None and adaptive_margin is not None,
                f"top-{ROOT_TOP_COUNT} has no hitting gate: {body}",
            )
            row = {
                "stratum": stratum,
                "body": body,
                "fixed_margin": fixed_margin,
                "adaptive_k": adaptive_k,
                "adaptive_margin": adaptive_margin,
                "threshold_first": root["threshold_first"],
                "tail_first": root["tail_first"],
                "controls": root["controls"],
            }
            rows.append(row)
            digest_rows.append(
                f"S={stratum};E={','.join(map(str, body))};"
                f"K12={ftext(fixed_margin)};K={adaptive_k};"
                f"margin={ftext(adaptive_margin)};"
                f"threshold={root['threshold_first']};"
                f"tail={root['tail_first']}\n"
            )
    return rows, digest_text(
        "LRC14/j6/adaptive-root-gate-sample/v1\n",
        digest_rows,
    )


def exact_triple_cap(
    body: tuple[int, ...],
    apex_rank: int,
    apex: int,
    root_good: list[tuple[F, F]],
) -> dict[str, object]:
    """Compute a global triple-union cap on one first-apex carrier."""

    pair_row = T.profile_apex(body, root_good, apex_rank, apex)
    carrier = T.THM2883.subtract_local(root_good, apex)
    mass = interval_mass(carrier)
    components = len(carrier)
    direct_good, direct_r, direct_m = T.CORE.good_norm((*body, apex))
    require(
        carrier == direct_good
        and components == direct_r
        and mass == direct_m == pair_row["m"],
        f"first-apex reconstruction changed: {body}, {apex}",
    )
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed != apex
    ]
    coverage_rows = T.THM2885.coverages_many(carrier, speeds)
    by_speed = {speed: value for value, speed in coverage_rows}
    ranked = sorted(coverage_rows, key=lambda item: (-item[0], item[1]))
    require(
        tuple(ranked[:4]) == pair_row["top4"],
        f"pair/triple scalar head changed: {body}, {apex}",
    )
    pair_cap = pair_row["global_cap"]
    require(
        pair_row["margin"] > 0,
        f"pair cap misses 5h/7 in battery: {body}, {apex}",
    )
    tail_single = mass / 7 + S2 * components / (
        7 * (PAIR_HORIZON + 1)
    )
    tail_cap = pair_cap + tail_single
    head_cap = F(0)
    maximizing_triple: tuple[int, int, int] | None = None
    paid = 0
    paid_hash = hashlib.sha256(
        b"LRC14/j6/exact-triple-paid/v1\n"
    )
    for first_index, (first_value, first) in enumerate(ranked[:-2]):
        if (
            first_value
            + ranked[first_index + 1][0]
            + ranked[first_index + 2][0]
            <= head_cap
        ):
            break
        after_first = T.THM2883.subtract_local(carrier, first)
        for second_index in range(first_index + 1, len(ranked) - 1):
            second_value, second = ranked[second_index]
            if (
                first_value
                + second_value
                + ranked[second_index + 1][0]
                <= head_cap
            ):
                break
            after_pair = T.THM2883.subtract_local(after_first, second)
            for third_value, third in ranked[second_index + 1 :]:
                if first_value + second_value + third_value <= head_cap:
                    break
                after_triple = T.THM2883.subtract_local(after_pair, third)
                union = mass - interval_mass(after_triple)
                triple = tuple(sorted((first, second, third)))
                paid_hash.update(
                    (
                        f"T={triple[0]},{triple[1]},{triple[2]};"
                        f"U={ftext(union)}\n"
                    ).encode()
                )
                paid += 1
                if union > head_cap:
                    head_cap = union
                    maximizing_triple = triple
    require(
        paid > 0 and maximizing_triple is not None,
        f"empty triple maximization: {body}, {apex}",
    )
    reverse = T.THM2883.subtract_local_multi(
        carrier,
        tuple(reversed(maximizing_triple)),
    )
    forward = T.THM2883.subtract_local_multi(
        carrier,
        maximizing_triple,
    )
    family = tuple(sorted((*body, apex, *maximizing_triple)))
    family_good, family_r, family_m = T.CORE.good_norm(family)
    require(
        forward == reverse == family_good
        and family_r == len(forward)
        and mass - family_m == head_cap,
        f"triple direct control failed: {family}",
    )
    global_cap = max(head_cap, tail_cap)
    margin = F(5, 7) * mass - global_cap
    require(margin > 0, f"triple cap misses 5h/7: {body}, {apex}")
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "carrier": carrier,
        "m": mass,
        "r": components,
        "q1": ranked[0][0],
        "pair_cap": pair_cap,
        "cheap_margin": F(5, 7) * mass - (pair_cap + ranked[0][0]),
        "head_cap": head_cap,
        "tail_cap": tail_cap,
        "global_cap": global_cap,
        "margin": margin,
        "maximizing_triple": maximizing_triple,
        "paid": paid,
        "paid_digest": paid_hash.hexdigest(),
        "by_speed": by_speed,
    }


def residual_pair_cap(
    body: tuple[int, ...],
    apex: int,
    triangle: tuple[int, int, int],
    carrier: list[tuple[F, F]],
) -> dict[str, object]:
    """Globally cap two-comb coverage of one literal triangle residual."""

    residual = T.THM2883.subtract_local_multi(carrier, triangle)
    mass = interval_mass(residual)
    require(mass > 0, f"heavy triangle already covers: {body}, {apex}, {triangle}")
    components = len(residual)
    family = tuple(sorted((*body, apex, *triangle)))
    direct_good, direct_r, direct_m = T.CORE.good_norm(family)
    require(
        residual == direct_good
        and components == direct_r
        and mass == direct_m,
        f"triangle residual reconstruction changed: {family}",
    )
    excluded = {apex, *triangle}
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed not in excluded
    ]
    coverage_rows = T.THM2885.coverages_many(residual, speeds)
    ranked = sorted(coverage_rows, key=lambda item: (-item[0], item[1]))
    q1 = ranked[0][0]
    tail_single = mass / 7 + S2 * components / (
        7 * (PAIR_HORIZON + 1)
    )
    require(
        tail_single <= q1,
        f"triangle rank-one tail not sealed: {family}",
    )
    pair = T.pair_maximum(
        residual,
        mass,
        ranked,
        family,
    )
    tail_pair = q1 + tail_single
    global_cap = max(pair["head_cap"], tail_pair)
    margin = mass - global_cap
    require(margin > 0, f"triangle residual pair cap fails: {family}")
    return {
        "m": mass,
        "r": components,
        "head_cap": pair["head_cap"],
        "tail_cap": tail_pair,
        "global_cap": global_cap,
        "margin": margin,
        "maximizing_pair": pair["maximizing_pair"],
        "paid": pair["paid"],
        "paid_digest": pair["paid_digest"],
    }


def heavy_triangle_battery(root: dict[str, object]) -> dict[str, object]:
    body = root["body"]
    top = root["top"]
    apex_value, apex = top[0]
    del apex_value
    triple = exact_triple_cap(body, 1, apex, root["good"])
    mass = triple["m"]
    theta = mass - triple["global_cap"]
    require(theta > 2 * mass / 7, f"triple theta too small: {body}, {apex}")
    level = theta / 2
    threshold = S2 * triple["r"] / (7 * (level - mass / 7))
    tail_first = max(FIRST_EXTERNAL, ceiling(threshold))
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    coverage_rows = T.THM2885.coverages_many(triple["carrier"], speeds)
    by_speed = {speed: value for value, speed in coverage_rows}
    require(
        mass / 7 + S2 * triple["r"] / (7 * tail_first) <= level,
        f"H tail did not seal: {body}, {apex}",
    )
    head = tuple(
        sorted(speed for speed, value in by_speed.items() if value >= level)
    )
    edges: dict[tuple[int, int], F] = {}
    for first, second in combinations(head, 2):
        if by_speed[first] + by_speed[second] < theta:
            continue
        after = T.THM2883.subtract_local_multi(
            triple["carrier"],
            (first, second),
        )
        union = mass - interval_mass(after)
        if union >= theta:
            edges[(first, second)] = union
    adjacency = {speed: set() for speed in head}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    triangles = tuple(
        (first, second, third)
        for first, second, third in combinations(head, 3)
        if (
            second in adjacency[first]
            and third in adjacency[first]
            and third in adjacency[second]
        )
    )
    residuals = [
        residual_pair_cap(
            body,
            apex,
            triangle,
            triple["carrier"],
        )
        for triangle in triangles
    ]
    minimum_residual = min(
        (
            row["margin"],
            triangle,
            row["maximizing_pair"],
        )
        for triangle, row in zip(triangles, residuals)
    ) if residuals else None
    return {
        "body": body,
        "adaptive_k": root["adaptive_k"],
        "apex": apex,
        "triple": triple,
        "theta": theta,
        "H_tail": tail_first,
        "H": head,
        "edges": edges,
        "triangles": triangles,
        "residuals": residuals,
        "minimum_residual": minimum_residual,
    }


def strong_pair_flag_battery(
    row: dict[str, object],
) -> dict[str, object]:
    """Test THM-2893's (k,s)=(5,3) route on the same carrier."""

    body = row["body"]
    apex = row["apex"]
    triple = row["triple"]
    carrier = triple["carrier"]
    mass = triple["m"]
    pair_cap = triple["pair_cap"]
    strong_margin = F(4, 7) * mass - pair_cap
    require(
        strong_margin > 0,
        f"strong pair flag misses 4h/7: {body}, {apex}",
    )
    theta = mass - pair_cap
    level = theta / 3
    require(level > mass / 7, f"three-core threshold not finite: {body}, {apex}")
    threshold = S2 * triple["r"] / (7 * (level - mass / 7))
    tail_first = max(FIRST_EXTERNAL, ceiling(threshold))
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    coverage_rows = T.THM2885.coverages_many(carrier, speeds)
    by_speed = {speed: value for value, speed in coverage_rows}
    require(
        mass / 7 + S2 * triple["r"] / (7 * tail_first) <= level,
        f"strong-pair H3 tail did not seal: {body}, {apex}",
    )
    head = tuple(
        sorted(speed for speed, value in by_speed.items() if value >= level)
    )
    heavy: list[tuple[int, int, int]] = []
    candidate_count = 0
    for first, second, third in combinations(head, 3):
        if (
            by_speed[first] + by_speed[second] + by_speed[third]
            < theta
        ):
            continue
        candidate_count += 1
        after = T.THM2883.subtract_local_multi(
            carrier,
            (first, second, third),
        )
        union = mass - interval_mass(after)
        if union >= theta:
            heavy.append((first, second, third))
    known = {
        triangle: residual
        for triangle, residual in zip(row["triangles"], row["residuals"])
    }
    residuals: list[dict[str, object]] = []
    reused = 0
    for triangle in heavy:
        residual = known.get(triangle)
        if residual is None:
            residual = residual_pair_cap(
                body,
                apex,
                triangle,
                carrier,
            )
        else:
            reused += 1
        residuals.append(residual)
    minimum_residual = min(
        (
            residual["margin"],
            triangle,
            residual["maximizing_pair"],
        )
        for triangle, residual in zip(heavy, residuals)
    ) if residuals else None
    return {
        "body": body,
        "apex": apex,
        "strong_margin": strong_margin,
        "theta": theta,
        "H_tail": tail_first,
        "H": head,
        "candidate_count": candidate_count,
        "heavy": tuple(heavy),
        "residuals": residuals,
        "reused": reused,
        "minimum_residual": minimum_residual,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    require(
        S2**2 > 2
        and file_sha256(ATLAS_PATH) == ATLAS_SHA256
        and file_sha256(ATLAS_OUTPUT_PATH) == ATLAS_OUTPUT_SHA256,
        "inherited exact kernel changed",
    )

    sample, sample_digest = root_gate_sample()
    sample_counts = tuple(
        (
            stratum,
            sum(row["stratum"] == stratum for row in sample),
            sum(
                row["stratum"] == stratum and row["fixed_margin"] > 0
                for row in sample
            ),
            max(
                row["adaptive_k"]
                for row in sample
                if row["stratum"] == stratum
            ),
        )
        for stratum in ("low", "one", "both")
    )
    if EXPECTED_SAMPLE_COUNTS is not None:
        require(
            sample_counts == EXPECTED_SAMPLE_COUNTS
            and sample_digest == EXPECTED_SAMPLE_DIGEST,
            "adaptive root sample changed",
        )

    root_by_body: dict[tuple[int, ...], dict[str, object]] = {}
    for body in BATTERY_ROOTS:
        root = global_top(body)
        adaptive_k = None
        adaptive_margin = None
        for k in range(0, ROOT_TOP_COUNT - 5):
            margin = root["m"] - sum(
                (value for value, _ in root["top"][k : k + 6]),
                F(0),
            )
            if margin > 0:
                adaptive_k = k
                adaptive_margin = margin
                break
        require(adaptive_k is not None, f"battery root gate absent: {body}")
        root["adaptive_k"] = adaptive_k
        root["adaptive_margin"] = adaptive_margin
        root_by_body[body] = root

    battery = [heavy_triangle_battery(root_by_body[body]) for body in BATTERY_ROOTS]
    strong_pair = [strong_pair_flag_battery(row) for row in battery]
    battery_lines: list[str] = []
    for row, strong in zip(battery, strong_pair):
        triple = row["triple"]
        battery_lines.append(
            f"E={','.join(map(str, row['body']))};K={row['adaptive_k']};"
            f"a={row['apex']};h={ftext(triple['m'])};"
            f"B2={ftext(triple['pair_cap'])};"
            f"cheap={ftext(triple['cheap_margin'])};"
            f"H3={ftext(triple['head_cap'])};"
            f"tail3={ftext(triple['tail_cap'])};"
            f"T3={ftext(triple['global_cap'])};"
            f"margin={ftext(triple['margin'])};"
            f"maxT={','.join(map(str, triple['maximizing_triple']))};"
            f"paidT={triple['paid']};shaT={triple['paid_digest']};"
            f"Htail={row['H_tail']};H={','.join(map(str, row['H']))};"
            f"edges={len(row['edges'])};triangles={len(row['triangles'])}\n"
        )
        battery_lines.extend(
            f" E={','.join(map(str, row['body']))};a={row['apex']};"
            f"T={','.join(map(str, triangle))};"
            f"L={ftext(residual['m'])};"
            f"cap={ftext(residual['global_cap'])};"
            f"margin={ftext(residual['margin'])};"
            f"maxP={','.join(map(str, residual['maximizing_pair']))};"
            f"paid={residual['paid']};sha={residual['paid_digest']}\n"
            for triangle, residual in zip(row["triangles"], row["residuals"])
        )
        battery_lines.append(
            f" S3E={','.join(map(str, strong['body']))};a={strong['apex']};"
            f"margin4={ftext(strong['strong_margin'])};"
            f"Htail={strong['H_tail']};"
            f"H={','.join(map(str, strong['H']))};"
            f"candidates={strong['candidate_count']};"
            f"heavy={len(strong['heavy'])};reused={strong['reused']}\n"
        )
        battery_lines.extend(
            f"  S3E={','.join(map(str, strong['body']))};a={strong['apex']};"
            f"T={','.join(map(str, triangle))};"
            f"L={ftext(residual['m'])};"
            f"cap={ftext(residual['global_cap'])};"
            f"margin={ftext(residual['margin'])};"
            f"maxP={','.join(map(str, residual['maximizing_pair']))};"
            f"paid={residual['paid']};sha={residual['paid_digest']}\n"
            for triangle, residual in zip(strong["heavy"], strong["residuals"])
        )
    battery_digest = digest_text(
        "LRC14/j6/triple-heavy-residual-pair-battery/v1\n",
        battery_lines,
    )
    battery_counts = (
        len(battery),
        sum(len(row["H"]) for row in battery),
        max(len(row["H"]) for row in battery),
        sum(len(row["edges"]) for row in battery),
        sum(len(row["triangles"]) for row in battery),
        sum(
            len(row["residuals"])
            for row in battery
        ),
        sum(
            residual["paid"]
            for row in battery
            for residual in row["residuals"]
        ),
        sum(len(row["H"]) for row in strong_pair),
        max(len(row["H"]) for row in strong_pair),
        sum(row["candidate_count"] for row in strong_pair),
        sum(len(row["heavy"]) for row in strong_pair),
        sum(row["reused"] for row in strong_pair),
        sum(
            residual["paid"]
            for row in strong_pair
            for residual in row["residuals"]
        ),
    )
    if EXPECTED_BATTERY_COUNTS is not None:
        require(
            battery_counts == EXPECTED_BATTERY_COUNTS
            and battery_digest == EXPECTED_BATTERY_DIGEST,
            "triple-heavy battery changed",
        )

    print("LRC14 J=6 ADAPTIVE GATE / TRIPLE-HEAVY SCOPED BATTERY")
    print("status=FINITE-EXACT-SCOPED-PROBE;NO-UNIFORM-THEOREM")
    print(
        "seven_body_strata=low:792,one_of_13_14:1848,both:792;"
        "THM885_j6_low_inherited=NO"
    )
    print(
        "root_sample="
        + "|".join(
            f"{name}:n={count},top12_pass={passed},max_adaptive_K={maximum}"
            for name, count, passed, maximum in sample_counts
        )
        + f";max_top30_threshold={max(row['threshold_first'] for row in sample)};"
        f"sealed_tail_first={max(row['tail_first'] for row in sample)}"
    )
    fixed_failures = [row for row in sample if row["fixed_margin"] <= 0]
    print(
        f"fixed_top12_failures={len(fixed_failures)}/36;"
        f"first={fixed_failures[0]['body']};"
        f"margin={ftext(fixed_failures[0]['fixed_margin'])};"
        f"adaptive_K={fixed_failures[0]['adaptive_k']}"
    )
    print(f"sample_digest_sha256={sample_digest}")
    for row, strong in zip(battery, strong_pair):
        triple = row["triple"]
        print(
            f"battery_E={row['body']};adaptive_K={row['adaptive_k']};"
            f"rank1_apex={row['apex']};"
            f"B2_margin={ftext(F(5,7)*triple['m']-triple['pair_cap'])};"
            f"cheap_T3_margin={ftext(triple['cheap_margin'])};"
            f"exact_T3_margin={ftext(triple['margin'])};"
            f"triple_paid={triple['paid']};"
            f"triple_max={triple['maximizing_triple']};"
            f"H={len(row['H'])};edges={len(row['edges'])};"
            f"triangles={len(row['triangles'])};"
            "all_triangle_residual_pair_caps=PASS;"
            "minimum_pair_margin="
            + (
                "NA"
                if row["minimum_residual"] is None
                else f"{ftext(row['minimum_residual'][0])};"
                f"triangle={row['minimum_residual'][1]};"
                f"maxpair={row['minimum_residual'][2]}"
            )
        )
        print(
            f"strong_pair_E={strong['body']};apex={strong['apex']};"
            f"B2_below_4h7_margin={ftext(strong['strong_margin'])};"
            f"H3={len(strong['H'])};"
            f"scalar_triple_candidates={strong['candidate_count']};"
            f"heavy_triples={len(strong['heavy'])};"
            f"reused_pair_residuals={strong['reused']};"
            "all_heavy_triple_residual_pair_caps=PASS;"
            "minimum_pair_margin="
            + (
                "NA"
                if strong["minimum_residual"] is None
                else f"{ftext(strong['minimum_residual'][0])};"
                f"triangle={strong['minimum_residual'][1]};"
                f"maxpair={strong['minimum_residual'][2]}"
            )
        )
    print(
        f"battery_counts=roots:{battery_counts[0]},"
        f"H:{battery_counts[1]},maxH:{battery_counts[2]},"
        f"edges:{battery_counts[3]},triangles:{battery_counts[4]},"
        f"pair_residuals:{battery_counts[5]},"
        f"pair_paid:{battery_counts[6]},"
        f"strong_H3:{battery_counts[7]},max_strong_H3:{battery_counts[8]},"
        f"strong_candidates:{battery_counts[9]},"
        f"strong_heavy_triples:{battery_counts[10]},"
        f"strong_reused:{battery_counts[11]},"
        f"strong_pair_paid_including_reuse:{battery_counts[12]}"
    )
    print(f"battery_digest_sha256={battery_digest}")
    print(
        "scope=36 deterministic roots plus three rank-one hostile apices;"
        "not all 3432 roots,not the seven-body rung,not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
