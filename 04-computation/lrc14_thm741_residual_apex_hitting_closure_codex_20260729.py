#!/usr/bin/env python3
r"""Exact ranked-apex hitting closure of the residual THM-741 chamber.

This is an independent closure companion to the ranked-head atlas.  It
reconstructs the 602 rank-impossible nine-speed roots, seals their global coverage
ranks through rank fourteen, and tests the top-ten apex reduction

    q_11 + q_12 + q_13 + q_14 < |G_E|.

Consequently every quadruple whose individual coverages can exhaust ``G_E``
contains one of the chosen ten global apices.

For each root-apex pair it then removes the apex comb and profiles the
three-slot carrier.  A carrier is:

* ``direct`` if its three largest individual coverages sum to less than its
  mass;
* ``rank2`` if direct closure fails but the residual after its two largest
  coverages is strictly above the limiting mass/7 level;
* ``failure`` otherwise.

The remaining nominal carrier failures are disposed of by a weighted hitting
gate: the four largest coverages in the complement of the closed apices have
sum below the root mass.  All comparisons are exact.  The periodic-tooth
primitive is implemented locally instead of imported from the ranked-head
atlas so that the two computations are independent.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FIRST_EXTERNAL = 15
BASE_HORIZON = 600
S2 = F(99, 70)
BODIES = tuple(combinations(range(1, 15), 9))

EXPECTED_ROOT_PARTITION = (584, 816, 602)
EXPECTED_IMPOSSIBLE_DIGEST = (
    "a4ad9f9b5a8ce16103450ac05684cd13ec33637e8e3737218a9186086cd639d4"
)
EXPECTED_MAXIMUM_ROOT_MASS = F(353, 1176)
EXPECTED_APEX_MINIMUM = (
    F(67759, 5045040),
    (1, 2, 3, 5, 7, 8, 9, 11, 13),
)
EXPECTED_RANK14_MAXIMUM = (989, (1, 2, 8, 9, 10, 11, 12, 13, 14))
EXPECTED_APEX_DIGEST = (
    "af026722b290a2e578ded6f105722a104836a2cab31dee2272d836a7a7c14ecf"
)
EXPECTED_CARRIER_COUNTS = (6020, 4272, 1657, 91)
EXPECTED_RANK3_MAXIMUM = (
    1166,
    (2, 3, 4, 8, 9, 10, 11, 12, 14),
    7,
    156,
)
EXPECTED_CLASSIFICATION_DIGEST = (
    "24a4b9b07424b23f6fe334f84c45ef45efa0e2c5ea695d402d1be88687bd87e0"
)
EXPECTED_FAILURE_DISTINCT_ROOTS = 64
EXPECTED_FAILURE_RANK_HISTOGRAM = {
    1: 0,
    2: 9,
    3: 14,
    4: 20,
    5: 14,
    6: 6,
    7: 10,
    8: 7,
    9: 7,
    10: 4,
}
EXPECTED_FAILURE_DIGEST = (
    "2f97c76ccbf5ae6f53c9c32b7a8366a3bd62a4b83957c63fb96283c9cbec3743"
)
EXPECTED_FAILURE_ROOT_DIGEST = (
    "43ac76ed0be64801affe402fdd3f6b0374896eb0f35e9b18671b52701edf2230"
)
EXPECTED_HITTING_MINIMUM = (
    F(57684167, 7467740280),
    (1, 3, 4, 5, 7, 9, 10, 11, 13),
    (23, 17, 53, 72),
)
EXPECTED_HITTING_DIGEST = (
    "f6c5f21c5a16651b69a681d2709a498b83cabf6bc845cce27dfd37ff0b8cc000"
)
EXPECTED_RANK2_COUNTS = (1657, 29622)
EXPECTED_K2_MAXIMUM = (
    5711,
    (1, 4, 5, 6, 7, 8, 9, 10, 13),
    6,
    23,
)
EXPECTED_RANK2_TAIL_MAXIMUM = (
    168888,
    (1, 4, 5, 6, 7, 8, 9, 10, 13),
    6,
    23,
)
EXPECTED_DANGEROUS_MAXIMUM = (
    5807,
    (1, 4, 5, 6, 7, 8, 9, 10, 13),
    6,
    23,
)
EXPECTED_LITERAL_MINIMUM = (
    F(183, 43792),
    (1, 3, 4, 7, 9, 10, 11, 12, 13),
    16,
    (17, 19, 23),
    6,
)
EXPECTED_MINIMUM_COMPONENTS = 4
EXPECTED_CANDIDATE_DIGEST = (
    "9df78b9841ba5fbac29b49965dd45498656d5093d377d16c8929699bcb308243"
)
EXPECTED_LITERAL_DIGEST = (
    "dbef39aa9b7d681706321835dfd34143c1303392c58ea5788a41c24aa05780a4"
)
EXPECTED_CONTROL_COUNTS = (6020, 12040, 1657, 1657)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(file_sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_apex_core", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def tooth_primitive(value: F) -> F:
    integer = value.numerator // value.denominator
    remainder = value - integer
    return (
        F(integer, 7)
        + min(remainder, F(1, 14))
        + max(F(0), remainder - F(13, 14))
    )


def coverage(good: list[tuple[F, F]], speed: int) -> F:
    return sum(
        (
            tooth_primitive(speed * right)
            - tooth_primitive(speed * left)
        )
        / speed
        for left, right in good
    )


def coverages_many(
    good: list[tuple[F, F]],
    speeds: list[int] | range,
) -> list[tuple[F, int]]:
    """Exact vectorized primitive evaluations on a finite speed list."""

    speed_list = list(speeds)
    if not speed_list:
        return []
    denominator = lcm(
        *(endpoint.denominator for interval in good for endpoint in interval)
    )
    endpoint_ints = np.array(
        [
            endpoint.numerator * (denominator // endpoint.denominator)
            for interval in good
            for endpoint in interval
        ],
        dtype=np.int64,
    )
    speed_array = np.array(speed_list, dtype=np.int64)
    carrier_m = sum((right - left for left, right in good), F(0))
    # Every profiled object here is an apex carrier.  Its exact mass is below
    # 2/7.  Thus the nonnegative final numerator is below 4*D*w, while each
    # scaled primitive is smaller still.  This guard protects both elementwise
    # arithmetic and the int64 axis reduction.
    require(carrier_m < F(2, 7), "vectorized object is not an apex carrier")
    require(
        4 * denominator * max(speed_list) < 2**63,
        "vector primitive would overflow int64",
    )
    scaled = endpoint_ints[:, None] * speed_array[None, :]
    quotient, remainder = np.divmod(scaled, denominator)
    primitive = (
        2 * denominator * quotient
        + np.minimum(14 * remainder, denominator)
        + np.maximum(0, 14 * remainder - 13 * denominator)
    )
    numerators = (primitive[1::2] - primitive[0::2]).sum(axis=0)
    return [
        (F(int(numerator), 14 * denominator * speed), speed)
        for numerator, speed in zip(numerators, speed_list)
    ]


def local_bad_pieces(left: F, right: F, speed: int) -> list[tuple[F, F]]:
    radius = F(1, 14 * speed)
    lower = speed * (left - radius)
    upper = speed * (right + radius)
    first = lower.numerator // lower.denominator - 1
    last = -((-upper.numerator) // upper.denominator) + 1
    pieces: list[tuple[F, F]] = []
    for tooth in range(first, last + 1):
        center = F(tooth, speed)
        lo = max(left, center - radius)
        hi = min(right, center + radius)
        if lo < hi:
            pieces.append((lo, hi))
    return pieces


def subtract_local(
    carrier: list[tuple[F, F]],
    speed: int,
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    for left, right in carrier:
        pieces = sorted(local_bad_pieces(left, right, speed))
        cursor = left
        for lo, hi in pieces:
            if hi <= cursor:
                continue
            if lo > cursor:
                out.append((cursor, lo))
            cursor = max(cursor, hi)
            if cursor >= right:
                break
        if cursor < right:
            out.append((cursor, right))
    return out


def subtract_local_multi(
    carrier: list[tuple[F, F]],
    speeds: tuple[int, ...],
) -> list[tuple[F, F]]:
    """Remove several combs in one endpoint-event merge."""

    out: list[tuple[F, F]] = []
    for left, right in carrier:
        pieces = sorted(
            piece
            for speed in speeds
            for piece in local_bad_pieces(left, right, speed)
        )
        cursor = left
        for lo, hi in pieces:
            if hi <= cursor:
                continue
            if lo > cursor:
                out.append((cursor, lo))
            cursor = max(cursor, hi)
            if cursor >= right:
                break
        if cursor < right:
            out.append((cursor, right))
    return out


def base_profile(body: tuple[int, ...]) -> dict[str, object]:
    good, root_r, root_m = CORE.good_norm(body)
    rows = [
        (coverage(good, speed), speed)
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
    ]
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    fourth = ranked[3][0]
    require(fourth > root_m / 7, f"fourth rank below limit: {body}")
    threshold4 = S2 * root_r / (7 * (fourth - root_m / 7))
    require(threshold4 < BASE_HORIZON + 1, f"top four not sealed: {body}")
    top4_margin = root_m - sum((value for value, _ in ranked[:4]), F(0))
    residual3 = root_m - sum((value for value, _ in ranked[:3]), F(0))
    return {
        "body": body,
        "good": good,
        "r": root_r,
        "m": root_m,
        "rows": rows,
        "top4_margin": top4_margin,
        "residual3": residual3,
    }


def residual_profiles(workers: int) -> list[dict[str, object]]:
    context = mp.get_context("spawn")
    if workers == 1:
        profiles = list(map(base_profile, BODIES))
    else:
        with context.Pool(workers) as pool:
            profiles = pool.map(base_profile, BODIES, chunksize=2)
    require(
        tuple(row["body"] for row in profiles) == BODIES,
        "parallel root order changed",
    )
    positive = [row for row in profiles if row["top4_margin"] > 0]
    nonpositive = [row for row in profiles if row["top4_margin"] <= 0]
    feasible = [
        row for row in nonpositive if row["residual3"] > row["m"] / 7
    ]
    impossible = [
        row
        for row in nonpositive
        if row["residual3"] <= row["m"] / 7
    ]
    require(
        (len(positive), len(feasible), len(impossible))
        == EXPECTED_ROOT_PARTITION[:3],
        "global root partition changed",
    )
    require(
        max(row["m"] for row in profiles) == EXPECTED_MAXIMUM_ROOT_MASS,
        "maximum root mass changed",
    )
    impossible_text = "\n".join(
        ",".join(map(str, row["body"])) for row in impossible
    )
    require(
        hashlib.sha256(impossible_text.encode()).hexdigest()
        == EXPECTED_IMPOSSIBLE_DIGEST,
        "rank-impossible root atlas changed",
    )
    return impossible


def seal_root_top14(row: dict[str, object]) -> dict[str, object]:
    body = row["body"]
    good = row["good"]
    root_r = row["r"]
    root_m = row["m"]
    rows = list(row["rows"])
    require(isinstance(body, tuple), "body type changed")
    require(isinstance(good, list), "carrier type changed")
    require(isinstance(root_r, int), "component count type changed")
    require(isinstance(root_m, F), "mass type changed")
    ranked600 = sorted(rows, key=lambda item: (-item[0], item[1]))
    q14 = ranked600[13][0]
    require(q14 > root_m / 7, f"rank fourteen misses limiting level: {body}")
    threshold14 = S2 * root_r / (7 * (q14 - root_m / 7))
    tail_first = ceiling(threshold14)
    for speed in range(BASE_HORIZON + 1, tail_first):
        rows.append((coverage(good, speed), speed))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q14_global = ranked[13][0]
    # The threshold made from the possibly smaller W=600 q14 is conservative.
    require(
        root_m / 7 + S2 * root_r / (7 * tail_first) <= q14,
        f"rank-fourteen tail seal failed: {body}",
    )
    outside4 = sum((value for value, _ in ranked[10:14]), F(0))
    margin = root_m - outside4
    return {
        "body": body,
        "good": good,
        "r": root_r,
        "m": root_m,
        "top10": tuple(ranked[:10]),
        "rank14": tuple(ranked[:14]),
        "q14": q14_global,
        "threshold14": threshold14,
        "tail_first14": tail_first,
        "apex_margin": margin,
        "apex_ok": margin > 0,
    }


def carrier_base_profile(
    task: tuple[dict[str, object], int, int],
) -> dict[str, object]:
    root, apex_rank, apex = task
    body = root["body"]
    root_good = root["good"]
    require(isinstance(body, tuple), "body type changed")
    require(isinstance(root_good, list), "root carrier type changed")
    carrier = subtract_local(root_good, apex)
    carrier_m = sum((right - left for left, right in carrier), F(0))
    carrier_r = len(carrier)
    require(carrier_m > 0 and carrier_r > 0, f"empty apex carrier: {body}, {apex}")
    direct_r, direct_m, direct_good = CORE.subtract(root_good, apex)
    require(
        direct_good == carrier
        and direct_r == carrier_r
        and direct_m == carrier_m,
        f"local/full apex-carrier mismatch: {body}, {apex}",
    )
    base_speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed != apex
    ]
    rows = coverages_many(carrier, base_speeds)
    for control_speed in (base_speeds[0], base_speeds[8]):
        require(
            dict((speed, value) for value, speed in rows)[control_speed]
            == coverage(carrier, control_speed),
            f"vector/scalar coverage mismatch: {body}, {apex}, {control_speed}",
        )
    ranked600 = sorted(rows, key=lambda item: (-item[0], item[1]))
    q3 = ranked600[2][0]
    require(q3 > carrier_m / 7, f"carrier rank three misses limit: {body}, {apex}")
    threshold3 = S2 * carrier_r / (7 * (q3 - carrier_m / 7))
    tail_first3 = ceiling(threshold3)
    tail_speeds = [
        speed
        for speed in range(BASE_HORIZON + 1, tail_first3)
        if speed != apex
    ]
    rows.extend(coverages_many(carrier, tail_speeds))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    direct_margin = carrier_m - sum((value for value, _ in ranked[:3]), F(0))
    residual2 = carrier_m - sum((value for value, _ in ranked[:2]), F(0))
    classification = (
        "direct"
        if direct_margin > 0
        else ("rank2" if residual2 > carrier_m / 7 else "failure")
    )
    return {
        "body": body,
        "apex": apex,
        "apex_rank": apex_rank,
        "r": carrier_r,
        "m": carrier_m,
        "top3": tuple(ranked[:3]),
        "direct_margin": direct_margin,
        "residual2": residual2,
        "classification": classification,
        "threshold3": threshold3,
        "tail_first3": tail_first3,
        "carrier": carrier if classification == "rank2" else None,
    }


def enrich_rank2(row: dict[str, object]) -> dict[str, object]:
    if row["classification"] != "rank2":
        return row
    carrier = row["carrier"]
    carrier_m = row["m"]
    carrier_r = row["r"]
    residual2 = row["residual2"]
    apex = row["apex"]
    require(isinstance(carrier, list), "carrier type changed")
    require(isinstance(carrier_m, F), "mass type changed")
    require(isinstance(carrier_r, int), "component count type changed")
    require(isinstance(residual2, F), "residual type changed")
    threshold = S2 * carrier_r / (7 * (residual2 - carrier_m / 7))
    tail_first = ceiling(threshold)
    head: list[tuple[F, int]] = []
    chunk = 10_000
    for start in range(FIRST_EXTERNAL, tail_first, chunk):
        speeds = [
            speed
            for speed in range(start, min(tail_first, start + chunk))
            if speed != apex
        ]
        head.extend(
            (value, speed)
            for value, speed in coverages_many(carrier, speeds)
            if value >= residual2
        )
    out = dict(row)
    out["rank2_threshold"] = threshold
    out["rank2_tail_first"] = tail_first
    out["rank2_head"] = tuple(sorted(head, key=lambda item: (-item[0], item[1])))
    out["K2"] = len(head)
    ranked = out["rank2_head"]
    dangerous: list[tuple[int, int, int]] = []
    for first in range(len(ranked) - 2):
        if (
            ranked[first][0]
            + ranked[first + 1][0]
            + ranked[first + 2][0]
            < carrier_m
        ):
            break
        for second in range(first + 1, len(ranked) - 1):
            prefix = ranked[first][0] + ranked[second][0]
            if prefix + ranked[second + 1][0] < carrier_m:
                break
            for third in range(second + 1, len(ranked)):
                if prefix + ranked[third][0] < carrier_m:
                    break
                dangerous.append(
                    tuple(
                        sorted(
                            (
                                ranked[first][1],
                                ranked[second][1],
                                ranked[third][1],
                            )
                        )
                    )
                )
    require(
        len(dangerous) == len(set(dangerous)),
        f"duplicate dangerous carrier triple: {row['body']}, {row['apex']}",
    )
    out["dangerous_triples"] = tuple(sorted(dangerous))
    out["dangerous_count"] = len(dangerous)
    return out


def close_rank2(row: dict[str, object]) -> dict[str, object]:
    carrier = row["carrier"]
    body = row["body"]
    apex = row["apex"]
    head = row["rank2_head"]
    candidates = row["dangerous_triples"]
    require(isinstance(carrier, list), "carrier type changed")
    require(isinstance(body, tuple), "body type changed")
    require(isinstance(apex, int), "apex type changed")
    require(isinstance(head, tuple), "head type changed")
    require(isinstance(candidates, tuple), "candidate type changed")
    rank_index = {speed: index for index, (_, speed) in enumerate(head)}
    pair_cache: dict[tuple[int, int], list[tuple[F, F]]] = {}
    literal_rows: list[str] = []
    minimum = None
    minimum_components = None
    control_done = False
    simultaneous_controls = 0
    direct_controls = 0
    for triple in candidates:
        rank_triple = tuple(sorted(triple, key=rank_index.__getitem__))
        pair = tuple(sorted(rank_triple[:2]))
        third = rank_triple[2]
        pair_carrier = pair_cache.get(pair)
        if pair_carrier is None:
            pair_carrier = subtract_local_multi(carrier, pair)
            pair_cache[pair] = pair_carrier
        nested = subtract_local(pair_carrier, third)
        survivor = sum((right - left for left, right in nested), F(0))
        components = len(nested)
        require(
            survivor > 0,
            f"nonpositive rank2 carrier survivor: {body}, {apex}, {triple}",
        )
        if not control_done:
            require(
                subtract_local_multi(carrier, triple) == nested,
                f"cached-pair/simultaneous mismatch: {body}, {apex}, {triple}",
            )
            simultaneous_controls += 1
            if max(triple) <= 2500:
                family = tuple(sorted((*body, apex, *triple)))
                direct_good, direct_r, direct_m = CORE.good_norm(family)
                require(
                    direct_good == nested
                    and direct_r == components
                    and direct_m == survivor,
                    f"direct thirteen-comb mismatch: {family}",
                )
                direct_controls += 1
            control_done = True
        record = (survivor, body, apex, triple, components)
        if minimum is None or record < minimum:
            minimum = record
        component_record = (components, survivor, body, apex, triple)
        if minimum_components is None or component_record < minimum_components:
            minimum_components = component_record
        literal_rows.append(
            "E="
            + ",".join(map(str, body))
            + f";rank={row['apex_rank']};a={apex};"
            + "T="
            + ",".join(map(str, triple))
            + f";L={ftext(survivor)};r={components}\n"
        )
    return {
        "body": body,
        "apex": apex,
        "apex_rank": row["apex_rank"],
        "dangerous_count": len(candidates),
        "minimum": minimum,
        "minimum_components": minimum_components,
        "simultaneous_controls": simultaneous_controls,
        "direct_controls": direct_controls,
        "literal_rows": tuple(literal_rows),
    }


def digest_rows(rows: list[str], header: str) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    parser.add_argument("--apex-only", action="store_true")
    parser.add_argument("--skip-rank2-heads", action="store_true")
    parser.add_argument("--literal-rank2", action="store_true")
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    print("THM-2883 RANKED-APEX HITTING CLOSURE")
    print("status=FINITE-EXACT+GLOBAL-TAIL-SEALED")
    print(
        "universe=all 602 rank-impossible roots after the disjoint "
        "584+816+602 THM-741 root partition"
    )
    residual = residual_profiles(args.workers)
    context = mp.get_context("spawn")
    if args.workers == 1:
        roots = list(map(seal_root_top14, residual))
    else:
        with context.Pool(args.workers) as pool:
            roots = pool.map(seal_root_top14, residual, chunksize=1)
    require(all(row["apex_ok"] for row in roots), "top-ten apex reduction failed")
    minimum_apex = min(
        (row["apex_margin"], row["body"], row["rank14"][10:14])
        for row in roots
    )
    maximum_rank14_horizon = max(
        (row["tail_first14"], row["body"], row["threshold14"])
        for row in roots
    )
    apex_rows = [
        "E="
        + ",".join(map(str, row["body"]))
        + ";top10="
        + ",".join(str(speed) for _, speed in row["top10"])
        + f";margin={ftext(row['apex_margin'])};tail={row['tail_first14']}\n"
        for row in roots
    ]
    apex_digest = digest_rows(
        apex_rows,
        "THM741/residual-top10-apex-reduction/v1\n",
    )
    require(
        minimum_apex[:2] == EXPECTED_APEX_MINIMUM,
        "minimum top-ten apex margin changed",
    )
    require(
        maximum_rank14_horizon[:2] == EXPECTED_RANK14_MAXIMUM,
        "maximum rank-fourteen horizon changed",
    )
    require(apex_digest == EXPECTED_APEX_DIGEST, "apex ledger changed")
    print(
        "rank_impossible_roots=602;"
        f"min_apex_margin={ftext(minimum_apex[0])};"
        f"min_body={minimum_apex[1]};"
        f"max_rank14_tail={maximum_rank14_horizon[0]};"
        f"max_tail_body={maximum_rank14_horizon[1]};"
        f"apex_digest={apex_digest}",
        flush=True,
    )
    if args.apex_only:
        return

    tasks = [
        (root, apex_rank, speed)
        for root in roots
        for apex_rank, (_, speed) in enumerate(root["top10"], start=1)
    ]
    require(len(tasks) == 6020, "root-apex task census changed")
    if args.workers == 1:
        carriers = list(map(carrier_base_profile, tasks))
    else:
        with context.Pool(args.workers) as pool:
            carriers = pool.map(carrier_base_profile, tasks, chunksize=1)
    counts = {
        name: sum(row["classification"] == name for row in carriers)
        for name in ("direct", "rank2", "failure")
    }
    classification_rows = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']};class={row['classification']}\n"
        for row in carriers
    ]
    classification_digest = digest_rows(
        classification_rows,
        "THM741/residual-apex-carrier-classification/v1\n",
    )
    failures = [row for row in carriers if row["classification"] == "failure"]
    failure_roots = sorted({row["body"] for row in failures})
    failure_rows = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']}\n"
        for row in failures
    ]
    failure_digest = digest_rows(
        failure_rows,
        "THM741/residual-apex-carrier-failures/v1\n",
    )
    failure_root_text = "\n".join(",".join(map(str, body)) for body in failure_roots)
    failure_root_digest = hashlib.sha256(failure_root_text.encode()).hexdigest()
    failure_rank_histogram = {
        rank: sum(row["apex_rank"] == rank for row in failures)
        for rank in range(1, 11)
    }
    failures_by_body: dict[tuple[int, ...], set[int]] = {}
    for row in failures:
        failures_by_body.setdefault(row["body"], set()).add(row["apex"])
    hitting_rows = []
    for root in roots:
        failed_speeds = failures_by_body.get(root["body"], set())
        allowed = [
            (value, speed)
            for index, (value, speed) in enumerate(root["rank14"])
            if index >= 10 or speed in failed_speeds
        ]
        require(len(allowed) >= 4, f"short hitting complement: {root['body']}")
        allowed_top4 = tuple(allowed[:4])
        margin = root["m"] - sum((value for value, _ in allowed_top4), F(0))
        hitting_rows.append(
            {
                "body": root["body"],
                "margin": margin,
                "allowed_top4": allowed_top4,
                "closed": margin > 0,
            }
        )
    hitting_failures = [row for row in hitting_rows if not row["closed"]]
    minimum_hitting = min(
        (row["margin"], row["body"], row["allowed_top4"])
        for row in hitting_rows
    )
    hitting_ledger = "THM741/residual-closed-apex-hitting-gate/v1\n" + "".join(
        "E="
        + ",".join(map(str, row["body"]))
        + ";margin="
        + ftext(row["margin"])
        + ";allowed="
        + ",".join(str(speed) for _, speed in row["allowed_top4"])
        + "\n"
        for row in hitting_rows
    )
    maximum_rank3 = max(
        (row["tail_first3"], row["body"], row["apex_rank"], row["apex"])
        for row in carriers
    )
    require(
        (
            len(carriers),
            counts["direct"],
            counts["rank2"],
            counts["failure"],
        )
        == EXPECTED_CARRIER_COUNTS,
        "apex-carrier classification census changed",
    )
    require(
        maximum_rank3 == EXPECTED_RANK3_MAXIMUM,
        "maximum carrier rank-three horizon changed",
    )
    require(
        classification_digest == EXPECTED_CLASSIFICATION_DIGEST,
        "carrier classification ledger changed",
    )
    require(
        len(failure_roots) == EXPECTED_FAILURE_DISTINCT_ROOTS
        and failure_rank_histogram == EXPECTED_FAILURE_RANK_HISTOGRAM
        and failure_digest == EXPECTED_FAILURE_DIGEST
        and failure_root_digest == EXPECTED_FAILURE_ROOT_DIGEST,
        "nominal failure atlas changed",
    )
    hitting_digest = hashlib.sha256(hitting_ledger.encode()).hexdigest()
    require(
        not hitting_failures
        and (
            minimum_hitting[0],
            minimum_hitting[1],
            tuple(speed for _, speed in minimum_hitting[2]),
        )
        == EXPECTED_HITTING_MINIMUM
        and hitting_digest == EXPECTED_HITTING_DIGEST,
        "closed-apex hitting gate changed",
    )
    print(
        f"carrier_tasks={len(carriers)};"
        f"direct={counts['direct']};rank2={counts['rank2']};"
        f"failure={counts['failure']};"
        f"max_rank3_tail={maximum_rank3[0]};"
        f"max_rank3_body={maximum_rank3[1]};"
        f"max_rank3_apex_rank={maximum_rank3[2]};"
        f"max_rank3_apex={maximum_rank3[3]};"
        f"full_carrier_controls={len(carriers)};"
        f"vector_scalar_controls={2*len(carriers)};"
        f"classification_digest={classification_digest}",
        flush=True,
    )
    print(
        f"failure_distinct_roots={len(failure_roots)};"
        f"potential_whole_roots_after_rank2={602-len(failure_roots)};"
        f"failure_rank_histogram={failure_rank_histogram};"
        f"failure_digest={failure_digest};"
        f"failure_root_digest={failure_root_digest}",
        flush=True,
    )
    print(
        f"closed_apex_hitting_roots={602-len(hitting_failures)};"
        f"hitting_residual={len(hitting_failures)};"
        f"minimum_hitting_margin={ftext(minimum_hitting[0])};"
        f"minimum_hitting_body={minimum_hitting[1]};"
        "minimum_hitting_allowed="
        + ",".join(str(speed) for _, speed in minimum_hitting[2])
        + ";hitting_digest="
        + hitting_digest,
        flush=True,
    )
    if args.skip_rank2_heads:
        return

    rank2_rows = [row for row in carriers if row["classification"] == "rank2"]
    if args.workers == 1:
        enriched = list(map(enrich_rank2, rank2_rows))
    else:
        with context.Pool(args.workers) as pool:
            enriched = pool.map(enrich_rank2, rank2_rows, chunksize=1)
    if enriched:
        maximum_k = max(
            (row["K2"], row["body"], row["apex_rank"], row["apex"])
            for row in enriched
        )
        maximum_horizon = max(
            (
                row["rank2_tail_first"],
                row["body"],
                row["apex_rank"],
                row["apex"],
            )
            for row in enriched
        )
        total_dangerous = sum(row["dangerous_count"] for row in enriched)
        maximum_dangerous = max(
            (
                row["dangerous_count"],
                row["body"],
                row["apex_rank"],
                row["apex"],
            )
            for row in enriched
        )
        require(
            (len(enriched), total_dangerous) == EXPECTED_RANK2_COUNTS
            and maximum_k == EXPECTED_K2_MAXIMUM
            and maximum_horizon == EXPECTED_RANK2_TAIL_MAXIMUM
            and maximum_dangerous == EXPECTED_DANGEROUS_MAXIMUM,
            "rank-two head census changed",
        )
        print(
            f"rank2_tasks={len(enriched)};"
            f"max_K2={maximum_k[0]};max_K2_body={maximum_k[1]};"
            f"max_K2_apex_rank={maximum_k[2]};max_K2_apex={maximum_k[3]};"
            f"max_rank2_tail={maximum_horizon[0]};"
            f"max_tail_body={maximum_horizon[1]};"
            f"max_tail_apex_rank={maximum_horizon[2]};"
            f"max_tail_apex={maximum_horizon[3]};"
            f"dangerous_triples={total_dangerous};"
            f"max_dangerous_per_carrier={maximum_dangerous[0]};"
            f"max_dangerous_body={maximum_dangerous[1]};"
            f"max_dangerous_apex_rank={maximum_dangerous[2]};"
            f"max_dangerous_apex={maximum_dangerous[3]}",
            flush=True,
        )
    if not args.literal_rank2:
        return

    if args.workers == 1:
        closures = list(map(close_rank2, enriched))
    else:
        with context.Pool(args.workers) as pool:
            closures = pool.map(close_rank2, enriched, chunksize=1)
    require(
        sum(row["dangerous_count"] for row in closures)
        == sum(row["dangerous_count"] for row in enriched),
        "rank2 literal candidate census changed",
    )
    candidate_text = "THM741/residual-apex-rank2-dangerous-triples/v1\n" + "".join(
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']};T="
        + ",".join(map(str, triple))
        + "\n"
        for row in enriched
        for triple in row["dangerous_triples"]
    )
    literal_text = "THM741/residual-apex-rank2-literal-endpoints/v1\n" + "".join(
        line for row in closures for line in row["literal_rows"]
    )
    minima = [row["minimum"] for row in closures if row["minimum"] is not None]
    component_minima = [
        row["minimum_components"]
        for row in closures
        if row["minimum_components"] is not None
    ]
    require(minima and component_minima, "rank2 literal controls disappeared")
    minimum = min(minima)
    minimum_components = min(component_minima)
    candidate_digest = hashlib.sha256(candidate_text.encode()).hexdigest()
    literal_digest = hashlib.sha256(literal_text.encode()).hexdigest()
    simultaneous_controls = sum(row["simultaneous_controls"] for row in closures)
    direct_controls = sum(row["direct_controls"] for row in closures)
    require(
        minimum == EXPECTED_LITERAL_MINIMUM
        and minimum_components[0] == EXPECTED_MINIMUM_COMPONENTS
        and candidate_digest == EXPECTED_CANDIDATE_DIGEST
        and literal_digest == EXPECTED_LITERAL_DIGEST
        and (
            len(carriers),
            2 * len(carriers),
            simultaneous_controls,
            direct_controls,
        )
        == EXPECTED_CONTROL_COUNTS,
        "literal rank-two closure ledger changed",
    )
    print(
        f"literal_rank2_triples={sum(row['dangerous_count'] for row in closures)};"
        "positive=all;tight=0;"
        f"minimum_survivor={ftext(minimum[0])};"
        f"minimum_body={minimum[1]};minimum_apex={minimum[2]};"
        f"minimum_triple={minimum[3]};minimum_r={minimum[4]};"
        f"minimum_components={minimum_components[0]};"
        f"simultaneous_controls={simultaneous_controls};"
        f"direct_controls={direct_controls};"
        f"candidate_digest={candidate_digest};"
        f"literal_digest={literal_digest}",
        flush=True,
    )
    print(
        "pure_tail_root_consequence=2002/2002;"
        "whole_root_composition=THM-738 small-speed chamber;"
        "THM-741=PROVED;general_LRC14=OPEN"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
