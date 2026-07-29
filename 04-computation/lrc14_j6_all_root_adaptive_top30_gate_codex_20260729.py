#!/usr/bin/env python3
r"""Exact adaptive six-cover gate atlas for all seven-body LRC roots.

For every

    E in C({1,...,14},7),

let ``G_E`` be its gap-``1/14`` good set, let ``m_E=|G_E|``, and put

    c_E(w)=|G_E intersect D_w|,                  w>=15.

Write ``q_1(E)>=q_2(E)>=...`` for the globally ranked values ``c_E(w)``,
breaking ties by speed.  If six external teeth cover ``G_E`` while avoiding
the first ``K`` globally ranked speeds, then necessarily

    m_E <= q_{K+1}(E)+...+q_{K+6}(E).

Consequently any strict reverse inequality is a valid hitting gate: every
putative six-cover meets the first ``K`` speeds.  This verifier finds the
least such ``K`` in ``{0,...,24}`` from an exact global top-thirty ledger.

The ledger is global, not merely a finite-head guess.  It scans through
``BASE_HORIZON`` and uses the strict THM-735(ii) discrepancy estimate

    c_E(w) < m_E/7 + (99/70) r_E/(7w)

to seal every omitted tail speed below the provisional thirtieth value.
When the resulting threshold exceeds the base horizon, the missing finite
head is scanned exactly before ranking.  Scalar evaluations independently
check four boundary/rank entries per root.

The universe is also split into roots using zero, one, or both of 13 and 14.
This distinction is bookkeeping only: the 792 roots inside {1,...,12} are
not inherited terminals, because the earlier j=6 sweep is unfinished.

This is a hitting-gate atlas only.  It neither closes the first-apex
carriers nor proves the seven-body/six-slot rung or LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2885_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
THM2885_SHA256 = (
    "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
)

FIRST_EXTERNAL = 15
BASE_HORIZON = 1600
TOP_COUNT = 30
COVER_SIZE = 6
MAX_GATE_SIZE = TOP_COUNT - COVER_SIZE
S2 = F(99, 70)

BODIES = tuple(combinations(range(1, 15), 7))

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_GATE_COUNTS: tuple[int, ...] | None = None
EXPECTED_K_DISTRIBUTION: tuple[tuple[str, tuple[tuple[int, int], ...]], ...] | None = (
    None
)
EXPECTED_EXTREMA: tuple[object, ...] | None = None
EXPECTED_STRATIFIED_EXTREMA: tuple[tuple[str, tuple[object, ...]], ...] | None = (
    None
)
EXPECTED_LEDGER_DIGEST: str | None = None
EXPECTED_STRATIFIED_DIGESTS: tuple[tuple[str, str], ...] | None = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_thm2885():
    require(file_sha256(THM2885_PATH) == THM2885_SHA256, "THM-2885 changed")
    spec = importlib.util.spec_from_file_location("j6_all_root_thm2885", THM2885_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-2885")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_thm2885()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def stratum(body: tuple[int, ...]) -> str:
    high = int(13 in body) + int(14 in body)
    return ("low", "one", "both")[high]


def profile_body(body: tuple[int, ...]) -> dict[str, object]:
    """Return one exact globally sealed top-thirty gate profile."""

    good, components, mass = T.CORE.good_norm(body)
    require(
        components == len(good) and mass > 0,
        f"bad seven-body carrier: {body}",
    )

    rows = T.coverages_many(
        good,
        range(FIRST_EXTERNAL, BASE_HORIZON + 1),
    )
    base_ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q30_base = base_ranked[TOP_COUNT - 1][0]
    require(
        q30_base > mass / 7,
        f"base rank {TOP_COUNT} misses the discrepancy limit: {body}",
    )

    threshold = S2 * components / (7 * (q30_base - mass / 7))
    threshold_first = ceiling(threshold)
    tail_first = max(BASE_HORIZON + 1, threshold_first)
    if tail_first > BASE_HORIZON + 1:
        rows.extend(
            T.coverages_many(
                good,
                range(BASE_HORIZON + 1, tail_first),
            )
        )

    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q30_base,
        f"rank-{TOP_COUNT} tail did not seal: {body}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top = tuple(ranked[:TOP_COUNT])
    require(
        len(top) == TOP_COUNT
        and len({speed for _, speed in top}) == TOP_COUNT,
        f"bad top-{TOP_COUNT} ledger: {body}",
    )

    adaptive_k: int | None = None
    adaptive_margin: F | None = None
    margins: list[F] = []
    for k in range(MAX_GATE_SIZE + 1):
        margin = mass - sum(
            (value for value, _ in top[k : k + COVER_SIZE]),
            F(0),
        )
        margins.append(margin)
        if adaptive_k is None and margin > 0:
            adaptive_k = k
            adaptive_margin = margin
    require(
        adaptive_k is not None and adaptive_margin is not None,
        f"top-{TOP_COUNT} contains no valid adaptive gate: {body}",
    )
    fixed_k12_margin = margins[12]

    by_speed = {speed: value for value, speed in rows}
    control_speeds = tuple(
        dict.fromkeys(
            (
                FIRST_EXTERNAL,
                top[0][1],
                top[-1][1],
                tail_first - 1,
            )
        )
    )
    for speed in control_speeds:
        require(
            speed in by_speed,
            f"scalar control escaped finite head: {body}, {speed}",
        )
        require(
            by_speed[speed] == T.coverage(good, speed),
            f"vector/scalar mismatch: {body}, {speed}",
        )

    return {
        "body": body,
        "stratum": stratum(body),
        "m": mass,
        "r": components,
        "top": top,
        "adaptive_k": adaptive_k,
        "adaptive_margin": adaptive_margin,
        "fixed_k12_margin": fixed_k12_margin,
        "threshold": threshold,
        "threshold_first": threshold_first,
        "tail_first": tail_first,
        "controls": len(control_speeds),
    }


def profile_body_star(body: tuple[int, ...]) -> dict[str, object]:
    """Pickle-stable worker entry point."""

    return profile_body(body)


def extrema_key_fraction_body(
    row: dict[str, object],
    field: str,
) -> tuple[F, tuple[int, ...]]:
    return row[field], row["body"]


def ledger_line(row: dict[str, object]) -> str:
    """Canonical complete-ledger serialization of one root."""

    return (
        f"S={row['stratum']};"
        f"E={','.join(map(str, row['body']))};"
        f"m={ftext(row['m'])};r={row['r']};"
        f"K={row['adaptive_k']};"
        f"margin={ftext(row['adaptive_margin'])};"
        f"K12={ftext(row['fixed_k12_margin'])};"
        f"threshold={ftext(row['threshold'])};"
        f"threshold_first={row['threshold_first']};"
        f"tail_first={row['tail_first']};"
        + "top="
        + ",".join(
            f"{rank}:{speed}:{ftext(value)}"
            for rank, (value, speed) in enumerate(row["top"], start=1)
        )
        + "\n"
    )


def extrema_for(rows: list[dict[str, object]]) -> tuple[object, ...]:
    """Return deterministic hostile extrema for one nonempty root set."""

    require(rows, "cannot extract extrema from an empty root set")
    minimum_adaptive = min(
        rows,
        key=lambda row: extrema_key_fraction_body(row, "adaptive_margin"),
    )
    minimum_fixed_k12 = min(
        rows,
        key=lambda row: extrema_key_fraction_body(row, "fixed_k12_margin"),
    )
    maximum_k = max(
        rows,
        key=lambda row: (row["adaptive_k"], tuple(-x for x in row["body"])),
    )
    maximum_threshold = max(
        rows,
        key=lambda row: (row["threshold"], tuple(-x for x in row["body"])),
    )
    maximum_tail = max(
        rows,
        key=lambda row: (row["tail_first"], tuple(-x for x in row["body"])),
    )
    maximum_retained = max(
        (
            speed,
            row["body"],
            rank,
            value,
        )
        for row in rows
        for rank, (value, speed) in enumerate(row["top"], start=1)
    )
    return (
        (
            ftext(minimum_adaptive["adaptive_margin"]),
            minimum_adaptive["body"],
            minimum_adaptive["adaptive_k"],
        ),
        (
            ftext(minimum_fixed_k12["fixed_k12_margin"]),
            minimum_fixed_k12["body"],
        ),
        (
            maximum_k["adaptive_k"],
            maximum_k["body"],
            ftext(maximum_k["adaptive_margin"]),
        ),
        (
            ftext(maximum_threshold["threshold"]),
            maximum_threshold["body"],
            maximum_threshold["threshold_first"],
        ),
        (
            maximum_tail["tail_first"],
            maximum_tail["body"],
            ftext(maximum_tail["threshold"]),
        ),
        (
            maximum_retained[0],
            maximum_retained[1],
            maximum_retained[2],
            ftext(maximum_retained[3]),
        ),
    )


def summarize(profiles: list[dict[str, object]]) -> dict[str, object]:
    """Build exact counts, extrema, and a complete ledger digest."""

    require(
        tuple(row["body"] for row in profiles) == BODIES,
        "worker output order/universe changed",
    )
    strata = ("low", "one", "both")
    stratum_counts = tuple(
        (name, sum(row["stratum"] == name for row in profiles))
        for name in strata
    )
    require(
        stratum_counts == (("low", 792), ("one", 1848), ("both", 792)),
        "seven-body stratum counts changed",
    )

    distributions = tuple(
        (
            name,
            tuple(
                sorted(
                    Counter(
                        row["adaptive_k"]
                        for row in profiles
                        if row["stratum"] == name
                    ).items()
                )
            ),
        )
        for name in strata
    )
    global_distribution = tuple(
        sorted(Counter(row["adaptive_k"] for row in profiles).items())
    )
    fixed_k12_success = sum(
        row["fixed_k12_margin"] > 0 for row in profiles
    )
    extrema = extrema_for(profiles)
    stratified_extrema = tuple(
        (
            name,
            extrema_for([row for row in profiles if row["stratum"] == name]),
        )
        for name in strata
    )
    total_controls = sum(row["controls"] for row in profiles)

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-root-adaptive-top30-gate/v1\n")
    stratum_hashes = {
        name: hashlib.sha256(
            f"LRC14/j6/all-root-adaptive-top30-gate/{name}/v1\n".encode()
        )
        for name in strata
    }
    for row in profiles:
        encoded = ledger_line(row).encode()
        digest.update(encoded)
        stratum_hashes[row["stratum"]].update(encoded)
    stratified_digests = tuple(
        (name, stratum_hashes[name].hexdigest()) for name in strata
    )

    gate_counts = (
        len(profiles),
        fixed_k12_success,
        len(profiles) - fixed_k12_success,
        extrema[2][0],
        extrema[4][0],
        extrema[5][0],
        total_controls,
    )
    return {
        "stratum_counts": stratum_counts,
        "global_distribution": global_distribution,
        "distributions": distributions,
        "gate_counts": gate_counts,
        "extrema": extrema,
        "stratified_extrema": stratified_extrema,
        "digest": digest.hexdigest(),
        "stratified_digests": stratified_digests,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(
        len(BODIES) == 3432 and len(set(BODIES)) == 3432,
        "seven-body universe changed",
    )

    if args.workers == 1:
        profiles = [profile_body(body) for body in BODIES]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            profiles = list(
                pool.imap(profile_body_star, BODIES, chunksize=8)
            )
    summary = summarize(profiles)

    if EXPECTED_GATE_COUNTS is not None:
        require(
            summary["gate_counts"] == EXPECTED_GATE_COUNTS,
            "gate counts changed",
        )
    if EXPECTED_K_DISTRIBUTION is not None:
        require(
            summary["distributions"] == EXPECTED_K_DISTRIBUTION,
            "stratified K distribution changed",
        )
    if EXPECTED_EXTREMA is not None:
        require(summary["extrema"] == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_STRATIFIED_EXTREMA is not None:
        require(
            summary["stratified_extrema"] == EXPECTED_STRATIFIED_EXTREMA,
            "stratified extrema changed",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            summary["digest"] == EXPECTED_LEDGER_DIGEST,
            "complete ledger digest changed",
        )
    if EXPECTED_STRATIFIED_DIGESTS is not None:
        require(
            summary["stratified_digests"] == EXPECTED_STRATIFIED_DIGESTS,
            "stratified ledger digests changed",
        )

    print("LRC14 j6 all-root adaptive top30 gate atlas")
    print(f"universe={len(profiles)}")
    print(f"stratum_counts={summary['stratum_counts']}")
    print(f"global_K_distribution={summary['global_distribution']}")
    print(f"stratified_K_distribution={summary['distributions']}")
    print(f"gate_counts={summary['gate_counts']}")
    print(f"extrema={summary['extrema']}")
    print(f"stratified_extrema={summary['stratified_extrema']}")
    print(f"complete_ledger_sha256={summary['digest']}")
    print(f"stratified_ledger_sha256={summary['stratified_digests']}")
    if EXPECTED_GATE_COUNTS is None:
        print("mode=DISCOVERY (expected constants are not yet locked)")
    else:
        print("mode=LOCKED")
    print("scope=hitting gate only; no apex carrier closure and no LRC14")


if __name__ == "__main__":
    main()
