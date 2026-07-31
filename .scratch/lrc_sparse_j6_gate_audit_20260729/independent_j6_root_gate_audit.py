#!/usr/bin/env python3
"""Independent exact audit of the scoped seven-body adaptive root gate.

This script intentionally does not import any repository LRC computation
module.  It rebuilds the good set as the complement of a rational interval
union and evaluates every single-comb coverage with a scalar periodic
antiderivative.  Its purpose is to cross-check the vectorized implementation
and the frozen 36-root transcript in
``lrc14_j6_adaptive_gate_triple_heavy_battery_codex_20260729.py``.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from itertools import combinations


FIRST_EXTERNAL = 15
BASE_HORIZON = 1600
TOP_COUNT = 30
COVER_SIZE = 6
SAMPLE_PER_STRATUM = 12
S2 = F(99, 70)
ONE = F(1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def merge(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    merged: list[list[F]] = []
    for left, right in sorted(intervals):
        require(F(0) <= left <= right <= F(1), "interval left the unit circle")
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1][1] = right
        else:
            merged.append([left, right])
    return [(left, right) for left, right in merged]


def danger(speed: int) -> list[tuple[F, F]]:
    """Return D_speed as disjoint closed rational intervals in [0,1]."""

    radius = F(1, 14 * speed)
    pieces: list[tuple[F, F]] = []
    for residue in range(speed):
        center = F(residue, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            pieces.append((F(0), right))
            pieces.append((left + 1, F(1)))
        elif right > 1:
            pieces.append((left, F(1)))
            pieces.append((F(0), right - 1))
        else:
            pieces.append((left, right))
    return merge(pieces)


def good_set(body: tuple[int, ...]) -> list[tuple[F, F]]:
    bad = merge([interval for speed in body for interval in danger(speed)])
    good: list[tuple[F, F]] = []
    cursor = F(0)
    for left, right in bad:
        if cursor < left:
            good.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        good.append((cursor, F(1)))
    return good


def tooth_antiderivative(value: F) -> F:
    """Integral from 0 to value of 1_{||x|| <= 1/14}."""

    integer = value.numerator // value.denominator
    residue = value - integer
    return (
        F(integer, 7)
        + min(residue, F(1, 14))
        + max(F(0), residue - F(13, 14))
    )


def coverage(good: list[tuple[F, F]], speed: int) -> F:
    return sum(
        (
            tooth_antiderivative(speed * right)
            - tooth_antiderivative(speed * left)
        )
        / speed
        for left, right in good
    )


def quantile_sample(
    rows: list[tuple[int, ...]],
    count: int,
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        rows[(len(rows) - 1) * index // (count - 1)]
        for index in range(count)
    )


def sample_universes() -> tuple[tuple[str, tuple[tuple[int, ...], ...]], ...]:
    low = list(combinations(range(1, 13), 7))
    one = [
        tuple(sorted((*base, special)))
        for base in combinations(range(1, 13), 6)
        for special in (13, 14)
    ]
    both = [(*base, 13, 14) for base in combinations(range(1, 13), 5)]
    require(
        (len(low), len(one), len(both)) == (792, 1848, 792),
        "seven-body strata changed",
    )
    return (
        ("low", quantile_sample(low, SAMPLE_PER_STRATUM)),
        ("one", quantile_sample(one, SAMPLE_PER_STRATUM)),
        ("both", quantile_sample(both, SAMPLE_PER_STRATUM)),
    )


def profile(body: tuple[int, ...]) -> dict[str, object]:
    good = good_set(body)
    mass = sum((right - left for left, right in good), F(0))
    components = len(good)
    require(mass > 0, f"empty root carrier: {body}")

    rows = [
        (coverage(good, speed), speed)
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
    ]
    ranked_base = sorted(rows, key=lambda item: (-item[0], item[1]))
    q30 = ranked_base[TOP_COUNT - 1][0]
    require(q30 > mass / 7, f"rank thirty misses discrepancy limit: {body}")
    threshold = S2 * components / (7 * (q30 - mass / 7))
    threshold_first = ceiling(threshold)
    tail_first = max(BASE_HORIZON + 1, threshold_first)
    rows.extend(
        (coverage(good, speed), speed)
        for speed in range(BASE_HORIZON + 1, tail_first)
    )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q30,
        f"tail seal failed: {body}",
    )
    top = tuple(sorted(rows, key=lambda item: (-item[0], item[1]))[:TOP_COUNT])

    fixed_margin = mass - sum(
        (value for value, _ in top[12:18]),
        F(0),
    )
    adaptive_k = None
    adaptive_margin = None
    for k in range(TOP_COUNT - COVER_SIZE + 1):
        margin = mass - sum(
            (value for value, _ in top[k : k + COVER_SIZE]),
            F(0),
        )
        if margin > 0:
            adaptive_k = k
            adaptive_margin = margin
            break
    require(
        adaptive_k is not None and adaptive_margin is not None,
        f"no top-thirty gate: {body}",
    )
    return {
        "body": body,
        "good": good,
        "m": mass,
        "r": components,
        "top": top,
        "fixed_margin": fixed_margin,
        "adaptive_k": adaptive_k,
        "adaptive_margin": adaptive_margin,
        "threshold_first": threshold_first,
        "tail_first": tail_first,
    }


def main() -> None:
    rows: list[tuple[str, dict[str, object]]] = []
    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/adaptive-root-gate-sample/v1\n")
    for stratum, bodies in sample_universes():
        for body in bodies:
            row = profile(body)
            rows.append((stratum, row))
            digest.update(
                (
                    f"S={stratum};E={','.join(map(str, body))};"
                    f"K12={ftext(row['fixed_margin'])};"
                    f"K={row['adaptive_k']};"
                    f"margin={ftext(row['adaptive_margin'])};"
                    f"threshold={row['threshold_first']};"
                    f"tail={row['tail_first']}\n"
                ).encode()
            )

    counts = tuple(
        (
            stratum,
            sum(name == stratum for name, _ in rows),
            sum(
                name == stratum and row["fixed_margin"] > 0
                for name, row in rows
            ),
            max(
                row["adaptive_k"]
                for name, row in rows
                if name == stratum
            ),
        )
        for stratum in ("low", "one", "both")
    )
    failures = [
        row for _, row in rows if row["fixed_margin"] <= 0
    ]
    first = failures[0]
    expected_counts = (
        ("low", 12, 10, 16),
        ("one", 12, 5, 20),
        ("both", 12, 5, 21),
    )
    expected_digest = (
        "ce3be5476ee287f7d972d1853f4e9fe52175583347832a75f103fd6e9172ba51"
    )
    require(counts == expected_counts, "sample counts disagree with transcript")
    require(len(failures) == 16, "fixed top-twelve failure count changed")
    require(
        first["body"] == (1, 4, 7, 8, 10, 11, 12)
        and first["fixed_margin"] == F(-238717, 34474440)
        and first["adaptive_k"] == 16,
        "first hostile sample changed",
    )
    require(digest.hexdigest() == expected_digest, "sample digest changed")

    print("INDEPENDENT J6 ROOT-GATE AUDIT")
    print("engine=rational interval union + scalar periodic antiderivative")
    print(f"sample_counts={counts}")
    print(f"fixed_top12_failures={len(failures)}/36")
    print(
        f"first_failure={first['body']};"
        f"mass={ftext(first['m'])};"
        f"margin={ftext(first['fixed_margin'])};"
        f"adaptive_K={first['adaptive_k']}"
    )
    print(
        "first_failure_ranks_13_to_18="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in first["top"][12:18]
        )
    )
    print(
        f"maximum_threshold={max(row['threshold_first'] for _, row in rows)};"
        f"maximum_tail_first={max(row['tail_first'] for _, row in rows)}"
    )
    print(f"sample_digest_sha256={digest.hexdigest()}")
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
