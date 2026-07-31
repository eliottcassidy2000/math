#!/usr/bin/env python3
"""Probe the naive Hunter-star recursion at the critical seven-slot seam.

This is discovery telemetry, not a theorem verifier.  It asks whether the
same pair-cap/Hunter mechanism that closes the seven-body/six-slot rung can
make a finite first-centre core one rung earlier, on

    E in C({1,...,14}, 6),  seven external labels.

At seven slots the critical obstruction is exact: if q_7 >= h/7 and
beta_2 >= 2h/7, then G_7(h/7)=h and G_7(a)<h for a<h/7, so the first hostile
level is exactly h/7.  The THM-735 discrepancy tail does not make that
superlevel finite.  This program measures how uniformly that obstruction
occurs across the 3,003 six-body roots.
"""

from __future__ import annotations

import argparse
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE1 = (
    ROOT
    / "04-computation/"
    "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage1.py"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_stage1():
    spec = importlib.util.spec_from_file_location("j7_critical_stage1", STAGE1)
    require(spec is not None and spec.loader is not None, "cannot load stage1")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_stage1()


def hunter_data(
    qs: tuple[F, ...],
    pair_cap: F,
    mass: F,
) -> tuple[F, F | None]:
    require(len(qs) == 7, "critical probe requires seven singleton ranks")
    require(
        all(qs[index] >= qs[index + 1] for index in range(6)),
        "singleton ranks are not decreasing",
    )
    upper = min(qs[0], pair_cap)

    def clipped(value: F) -> F:
        return min(max(value, F(0)), upper)

    breakpoints = {F(0), upper, clipped(pair_cap / 2)}
    for singleton in qs[1:]:
        breakpoints.add(clipped(singleton))
        breakpoints.add(clipped(pair_cap - singleton))
    ordered = sorted(breakpoints)

    def objective(center: F) -> F:
        return center + sum(
            (
                min(center, singleton, pair_cap - center)
                for singleton in qs[1:]
            ),
            F(0),
        )

    envelope = max(objective(point) for point in ordered)
    if envelope < mass:
        return envelope, None
    for left, right in zip(ordered, ordered[1:]):
        left_value = objective(left)
        right_value = objective(right)
        if left_value >= mass:
            return envelope, left
        if left_value < mass <= right_value:
            threshold = (
                left
                + (mass - left_value)
                * (right - left)
                / (right_value - left_value)
            )
            require(objective(threshold) == mass, "threshold solve failed")
            return envelope, threshold
    require(objective(ordered[-1]) >= mass, "hostile level disappeared")
    return envelope, ordered[-1]


def analyze(body: tuple[int, ...]) -> dict[str, object]:
    carrier, components, mass = S.direct_carrier(body)
    top7 = S.exact_top_k(
        carrier,
        frozenset(),
        rank=7,
        initial_horizon=S.INITIAL_TOP3_HORIZON,
    )
    pair = S.exact_pair_cap(
        carrier,
        frozenset(),
        initial_horizon=S.INITIAL_PAIR_HORIZON,
    )
    qs = tuple(value for value, _ in top7["top"])
    envelope, threshold = hunter_data(qs, pair["head"], mass)
    return {
        "body": body,
        "components": components,
        "mass": mass,
        "q7_gap": qs[-1] - mass / 7,
        "pair_gap": pair["head"] - 2 * mass / 7,
        "envelope_gap": envelope - mass,
        "lambda_gap": None if threshold is None else threshold - mass / 7,
        "pair_horizon": pair["horizon"],
        "top7_horizon": top7["horizon"],
    }


def ftext(value: F | None) -> str:
    if value is None:
        return "None"
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    bodies = list(combinations(range(1, 15), 6))
    if args.limit is not None:
        bodies = bodies[: args.limit]
    if args.workers == 1:
        rows = list(map(analyze, bodies))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = pool.map(analyze, bodies, chunksize=1)
    rows.sort(key=lambda row: row["body"])

    classes = Counter(
        "closed"
        if row["lambda_gap"] is None
        else "positive"
        if row["lambda_gap"] > 0
        else "zero"
        if row["lambda_gap"] == 0
        else "negative"
        for row in rows
    )
    print(f"roots={len(rows)}")
    print(f"lambda_classes={tuple(sorted(classes.items()))}")
    for field in ("q7_gap", "pair_gap", "envelope_gap"):
        least = min(rows, key=lambda row: (row[field], row["body"]))
        greatest = max(rows, key=lambda row: (row[field], row["body"]))
        print(
            f"{field}_min={ftext(least[field])};body={least['body']}"
        )
        print(
            f"{field}_max={ftext(greatest[field])};body={greatest['body']}"
        )
    zero_under_critical_hypotheses = sum(
        row["q7_gap"] >= 0
        and row["pair_gap"] >= 0
        and row["lambda_gap"] == 0
        for row in rows
    )
    print(
        "critical_hypotheses_and_zero_lambda="
        f"{zero_under_critical_hypotheses}"
    )
    print(
        "maximum_pair_horizon="
        f"{max(row['pair_horizon'] for row in rows)}"
    )
    print(
        "maximum_top7_horizon="
        f"{max(row['top7_horizon'] for row in rows)}"
    )


if __name__ == "__main__":
    main()
