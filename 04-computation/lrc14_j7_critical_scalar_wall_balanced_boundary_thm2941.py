#!/usr/bin/env python3
"""Exact primary verifier for the critical seven-slot scalar wall.

The universe is

    E in C({1,...,14}, 6),          seven external labels w >= 15.

For the literal carrier G_E, the verifier seals the global top seven
singleton coverages and the global pair-union cap with the strict THM-735
tail.  It then evaluates the seven-slot pair-Hunter envelope exactly.  The
same carrier reconstruction records the longest-component obstruction to an
exact balanced partition.

This verifier intentionally reuses the independently locked interval engine
from THM-2923.  THM-2941 also has a separate integer-ruler implementation
which does not import that engine.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALIGNED_SAFE_SURPLUS = (
    (2, F(1, 91)),
    (3, F(3, 91)),
    (4, F(51, 1183)),
    (5, F(88, 1365)),
    (6, F(22, 273)),
)
STAGE1 = (
    ROOT
    / "04-computation"
    / "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage1.py"
)
EXPECTED_STAGE1_SHA256 = (
    "32751b2a5beb789b1657f06d3964cbf24e634251e57f57f299b8f50c647f2103"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge/results"
    / "lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.out"
)
EXPECTED_SEMANTIC_SHA256: str | None = (
    "a1c857a75356401a6e268ef79055b4806ccc7559620f28a592803e31273f5528"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load_stage1():
    require(
        file_sha256(STAGE1) == EXPECTED_STAGE1_SHA256,
        "THM-2923 stage-one verifier changed",
    )
    spec = importlib.util.spec_from_file_location("thm2941_stage1", STAGE1)
    require(spec is not None and spec.loader is not None, "cannot load stage one")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_stage1()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return -(-value.numerator // value.denominator)


def hunter_data(
    ranks: tuple[F, ...],
    pair_cap: F,
    mass: F,
) -> tuple[F, F]:
    require(len(ranks) == 7, "seven-slot envelope has the wrong arity")
    require(
        all(ranks[index] >= ranks[index + 1] for index in range(6)),
        "singleton ranks are not nonincreasing",
    )
    upper = min(ranks[0], pair_cap)

    def clipped(value: F) -> F:
        return min(max(value, F(0)), upper)

    breakpoints = {F(0), upper, clipped(pair_cap / 2)}
    for singleton in ranks[1:]:
        breakpoints.add(clipped(singleton))
        breakpoints.add(clipped(pair_cap - singleton))
    ordered = tuple(sorted(breakpoints))

    def objective(center: F) -> F:
        return center + sum(
            (
                min(center, singleton, pair_cap - center)
                for singleton in ranks[1:]
            ),
            F(0),
        )

    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        require(
            2 * objective(midpoint) == objective(left) + objective(right),
            "Hunter breakpoint list missed an affine change",
        )
    envelope = max(objective(point) for point in ordered)
    require(envelope >= mass, "critical envelope unexpectedly closes")

    threshold: F | None = None
    if objective(ordered[0]) >= mass:
        threshold = ordered[0]
    else:
        for left, right in zip(ordered, ordered[1:]):
            left_value = objective(left)
            right_value = objective(right)
            if left_value < mass <= right_value:
                threshold = (
                    left
                    + (mass - left_value)
                    * (right - left)
                    / (right_value - left_value)
                )
                require(
                    objective(threshold) == mass,
                    "affine first-crossing interpolation failed",
                )
                break
    require(threshold is not None, "critical first crossing disappeared")
    return envelope, threshold


def profile(body: tuple[int, ...]) -> dict[str, object]:
    require(len(body) == 6 and len(set(body)) == 6, "body arity changed")
    carrier, components, mass = S.direct_carrier(body)
    require(carrier and components == len(carrier) and mass > 0, "empty carrier")

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
    ranks = tuple(value for value, _ in top7["top"])
    q7_gap = ranks[-1] - mass / 7
    pair_gap = pair["head"] - 2 * mass / 7
    require(q7_gap > 0, f"{body}: q7 did not clear h/7")
    require(pair_gap > 0, f"{body}: pair cap did not clear 2h/7")
    # The pair seal must also control a pair with two omitted endpoints.
    # The finite winner dominates q1+tau; subadditivity bounds it by 2q1.
    # Hence tau<=q1 and 2tau<=q1+tau<=B2.
    q1 = ranks[0]
    tail_singleton = pair["singleton_tail"]
    require(
        pair["head"] >= q1
        and pair["head"] >= q1 + tail_singleton
        and pair["head"] <= 2 * q1
        and tail_singleton <= q1
        and 2 * tail_singleton <= pair["head"],
        f"{body}: two-omitted-label pair seal failed",
    )

    envelope, threshold = hunter_data(ranks, pair["head"], mass)
    require(
        threshold == mass / 7,
        f"{body}: critical first crossing is not h/7",
    )

    longest = max(right - left for left, right in carrier)
    owner_bound = ceiling(F(1, 1) / (7 * longest)) - 1
    require(owner_bound >= 0, f"{body}: negative component-owner bound")
    aligned_first_drift_bounds = tuple(
        (
            k,
            F(6 * (7 - k) * components, 49) / (surplus * mass),
        )
        for k, surplus in ALIGNED_SAFE_SURPLUS
    )

    return {
        "body": body,
        "mass": mass,
        "components": components,
        "components_per_mass": F(components, 1) / mass,
        "top7": top7["top"],
        "top7_horizon": top7["horizon"],
        "top7_tail_gap": top7["tail_gap"],
        "q7_gap": q7_gap,
        "pair_cap": pair["head"],
        "pair_witness": pair["witness"],
        "pair_horizon": pair["horizon"],
        "pair_tail_gap": pair["tail_gap"],
        "pair_gap": pair_gap,
        "envelope": envelope,
        "envelope_gap": envelope - mass,
        "threshold": threshold,
        "longest": longest,
        "owner_bound": owner_bound,
        "aligned_first_drift_bounds": aligned_first_drift_bounds,
    }


def extreme(
    rows: list[dict[str, object]],
    field: str,
    largest: bool = False,
) -> dict[str, object]:
    if largest:
        return min(rows, key=lambda row: (-row[field], row["body"]))
    return min(rows, key=lambda row: (row[field], row["body"]))


def extremum_text(row: dict[str, object], field: str) -> str:
    value = row[field]
    rendered = ftext(value) if isinstance(value, F) else str(value)
    return f"{rendered};body={row['body']}"


def row_digest(rows: list[dict[str, object]]) -> str:
    payload = tuple(
        (
            row["body"],
            row["mass"],
            row["components"],
            row["components_per_mass"],
            row["top7"],
            row["top7_horizon"],
            row["pair_cap"],
            row["pair_witness"],
            row["pair_horizon"],
            row["envelope"],
            row["threshold"],
            row["longest"],
            row["owner_bound"],
            row["aligned_first_drift_bounds"],
        )
        for row in rows
    )
    return hashlib.sha256(
        b"LRC14/THM2941/critical-seven-slot/v1\n" + repr(payload).encode()
    ).hexdigest()


def render(rows: list[dict[str, object]]) -> str:
    require(len(rows) == math.comb(14, 6) == 3_003, "root universe changed")
    require(
        all(row["q7_gap"] > 0 for row in rows)
        and all(row["pair_gap"] > 0 for row in rows)
        and all(row["threshold"] == row["mass"] / 7 for row in rows),
        "critical all-root classification failed",
    )
    owner_histogram = tuple(sorted(Counter(row["owner_bound"] for row in rows).items()))
    require(
        owner_histogram == ((1, 39), (2, 1601), (3, 1049), (4, 283), (5, 28), (6, 3)),
        "balanced-boundary owner histogram changed",
    )
    minimum_longest = extreme(rows, "longest")
    aligned_apex = extreme(
        rows,
        "components_per_mass",
        largest=True,
    )
    require(
        minimum_longest["longest"] == F(23, 1092)
        and minimum_longest["body"] == (1, 6, 7, 8, 10, 13)
        and max(row["owner_bound"] for row in rows) == 6,
        "uniform balanced-boundary obstruction changed",
    )
    require(
        aligned_apex["components_per_mass"] == F(3_993_990, 32_029)
        and aligned_apex["body"] == (1, 10, 11, 12, 13, 14)
        and aligned_apex["mass"] == F(32_029, 105_105)
        and aligned_apex["components"] == 38
        and aligned_apex["aligned_first_drift_bounds"]
        == (
            (2, F(222_522_300, 32_029)),
            (3, F(59_339_280, 32_029)),
            (4, F(578_557_980, 544_493)),
            (5, F(15_171_975, 32_029)),
            (6, F(6_068_790, 32_029)),
        )
        and tuple(
            bound.numerator // bound.denominator
            for _, bound in aligned_apex["aligned_first_drift_bounds"]
        )
        == (6_947, 1_852, 1_062, 473, 189),
        "aligned-sector first-drift caps changed",
    )

    lines = [
        "LRC14 THM2941 critical seven-slot scalar wall",
        "universe=(six_body_roots=3003,body_labels=1..14,external_labels>=15)",
        "counts=(q7_gt_h_over_7=3003,B2_gt_2h_over_7=3003,"
        "lambda_eq_h_over_7=3003,scalar_closures=0)",
    ]
    for field in ("q7_gap", "pair_gap", "envelope_gap"):
        lines.append(
            f"{field}_min={extremum_text(extreme(rows, field), field)}"
        )
        lines.append(
            f"{field}_max={extremum_text(extreme(rows, field, True), field)}"
        )
    lines.extend(
        [
            "pair_tail_gap_min="
            + extremum_text(extreme(rows, "pair_tail_gap"), "pair_tail_gap"),
            f"maximum_pair_horizon={max(row['pair_horizon'] for row in rows)}",
            f"maximum_top7_scan_horizon={max(row['top7_horizon'] for row in rows)}",
            "minimum_longest_component="
            + extremum_text(minimum_longest, "longest"),
            f"critical_owner_bound_histogram={owner_histogram}",
            "aligned_safe_surplus_first_drift_bounds="
            + repr(
                tuple(
                    (k, ftext(bound))
                    for k, bound in aligned_apex["aligned_first_drift_bounds"]
                )
            )
            + ";integral_caps=(6947,1852,1062,473,189);"
            + f"max_root={aligned_apex['body']}",
            "balanced_pointwise_cover_boundary=EMPTY;zero excess assigns each "
            "closed positive carrier component to one open tooth, forcing "
            "an owner<=6<15",
            "hypothetical_balanced_partition_gram=(h/49)(7I-J);"
            "rank=6;kernel=span(all_ones)",
            "next_object=multi-slope restricted Gram plus component owner/"
            "address-transition ledger; scalar pair cap is insufficient",
            f"row_semantic_sha256={row_digest(rows)}",
            "scope=complete scalar-wall and balanced-boundary audit on the "
            "six-body/seven-tail rung; arbitrary positive-excess covers OPEN",
            "not_unrestricted_LRC14=TRUE",
        ]
    )
    body = "\n".join(lines) + "\n"
    semantic = hashlib.sha256(
        b"LRC14/THM2941/summary/v1\n" + body.encode()
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic == EXPECTED_SEMANTIC_SHA256,
            "THM2941 summary semantic digest changed",
        )
    return body + f"semantic_sha256={semantic}\nall_exact_controls=PASS\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)

    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        rows = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])
    output = render(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
