#!/usr/bin/env python3
"""Exact kernel for the THM-2893 (k,s,ell)=(5,4,2) j=6 route.

For each selected first-apex carrier C, let q1 be the globally maximal
allowed singleton coverage.  When q1 < 3|C|/7, THM-2893 gives the finite
core

    H4 = {w : |C intersect D_w| >= (|C|-q1)/4}.

Every hypothetical five-cover contains at least two H4 labels.  For every
H4 pair this script constructs the literal residual R and tests the cheap
global triple-cover obstruction

    B2(R) + q1(R) < |R|,

where B2(R) is a globally tail-sealed exact pair-union cap.  The three
carriers are the hostile rank-one cases from the scoped j=6 battery.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SINGLE_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
SINGLE_SHA = "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
SUBTRACT_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
SUBTRACT_SHA = "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"

FIRST_EXTERNAL = 15
HORIZON = 2500
S2 = F(99, 70)
CASES = (
    ((2, 8, 9, 10, 11, 13, 14), 19),
    ((1, 3, 9, 10, 11, 12, 14), 39),
    ((2, 5, 9, 11, 12, 13, 14), 16),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected: str, name: str):
    require(file_sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load(SINGLE_PATH, SINGLE_SHA, "j6_h4_single")
R = load(SUBTRACT_PATH, SUBTRACT_SHA, "j6_h4_subtract")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def globally_ranked(
    carrier: list[tuple[F, F]],
    excluded: set[int],
) -> tuple[list[tuple[F, int]], F]:
    """Return the exact finite head and a strict omitted singleton majorant."""

    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, HORIZON + 1)
        if speed not in excluded
    ]
    rows = T.coverages_many(carrier, speeds)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    by_speed = {speed: value for value, speed in rows}
    carrier_mass = mass(carrier)
    tail = carrier_mass / 7 + S2 * len(carrier) / (7 * (HORIZON + 1))
    require(
        ranked and tail <= ranked[0][0],
        "rank-one finite head did not seal",
    )
    for speed in (speeds[0], speeds[-1], ranked[0][1]):
        require(
            by_speed[speed] == T.coverage(carrier, speed),
            f"vector/scalar disagreement at {speed}",
        )
    return ranked, tail


def pair_cap(
    carrier: list[tuple[F, F]],
    excluded: set[int],
) -> dict[str, object]:
    """Exact finite-head pair maximum plus a strict infinite-tail cap."""

    ranked, tail_single = globally_ranked(carrier, excluded)
    carrier_mass = mass(carrier)
    q1 = ranked[0][0]
    head = F(0)
    maximizer = None
    paid = 0
    paid_hash = hashlib.sha256(b"LRC14/j6/H4/residual-paid-pairs/v1\n")
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= head:
            break
        after_first = R.subtract_local(carrier, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= head:
                break
            survivor = R.subtract_local(after_first, second)
            survivor_mass = mass(survivor)
            union = carrier_mass - survivor_mass
            paid += 1
            paid_hash.update(
                (
                    f"P={min(first, second)},{max(first, second)};"
                    f"U={ftext(union)};L={ftext(survivor_mass)}\n"
                ).encode()
            )
            if union > head:
                head = union
                maximizer = tuple(sorted((first, second)))
    require(paid > 0 and maximizer is not None, "empty pair maximization")
    tail_pair = q1 + tail_single
    return {
        "q1": q1,
        "ranked": ranked,
        "tail_single": tail_single,
        "head": head,
        "tail": tail_pair,
        "cap": max(head, tail_pair),
        "maximizer": maximizer,
        "paid": paid,
        "paid_digest": paid_hash.hexdigest(),
        "top3": tuple(ranked[:3]),
    }


def exact_triple_cap(
    carrier: list[tuple[F, F]],
    pair_data: dict[str, object],
) -> dict[str, object]:
    """Exact finite triple maximum, sealed by B2 plus one tail singleton."""

    ranked = pair_data["ranked"]
    carrier_mass = mass(carrier)
    head = F(0)
    maximizer = None
    paid = 0
    for first_index, (first_value, first) in enumerate(ranked[:-2]):
        if (
            first_value
            + ranked[first_index + 1][0]
            + ranked[first_index + 2][0]
            <= head
        ):
            break
        after_first = R.subtract_local(carrier, first)
        for second_index in range(first_index + 1, len(ranked) - 1):
            second_value, second = ranked[second_index]
            if (
                first_value
                + second_value
                + ranked[second_index + 1][0]
                <= head
            ):
                break
            after_second = R.subtract_local(after_first, second)
            for third_value, third in ranked[second_index + 1 :]:
                if first_value + second_value + third_value <= head:
                    break
                survivor = R.subtract_local(after_second, third)
                union = carrier_mass - mass(survivor)
                paid += 1
                if union > head:
                    head = union
                    maximizer = tuple(sorted((first, second, third)))
    require(paid > 0 and maximizer is not None, "empty triple maximization")
    tail = pair_data["cap"] + pair_data["tail_single"]
    cap = max(head, tail)
    return {
        "head": head,
        "tail": tail,
        "cap": cap,
        "margin": carrier_mass - cap,
        "maximizer": maximizer,
        "paid": paid,
    }


def h4_core(
    body: tuple[int, ...],
    apex: int,
) -> dict[str, object]:
    carrier, components, carrier_mass = T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    excluded = {apex}
    ranked, tail_single = globally_ranked(carrier, excluded)
    q1 = ranked[0][0]
    singleton_margin = F(3, 7) * carrier_mass - q1
    require(singleton_margin > 0, f"B1 reaches 3h/7: {body}, {apex}")
    level = (carrier_mass - q1) / 4
    require(level > carrier_mass / 7, f"H4 threshold is not finite: {body}")
    threshold = S2 * components / (7 * (level - carrier_mass / 7))
    tail_first = max(FIRST_EXTERNAL, ceiling(threshold))
    if tail_first <= HORIZON + 1:
        by_speed = {speed: value for value, speed in ranked}
    else:
        extra = [
            speed
            for speed in range(HORIZON + 1, tail_first)
            if speed not in excluded
        ]
        by_speed = {speed: value for value, speed in ranked}
        by_speed.update(
            {
                speed: value
                for value, speed in T.coverages_many(carrier, extra)
            }
        )
    require(
        carrier_mass / 7 + S2 * components / (7 * tail_first) <= level,
        f"H4 tail did not seal: {body}, {apex}",
    )
    core = tuple(
        sorted(
            speed
            for speed, value in by_speed.items()
            if speed not in excluded and value >= level
        )
    )
    return {
        "body": body,
        "apex": apex,
        "carrier": carrier,
        "m": carrier_mass,
        "r": components,
        "q1": q1,
        "tail_single": tail_single,
        "singleton_margin": singleton_margin,
        "level": level,
        "tail_first": tail_first,
        "H": core,
    }


def probe_pair(
    root: dict[str, object],
    pair: tuple[int, int],
) -> dict[str, object]:
    residual = R.subtract_local_multi(root["carrier"], pair)
    residual_mass = mass(residual)
    require(residual_mass > 0, f"H4 pair covers carrier: {pair}")
    family = tuple(sorted((*root["body"], root["apex"], *pair)))
    direct, direct_r, direct_m = T.CORE.good_norm(family)
    require(
        residual == direct
        and len(residual) == direct_r
        and residual_mass == direct_m,
        f"literal/direct pair residual disagreement: {family}",
    )
    cap = pair_cap(residual, {root["apex"], *pair})
    q3_tail_margin = cap["top3"][2][0] - cap["tail_single"]
    require(
        q3_tail_margin >= 0,
        f"finite top three did not seal: {root['body']}, {pair}",
    )
    direct3_margin = residual_mass - sum(
        (value for value, _ in cap["top3"]),
        F(0),
    )
    cheap_margin = residual_mass - cap["cap"] - cap["q1"]
    return {
        "pair": pair,
        "m": residual_mass,
        "r": len(residual),
        "direct3_margin": direct3_margin,
        "q3_tail_margin": q3_tail_margin,
        "cheap_margin": cheap_margin,
        "q1": cap["q1"],
        "B2": cap["cap"],
        "maxpair": cap["maximizer"],
        "paid": cap["paid"],
        "residual": residual,
        "pair_data": cap,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-only", action="store_true")
    parser.add_argument("--exact-failures", action="store_true")
    args = parser.parse_args()

    roots = [h4_core(body, apex) for body, apex in CASES]
    print("J6 (k,s,ell)=(5,4,2) H4 PAIR-RESIDUAL PROBE")
    for root in roots:
        print(
            f"E={root['body']};apex={root['apex']};"
            f"h={ftext(root['m'])};r={root['r']};"
            f"q1={ftext(root['q1'])};"
            f"3h7-q1={ftext(root['singleton_margin'])};"
            f"level={ftext(root['level'])};"
            f"Htail={root['tail_first']};H={len(root['H'])};"
            f"Hlabels={root['H']}"
        )
    require(
        tuple(len(root["H"]) for root in roots) == (16, 12, 16),
        "reported H4 core sizes changed",
    )
    if args.h_only:
        print("scope=H4 only")
        print("all_exact_controls=PASS")
        return

    total_pairs = 0
    total_direct = 0
    total_cheap = 0
    total_paid = 0
    all_rows: list[dict[str, object]] = []
    failures: list[tuple[tuple[int, ...], int, dict[str, object]]] = []
    adaptive_failures: list[
        tuple[tuple[int, ...], int, dict[str, object]]
    ] = []
    sharp = None
    for root in roots:
        rows = [
            probe_pair(root, pair)
            for pair in combinations(root["H"], 2)
        ]
        all_rows.extend(rows)
        total_pairs += len(rows)
        total_direct += sum(row["direct3_margin"] > 0 for row in rows)
        total_cheap += sum(row["cheap_margin"] > 0 for row in rows)
        total_paid += sum(row["paid"] for row in rows)
        failures.extend(
            (root["body"], root["apex"], row)
            for row in rows
            if row["cheap_margin"] <= 0
        )
        adaptive_failures.extend(
            (root["body"], root["apex"], row)
            for row in rows
            if (
                row["direct3_margin"] <= 0
                and row["cheap_margin"] <= 0
            )
        )
        local_sharp = min(
            rows,
            key=lambda row: (row["cheap_margin"], row["pair"]),
        )
        candidate = (
            local_sharp["cheap_margin"],
            root["body"],
            root["apex"],
            local_sharp,
        )
        if sharp is None or candidate[:3] < sharp[:3]:
            sharp = candidate
        print(
            f"E={root['body']};apex={root['apex']};"
            f"pairs={len(rows)};"
            f"direct_top3={sum(row['direct3_margin'] > 0 for row in rows)};"
            f"cheap_B2_plus_q1={sum(row['cheap_margin'] > 0 for row in rows)};"
            f"minimum_cheap_margin={ftext(local_sharp['cheap_margin'])};"
            f"pair={local_sharp['pair']};"
            f"B2max={local_sharp['maxpair']}"
        )

    require(sharp is not None, "empty H4 pair universe")
    print(
        f"totals=pairs:{total_pairs},direct_top3:{total_direct},"
        f"cheap_B2_plus_q1:{total_cheap},"
        f"adaptive_union:{total_pairs-len(adaptive_failures)},"
        f"cheap_failures:{len(failures)},"
        f"adaptive_failures:{len(adaptive_failures)},"
        f"paid_pairs:{total_paid}"
    )
    print(
        "minimum_q3_tail_margin="
        + ftext(min(row["q3_tail_margin"] for row in all_rows))
    )
    print(
        f"global_minimum_cheap_margin={ftext(sharp[0])};"
        f"E={sharp[1]};apex={sharp[2]};"
        f"Hpair={sharp[3]['pair']};B2max={sharp[3]['maxpair']}"
    )
    if failures:
        for body, apex, row in failures[:12]:
            print(
                f"cheap_failure_E={body};apex={apex};Hpair={row['pair']};"
                f"margin={ftext(row['cheap_margin'])};"
                f"B2max={row['maxpair']}"
            )
    if args.exact_failures:
        exact_rows = [
            (
                body,
                apex,
                row,
                exact_triple_cap(row["residual"], row["pair_data"]),
            )
            for body, apex, row in adaptive_failures
        ]
        print(
            f"exact_failure_triples={len(exact_rows)};"
            f"closed={sum(exact['margin'] > 0 for _, _, _, exact in exact_rows)}"
        )
        for body, apex, row, exact in exact_rows:
            print(
                f"exact_E={body};apex={apex};Hpair={row['pair']};"
                f"margin={ftext(exact['margin'])};"
                f"maxtriple={exact['maximizer']};paid={exact['paid']};"
                f"tail_active={int(exact['tail'] >= exact['head'])}"
            )
    print("scope=three hostile rank-one carriers only;not uniform")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
