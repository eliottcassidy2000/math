#!/usr/bin/env python3
"""Third-level exact recursion on failed THM-2915 grandchildren."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = Path(__file__).with_name(
    "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage1.py"
)
EXPECTED_BASE_SHA256 = (
    "74e5eb2d0b23dda6366d115376890726ac32dc834a09a71479ea875f05e7615e"
)
FIRST_EXTERNAL = 15
INITIAL_HORIZON = 2_000
EXPECTED_GRANDCHILD_DIGEST = (
    "0b70905f638bb2a548e74a265df76506769f33424989c1385aa170b73c8a302a"
)
EXPECTED_THIRD_DIGEST = (
    "55c2e237e4b22475d2daa1dd08f8d8df25653952387c62813dd33edc2f66ae67"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    require(file_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "stage-one source changed")
    spec = importlib.util.spec_from_file_location("third_level_base", BASE_PATH)
    require(
        spec is not None and spec.loader is not None,
        "cannot load stage-one verifier",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = load_base()
H = G.H
T = G.T


def ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def fields(line: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in line.rstrip().split(";")[1:])


def first_identity(data: dict[str, str]) -> tuple[object, ...]:
    return (
        ints(data["E"]),
        int(data["rank"]),
        int(data["a"]),
        ints(data["P"]),
        int(data["x"]),
        ints(data["earlier"]),
    )


def parse_pair_rows(path: Path) -> dict[tuple[object, ...], dict[str, object]]:
    rows: dict[tuple[object, ...], dict[str, object]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        data = fields(line)
        key = first_identity(data)
        require(key not in rows, "duplicate pair identity")
        rows[key] = {
            "body": key[0],
            "rank": key[1],
            "apex": key[2],
            "prefix": key[3],
            "center": key[4],
            "earlier": key[5],
            "child_mass": F(data["h"]),
            "child_components": int(data["r"]),
            "ledger_horizon": int(data["M"]),
            "route": data["route"],
        }
    require(len(rows) == 4_866, "pair universe changed")
    return rows


def parse_failed_grandchildren(
    path: Path,
    first_rows: dict[tuple[object, ...], dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    dict[tuple[object, ...], set[tuple[object, ...]]],
]:
    failed: list[dict[str, object]] = []
    second_ids_by_first: dict[
        tuple[object, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    recursive_firsts: set[tuple[object, ...]] = set()
    for line in path.read_text().splitlines():
        if line.startswith("RECURSIVE;"):
            data = fields(line)
            key = first_identity(data)
            require(key in first_rows, "recursive first identity missing")
            recursive_firsts.add(key)
        elif line.startswith("SECOND;"):
            data = fields(line)
            first = first_identity(data)
            require(first in first_rows, "second row first identity missing")
            second = (
                *first,
                int(data["y"]),
                ints(data["yearlier"]),
            )
            require(second not in second_ids_by_first[first], "second identity collision")
            second_ids_by_first[first].add(second)
            margin = F(data["margin"])
            derived_closed = margin > 0
            require(
                (data["closed"] == "1") == derived_closed,
                f"{second}: second closed bit disagrees with exact margin",
            )
            if derived_closed:
                continue
            failed.append(
                {
                    **first_rows[first],
                    "second_center": second[-2],
                    "second_earlier": second[-1],
                    "grandchild_mass": F(data["h"]),
                    "grandchild_components": int(data["r"]),
                    "ledger_top3": G.parse_ranked(data["top3"]),
                    "ledger_margin": margin,
                    "second_identity": second,
                }
            )
    require(
        len(recursive_firsts) == 2_079
        and len(second_ids_by_first) == 2_079
        and sum(map(len, second_ids_by_first.values())) == 6_172
        and len(failed) == 228,
        "recursive/failed-grandchild universe changed",
    )
    failed.sort(key=lambda row: row["second_identity"])
    return failed, second_ids_by_first


def reconstruct_grandchild(
    row: dict[str, object],
) -> tuple[list[tuple[F, F]], frozenset[int]]:
    child, first_forbidden = G.reconstruct_child(row)
    grandchild = H.R.subtract_local(child, row["second_center"])
    grandchild_family = tuple(
        (
            *row["body"],
            row["apex"],
            row["center"],
            row["second_center"],
        )
    )
    require(
        len(grandchild_family) == 10,
        "grandchild arity is not 10 fixed + 3 remaining",
    )
    direct, components, mass = G.direct_carrier(grandchild_family)
    require(
        grandchild == direct
        and len(grandchild) == components == row["grandchild_components"]
        and G.interval_mass(grandchild) == mass == row["grandchild_mass"],
        f"{row['second_identity']}: grandchild reconstruction failed",
    )
    forbidden = frozenset(
        (
            *first_forbidden,
            row["second_center"],
            *row["second_earlier"],
        )
    )
    require(
        len(forbidden)
        == len(first_forbidden) + 1 + len(row["second_earlier"]),
        f"{row['second_identity']}: second forbidden prefix collided",
    )
    return grandchild, forbidden


def analyze_grandchild(row: dict[str, object]) -> dict[str, object]:
    grandchild, forbidden = reconstruct_grandchild(row)
    pair = G.exact_pair_cap(
        grandchild,
        forbidden,
        initial_horizon=INITIAL_HORIZON,
    )
    pair_residual = G.subtract_multi(grandchild, pair["witness"])
    pair_family = tuple(
        (
            *row["body"],
            row["apex"],
            row["center"],
            row["second_center"],
            *pair["witness"],
        )
    )
    require(
        len(pair_family) == 12,
        "grandchild B2 did not consume exactly two of three slots",
    )
    pair_direct, pair_components, pair_mass = G.direct_carrier(pair_family)
    require(
        pair_residual == pair_direct
        and len(pair_residual) == pair_components
        and G.interval_mass(pair_residual) == pair_mass
        and row["grandchild_mass"] - pair_mass == pair["head"],
        f"{row['second_identity']}: grandchild pair winner failed full-family reconstruction",
    )
    require(
        tuple(value for value, _ in pair["top4"][:3])
        == tuple(value for value, _ in row["ledger_top3"]),
        f"{row['second_identity']}: global top3 values changed",
    )
    qs = tuple(value for value, _ in pair["top4"][:3])
    require(
        row["grandchild_mass"] - sum(qs, F(0)) == row["ledger_margin"],
        f"{row['second_identity']}: top3 margin changed",
    )
    partition_margin = row["grandchild_mass"] - (
        pair["head"] + qs[0]
    )
    hunter, threshold = T.hunter_data(qs, pair["head"], row["grandchild_mass"])
    hunter_margin = row["grandchild_mass"] - hunter
    route = (
        "pair_single"
        if partition_margin > 0
        else "hunter3"
        if hunter_margin > 0
        else "third_recursive"
    )
    return {
        **row,
        "pair": pair,
        "partition_margin": partition_margin,
        "hunter": hunter,
        "hunter_margin": hunter_margin,
        "threshold": threshold,
        "route3": route,
    }


def third_pivot(
    row: dict[str, object],
    grandchild: list[tuple[F, F]],
    forbidden: frozenset[int],
    core: tuple[tuple[F, int], ...],
    index: int,
) -> dict[str, object]:
    coverage, center = core[index]
    earlier = tuple(label for _, label in core[:index])
    third_forbidden = frozenset((*forbidden, center, *earlier))
    require(
        center not in forbidden
        and not (set(earlier) & forbidden)
        and center not in earlier
        and len(third_forbidden) == len(forbidden) + 1 + len(earlier),
        f"{row['second_identity']};z={center}: third-centre forbidden sidecar collided",
    )
    greatgrandchild = H.R.subtract_local(grandchild, center)
    greatgrandchild_family = tuple(
        (
            *row["body"],
            row["apex"],
            row["center"],
            row["second_center"],
            center,
        )
    )
    require(
        len(greatgrandchild_family) == 11,
        "third pivot is not 11 fixed + 2 remaining",
    )
    direct, components, mass = G.direct_carrier(greatgrandchild_family)
    require(
        greatgrandchild == direct
        and len(greatgrandchild) == components
        and G.interval_mass(greatgrandchild) == mass,
        f"{row['second_identity']};z={center}: great-grandchild reconstruction failed",
    )
    pair = G.exact_pair_cap(
        greatgrandchild,
        third_forbidden,
        initial_horizon=INITIAL_HORIZON,
    )
    pair_residual = G.subtract_multi(greatgrandchild, pair["witness"])
    terminal_family = tuple(
        (
            *row["body"],
            row["apex"],
            row["center"],
            row["second_center"],
            center,
            *pair["witness"],
        )
    )
    require(
        len(terminal_family) == 13,
        "terminal B2 did not complete exactly 13 speeds",
    )
    pair_direct, pair_components, pair_mass = G.direct_carrier(terminal_family)
    require(
        pair_residual == pair_direct
        and len(pair_residual) == pair_components
        and G.interval_mass(pair_residual) == pair_mass
        and mass - pair_mass == pair["head"],
        f"{row['second_identity']};z={center}: terminal pair winner failed full-family reconstruction",
    )
    margin = mass - pair["head"]
    return {
        "center": center,
        "coverage": coverage,
        "earlier": earlier,
        "mass": mass,
        "components": components,
        "pair": pair,
        "margin": margin,
        "closed": margin > 0,
    }


def analyze_third_recursive(row: dict[str, object]) -> dict[str, object]:
    require(row["route3"] == "third_recursive", "wrong route in third recursion")
    grandchild, forbidden = reconstruct_grandchild(row)
    mass = row["grandchild_mass"]
    threshold = row["threshold"]
    require(threshold is not None, "Hunter3 failure has no threshold")
    delta = threshold - mass / 7
    require(delta > 0, "third-centre threshold is not discrepancy-finite")
    gamma = H.S2 * len(grandchild) / 7
    cutoff = G.ceiling(gamma / delta) - 1
    require(cutoff >= FIRST_EXTERNAL, "third-centre core is empty by cutoff")
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in forbidden
    ]
    coverages = H.exact_coverages(grandchild, labels)
    core = tuple(
        sorted(
            (
                (value, label)
                for value, label in coverages
                if value >= threshold
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    require(core, "third-centre hostile core is empty")
    require(
        mass / 7 + gamma / (cutoff + 1) <= threshold,
        "third-centre tail failed to seal",
    )
    pivots = tuple(
        third_pivot(row, grandchild, forbidden, core, index)
        for index in range(len(core))
    )
    return {
        **row,
        "delta3": delta,
        "cutoff3": cutoff,
        "core3": core,
        "pivots3": pivots,
        "third_closed": all(pivot["closed"] for pivot in pivots),
    }


def first_key_text(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"x={row['center']};earlier={','.join(map(str, row['earlier']))}"
    )


def second_key_text(row: dict[str, object]) -> str:
    return (
        first_key_text(row)
        + f";y={row['second_center']};"
        f"yearlier={','.join(map(str, row['second_earlier']))}"
    )


def grandchild_line(row: dict[str, object]) -> str:
    pair = row["pair"]
    return (
        "GRANDCHILD;"
        + second_key_text(row)
        + f";h={G.ftext(row['grandchild_mass'])};r={row['grandchild_components']};"
        + f"M={pair['horizon']};B2={G.ftext(pair['head'])};"
        + f"B2w={pair['witness'][0]},{pair['witness'][1]};"
        + f"B2tail={G.ftext(pair['tail'])};B2gap={G.ftext(pair['tail_gap'])};"
        + f"paid={pair['paid']};"
        + f"partitionmargin={G.ftext(row['partition_margin'])};"
        + f"G3={G.ftext(row['hunter'])};G3margin={G.ftext(row['hunter_margin'])};"
        + f"lambda3={None if row['threshold'] is None else G.ftext(row['threshold'])};"
        + f"route3={row['route3']}"
    )


def third_lines(row: dict[str, object]) -> list[str]:
    lines = [
        (
            "THIRD_RECURSIVE;"
            + second_key_text(row)
            + f";lambda3={G.ftext(row['threshold'])};"
            + f"delta3={G.ftext(row['delta3'])};N3={row['cutoff3']};"
            + "H1="
            + ",".join(
                f"{label}:{G.ftext(value)}" for value, label in row["core3"]
            )
            + f";closed={int(row['third_closed'])}"
        )
    ]
    for pivot in row["pivots3"]:
        pair = pivot["pair"]
        lines.append(
            "THIRD;"
            + second_key_text(row)
            + f";z={pivot['center']};"
            + f"zearlier={','.join(map(str, pivot['earlier']))};"
            + f"h={G.ftext(pivot['mass'])};r={pivot['components']};"
            + f"M={pair['horizon']};B2={G.ftext(pair['head'])};"
            + f"B2w={pair['witness'][0]},{pair['witness'][1]};"
            + f"B2tail={G.ftext(pair['tail'])};B2gap={G.ftext(pair['tail_gap'])};"
            + f"paid={pair['paid']};margin={G.ftext(pivot['margin'])};"
            + f"closed={int(pivot['closed'])}"
        )
    return lines


def semantic_digest(domain: bytes, lines: list[str]) -> str:
    digest = hashlib.sha256(domain + b"\n")
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-ledger", type=Path, required=True)
    parser.add_argument("--recursive-ledger", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--grandchild-ledger", type=Path, required=True)
    parser.add_argument("--third-ledger", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")

    first_rows = parse_pair_rows(args.pair_ledger)
    failed, second_ids_by_first = parse_failed_grandchildren(
        args.recursive_ledger, first_rows
    )
    if args.workers == 1:
        grandchildren = list(map(analyze_grandchild, failed))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            grandchildren = pool.map(analyze_grandchild, failed, chunksize=1)
    grandchildren.sort(key=lambda row: row["second_identity"])

    grandchild_lines = [grandchild_line(row) for row in grandchildren]
    grandchild_digest = semantic_digest(
        b"LRC14/THM2915-failures/third-level-pair-Hunter/v1",
        grandchild_lines,
    )
    args.grandchild_ledger.write_text(
        "LRC14 THM2915 failed-grandchild pair/Hunter ledger\n"
        + "\n".join(grandchild_lines)
        + "\n"
        + f"semantic_sha256={grandchild_digest}\n"
        + "scope=228 failed second pivots;seven-body/six-external rung;"
        + "not unrestricted LRC14\n",
        encoding="utf-8",
        newline="\n",
    )

    third_inputs = [
        row for row in grandchildren if row["route3"] == "third_recursive"
    ]
    if args.workers == 1:
        third_rows = list(map(analyze_third_recursive, third_inputs))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            third_rows = pool.map(analyze_third_recursive, third_inputs, chunksize=1)
    third_rows.sort(key=lambda row: row["second_identity"])
    third_lines_all = [
        line for row in third_rows for line in third_lines(row)
    ]
    third_digest = semantic_digest(
        b"LRC14/THM2915-failures/ordered-third-centre/v1",
        third_lines_all,
    )
    args.third_ledger.write_text(
        "LRC14 THM2915 ordered-third-centre ledger\n"
        + "\n".join(third_lines_all)
        + "\n"
        + f"semantic_sha256={third_digest}\n"
        + f"scope={len(third_rows)} failed-grandchild recursive rows;"
        + "seven-body/six-external rung;not unrestricted LRC14\n",
        encoding="utf-8",
        newline="\n",
    )

    third_by_id = {
        row["second_identity"]: row for row in third_rows
    }
    upgraded_second = {
        row["second_identity"]
        for row in grandchildren
        if row["route3"] in {"pair_single", "hunter3"}
        or (
            row["second_identity"] in third_by_id
            and third_by_id[row["second_identity"]]["third_closed"]
        )
    }
    failed_second_ids = {row["second_identity"] for row in failed}
    remaining_second = failed_second_ids - upgraded_second
    upgraded_first = {
        first
        for first, seconds in second_ids_by_first.items()
        if not (seconds & remaining_second)
    }
    initially_open_first = {
        first for first, seconds in second_ids_by_first.items()
        if seconds & failed_second_ids
    }
    require(len(initially_open_first) == 195, "open first-child count changed")
    newly_closed_first = initially_open_first & upgraded_first
    remaining_first = initially_open_first - newly_closed_first

    partition_positive = [
        row for row in grandchildren if row["partition_margin"] > 0
    ]
    partition_failed = [
        row for row in grandchildren if row["partition_margin"] <= 0
    ]
    hunter_positive = [
        row for row in partition_failed if row["hunter_margin"] > 0
    ]
    hunter_failed = [
        row for row in partition_failed if row["hunter_margin"] <= 0
    ]
    third_pivots = [
        (row, pivot)
        for row in third_rows
        for pivot in row["pivots3"]
    ]
    require(
        grandchild_digest == EXPECTED_GRANDCHILD_DIGEST,
        "grandchild semantic digest changed",
    )
    require(third_digest == EXPECTED_THIRD_DIGEST, "third semantic digest changed")
    require(
        (
            len(partition_positive),
            len(partition_failed),
            len(hunter_positive),
            len(hunter_failed),
            len(third_rows),
            sum(row["third_closed"] for row in third_rows),
            sum(not row["third_closed"] for row in third_rows),
            len(third_pivots),
            sum(pivot["closed"] for _, pivot in third_pivots),
            sum(not pivot["closed"] for _, pivot in third_pivots),
            len(upgraded_second),
            len(remaining_second),
            len(newly_closed_first),
            len(remaining_first),
            sum(row["pair"]["paid"] for row in grandchildren),
            sum(pivot["pair"]["paid"] for _, pivot in third_pivots),
        )
        == (
            65,
            163,
            114,
            49,
            49,
            48,
            1,
            117,
            116,
            1,
            227,
            1,
            194,
            1,
            2_153,
            1_722,
        ),
        "third-level census changed",
    )
    summary: dict[str, object] = {
        "failed_grandchildren": len(grandchildren),
        "pair_single_positive": len(partition_positive),
        "pair_single_nonpositive": len(partition_failed),
        "pair_single_equal": sum(
            row["partition_margin"] == 0 for row in grandchildren
        ),
        "hunter3_positive": len(hunter_positive),
        "hunter3_nonpositive": len(hunter_failed),
        "hunter3_equal": sum(
            row["hunter_margin"] == 0 for row in partition_failed
        ),
        "third_recursive_rows": len(third_rows),
        "third_recursive_closed": sum(row["third_closed"] for row in third_rows),
        "third_recursive_open": sum(not row["third_closed"] for row in third_rows),
        "third_pivots": len(third_pivots),
        "third_pivots_closed": sum(pivot["closed"] for _, pivot in third_pivots),
        "third_pivots_open": sum(not pivot["closed"] for _, pivot in third_pivots),
        "maximum_third_core": max((len(row["core3"]) for row in third_rows), default=0),
        "maximum_third_cutoff": max((row["cutoff3"] for row in third_rows), default=0),
        "maximum_pair_horizon": max(
            [
                row["pair"]["horizon"] for row in grandchildren
            ]
            + [
                pivot["pair"]["horizon"] for _, pivot in third_pivots
            ]
        ),
        "minimum_pair_tail_gap": min(
            [
                row["pair"]["tail_gap"] for row in grandchildren
            ]
            + [
                pivot["pair"]["tail_gap"] for _, pivot in third_pivots
            ]
        ),
        "upgraded_second_pivots": len(upgraded_second),
        "remaining_second_pivots": len(remaining_second),
        "newly_closed_first_children": len(newly_closed_first),
        "remaining_first_children": len(remaining_first),
        "slot_chain": ((10, 3), (11, 2), (13, 0)),
        "route3_histogram": tuple(
            sorted(Counter(row["route3"] for row in grandchildren).items())
        ),
        "grandchild_semantic_sha256": grandchild_digest,
        "third_semantic_sha256": third_digest,
    }
    if partition_positive:
        row = min(
            partition_positive,
            key=lambda row: (row["partition_margin"], row["second_identity"]),
        )
        summary["closest_pair_single_positive"] = (
            second_key_text(row)
            + f";margin={ftext(row['partition_margin'])}"
        )
    if partition_failed:
        row = max(
            partition_failed,
            key=lambda row: (row["partition_margin"], row["second_identity"]),
        )
        summary["closest_pair_single_failure"] = (
            second_key_text(row)
            + f";margin={ftext(row['partition_margin'])}"
        )
    if hunter_positive:
        row = min(
            hunter_positive,
            key=lambda row: (row["hunter_margin"], row["second_identity"]),
        )
        summary["closest_hunter3_positive"] = (
            second_key_text(row)
            + f";margin={ftext(row['hunter_margin'])}"
        )
    if hunter_failed:
        row = max(
            hunter_failed,
            key=lambda row: (row["hunter_margin"], row["second_identity"]),
        )
        summary["closest_hunter3_failure"] = (
            second_key_text(row)
            + f";margin={ftext(row['hunter_margin'])}"
        )
    positive_third = [
        pair for pair in third_pivots if pair[1]["margin"] > 0
    ]
    failed_third = [
        pair for pair in third_pivots if pair[1]["margin"] <= 0
    ]
    if positive_third:
        row, pivot = min(
            positive_third,
            key=lambda pair: (
                pair[1]["margin"],
                pair[0]["second_identity"],
                pair[1]["center"],
            ),
        )
        summary["closest_third_positive"] = (
            second_key_text(row)
            + f";z={pivot['center']};margin={ftext(pivot['margin'])}"
        )
    if failed_third:
        row, pivot = max(
            failed_third,
            key=lambda pair: (
                pair[1]["margin"],
                pair[0]["second_identity"],
                pair[1]["center"],
            ),
        )
        summary["closest_third_failure"] = (
            second_key_text(row)
            + f";z={pivot['center']};margin={ftext(pivot['margin'])}"
        )

    rendered = "\n".join(f"{key}={value}" for key, value in summary.items()) + "\n"
    if args.summary is not None:
        args.summary.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
