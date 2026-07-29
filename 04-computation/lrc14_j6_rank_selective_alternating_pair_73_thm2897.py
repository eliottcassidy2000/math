#!/usr/bin/env python3
"""Locked four-root census of THM-2897's q5 + 2 B2 < h certificate."""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.py"
)
SUFFIX_SHA = "6434f020c5aa4000ac81fa081881d93ac0b4190516f854fbd9d8493475baf539"
PAIR_PATH = (
    ROOT
    / "04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py"
)
PAIR_SHA = "b82f318bf89ffd3ab4c918c87736461d068e03f25941aa25a0961d0f74b4d70a"
EXPECTED_COUNTS: tuple[int, ...] | None = (73, 48, 50, 51, 3, 1, 22, 1_388)
EXPECTED_ROOT_COUNTS: tuple[tuple[object, ...], ...] | None = (
    ((2, 8, 9, 10, 11, 13, 14), 19, 13, 13, 14, 1),
    ((1, 3, 9, 10, 11, 12, 14), 20, 15, 17, 17, 2),
    ((2, 5, 9, 11, 12, 13, 14), 21, 13, 13, 13, 0),
    ((2, 3, 4, 5, 6, 7, 8), 13, 7, 7, 7, 0),
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "97e4f5c0b8194f98b568245e88517e2c239875fb0512bb1cb8b40e0e6f90161b"
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


S = load(SUFFIX_PATH, SUFFIX_SHA, "j6_alt_pair_suffix")
H = load(PAIR_PATH, PAIR_SHA, "j6_alt_pair_kernel")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def main() -> None:
    rows: list[dict[str, object]] = []
    ledger = hashlib.sha256(b"LRC14/j6/alternating-pair-73/v1\n")
    for body in S.BATTERY_ROOTS:
        root = S.global_root(body)
        for rank in range(1, root["K"] + 1):
            branch = S.suffix_branch(root, rank)
            carrier = H.R.subtract_local(root["good"], branch["apex"])
            require(
                mass(carrier) == branch["m"]
                and len(carrier) == branch["r"],
                f"carrier mismatch: {body}, rank {rank}",
            )
            excluded = set(branch["excluded_prefix"])
            pair = H.pair_cap(carrier, excluded)
            q1 = pair["q1"]
            q5 = branch["top5"][4][0]
            require(
                pair["ranked"][4][0] == q5,
                f"global q5 seal mismatch: {body}, rank {rank}",
            )
            b2 = pair["cap"]
            margin = branch["m"] - q5 - 2 * b2
            row = {
                "body": body,
                "rank": rank,
                "apex": branch["apex"],
                "prefix": branch["excluded_prefix"],
                "h": branch["m"],
                "scalar_closed": branch["closed"],
                "scalar_margin": branch["margin"],
                "q1": q1,
                "q5": q5,
                "b2": b2,
                "head": pair["head"],
                "tail": pair["tail"],
                "maximizer": pair["maximizer"],
                "paid": pair["paid"],
                "margin": margin,
                "closed": margin > 0,
            }
            rows.append(row)
            ledger.update(
                (
                    f"E={body};rank={rank};apex={branch['apex']};"
                    f"P={branch['excluded_prefix']};h={ftext(branch['m'])};"
                    f"scalar={int(branch['closed'])}:{ftext(branch['margin'])};"
                    f"q1={ftext(q1)};q5={ftext(q5)};B2={ftext(b2)};"
                    f"head={ftext(pair['head'])}:{pair['maximizer']};"
                    f"tail={ftext(pair['tail'])};paid={pair['paid']};"
                    f"paid_digest={pair['paid_digest']};"
                    f"margin={ftext(margin)};closed={int(margin > 0)}\n"
                ).encode()
            )

    require(len(rows) == 73, "suffix universe changed")
    root_counts = tuple(
        (
            body,
            sum(row["body"] == body for row in rows),
            sum(
                row["body"] == body and row["scalar_closed"]
                for row in rows
            ),
            sum(row["body"] == body and row["closed"] for row in rows),
            sum(
                row["body"] == body
                and (row["scalar_closed"] or row["closed"])
                for row in rows
            ),
            sum(
                row["body"] == body
                and not row["scalar_closed"]
                and row["closed"]
                for row in rows
            ),
        )
        for body in S.BATTERY_ROOTS
    )
    counts = (
        len(rows),
        sum(row["scalar_closed"] for row in rows),
        sum(row["closed"] for row in rows),
        sum(row["scalar_closed"] or row["closed"] for row in rows),
        sum(not row["scalar_closed"] and row["closed"] for row in rows),
        sum(row["scalar_closed"] and not row["closed"] for row in rows),
        sum(
            not row["scalar_closed"] and not row["closed"]
            for row in rows
        ),
        sum(row["paid"] for row in rows),
    )
    failures = sorted(
        (row for row in rows if not row["closed"]),
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    newly_closed = sorted(
        (
            row
            for row in rows
            if not row["scalar_closed"] and row["closed"]
        ),
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    min_row = min(rows, key=lambda row: row["margin"])
    max_paid = max(rows, key=lambda row: row["paid"])
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, f"aggregate counts changed: {counts}")
    if EXPECTED_ROOT_COUNTS is not None:
        require(
            root_counts == EXPECTED_ROOT_COUNTS,
            f"root counts changed: {root_counts}",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger.hexdigest() == EXPECTED_LEDGER_DIGEST,
            "canonical ledger changed",
        )

    print("LRC14 j6 THM-2897 rank-selective 73-branch verifier")
    print(f"counts={counts}")
    print(f"root_counts={root_counts}")
    print(
        "minimum_margin="
        f"{ftext(min_row['margin'])};E={min_row['body']};"
        f"rank={min_row['rank']};apex={min_row['apex']};"
        f"q1={ftext(min_row['q1'])};q5={ftext(min_row['q5'])};"
        f"B2={ftext(min_row['b2'])};"
        f"maximizer={min_row['maximizer']}"
    )
    print(
        "maximum_paid="
        f"{max_paid['paid']};E={max_paid['body']};"
        f"rank={max_paid['rank']};apex={max_paid['apex']}"
    )
    print(f"newly_closed_count={len(newly_closed)}")
    for row in newly_closed:
        print(
            f"NEW E={row['body']};rank={row['rank']};"
            f"apex={row['apex']};margin={ftext(row['margin'])};"
            f"q1={ftext(row['q1'])};q5={ftext(row['q5'])};"
            f"B2={ftext(row['b2'])};"
            f"maximizer={row['maximizer']}"
        )
    print(f"failure_count={len(failures)}")
    for row in failures:
        print(
            f"FAIL E={row['body']};rank={row['rank']};"
            f"apex={row['apex']};prefix={row['prefix']};"
            f"scalar_closed={int(row['scalar_closed'])};"
            f"margin={ftext(row['margin'])};"
            f"q1={ftext(row['q1'])};q5={ftext(row['q5'])};"
            f"B2={ftext(row['b2'])};"
            f"head={ftext(row['head'])};tail={ftext(row['tail'])};"
            f"maximizer={row['maximizer']};paid={row['paid']}"
        )
    print(f"canonical_ledger_sha256={ledger.hexdigest()}")
    print("scope=four roots;73 actual marked suffixes;not uniform LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
