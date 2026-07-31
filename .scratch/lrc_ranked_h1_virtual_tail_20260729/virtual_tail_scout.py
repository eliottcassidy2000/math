#!/usr/bin/env python3
"""Exact virtual-tail certificate for deferred ranked-r4/H1 branches.

For a literal carrier C of mass h with r interval components, THM-735 gives

    c_C(w) < h/7 + (99/70) r / (7w).

After scanning labels through M, every unscanned label has coverage strictly
below

    U_M = h/7 + (99/70) r / (7(M+1)).

This script uses four virtual children of weight U_M in the ordered
literal-child certificate.  Thus it reasons about the complete infinite tail
without enumerating labels to the discrepancy cutoff.

Discovery-only scratch code.  No theorem is claimed by this file.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py"
)
ENGINE_SHA = "e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9"
BASE_CUTOFF = 15_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_engine():
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA,
        "ranked-suffix engine changed",
    )
    spec = importlib.util.spec_from_file_location("virtual_tail_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_engine()


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def cutoff_for(row: dict[str, object]) -> tuple[F, F, int]:
    r4 = sum((value for value, _ in row["top5"][:4]), F(0))
    epsilon = row["m"] - r4 - row["m"] / 7
    require(epsilon > 0, "cutoff requested outside ranked-r4 gate")
    cutoff = S.ceiling(S.S2 * row["r"] / (7 * epsilon)) - 1
    return r4, epsilon, cutoff


def tail_certificate(
    carrier: list[tuple[F, F]],
    h: F,
    components: int,
    prefix: tuple[int, ...],
    level: F,
    head_max: int,
) -> dict[str, object]:
    excluded = set(prefix)
    labels = [
        w
        for w in range(S.FIRST_EXTERNAL, head_max + 1)
        if w not in excluded
    ]
    coverage_rows = S.G.T.coverages_many(carrier, labels)
    head = tuple(
        sorted(
            (
                (value, speed)
                for value, speed in coverage_rows
                if value >= level
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    tail_upper = h / 7 + S.S2 * components / (7 * (head_max + 1))
    all_tail_margin = h - 5 * tail_upper
    if all_tail_margin < 0:
        return {
            "closed": False,
            "reason": "ALL_TAIL",
            "head": head,
            "tail_upper": tail_upper,
            "all_tail_margin": all_tail_margin,
            "parent": 0,
            "local": 0,
            "open": (),
            "min_margin": None,
        }

    parent_count = 0
    local_count = 0
    opens: list[tuple[int, int, F, F, F]] = []
    margins: list[F] = [all_tail_margin]
    for index, (_, speed) in enumerate(head):
        residual = S.R.subtract_local(carrier, speed)
        residual_mass = mass(residual)
        if not residual:
            opens.append((index, speed, residual_mass, F(0), F(0)))
            continue

        later = head[index + 1 :]
        inherited_values = [value for value, _ in later]
        inherited_values.extend([tail_upper] * 4)
        inherited_values.sort(reverse=True)
        inherited_upper = sum(inherited_values[:4], F(0))
        inherited_margin = residual_mass - inherited_upper
        if inherited_margin > 0:
            parent_count += 1
            margins.append(inherited_margin)
            continue

        child_rows = S.G.T.coverages_many(
            residual,
            [other for _, other in later],
        )
        local_values = [value for value, _ in child_rows]
        local_values.extend([tail_upper] * 4)
        local_values.sort(reverse=True)
        local_upper = sum(local_values[:4], F(0))
        local_margin = residual_mass - local_upper
        if local_margin > 0:
            local_count += 1
            margins.append(local_margin)
        else:
            opens.append(
                (
                    index,
                    speed,
                    residual_mass,
                    inherited_upper,
                    local_upper,
                )
            )

    return {
        "closed": not opens,
        "reason": "CLOSED" if not opens else "PIVOT",
        "head": head,
        "tail_upper": tail_upper,
        "all_tail_margin": all_tail_margin,
        "parent": parent_count,
        "local": local_count,
        "open": tuple(opens),
        "min_margin": min(margins),
    }


def profile_body(task: tuple[tuple[int, ...], int]) -> dict[str, object]:
    body, head_max = task
    root = S.G.profile_body(body)
    good, components, h0 = S.G.T.CORE.good_norm(body)
    require(
        components == root["r"] and h0 == root["m"],
        f"root reconstruction changed: {body}",
    )
    root["good"] = good
    hard = [
        row
        for rank in range(1, root["adaptive_k"] + 1)
        for row in (S.profile_branch(root, rank),)
        if not row["closed"]
    ]
    records = []
    for row in hard:
        margins = S.ranked_core_margins(row)
        if margins[3] <= 0:
            continue
        r4, epsilon, cutoff = cutoff_for(row)
        if cutoff <= BASE_CUTOFF:
            continue
        carrier = S.R.subtract_local(good, row["apex"])
        direct, direct_r, direct_h = S.G.T.CORE.good_norm(
            tuple(sorted((*body, row["apex"])))
        )
        require(
            carrier == direct
            and len(carrier) == direct_r == row["r"]
            and mass(carrier) == direct_h == row["m"],
            f"literal carrier changed: {body}, rank {row['rank']}",
        )
        level = row["m"] - r4
        cert = tail_certificate(
            carrier,
            row["m"],
            row["r"],
            row["prefix"],
            level,
            head_max,
        )
        records.append(
            {
                "body": body,
                "K": root["adaptive_k"],
                "rank": row["rank"],
                "apex": row["apex"],
                "prefix": row["prefix"],
                "h": row["m"],
                "r": row["r"],
                "r4": r4,
                "epsilon": epsilon,
                "cutoff": cutoff,
                **cert,
            }
        )
    return {"body": body, "records": tuple(records)}


def record_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};"
        f"a={row['apex']};P={row['prefix']};h={ftext(row['h'])};"
        f"r={row['r']};R4={ftext(row['r4'])};"
        f"epsilon={ftext(row['epsilon'])};N={row['cutoff']};"
        f"M={row['head_max']};Hhead={len(row['head'])};"
        f"U={ftext(row['tail_upper'])};"
        f"tailmargin={ftext(row['all_tail_margin'])};"
        f"parent={row['parent']};local={row['local']};"
        f"minmargin={ftext(row['min_margin'])};"
        f"closed={row['closed']};reason={row['reason']};"
        f"open={row['open']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--head-max", type=int, default=BASE_CUTOFF)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 2),
    )
    args = parser.parse_args()
    require(args.head_max >= S.FIRST_EXTERNAL, "head maximum too small")
    require(args.workers >= 1, "workers must be positive")

    tasks = tuple((body, args.head_max) for body in S.G.BODIES)
    if args.workers == 1:
        roots = [profile_body(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            roots = list(pool.imap(profile_body, tasks, chunksize=1))
    records = [row for root in roots for row in root["records"]]
    for row in records:
        row["head_max"] = args.head_max
    require(len(records) == 181, "deferred-row universe changed")

    status = Counter(row["reason"] for row in records)
    closed_roots = tuple(
        root["body"]
        for root in roots
        if root["records"]
        and all(row["closed"] for row in root["records"])
    )
    digest = hashlib.sha256(b"LRC14/j6/ranked-H1-virtual-tail/v1\n")
    for row in records:
        digest.update(record_line(row).encode())

    print("LRC14 j6 ranked-r4/H1 virtual-tail scout")
    print(
        f"parameters=head_max:{args.head_max},workers:{args.workers},"
        f"base_cutoff:{BASE_CUTOFF}"
    )
    print(
        "counts="
        + repr(
            (
                len(S.G.BODIES),
                len(records),
                sum(row["closed"] for row in records),
                len(closed_roots),
                sum(row["parent"] for row in records),
                sum(row["local"] for row in records),
            )
        )
    )
    print(f"statuses={tuple(sorted(status.items()))}")
    print(
        "head_sizes="
        + repr(
            (
                min(len(row["head"]) for row in records),
                max(len(row["head"]) for row in records),
            )
        )
    )
    print(
        "cutoffs="
        + repr(
            (
                min(row["cutoff"] for row in records),
                max(row["cutoff"] for row in records),
            )
        )
    )
    print(
        "minimum_tail_margin="
        + ftext(min(row["all_tail_margin"] for row in records))
    )
    print(
        "minimum_certificate_margin="
        + ftext(
            min(
                row["min_margin"]
                for row in records
                if row["min_margin"] is not None
            )
        )
    )
    print(f"closed_roots={closed_roots}")
    print("OPEN")
    for row in records:
        if not row["closed"]:
            print(record_line(row).rstrip())
    print("LARGEST_HEADS")
    for row in sorted(
        records,
        key=lambda item: (
            len(item["head"]),
            item["cutoff"],
            item["body"],
            item["rank"],
        ),
        reverse=True,
    )[:20]:
        print(record_line(row).rstrip())
    print(f"ledger_sha256={digest.hexdigest()}")
    print("scope=181 cutoff-deferred ranked-r4/H1 rows;discovery-only")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
