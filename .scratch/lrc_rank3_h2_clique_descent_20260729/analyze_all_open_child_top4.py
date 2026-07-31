#!/usr/bin/env python3
"""Recompose the full all-open-centre child-top-four discovery ledger."""

from __future__ import annotations

import ast
import hashlib
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LEDGER = (
    ROOT
    / ".scratch/lrc_rank3_h2_clique_descent_20260729/"
    "all_open_child_top4_ALL.dynamic.ledger.out"
)
THM2912_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out"
)
THM2912_UNION_DIGEST = (
    "7f7356072167bd09c7eee35560c3a5076782c1f8704d06b4a86b64e273f4b32a"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def body_digest(bodies: set[tuple[int, ...]]) -> str:
    payload = "".join(
        ",".join(map(str, body)) + "\n"
        for body in sorted(bodies)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def output_literal(path: Path, prefix: str):
    rows = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(rows) == 1, f"expected one {prefix!r} line")
    return ast.literal_eval(rows[0])


def parse_children() -> list[dict[str, object]]:
    rows = []
    for line in LEDGER.read_text().splitlines():
        if not line.startswith("CHILD;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        rows.append(
            {
                "body": tuple(map(int, fields["E"].split(","))),
                "rank": int(fields["rank"]),
                "apex": int(fields["a"]),
                "center": int(fields["x"]),
                "margin": F(fields["margin"]),
                "tail_gap": F(fields["tailgap"]),
                "head_limit": int(fields["M"]),
                "required_limit": (
                    None if fields["N"] == "None" else int(fields["N"])
                ),
                "sealed": bool(int(fields["sealed"])),
                "closed": bool(int(fields["closed"])),
            }
        )
    require(len(rows) == 51_222, "full child count changed")
    return rows


def main() -> None:
    children = parse_children()
    by_parent = defaultdict(list)
    by_body = defaultdict(list)
    for row in children:
        by_parent[(row["body"], row["rank"], row["apex"])].append(row)
        by_body[row["body"]].append(row)
    require(len(by_parent) == 11_563, "open parent-row count changed")

    closed_parents = {
        parent
        for parent, rows in by_parent.items()
        if all(row["closed"] for row in rows)
    }
    route_roots = {
        body
        for body, rows in by_body.items()
        if all(row["closed"] for row in rows)
    }
    failed_roots = set(by_body) - route_roots
    baseline = set(output_literal(THM2912_OUT, "proved_union_roots="))
    require(
        len(baseline) == 314
        and body_digest(baseline) == THM2912_UNION_DIGEST,
        "THM2912 proved union changed",
    )
    overlap = route_roots & baseline
    additive = route_roots - baseline
    new_union = route_roots | baseline
    failed_children = [row for row in children if not row["closed"]]

    print("LRC14 all-open-centre child-top-four recomposition")
    print(
        f"children={len(children)};"
        f"closed={sum(row['closed'] for row in children)};"
        f"failed={len(failed_children)};"
        f"sealed={sum(row['sealed'] for row in children)};"
        f"unsealed={sum(not row['sealed'] for row in children)}"
    )
    print(
        f"open_parent_rows={len(by_parent)};"
        f"closed_parent_rows={len(closed_parents)};"
        f"failed_parent_rows={len(by_parent)-len(closed_parents)}"
    )
    print(
        f"route_bodies={len(by_body)};route_closed_roots={len(route_roots)};"
        f"route_failed_roots={len(failed_roots)}"
    )
    print(
        f"baseline314={len(baseline)};overlap={len(overlap)};"
        f"additive={len(additive)};new_union={len(new_union)};"
        f"official_residual={3432-len(new_union)}"
    )
    print(
        f"extended={sum(row['head_limit']>2000 for row in children)};"
        f"max_head={max(row['head_limit'] for row in children)};"
        f"max_required={max((row['required_limit'] or 0) for row in children)};"
        f"minimum_tail_gap={min(row['tail_gap'] for row in children)}"
    )
    print(
        "failed_children_per_root_histogram="
        + repr(
            tuple(
                sorted(
                    Counter(
                        sum(not row["closed"] for row in rows)
                        for rows in by_body.values()
                        if any(not row["closed"] for row in rows)
                    ).items()
                )
            )
        )
    )
    print(f"route_root_digest={body_digest(route_roots)}")
    print(f"overlap_digest={body_digest(overlap)}")
    print(f"additive_digest={body_digest(additive)}")
    print(f"new_union_digest={body_digest(new_union)}")
    print(f"route_roots={tuple(sorted(route_roots))}")
    print(f"additive_roots={tuple(sorted(additive))}")
    print(f"failed_roots={tuple(sorted(failed_roots))}")
    if failed_children:
        closest = max(
            failed_children,
            key=lambda row: (
                row["margin"],
                tuple(-x for x in row["body"]),
                -row["rank"],
                -row["apex"],
                -row["center"],
            ),
        )
        print(
            "closest_failure="
            f"{closest['body']};rank={closest['rank']};a={closest['apex']};"
            f"x={closest['center']};margin={closest['margin']}"
        )
    print("all_recomposition_controls=PASS")


if __name__ == "__main__":
    main()
