#!/usr/bin/env python3
"""Exact child-top-four scout on every open THM-2904 centre.

This is a discovery generalization of the THM-2912 one-H3-row mechanism.
It reads the locked THM-2904 hostile-centre ledger, retains every pivot
which the parent Hunter bound did not close, reconstructs its literal
child in both sequential and direct-denominator forms, and applies the
same exact head-2000 plus discrepancy-tail seal.

The default is a deterministic evenly spaced sample of open parent rows.
Pass ``--all`` for the full 51,222-centre universe.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PARENT_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py"
)
PARENT_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_ranked_h1_thm2911/"
    "thm2904_hostile_centre.ledger.out"
)
PARENT_SHA256 = (
    "99f1938f264d90c2b34ec3c64566605cc8fd12520424ad2f5cd0957342202ba0"
)
PARENT_LEDGER_SHA256 = (
    "bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a"
)
PARENT_SEMANTIC_SHA256 = (
    "ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab"
)

FIRST_EXTERNAL = 15
HEAD_LIMIT = 2_000
MAX_HEAD_LIMIT = 100_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(text: str) -> F:
    return F(text)


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def load_parent():
    require(file_sha256(PARENT_PATH) == PARENT_SHA256, "THM2904 source changed")
    spec = importlib.util.spec_from_file_location(
        "all_open_child_top4_parent",
        PARENT_PATH,
    )
    require(spec is not None and spec.loader is not None, "cannot load THM2904")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_parent()


def parse_ranked_pairs(
    text: str,
    field_count: int,
) -> tuple[tuple[object, ...], ...]:
    if not text:
        return ()
    out = []
    for item in text.split(","):
        fields = item.split(":")
        require(len(fields) == field_count, "ranked field width changed")
        out.append(tuple(fields))
    return tuple(out)


def parse_parent_rows() -> list[dict[str, object]]:
    require(
        file_sha256(PARENT_LEDGER) == PARENT_LEDGER_SHA256,
        "locked THM2904 hostile-centre ledger changed",
    )
    lines = PARENT_LEDGER.read_text().splitlines()
    require(
        lines[-2] == f"ledger_sha256={PARENT_SEMANTIC_SHA256}"
        and lines[-1].endswith("not LRC14"),
        "THM2904 ledger footer changed",
    )
    rows = []
    for line in lines:
        if not line.startswith("H1;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        body = tuple(map(int, fields["E"].split(",")))
        prefix = tuple(map(int, fields["P"].split(",")))
        core_rows = tuple(
            (parse_fraction(value), int(label))
            for label, value in parse_ranked_pairs(fields["H1"], 2)
        )
        pivots = tuple(
            (
                int(label),
                parse_fraction(bound),
                parse_fraction(margin),
                bool(int(closed)),
                bool(int(repaired)),
            )
            for label, bound, margin, closed, repaired in parse_ranked_pairs(
                fields["pivot"],
                5,
            )
        )
        require(
            tuple(label for _, label in core_rows)
            == tuple(label for label, *_ in pivots),
            "core/pivot order changed",
        )
        rows.append(
            {
                "body": body,
                "rank": int(fields["rank"]),
                "apex": int(fields["a"]),
                "prefix": prefix,
                "mass": parse_fraction(fields["h"]),
                "threshold": parse_fraction(fields["lambda"]),
                "core_rows": core_rows,
                "pivot_rows": pivots,
            }
        )
    require(len(rows) == 11_842, "THM2904 row count changed")
    require(
        sum(len(row["core_rows"]) for row in rows) == 55_293,
        "THM2904 centre count changed",
    )
    require(
        sum(
            not pivot[3]
            for row in rows
            for pivot in row["pivot_rows"]
        )
        == 51_222,
        "THM2904 open-centre count changed",
    )
    require(
        sum(
            any(not pivot[3] for pivot in row["pivot_rows"])
            for row in rows
        )
        == 11_563,
        "THM2904 open-row count changed",
    )
    return rows


def child_row(
    row: dict[str, object],
    index: int,
    carrier: list[tuple[F, F]],
    components: int,
    max_head_limit: int,
) -> dict[str, object]:
    coverage, center = row["core_rows"][index]
    _, _, parent_margin, parent_closed, _ = row["pivot_rows"][index]
    require(not parent_closed, "closed parent pivot entered child census")

    child = H.R.subtract_local(carrier, center)
    child_mass = sum((right - left for left, right in child), F(0))
    require(
        child_mass == row["mass"] - coverage and child_mass > 0,
        "child mass changed",
    )
    direct_child, direct_components, direct_mass = H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"], center)))
    )
    require(
        child == direct_child
        and len(child) == direct_components
        and child_mass == direct_mass,
        "subtracted child differs from direct denominator reconstruction",
    )

    earlier = tuple(label for _, label in row["core_rows"][:index])
    forbidden = frozenset((*row["prefix"], center, *earlier))
    require(
        len(forbidden) == len(row["prefix"]) + 1 + len(earlier)
        and not (set(earlier) & set(row["prefix"]))
        and center not in row["prefix"]
        and center not in earlier,
        "prefix/centre/earlier-centre exclusions collided",
    )
    require(max(forbidden) <= HEAD_LIMIT, "forbidden label escaped sealed head")
    labels = [
        label
        for label in range(FIRST_EXTERNAL, HEAD_LIMIT + 1)
        if label not in forbidden
    ]
    head = H.exact_coverages(child, labels)
    require(len(head) >= 4, "child head has fewer than four labels")
    top_four = tuple(sorted(head, reverse=True)[:4])
    require(
        len({label for _, label in top_four}) == 4
        and all(label not in forbidden for _, label in top_four),
        "child top four violated distinctness or exclusion",
    )
    gamma = H.S2 * len(child) / 7
    excess = top_four[-1][0] - child_mass / 7
    initial_required_limit = (
        None
        if excess <= 0
        else max(HEAD_LIMIT, ceiling(gamma / excess) - 1)
    )
    evaluated_limit = HEAD_LIMIT
    if (
        initial_required_limit is not None
        and initial_required_limit > HEAD_LIMIT
    ):
        evaluated_limit = min(initial_required_limit, max_head_limit)
        extra_labels = [
            label
            for label in range(HEAD_LIMIT + 1, evaluated_limit + 1)
            if label not in forbidden
        ]
        head.extend(H.exact_coverages(child, extra_labels))
        top_four = tuple(sorted(head, reverse=True)[:4])
    final_excess = top_four[-1][0] - child_mass / 7
    required_limit = (
        None
        if final_excess <= 0
        else max(HEAD_LIMIT, ceiling(gamma / final_excess) - 1)
    )
    tail_bound = child_mass / 7 + gamma / (evaluated_limit + 1)
    tail_gap = top_four[-1][0] - tail_bound
    sealed = tail_gap >= 0
    require(
        sealed
        == (required_limit is not None and required_limit <= evaluated_limit),
        "dynamic child-tail seal logic changed",
    )
    margin = child_mass - sum((value for value, _ in top_four), F(0))
    return {
        "body": row["body"],
        "rank": row["rank"],
        "apex": row["apex"],
        "prefix": row["prefix"],
        "parent_components": components,
        "center": center,
        "earlier": earlier,
        "parent_margin": parent_margin,
        "child_mass": child_mass,
        "child_components": len(child),
        "head_limit": evaluated_limit,
        "head_scanned": len(head),
        "required_limit": required_limit,
        "tail_bound": tail_bound,
        "tail_gap": tail_gap,
        "sealed": sealed,
        "top_four": top_four,
        "margin": margin,
        "closed": sealed and margin > 0,
    }


def profile_parent(
    task: tuple[dict[str, object], int],
) -> list[dict[str, object]]:
    row, max_head_limit = task
    carrier, components, mass = H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(mass == row["mass"] and mass > 0, "parent reconstruction changed")
    return [
        child_row(row, index, carrier, components, max_head_limit)
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]


def evenly_spaced(rows: list[dict[str, object]], count: int) -> list[dict[str, object]]:
    require(count >= 1, "sample size must be positive")
    if count >= len(rows):
        return rows
    if count == 1:
        return [rows[len(rows) // 2]]
    indices = tuple(
        (index * (len(rows) - 1)) // (count - 1)
        for index in range(count)
    )
    require(len(set(indices)) == count, "sample indices collided")
    return [rows[index] for index in indices]


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"x={row['center']};earlier={','.join(map(str, row['earlier']))};"
        f"pmargin={ftext(row['parent_margin'])};"
        f"hR={ftext(row['child_mass'])};rR={row['child_components']};"
        f"M={row['head_limit']};scan={row['head_scanned']};"
        f"N={row['required_limit']};sealed={int(row['sealed'])};"
        f"tail={ftext(row['tail_bound'])};tailgap={ftext(row['tail_gap'])};"
        + "top4="
        + ",".join(
            f"{label}:{ftext(value)}"
            for value, label in row["top_four"]
        )
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--sample-rows", type=int, default=256)
    parser.add_argument("--max-head", type=int, default=MAX_HEAD_LIMIT)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--emit-ledger", action="store_true")
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.max_head >= HEAD_LIMIT, "max-head is below initial head")

    rows = parse_parent_rows()
    open_rows = [
        row
        for row in rows
        if any(not pivot[3] for pivot in row["pivot_rows"])
    ]
    selected = open_rows if args.all else evenly_spaced(open_rows, args.sample_rows)
    tasks = [(row, args.max_head) for row in selected]
    if args.workers == 1:
        nested = list(map(profile_parent, tasks))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            nested = pool.map(profile_parent, tasks, chunksize=1)
    children = [child for row in nested for child in row]
    children.sort(
        key=lambda row: (
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        )
    )
    require(children, "empty child sample")

    digest = hashlib.sha256(b"LRC14/j6/all-open-centres/child-top4/scout/v1\n")
    for row in children:
        line = ledger_line(row)
        digest.update(line.encode())
        if args.emit_ledger:
            print("CHILD;" + line.rstrip())

    by_parent = defaultdict(list)
    by_body = defaultdict(list)
    for row in children:
        parent = (row["body"], row["rank"], row["apex"])
        by_parent[parent].append(row)
        by_body[row["body"]].append(row)
    closed_parents = sum(
        all(row["closed"] for row in group)
        for group in by_parent.values()
    )
    closed_sample_bodies = sum(
        all(row["closed"] for row in group)
        for group in by_body.values()
    )
    closed = [row for row in children if row["closed"]]
    failed = [row for row in children if not row["closed"]]
    closest_positive = min(
        closed,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        ),
    ) if closed else None
    closest_failure = max(
        failed,
        key=lambda row: (
            row["margin"],
            tuple(-value for value in row["body"]),
            -row["rank"],
            -row["apex"],
            -row["center"],
        ),
    ) if failed else None

    print("LRC14 all-open-centre exact child-top-four scout")
    print(
        f"mode={'ALL' if args.all else 'EVEN_SAMPLE'};"
        f"selected_parent_rows={len(selected)};"
        f"selected_children={len(children)}"
    )
    print(
        f"closed_children={len(closed)};failed_children={len(failed)};"
        f"sealed_children={sum(row['sealed'] for row in children)};"
        f"unsealed_children={sum(not row['sealed'] for row in children)};"
        f"closed_parent_rows={closed_parents}/{len(by_parent)};"
        f"closed_sample_bodies={closed_sample_bodies}/{len(by_body)}"
    )
    print(
        f"extended_children={sum(row['head_limit'] > HEAD_LIMIT for row in children)};"
        f"maximum_head_limit={max(row['head_limit'] for row in children)};"
        f"maximum_required_limit="
        f"{max((row['required_limit'] or 0) for row in children)}"
    )
    print(
        "child_component_histogram="
        + repr(tuple(sorted(Counter(row["child_components"] for row in children).items())))
    )
    print(
        "open_centres_per_selected_parent="
        + repr(tuple(sorted(Counter(len(group) for group in by_parent.values()).items())))
    )
    print(
        "closest_positive="
        + (
            "None"
            if closest_positive is None
            else (
                f"{closest_positive['body']};rank={closest_positive['rank']};"
                f"a={closest_positive['apex']};x={closest_positive['center']};"
                f"margin={ftext(closest_positive['margin'])}"
            )
        )
    )
    print(
        "closest_failure="
        + (
            "None"
            if closest_failure is None
            else (
                f"{closest_failure['body']};rank={closest_failure['rank']};"
                f"a={closest_failure['apex']};x={closest_failure['center']};"
                f"margin={ftext(closest_failure['margin'])}"
            )
        )
    )
    print(f"minimum_tail_gap={ftext(min(row['tail_gap'] for row in children))}")
    print(f"semantic_digest={digest.hexdigest()}")
    print(
        "scope=G5-surviving THM2904 rows;every inherited open ordered centre;"
        "exact child top4;sample unless mode=ALL;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
