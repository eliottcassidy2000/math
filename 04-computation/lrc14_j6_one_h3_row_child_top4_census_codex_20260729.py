#!/usr/bin/env python3
r"""Exact child top-four census on the current one-H3-row LRC14 roots.

THM-2904 allocates every hypothetical five-cover to the earliest
maximum-coverage centre in its finite hostile H1 core.  Once a centre is
fixed, the residual carrier must be covered by four distinct allowed
labels, excluding the marked prefix, the centre, and every earlier core
label.

This verifier reconstructs the full THM-2904 parent universe, selects the
210 current residual bodies having exactly one ordinary H3 row after
THM-2907, and computes the exact four largest singleton coverages on every
open centre child.  Labels through 2,000 are evaluated literally.  The
THM-735 discrepancy bound seals every omitted label below the retained
fourth rank, so a strict top-four sum deficit excludes that centre.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py"
)
PARENT_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.out"
)
H4_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.out"
)
H4_ENDPOINT_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.out"
)
PARENT_SHA256 = (
    "99f1938f264d90c2b34ec3c64566605cc8fd12520424ad2f5cd0957342202ba0"
)
PARENT_OUT_SHA256 = (
    "0933c67a108b6d588e36737fb2b17b325ca36146976cfb035bebe036a6234036"
)
H4_OUT_SHA256 = (
    "1df929e106cd16c094886d59f3702ba9bafa395ee906527fe4592a1552e9b458"
)
H4_ENDPOINT_OUT_SHA256 = (
    "a93a4a724dac6c55806f3358c2f5ab25de8f0261c92906a0161414781a717d20"
)

FIRST_EXTERNAL = 15
HEAD_LIMIT = 2_000

EXPECTED_COUNTS: tuple[object, ...] | None = (
    11_842,
    52,
    11_511,
    76,
    11_435,
    3_344,
    (
        (1, 210),
        (2, 726),
        (3, 943),
        (4, 754),
        (5, 449),
        (6, 181),
        (7, 60),
        (8, 13),
        (9, 7),
        (10, 1),
    ),
    210,
    807,
    1_599_450,
    20,
    56,
    807,
    765,
    42,
    172,
    38,
    260,
    3_172,
)
EXPECTED_CLOSED_ROOT_DIGEST: str | None = (
    "fcb9e88dd0b5e8aa855e0164636a17f42b4f1620407b9a0a780e50fd48336217"
)
EXPECTED_LEDGER_SHA256: str | None = (
    "39f3a9fb8ec6447baf96512bee3ee174e2390639ab3bbf0a6a36dcb5cdf0274e"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_output_text(path: Path, prefix: str) -> str:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return matches[0]


def literal_output_value(path: Path, prefix: str) -> object:
    return ast.literal_eval(unique_output_text(path, prefix))


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def load_parent():
    require(file_sha256(PARENT_PATH) == PARENT_SHA256, "THM2904 source changed")
    spec = importlib.util.spec_from_file_location("thm2904_child_top4", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM2904")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_parent()


def check_dependency_outputs() -> None:
    require(
        file_sha256(PARENT_OUT) == PARENT_OUT_SHA256,
        "THM2904 output changed",
    )
    require(file_sha256(H4_OUT) == H4_OUT_SHA256, "THM2907 output changed")
    require(
        file_sha256(H4_ENDPOINT_OUT) == H4_ENDPOINT_OUT_SHA256,
        "THM2907 endpoint output changed",
    )
    parent_counts = literal_output_value(PARENT_OUT, "counts=")
    require(
        parent_counts[0] == 11_842
        and parent_counts[1] == 52
        and parent_counts[11:26]
        == (
            55_293,
            4_071,
            51_222,
            4,
            4,
            0,
            279,
            0,
            11_563,
            3_411,
            10,
            4,
            6,
            88,
            3_344,
        ),
        "THM2904 root/pivot controls changed",
    )
    require(
        unique_output_text(PARENT_OUT, "ledger_sha256=")
        == "ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab"
        and unique_output_text(PARENT_OUT, "mode=") == "LOCKED"
        and unique_output_text(PARENT_OUT, "all_exact_controls=") == "PASS",
        "THM2904 replay controls changed",
    )
    h4_counts = literal_output_value(H4_OUT, "counts=")
    require(
        h4_counts[:10] == (18_290, 16_357, 1_933, 179, 529, 1_225, 18_280, 10, 51, 1)
        and h4_counts[24:26] == (52, 52),
        "THM2907 branch controls changed",
    )
    require(
        unique_output_text(H4_OUT, "semantic_digest=")
        == "e8d4119f101ae3ac310fe5ca8a056607390ff4d7aa166cb90168983ea7069356"
        and unique_output_text(H4_OUT, "mode=") == "LOCKED"
        and unique_output_text(H4_OUT, "all_exact_controls=") == "PASS"
        and unique_output_text(H4_ENDPOINT_OUT, "semantic_digest=")
        == "b521e84294bc8b7ac9f431ec7cc6841bbc326bddf26ce60fcd546c55c9479c05"
        and unique_output_text(H4_ENDPOINT_OUT, "mode=") == "LOCKED"
        and unique_output_text(H4_ENDPOINT_OUT, "all_exact_controls=") == "PASS",
        "THM2907 endpoint controls changed",
    )


def compute_parent(row: dict[str, object]) -> dict[str, object]:
    return H.actual_core(row)


def current_proved_union() -> set[tuple[int, ...]]:
    through_2905 = H.current_proved_union()
    require(len(through_2905) == 82, "THM2905 union changed")
    additive_2904 = set(H.EXPECTED_ADDITIVE_ROOTS)
    require(len(additive_2904) == 6, "THM2904 additive roots changed")
    current = through_2905 | additive_2904
    require(len(current) == 88, "proved union through THM2904 changed")
    return current


def residual_one_rows(
    rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], tuple[tuple[int, int], ...]]:
    proved = current_proved_union()
    ordinary_all = [
        row
        for row in rows
        if not row["pair_exception"] and not row["branch_closed"]
    ]
    residual = [row for row in ordinary_all if row["body"] not in proved]
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in residual:
        by_body[row["body"]].append(row)
    histogram = tuple(sorted(Counter(map(len, by_body.values())).items()))
    require(
        len(ordinary_all) == 11_511
        and len(ordinary_all) - len(residual) == 76
        and len(residual) == 11_435
        and len(by_body) == 3_344
        and histogram
        == (
            (1, 210),
            (2, 726),
            (3, 943),
            (4, 754),
            (5, 449),
            (6, 181),
            (7, 60),
            (8, 13),
            (9, 7),
            (10, 1),
        ),
        "current residual H3 universe changed",
    )
    one_rows = [body_rows[0] for body_rows in by_body.values() if len(body_rows) == 1]
    one_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    require(len(one_rows) == 210, "one-H3-row body count changed")
    return one_rows, histogram


def child_task(
    item: tuple[dict[str, object], int],
) -> dict[str, object]:
    row, index = item
    coverage, center = row["core_rows"][index]
    _, _, _, pivot_closed, _ = row["pivot_rows"][index]
    require(not pivot_closed, "closed parent pivot entered child census")
    carrier, components, mass = H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(
        components == row["components"] and mass == row["mass"],
        "parent carrier reconstruction changed",
    )
    child = H.R.subtract_local(carrier, center)
    child_mass = sum((right - left for left, right in child), F(0))
    require(child_mass == mass - coverage and child_mass > 0, "child mass changed")

    earlier = tuple(label for _, label in row["core_rows"][:index])
    forbidden = frozenset((*row["prefix"], center, *earlier))
    labels = [
        label
        for label in range(FIRST_EXTERNAL, HEAD_LIMIT + 1)
        if label not in forbidden
    ]
    head = H.exact_coverages(child, labels)
    require(len(head) >= 4, "child head has fewer than four labels")
    top_four = tuple(sorted(head, reverse=True)[:4])
    child_components = len(child)
    tail_bound = (
        child_mass / 7
        + H.S2 * child_components / (7 * (HEAD_LIMIT + 1))
    )
    tail_gap = top_four[-1][0] - tail_bound
    require(
        tail_gap >= 0,
        "head limit does not seal the fourth child singleton rank",
    )
    top_four_sum = sum((value for value, _ in top_four), F(0))
    margin = child_mass - top_four_sum
    return {
        "body": row["body"],
        "rank": row["rank"],
        "apex": row["apex"],
        "prefix": row["prefix"],
        "center": center,
        "earlier": earlier,
        "child_mass": child_mass,
        "child_components": child_components,
        "head_scanned": len(labels),
        "tail_bound": tail_bound,
        "tail_gap": tail_gap,
        "top_four": top_four,
        "margin": margin,
        "closed": margin > 0,
    }


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"x={row['center']};earlier={','.join(map(str, row['earlier']))};"
        f"hR={ftext(row['child_mass'])};rR={row['child_components']};"
        f"scan={row['head_scanned']};tail={ftext(row['tail_bound'])};"
        f"tailgap={ftext(row['tail_gap'])};"
        + "top4="
        + ",".join(f"{label}:{ftext(value)}" for value, label in row["top_four"])
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    check_dependency_outputs()

    inputs = H.survivor_inputs()
    context = mp.get_context("spawn")
    if args.workers == 1:
        parent_rows = list(map(compute_parent, inputs))
    else:
        with context.Pool(args.workers) as pool:
            parent_rows = pool.map(compute_parent, inputs)
    parent_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    one_rows, row_histogram = residual_one_rows(parent_rows)

    tasks = [
        (row, index)
        for row in one_rows
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]
    require(len(tasks) == 807, "one-row open-pivot count changed")
    if args.workers == 1:
        children = list(map(child_task, tasks))
    else:
        with context.Pool(args.workers) as pool:
            children = pool.map(child_task, tasks)
    children.sort(
        key=lambda row: (row["body"], row["rank"], row["apex"], row["center"])
    )

    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in children:
        by_body[row["body"]].append(row)
    closed_roots = tuple(
        sorted(
            body
            for body, body_rows in by_body.items()
            if all(row["closed"] for row in body_rows)
        )
    )
    failed_roots = tuple(sorted(set(by_body) - set(closed_roots)))
    prior = current_proved_union()
    require(not (set(closed_roots) & prior), "new child roots overlap prior union")
    new_union = prior | set(closed_roots)

    counts = (
        len(parent_rows),
        sum(row["pair_exception"] for row in parent_rows),
        11_511,
        76,
        11_435,
        3_344,
        row_histogram,
        len(one_rows),
        len(children),
        sum(row["head_scanned"] for row in children),
        min(row["child_components"] for row in children),
        max(row["child_components"] for row in children),
        sum(row["tail_gap"] >= 0 for row in children),
        sum(row["closed"] for row in children),
        sum(not row["closed"] for row in children),
        len(closed_roots),
        len(failed_roots),
        len(new_union),
        3_432 - len(new_union),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "one-row child counts changed")

    root_digest = hashlib.sha256(repr(closed_roots).encode()).hexdigest()
    if EXPECTED_CLOSED_ROOT_DIGEST is not None:
        require(
            root_digest == EXPECTED_CLOSED_ROOT_DIGEST,
            "one-row closed-root list changed",
        )
    ledger = hashlib.sha256()
    ledger.update(b"LRC14/j6/one-H3-row/child-top4/v1\n")
    for row in children:
        ledger.update(ledger_line(row).encode())
    ledger_sha256 = ledger.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "child ledger changed")

    closed_children = [row for row in children if row["closed"]]
    failed_children = [row for row in children if not row["closed"]]
    closest_positive = min(
        closed_children,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["rank"],
            row["center"],
        ),
    )
    closest_failure = max(
        failed_children,
        key=lambda row: (
            row["margin"],
            tuple(-value for value in row["body"]),
            -row["rank"],
            -row["center"],
        ),
    )

    print("LRC14 one-H3-row exact child top-four census")
    print(f"head_limit={HEAD_LIMIT}")
    print(f"counts={counts}")
    print(f"closed_roots={closed_roots}")
    print(f"failed_roots={failed_roots}")
    print(f"closed_root_digest={root_digest}")
    print(
        "closest_positive="
        f"{closest_positive['body']};rank={closest_positive['rank']};"
        f"a={closest_positive['apex']};x={closest_positive['center']};"
        f"margin={ftext(closest_positive['margin'])}"
    )
    print(
        "closest_failure="
        f"{closest_failure['body']};rank={closest_failure['rank']};"
        f"a={closest_failure['apex']};x={closest_failure['center']};"
        f"margin={ftext(closest_failure['margin'])}"
    )
    print(
        f"minimum_tail_gap={ftext(min(row['tail_gap'] for row in children))}"
    )
    print(f"ledger_sha256={ledger_sha256}")
    print(
        "mode=DISCOVERY"
        if any(
            value is None
            for value in (
                EXPECTED_COUNTS,
                EXPECTED_CLOSED_ROOT_DIGEST,
                EXPECTED_LEDGER_SHA256,
            )
        )
        else "mode=LOCKED"
    )
    print(
        "scope=210 current residual bodies with one ordinary H3 row;"
        "807 ordered centre children;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
