#!/usr/bin/env python3
r"""Exact ordered-child top-four census on the post-THM2913 two-H3 stratum.

THM-2904 allocates a hypothetical five-cover to the earliest one of its
maximum-coverage labels in a finite hostile-centre core.  After deleting
that centre, the four remaining labels must avoid the inherited marked
prefix, the centre, and all earlier hostile centres.

The post-THM2913 residual has no one-H3-row body.  This program isolates
its exact next stratum, the bodies with two surviving ordinary H3 rows,
and tests every open ordered-centre child on its literal residual carrier.
Instead of imposing one uniform scan horizon, it doubles a finite head
until the THM-735 discrepancy tail lies below the retained fourth
singleton rank.  The resulting top four are therefore global.
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
BASE_PATH = (
    ROOT
    / "04-computation/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py"
)
BASE_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out"
)
THM2913_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out"
)

FIRST_EXTERNAL = 15
INITIAL_HORIZON = 127
MAX_HORIZON = 4_095

# Discovery locks.  Fill these only after an independent ordinary/optimized
# replay and exact root recomposition audit.
EXPECTED_COUNTS: tuple[object, ...] | None = None
EXPECTED_CLOSED_ROOT_DIGEST: str | None = None
EXPECTED_ADDITIVE_ROOT_DIGEST: str | None = None
EXPECTED_LEDGER_SHA256: str | None = None


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


def load_base():
    spec = importlib.util.spec_from_file_location("thm2912_two_h3", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM2912")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = load_base()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def interval_mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def literal_output_value(path: Path, prefix: str) -> object:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return ast.literal_eval(matches[0])


def compute_parent(row: dict[str, object]) -> dict[str, object]:
    return C.compute_parent(row)


def current_union_351() -> set[tuple[int, ...]]:
    baseline = C.proved_union_through_2911()
    top4 = set(literal_output_value(BASE_OUT, "closed_roots="))
    repaired = set(literal_output_value(THM2913_OUT, "closed_roots="))
    require(
        len(baseline) == 181
        and len(top4) == 172
        and len(repaired) == 38
        and len(baseline | top4) == 314
        and len(baseline | top4 | repaired) == 351,
        "THM2911--2913 root union changed",
    )
    return baseline | top4 | repaired


def residual_two_rows(
    rows: list[dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
]:
    current = current_union_351()
    ordinary = [
        row
        for row in rows
        if not row["pair_exception"] and not row["branch_closed"]
    ]
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in ordinary:
        by_body[row["body"]].append(row)
    residual = {
        body: body_rows
        for body, body_rows in by_body.items()
        if body not in current
    }
    all_histogram = tuple(sorted(Counter(map(len, by_body.values())).items()))
    residual_histogram = tuple(
        sorted(Counter(map(len, residual.values())).items())
    )
    require(
        len(rows) == 11_842
        and len(ordinary) == 11_511
        and len(by_body) == 3_401
        and len(current) == 351
        and len(residual) == 3_081
        and all_histogram
        == (
            (1, 262),
            (2, 726),
            (3, 944),
            (4, 756),
            (5, 450),
            (6, 181),
            (7, 60),
            (8, 14),
            (9, 7),
            (10, 1),
        )
        and residual_histogram
        == (
            (2, 690),
            (3, 929),
            (4, 751),
            (5, 449),
            (6, 181),
            (7, 60),
            (8, 13),
            (9, 7),
            (10, 1),
        ),
        "post-THM2913 residual stratification changed",
    )
    selected = [
        row
        for body_rows in residual.values()
        if len(body_rows) == 2
        for row in body_rows
    ]
    selected.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    require(
        len(selected) == 1_380
        and len({row["body"] for row in selected}) == 690,
        "two-H3-row stratum changed",
    )
    return selected, all_histogram, residual_histogram


def dynamic_ranked_head(
    carrier: list[tuple[F, F]],
    forbidden: frozenset[int],
    rank: int,
) -> tuple[tuple[tuple[F, int], ...], F, int, int]:
    mass = interval_mass(carrier)
    gamma = C.H.S2 * len(carrier) / 7
    scanned_rows: list[tuple[F, int]] = []
    first = FIRST_EXTERNAL
    horizon = INITIAL_HORIZON
    while True:
        labels = [
            label
            for label in range(first, horizon + 1)
            if label not in forbidden
        ]
        scanned_rows.extend(C.H.exact_coverages(carrier, labels))
        require(len(scanned_rows) >= rank, "finite head has too few labels")
        ranked = tuple(
            sorted(scanned_rows, key=lambda item: (-item[0], item[1]))
        )
        tail = mass / 7 + gamma / (horizon + 1)
        if ranked[rank - 1][0] >= tail:
            return ranked[:rank], tail, horizon, len(scanned_rows)
        require(
            horizon < MAX_HORIZON,
            f"rank {rank} not sealed by horizon {MAX_HORIZON}",
        )
        first = horizon + 1
        horizon = min(2 * horizon + 1, MAX_HORIZON)


def child_task(
    item: tuple[dict[str, object], int],
) -> dict[str, object]:
    row, index = item
    coverage, center = row["core_rows"][index]
    _, _, _, pivot_closed, _ = row["pivot_rows"][index]
    require(not pivot_closed, "closed parent pivot entered child census")
    parent, parent_components, parent_mass = C.H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(
        parent_components == row["components"] and parent_mass == row["mass"],
        "parent carrier reconstruction changed",
    )
    child = C.H.R.subtract_local(parent, center)
    child_mass = interval_mass(child)
    direct, direct_components, direct_mass = C.H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"], center)))
    )
    require(
        child == direct
        and len(child) == direct_components
        and child_mass == direct_mass
        and child_mass == parent_mass - coverage
        and child_mass > 0,
        "literal/direct child reconstruction changed",
    )

    earlier = tuple(label for _, label in row["core_rows"][:index])
    forbidden = frozenset((*row["prefix"], center, *earlier))
    top_four, tail, horizon, scanned = dynamic_ranked_head(
        child,
        forbidden,
        4,
    )
    tail_gap = top_four[-1][0] - tail
    margin = child_mass - sum((value for value, _ in top_four), F(0))
    return {
        "body": row["body"],
        "rank": row["rank"],
        "apex": row["apex"],
        "prefix": row["prefix"],
        "center": center,
        "earlier": earlier,
        "child_mass": child_mass,
        "child_components": len(child),
        "horizon": horizon,
        "scanned": scanned,
        "tail": tail,
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
        f"h={ftext(row['child_mass'])};r={row['child_components']};"
        f"N={row['horizon']};scan={row['scanned']};"
        f"tail={ftext(row['tail'])};tailgap={ftext(row['tail_gap'])};"
        + "top4="
        + ",".join(
            f"{label}:{ftext(value)}" for value, label in row["top_four"]
        )
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}\n"
    )


def boundary_text(row: dict[str, object]) -> str:
    return (
        f"{row['body']};rank={row['rank']};a={row['apex']};"
        f"x={row['center']};margin={ftext(row['margin'])}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    C.check_dependency_outputs()

    context = mp.get_context("spawn")
    inputs = C.H.survivor_inputs()
    if args.workers == 1:
        parent_rows = list(map(compute_parent, inputs))
    else:
        with context.Pool(args.workers) as pool:
            parent_rows = pool.map(compute_parent, inputs)
    parent_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    selected_rows, all_histogram, residual_histogram = residual_two_rows(
        parent_rows
    )
    tasks = [
        (row, index)
        for row in selected_rows
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]
    require(len(tasks) == 5_618, "two-H3-row open-pivot count changed")
    if args.workers == 1:
        children = list(map(child_task, tasks))
    else:
        with context.Pool(args.workers) as pool:
            children = pool.map(child_task, tasks)
    children.sort(
        key=lambda row: (
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        )
    )

    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in children:
        by_body[row["body"]].append(row)
    require(len(by_body) == 690, "two-H3 child body count changed")
    closed_roots = tuple(
        sorted(
            body
            for body, body_rows in by_body.items()
            if all(row["closed"] for row in body_rows)
        )
    )
    failed_roots = tuple(sorted(set(by_body) - set(closed_roots)))
    baseline = current_union_351()
    overlap = tuple(sorted(set(closed_roots) & baseline))
    additive_roots = tuple(sorted(set(closed_roots) - baseline))
    require(not overlap, "residual two-H3 root overlapped current union")
    new_union = baseline | set(closed_roots)

    horizon_histogram = tuple(
        sorted(Counter(row["horizon"] for row in children).items())
    )
    counts = (
        len(parent_rows),
        11_511,
        3_401,
        len(baseline),
        3_432 - len(baseline),
        all_histogram,
        residual_histogram,
        len(selected_rows),
        len(by_body),
        len(children),
        sum(row["scanned"] for row in children),
        horizon_histogram,
        min(row["child_components"] for row in children),
        max(row["child_components"] for row in children),
        sum(row["tail_gap"] >= 0 for row in children),
        sum(row["margin"] > 0 for row in children),
        sum(row["margin"] == 0 for row in children),
        sum(row["margin"] < 0 for row in children),
        len(closed_roots),
        len(failed_roots),
        len(additive_roots),
        len(new_union),
        3_432 - len(new_union),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "two-H3 child counts changed")

    closed_digest = hashlib.sha256(repr(closed_roots).encode()).hexdigest()
    additive_digest = hashlib.sha256(repr(additive_roots).encode()).hexdigest()
    if EXPECTED_CLOSED_ROOT_DIGEST is not None:
        require(
            closed_digest == EXPECTED_CLOSED_ROOT_DIGEST,
            "two-H3 closed roots changed",
        )
    if EXPECTED_ADDITIVE_ROOT_DIGEST is not None:
        require(
            additive_digest == EXPECTED_ADDITIVE_ROOT_DIGEST,
            "two-H3 additive roots changed",
        )
    ledger_lines = [ledger_line(row) for row in children]
    ledger_hash = hashlib.sha256()
    ledger_hash.update(b"LRC14/j6/two-H3-row/child-top4/dynamic-tail/v1\n")
    for line in ledger_lines:
        ledger_hash.update(line.encode())
    ledger_sha256 = ledger_hash.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "child ledger changed")
    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 two-H3-row dynamic-tail child top-four ledger\n"
            + "".join(ledger_lines)
            + f"ledger_sha256={ledger_sha256}\n"
            + "scope=690 post-THM2913 two-H3 roots;not LRC14\n"
        )

    positive = [row for row in children if row["closed"]]
    nonpositive = [row for row in children if not row["closed"]]
    closest_positive = min(
        positive,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["rank"],
            row["center"],
        ),
    )
    closest_failure = max(
        nonpositive,
        key=lambda row: (
            row["margin"],
            tuple(-value for value in row["body"]),
            -row["rank"],
            -row["center"],
        ),
    )
    lines = [
        "LRC14 post-THM2913 two-H3-row dynamic-tail child top-four census",
        (
            f"initial_horizon={INITIAL_HORIZON};"
            f"maximum_horizon={MAX_HORIZON}"
        ),
        f"counts={counts}",
        f"closed_roots={closed_roots}",
        f"failed_roots={failed_roots}",
        f"additive_roots={additive_roots}",
        f"closed_root_digest={closed_digest}",
        f"additive_root_digest={additive_digest}",
        f"closest_positive={boundary_text(closest_positive)}",
        f"closest_failure={boundary_text(closest_failure)}",
        (
            "minimum_tail_gap="
            f"{ftext(min(row['tail_gap'] for row in children))}"
        ),
        f"ledger_sha256={ledger_sha256}",
        (
            "mode=DISCOVERY"
            if any(
                value is None
                for value in (
                    EXPECTED_COUNTS,
                    EXPECTED_CLOSED_ROOT_DIGEST,
                    EXPECTED_ADDITIVE_ROOT_DIGEST,
                    EXPECTED_LEDGER_SHA256,
                )
            )
            else "mode=LOCKED"
        ),
        (
            "scope=690 post-THM2913 roots with exactly two ordinary H3 rows;"
            "all open ordered-centre children;not LRC14"
        ),
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
