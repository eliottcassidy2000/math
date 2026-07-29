#!/usr/bin/env python3
r"""Pair/Hunter/toothpick closure of the 367 THM-2916-open children.

The first repair stage applies the exact THM-2913 child pair cap and
four-slot Hunter envelope.  A failed Hunter row has a discrepancy-finite
second-centre core; after allocating the earliest possible second centre,
the remaining three labels are tested by their globally sealed top three.

Five ordered children retain one open second-centre pivot.  The residual
after that pivot still has THREE slots, so a top-two singleton sum is not
a lawful terminal.  This verifier instead computes its exact global pair
cap and the three-slot Hunter envelope.  All five close strictly by G3.
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
    / "04-computation/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.py"
)
BASE_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.out"
)
BASE_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_two_h3_row_child_top4_census_codex_20260729.ledger.out"
)
THM2913_PATH = (
    ROOT
    / "04-computation/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.py"
)
THM2913_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out"
)

BASE_SHA256 = "7d3c4e82abb0ce3af13c43c1e03f09d4be1ee3dbfd496ca05273672cd2a462a6"
BASE_OUT_SHA256 = (
    "6b31abaadd4f089a9f98a9eea49845c0ed0123a810bb4bf4f3c0309c519e7df9"
)
BASE_LEDGER_SHA256 = (
    "a202807a1da472ed6c19109425478922b2218c48f2afdb39c67a0b9f2f5fb7c4"
)
THM2913_SHA256 = (
    "14e56e124197cd1bdae841efa195a1e7c282e7ea368a610e5f4d56509431858b"
)
THM2913_OUT_SHA256 = (
    "3604644a9691b13e7fa245249b68c9586ec2775996834f7761f32eb0b89f3e64"
)

EXPECTED_COUNTS: tuple[object, ...] | None = (
    11_842,
    1_380,
    5_618,
    5_251,
    367,
    149,
    218,
    111,
    107,
    102,
    5,
    283,
    278,
    5,
    3,
    2,
    102,
    195,
    292,
    296,
    745,
    1_041,
    2_391,
    911_490,
    3_011,
    29_548,
    702_330,
    86,
    ((1, 17), (2, 33), (3, 32), (4, 21), (5, 4)),
    0,
    0,
    0,
    0,
    0,
)
EXPECTED_CLOSED_ROOT_DIGEST: str | None = (
    "e3045198e08804c78025bd532111377309882911e08bc50604aa7119ac266c71"
)
EXPECTED_FINAL_TWO_ROW_DIGEST: str | None = (
    "772f67e5711cb009012a9c4abeb1b9a288195126f382244231f1f93362b63efc"
)
EXPECTED_LEDGER_SHA256: str | None = (
    "e1434329c46118991b1fe357be5f87a01a22b81a814d760166ba3582de79e83a"
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


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_module("thm2916_recursive_repair", BASE_PATH)
T = load_module("thm2913_recursive_repair", THM2913_PATH)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def literal_output_value(path: Path, prefix: str) -> object:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return ast.literal_eval(matches[0])


def dependency_controls() -> tuple[set[tuple[int, ...]], set[tuple[int, ...]]]:
    require(file_sha256(BASE_PATH) == BASE_SHA256, "THM2916 source changed")
    require(file_sha256(BASE_OUT) == BASE_OUT_SHA256, "THM2916 output changed")
    require(
        file_sha256(BASE_LEDGER) == BASE_LEDGER_SHA256,
        "THM2916 ledger changed",
    )
    require(file_sha256(THM2913_PATH) == THM2913_SHA256, "THM2913 source changed")
    require(
        file_sha256(THM2913_OUT) == THM2913_OUT_SHA256,
        "THM2913 output changed",
    )
    base_text = BASE_OUT.read_text()
    thm2913_text = THM2913_OUT.read_text()
    require(
        "mode=LOCKED" in base_text
        and "all_exact_controls=PASS" in base_text
        and "mode=LOCKED" in thm2913_text
        and "all_exact_controls=PASS" in thm2913_text,
        "THM2913/2916 replay controls changed",
    )
    base_closed = set(literal_output_value(BASE_OUT, "closed_roots="))
    base_failed = set(literal_output_value(BASE_OUT, "failed_roots="))
    require(
        len(base_closed) == 394
        and len(base_failed) == 296
        and not (base_closed & base_failed),
        "THM2916 root partition changed",
    )
    return base_closed, base_failed


def compute_parent(row: dict[str, object]) -> dict[str, object]:
    return M.compute_parent(row)


def compute_child(
    item: tuple[dict[str, object], int],
) -> dict[str, object]:
    return M.child_task(item)


def analyze_first_repair(child_row: dict[str, object]) -> dict[str, object]:
    parent, _, _ = T.direct_carrier(
        tuple((*child_row["body"], child_row["apex"]))
    )
    child = T.C.H.R.subtract_local(parent, child_row["center"])
    mass = T.interval_mass(child)
    require(
        mass == child_row["child_mass"]
        and len(child) == child_row["child_components"],
        "THM2916 child reconstruction changed",
    )
    direct, components, direct_mass = T.direct_carrier(
        tuple((*child_row["body"], child_row["apex"], child_row["center"]))
    )
    require(
        child == direct and len(child) == components and mass == direct_mass,
        "literal/direct child disagreement",
    )
    forbidden = frozenset(
        (*child_row["prefix"], child_row["center"], *child_row["earlier"])
    )
    pair = T.exact_pair_cap(child, forbidden)
    require(
        pair["top4"] == child_row["top_four"],
        "THM2916 and THM2913 top-four ranks disagree",
    )
    pair_margin = mass - 2 * pair["head"]
    qs = tuple(value for value, _ in pair["top4"])
    envelope, threshold = T.hunter_data(qs, pair["head"], mass)
    hunter_margin = mass - envelope
    second_core: tuple[tuple[F, int], ...] = ()
    second_rows: tuple[dict[str, object], ...] = ()
    cutoff = None
    delta = None
    core_scan = 0
    if threshold is not None:
        delta = threshold - mass / 7
        require(delta > 0, "second-centre threshold is not discrepancy-finite")
        gamma = T.C.H.S2 * len(child) / 7
        cutoff = T.C.H.ceiling(gamma / delta) - 1
        require(cutoff >= T.FIRST_EXTERNAL, "second-centre cutoff is empty")
        labels = [
            label
            for label in range(T.FIRST_EXTERNAL, cutoff + 1)
            if label not in forbidden
        ]
        coverages = T.C.H.exact_coverages(child, labels)
        core_scan = len(labels)
        second_core = tuple(
            sorted(
                (
                    (value, label)
                    for value, label in coverages
                    if value >= threshold
                ),
                key=lambda item: (-item[0], item[1]),
            )
        )
        require(second_core, "hostile second-centre core is empty")
        require(
            mass / 7 + gamma / (cutoff + 1) <= threshold,
            "second-centre discrepancy tail did not seal",
        )
        second_rows = tuple(
            T.second_pivot((child_row, child, forbidden, second_core, index))
            for index in range(len(second_core))
        )
    if pair_margin > 0:
        route = "pair"
    elif hunter_margin > 0:
        route = "hunter"
    elif second_rows and all(row["closed"] for row in second_rows):
        route = "toothpick"
    else:
        require(
            threshold is not None
            and second_rows
            and sum(not row["closed"] for row in second_rows) == 1,
            "deep route does not have exactly one open second pivot",
        )
        route = "deep"
    return {
        **child_row,
        "pair": pair,
        "pair_margin": pair_margin,
        "hunter": envelope,
        "hunter_margin": hunter_margin,
        "threshold": threshold,
        "delta": delta,
        "cutoff": cutoff,
        "core_scan": core_scan,
        "second_core": second_core,
        "second_rows": second_rows,
        "route": route,
        "closed": route != "deep",
    }


def analyze_deep(row: dict[str, object]) -> dict[str, object]:
    open_rows = [second for second in row["second_rows"] if not second["closed"]]
    require(len(open_rows) == 1, "deep child lost its unique open second pivot")
    second = open_rows[0]
    second_forbidden = frozenset(
        (
            *row["prefix"],
            row["center"],
            *row["earlier"],
            second["center"],
            *second["earlier"],
        )
    )
    grandchild, components, mass = T.direct_carrier(
        tuple(
            (
                *row["body"],
                row["apex"],
                row["center"],
                second["center"],
            )
        )
    )
    require(
        components == second["components"]
        and mass == second["mass"],
        "deep grandchild direct reconstruction changed",
    )
    parent, _, _ = T.direct_carrier(
        tuple((*row["body"], row["apex"], row["center"]))
    )
    literal = T.C.H.R.subtract_local(parent, second["center"])
    require(literal == grandchild, "deep grandchild literal/direct mismatch")

    ranked, singleton_tail, scanned = T.exact_ranked_head(
        grandchild,
        second_forbidden,
        T.TOP3_HORIZON,
        3,
    )
    top3 = tuple(ranked[:3])
    require(
        top3 == second["top3"]
        and mass - sum((value for value, _ in top3), F(0))
        == second["margin"],
        "deep top-three identity changed",
    )
    pair = T.exact_pair_cap(grandchild, second_forbidden)
    require(pair["top4"][:3] == top3, "deep pair/top-three ranks disagree")
    qs = tuple(value for value, _ in top3)
    pair_partition_margin = mass - (qs[0] + pair["head"])
    envelope, threshold = T.hunter_data(qs, pair["head"], mass)
    hunter_margin = mass - envelope
    require(
        hunter_margin > 0 and threshold is None,
        "three-slot Hunter repair did not close",
    )
    return {
        "second": second,
        "forbidden": tuple(sorted(second_forbidden)),
        "mass": mass,
        "components": components,
        "scanned": scanned,
        "singleton_tail": singleton_tail,
        "top3": top3,
        "pair": pair,
        "pair_partition_margin": pair_partition_margin,
        "hunter": envelope,
        "hunter_margin": hunter_margin,
        "route": "pair3" if pair_partition_margin > 0 else "hunter3",
        "closed": True,
    }


def deep_line(parent: dict[str, object], row: dict[str, object]) -> str:
    second = row["second"]
    pair = row["pair"]
    return (
        f"DEEP;E={','.join(map(str, parent['body']))};rank={parent['rank']};"
        f"a={parent['apex']};x={parent['center']};y={second['center']};"
        f"F={','.join(map(str, row['forbidden']))};"
        f"h={ftext(row['mass'])};r={row['components']};scan={row['scanned']};"
        f"Stail={ftext(row['singleton_tail'])};"
        + "top3="
        + ",".join(f"{label}:{ftext(value)}" for value, label in row["top3"])
        + f";B2={ftext(pair['head'])};B2w={pair['witness']};"
        f"B2tail={ftext(pair['tail'])};B2gap={ftext(pair['tail_gap'])};"
        f"pair3margin={ftext(row['pair_partition_margin'])};"
        f"G3={ftext(row['hunter'])};G3margin={ftext(row['hunter_margin'])};"
        f"route={row['route']};closed=1\n"
    )


def boundary_text(row: dict[str, object], key: str) -> str:
    return (
        f"{row['body']};rank={row['rank']};a={row['apex']};"
        f"x={row['center']};margin={ftext(row[key])}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    base_closed_artifact, base_failed_artifact = dependency_controls()
    M.C.check_dependency_outputs()

    context = mp.get_context("spawn")
    inputs = M.C.H.survivor_inputs()
    if args.workers == 1:
        parent_rows = list(map(compute_parent, inputs))
    else:
        with context.Pool(args.workers) as pool:
            parent_rows = pool.map(compute_parent, inputs)
    parent_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    selected_rows, _, _ = M.residual_two_rows(parent_rows)
    tasks = [
        (row, index)
        for row in selected_rows
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]
    require(len(tasks) == 5_618, "THM2916 child task count changed")
    if args.workers == 1:
        children = list(map(compute_child, tasks))
    else:
        with context.Pool(args.workers) as pool:
            children = pool.map(compute_child, tasks)
    children.sort(
        key=lambda row: (
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        )
    )
    base_closed = [row for row in children if row["closed"]]
    base_failed = [row for row in children if not row["closed"]]
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in children:
        by_body[row["body"]].append(row)
    computed_closed = {
        body
        for body, body_rows in by_body.items()
        if all(row["closed"] for row in body_rows)
    }
    require(
        len(base_closed) == 5_251
        and len(base_failed) == 367
        and computed_closed == base_closed_artifact
        and set(by_body) - computed_closed == base_failed_artifact,
        "THM2916 child/root partition changed",
    )

    if args.workers == 1:
        repaired = list(map(analyze_first_repair, base_failed))
    else:
        with context.Pool(args.workers) as pool:
            repaired = pool.map(analyze_first_repair, base_failed)
    repaired.sort(
        key=lambda row: (
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        )
    )
    deep_parents = [row for row in repaired if row["route"] == "deep"]
    if args.workers == 1:
        deep_rows = list(map(analyze_deep, deep_parents))
    else:
        with context.Pool(args.workers) as pool:
            deep_rows = pool.map(analyze_deep, deep_parents)
    deep_by_key = {
        (parent["body"], parent["rank"], parent["apex"], parent["center"]): row
        for parent, row in zip(deep_parents, deep_rows)
    }
    require(len(deep_by_key) == len(deep_rows), "deep identity collision")

    routes_by_body: dict[tuple[int, ...], list[str]] = defaultdict(list)
    for row in repaired:
        routes_by_body[row["body"]].append(row["route"])
    pair_roots = {
        body
        for body, routes in routes_by_body.items()
        if all(route == "pair" for route in routes)
    }
    hunter_roots = {
        body
        for body, routes in routes_by_body.items()
        if all(route in {"pair", "hunter"} for route in routes)
    }
    toothpick_roots = {
        body
        for body, routes in routes_by_body.items()
        if all(route in {"pair", "hunter", "toothpick"} for route in routes)
    }
    final_roots = set(routes_by_body)
    require(
        final_roots == base_failed_artifact
        and len(final_roots) == 296
        and len(deep_rows) == 5
        and all(row["closed"] for row in deep_rows),
        "two-H3 repair root composition changed",
    )
    baseline745 = M.current_union_351() | base_closed_artifact
    final_union = baseline745 | final_roots
    require(
        len(baseline745) == 745
        and not (baseline745 & final_roots)
        and len(final_union) == 1_041
        and 3_432 - len(final_union) == 2_391,
        "final root union changed",
    )

    second_rows = [
        (parent, second)
        for parent in repaired
        for second in parent["second_rows"]
    ]
    counts = (
        len(parent_rows),
        len(selected_rows),
        len(children),
        len(base_closed),
        len(base_failed),
        sum(row["route"] == "pair" for row in repaired),
        sum(row["pair_margin"] <= 0 for row in repaired),
        sum(row["route"] == "hunter" for row in repaired),
        sum(row["route"] in {"toothpick", "deep"} for row in repaired),
        sum(row["route"] == "toothpick" for row in repaired),
        sum(row["route"] == "deep" for row in repaired),
        len(second_rows),
        sum(second["closed"] for _, second in second_rows),
        sum(not second["closed"] for _, second in second_rows),
        sum(row["route"] == "pair3" for row in deep_rows),
        sum(row["route"] == "hunter3" for row in deep_rows),
        len(pair_roots),
        len(hunter_roots),
        len(toothpick_roots),
        len(final_roots),
        len(baseline745),
        len(final_union),
        3_432 - len(final_union),
        sum(row["pair"]["scanned"] for row in repaired),
        sum(row["pair"]["paid"] for row in repaired),
        sum(row["core_scan"] for row in repaired),
        sum(second["scanned"] for _, second in second_rows),
        sum(row["pair"]["paid"] for row in deep_rows),
        tuple(sorted(Counter(len(row["second_core"]) for row in repaired if row["second_core"]).items())),
        sum(row["pair_margin"] == 0 for row in repaired),
        sum(row["hunter_margin"] == 0 for row in repaired),
        sum(second["margin"] == 0 for _, second in second_rows),
        sum(row["pair_partition_margin"] == 0 for row in deep_rows),
        sum(row["hunter_margin"] == 0 for row in deep_rows),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "THM2920 counts changed")

    closed_roots = tuple(sorted(final_roots))
    final_two_row_roots = tuple(sorted(base_closed_artifact | final_roots))
    closed_digest = hashlib.sha256(repr(closed_roots).encode()).hexdigest()
    final_two_row_digest = hashlib.sha256(
        repr(final_two_row_roots).encode()
    ).hexdigest()
    if EXPECTED_CLOSED_ROOT_DIGEST is not None:
        require(
            closed_digest == EXPECTED_CLOSED_ROOT_DIGEST,
            "THM2920 additive roots changed",
        )
    if EXPECTED_FINAL_TWO_ROW_DIGEST is not None:
        require(
            final_two_row_digest == EXPECTED_FINAL_TWO_ROW_DIGEST,
            "complete two-H3 stratum changed",
        )

    ledger_lines: list[str] = []
    for row in repaired:
        ledger_lines.append(T.child_line(row))
        ledger_lines.extend(
            T.second_line(row, second) for second in row["second_rows"]
        )
        key = (row["body"], row["rank"], row["apex"], row["center"])
        if key in deep_by_key:
            ledger_lines.append(deep_line(row, deep_by_key[key]))
    ledger_hash = hashlib.sha256()
    ledger_hash.update(b"LRC14/j6/two-H3/pair-Hunter-recursive-toothpick/v1\n")
    for line in ledger_lines:
        ledger_hash.update(line.encode())
    ledger_sha256 = ledger_hash.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "THM2920 ledger changed")
    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 two-H3 pair-Hunter recursive-toothpick ledger\n"
            + "".join(ledger_lines)
            + f"ledger_sha256={ledger_sha256}\n"
            + "scope=367 THM2916-open children;296 roots;not LRC14\n"
        )

    pair_positive = [row for row in repaired if row["pair_margin"] > 0]
    hunter_positive = [
        row
        for row in repaired
        if row["pair_margin"] <= 0 and row["hunter_margin"] > 0
    ]
    first_toothpick = [
        row for row in repaired if row["route"] == "toothpick"
    ]
    deep_pair_positive = [
        row for row in deep_rows if row["pair_partition_margin"] > 0
    ]
    deep_hunter_positive = [
        row for row in deep_rows if row["route"] == "hunter3"
    ]
    closest_pair = min(pair_positive, key=lambda row: row["pair_margin"])
    closest_hunter = min(
        hunter_positive,
        key=lambda row: row["hunter_margin"],
    )
    closest_second = min(
        (
            (parent, second)
            for parent in first_toothpick
            for second in parent["second_rows"]
        ),
        key=lambda item: item[1]["margin"],
    )
    closest_pair3 = min(
        deep_pair_positive,
        key=lambda row: row["pair_partition_margin"],
    )
    closest_hunter3 = min(
        deep_hunter_positive,
        key=lambda row: row["hunter_margin"],
    )

    lines = [
        "LRC14 two-H3-row pair-Hunter recursive-toothpick closure",
        f"pair_horizon={T.PAIR_HORIZON};top3_horizon={T.TOP3_HORIZON}",
        f"counts={counts}",
        f"closed_roots={closed_roots}",
        f"complete_two_h3_roots={final_two_row_roots}",
        f"closed_root_digest={closed_digest}",
        f"complete_two_h3_digest={final_two_row_digest}",
        f"closest_pair_positive={boundary_text(closest_pair, 'pair_margin')}",
        f"closest_hunter_positive={boundary_text(closest_hunter, 'hunter_margin')}",
        (
            "closest_second_positive="
            f"{closest_second[0]['body']};rank={closest_second[0]['rank']};"
            f"a={closest_second[0]['apex']};x={closest_second[0]['center']};"
            f"y={closest_second[1]['center']};"
            f"margin={ftext(closest_second[1]['margin'])}"
        ),
        (
            "closest_pair3_positive="
            f"{deep_parents[deep_rows.index(closest_pair3)]['body']};"
            f"x={deep_parents[deep_rows.index(closest_pair3)]['center']};"
            f"y={closest_pair3['second']['center']};"
            f"margin={ftext(closest_pair3['pair_partition_margin'])}"
        ),
        (
            "closest_hunter3_positive="
            f"{deep_parents[deep_rows.index(closest_hunter3)]['body']};"
            f"x={deep_parents[deep_rows.index(closest_hunter3)]['center']};"
            f"y={closest_hunter3['second']['center']};"
            f"margin={ftext(closest_hunter3['hunter_margin'])}"
        ),
        (
            "minimum_pair_tail_gap="
            f"{ftext(min(row['pair']['tail_gap'] for row in repaired))}"
        ),
        (
            "minimum_second_tail_gap="
            f"{ftext(min(second['tail_gap'] for _, second in second_rows))}"
        ),
        (
            "minimum_deep_pair_tail_gap="
            f"{ftext(min(row['pair']['tail_gap'] for row in deep_rows))}"
        ),
        f"deep_records={tuple((parent['body'], parent['rank'], parent['apex'], parent['center'], row['second']['center'], row['route']) for parent, row in zip(deep_parents, deep_rows))}",
        f"ledger_sha256={ledger_sha256}",
        (
            "mode=DISCOVERY"
            if any(
                value is None
                for value in (
                    EXPECTED_COUNTS,
                    EXPECTED_CLOSED_ROOT_DIGEST,
                    EXPECTED_FINAL_TWO_ROW_DIGEST,
                    EXPECTED_LEDGER_SHA256,
                )
            )
            else "mode=LOCKED"
        ),
        (
            "scope=367 THM2916-open ordered children;296 additive roots;"
            "complete 690-root two-H3 stratum;not LRC14"
        ),
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
