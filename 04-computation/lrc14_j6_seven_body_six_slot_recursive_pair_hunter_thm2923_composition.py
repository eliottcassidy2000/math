#!/usr/bin/env python3
"""Compose the recursive closures with the seven-body THM-2915 routes."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_THM2915 = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_open_centre_child_top4_closure_thm2915.py"
)
DEFAULT_CHILD_LEDGER = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_all_open_centre_child_top4_closure_thm2915.ledger.out"
)
THM2913_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out"
)
EXPECTED_THM2915_SOURCE_SHA256 = (
    "9d2e6227a8cbda763fbd73f21dc4d162949e5d5fcd147abd6e8ea37513775215"
)
EXPECTED_THM2915_LEDGER_SHA256 = (
    "798cd660ab60e2021b28074a1390af3f6b1367c99f2d0ab63a581513f7871071"
)
EXPECTED_THM2913_OUT_SHA256 = (
    "3604644a9691b13e7fa245249b68c9586ec2775996834f7761f32eb0b89f3e64"
)
EXPECTED_ROOT_SEMANTIC_SHA256 = (
    "75c4e9f27ed81fb08250b6cf808ed8a0db762031246db9358dc287a081a51aa3"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def ranked_labels(text: str) -> tuple[int, ...]:
    return tuple(int(item.split(":", 1)[0]) for item in text.split(","))


def load(path: Path):
    spec = importlib.util.spec_from_file_location("thm2923_thm2915_composer", path)
    require(spec is not None and spec.loader is not None, "cannot load THM-2915")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fields(line: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in line.rstrip().split(";")[1:])


def parent_key(row: dict[str, object]) -> tuple[object, ...]:
    return row["body"], row["rank"], row["apex"], row["prefix"]


def line_parent_key(data: dict[str, str]) -> tuple[object, ...]:
    return (
        ints(data["E"]),
        int(data["rank"]),
        int(data["a"]),
        ints(data["P"]),
    )


def line_child_key(data: dict[str, str]) -> tuple[object, ...]:
    return (
        ints(data["E"]),
        int(data["rank"]),
        int(data["a"]),
        ints(data["P"]),
        int(data["x"]),
        ints(data["earlier"]),
    )


def line_second_key(data: dict[str, str]) -> tuple[object, ...]:
    return (
        *line_child_key(data),
        int(data["y"]),
        ints(data["yearlier"]),
    )


def digest_bodies(name: str, bodies: set[tuple[int, ...]]) -> str:
    digest = hashlib.sha256((name + "\n").encode())
    for body in sorted(bodies):
        digest.update((",".join(map(str, body)) + "\n").encode())
    return digest.hexdigest()


def root_projection(
    rows_by_body: dict[tuple[int, ...], set[tuple[object, ...]]],
    routes: set[tuple[object, ...]],
) -> set[tuple[int, ...]]:
    return {
        body for body, keys in rows_by_body.items() if keys <= routes
    }


def literal_line(path: Path, prefix: str) -> object:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix}")
    return ast.literal_eval(matches[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--thm2915-source", type=Path, default=DEFAULT_THM2915)
    parser.add_argument("--child-ledger", type=Path, default=DEFAULT_CHILD_LEDGER)
    parser.add_argument("--pair-ledger", type=Path, required=True)
    parser.add_argument("--recursive-ledger", type=Path, required=True)
    parser.add_argument("--grandchild-ledger", type=Path)
    parser.add_argument("--third-ledger", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--root-ledger", type=Path)
    args = parser.parse_args()

    require(
        file_sha256(args.thm2915_source) == EXPECTED_THM2915_SOURCE_SHA256,
        "THM-2915 source changed",
    )
    require(
        file_sha256(args.child_ledger) == EXPECTED_THM2915_LEDGER_SHA256,
        "THM-2915 child ledger changed",
    )
    require(
        file_sha256(THM2913_OUT) == EXPECTED_THM2913_OUT_SHA256,
        "THM-2913 output changed",
    )
    X = load(args.thm2915_source)
    parents = X.joined_parent_rows()
    require(len(parents) == 11_842, "parent universe changed")
    by_key = {parent_key(row): row for row in parents}
    require(len(by_key) == len(parents), "parent key collision")
    branch_closed_keys = {
        parent_key(row) for row in parents if row["branch_closed"]
    }
    require(len(branch_closed_keys) == 279, "literal parent branch changed")
    rows_by_body: dict[
        tuple[int, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    for row in parents:
        rows_by_body[row["body"]].add(parent_key(row))

    child_groups: dict[
        tuple[object, ...], list[tuple[tuple[object, ...], bool]]
    ] = defaultdict(list)
    failure_by_body: dict[
        tuple[int, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    all_failures: set[tuple[object, ...]] = set()
    all_children: set[tuple[object, ...]] = set()
    for line in args.child_ledger.read_text().splitlines():
        if not line.startswith("CHILD;"):
            continue
        data = fields(line)
        pkey = line_parent_key(data)
        ckey = line_child_key(data)
        require(pkey in by_key, "child parent is absent from the 11,842-row universe")
        require(pkey not in branch_closed_keys, "child emitted for literal-closed parent")
        require(ckey not in all_children, "duplicate THM-2915 child identity")
        all_children.add(ckey)
        margin = F(data["margin"])
        closed = margin > 0
        require(
            (data["closed"] == "1") == closed,
            f"{ckey}: THM-2915 closed bit disagrees with exact margin",
        )
        child_groups[pkey].append((ckey, closed))
        if not closed:
            all_failures.add(ckey)
            failure_by_body[ckey[0]].add(ckey)
    require(
        len(child_groups) == 11_563
        and sum(map(len, child_groups.values())) == 51_222
        and len(all_children) == 51_222
        and len(all_failures) == 4_866,
        "THM-2915 child universe changed",
    )
    require(
        set(child_groups) == set(by_key) - branch_closed_keys,
        "THM-2915 child groups do not exactly join the non-literal parent rows",
    )
    T0 = {
        key for key, group in child_groups.items() if all(closed for _, closed in group)
    }
    require(len(T0) == 8_112, "THM-2915 top-four T route changed")

    pair_closed: set[tuple[object, ...]] = set()
    recursive_ids: set[tuple[object, ...]] = set()
    pair_rows_seen: set[tuple[object, ...]] = set()
    for line in args.pair_ledger.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        data = fields(line)
        key = line_child_key(data)
        require(key in all_failures, "pair row not in failure universe")
        require(key not in pair_rows_seen, "duplicate pair-ledger child identity")
        pair_rows_seen.add(key)
        pair_margin = F(data["pairmargin"])
        hunter_margin = F(data["G4margin"])
        route = (
            "pair"
            if pair_margin > 0
            else "hunter"
            if hunter_margin > 0
            else "recursive"
        )
        require(data["route"] == route, f"{key}: pair route disagrees with margins")
        if route == "recursive":
            recursive_ids.add(key)
        else:
            pair_closed.add(key)
    require(
        len(pair_closed) + len(recursive_ids) == len(all_failures)
        and pair_rows_seen == all_failures
        and not (pair_closed & recursive_ids),
        "pair ledger is not a partition",
    )

    recursive_rows_seen: set[tuple[object, ...]] = set()
    second_rows_seen: set[tuple[object, ...]] = set()
    recursive_core_by_first: dict[
        tuple[object, ...], tuple[int, ...]
    ] = {}
    recursive_claim_by_first: dict[tuple[object, ...], bool] = {}
    second_ids_by_first_all: dict[
        tuple[object, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    second_closed_by_first: dict[
        tuple[object, ...], list[bool]
    ] = defaultdict(list)
    failed_seconds_by_first: dict[
        tuple[object, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    for line in args.recursive_ledger.read_text().splitlines():
        if line.startswith("RECURSIVE;"):
            data = fields(line)
            key = line_child_key(data)
            require(key in recursive_ids, "recursive row not in pair-Hunter residual")
            require(key not in recursive_rows_seen, "duplicate recursive identity")
            recursive_rows_seen.add(key)
            recursive_core_by_first[key] = ranked_labels(data["H1"])
            recursive_claim_by_first[key] = data["closed"] == "1"
        elif line.startswith("SECOND;"):
            data = fields(line)
            first_key = line_child_key(data)
            require(first_key in recursive_ids, "second pivot escaped recursive universe")
            second_key = line_second_key(data)
            require(second_key not in second_rows_seen, "duplicate second-pivot identity")
            second_rows_seen.add(second_key)
            second_ids_by_first_all[first_key].add(second_key)
            derived_closed = F(data["margin"]) > 0
            require(
                (data["closed"] == "1") == derived_closed,
                f"{second_key}: SECOND closed bit disagrees with exact margin",
            )
            second_closed_by_first[first_key].append(derived_closed)
            if not derived_closed:
                failed_seconds_by_first[first_key].add(second_key)
    require(
        recursive_rows_seen == recursive_ids
        and len(second_rows_seen) == 6_172
        and len(set().union(*failed_seconds_by_first.values())) == 228,
        "recursive ledger is not a partition",
    )
    for first, core in recursive_core_by_first.items():
        expected_seconds = {
            (*first, center, tuple(core[:index]))
            for index, center in enumerate(core)
        }
        require(
            second_ids_by_first_all[first] == expected_seconds,
            f"{first}: SECOND pivots do not exactly realize the ordered H1 core",
        )
        require(
            recursive_claim_by_first[first]
            == all(second_closed_by_first[first]),
            f"{first}: recursive closed bit disagrees with SECOND pivots",
        )
    recursive_closed = {
        first
        for first in recursive_ids
        if all(second_closed_by_first[first])
    }
    recursive_open = recursive_ids - recursive_closed
    require(
        len(recursive_closed) == 1_884 and len(recursive_open) == 195,
        "derived recursive partition changed",
    )
    deep_closed_first: set[tuple[object, ...]] = set()
    remaining_second: set[tuple[object, ...]] = set().union(
        *failed_seconds_by_first.values()
    )
    if args.grandchild_ledger is not None or args.third_ledger is not None:
        require(
            args.grandchild_ledger is not None and args.third_ledger is not None,
            "both deeper ledgers are required",
        )
        upgraded_second: set[tuple[object, ...]] = set()
        third_recursive_second: set[tuple[object, ...]] = set()
        grandchild_rows_seen: set[tuple[object, ...]] = set()
        for line in args.grandchild_ledger.read_text().splitlines():
            if not line.startswith("GRANDCHILD;"):
                continue
            data = fields(line)
            key = line_second_key(data)
            require(key in remaining_second, "grandchild not in failed-second universe")
            require(key not in grandchild_rows_seen, "duplicate grandchild identity")
            grandchild_rows_seen.add(key)
            partition_margin = F(data["partitionmargin"])
            hunter3_margin = F(data["G3margin"])
            route3 = (
                "pair_single"
                if partition_margin > 0
                else "hunter3"
                if hunter3_margin > 0
                else "third_recursive"
            )
            require(
                data["route3"] == route3,
                f"{key}: grandchild route disagrees with exact margins",
            )
            if route3 == "third_recursive":
                third_recursive_second.add(key)
            else:
                upgraded_second.add(key)
        require(
            grandchild_rows_seen == remaining_second,
            "grandchild ledger does not exactly cover the 228 failed pivots",
        )
        third_rows_seen: set[tuple[object, ...]] = set()
        third_core_by_second: dict[
            tuple[object, ...], tuple[int, ...]
        ] = {}
        third_claim_by_second: dict[tuple[object, ...], bool] = {}
        third_pivots_by_second: dict[
            tuple[object, ...], set[tuple[object, ...]]
        ] = defaultdict(set)
        third_pivot_closed_by_second: dict[
            tuple[object, ...], list[bool]
        ] = defaultdict(list)
        for line in args.third_ledger.read_text().splitlines():
            if line.startswith("THIRD_RECURSIVE;"):
                data = fields(line)
                key = line_second_key(data)
                require(key in third_recursive_second, "third row not in recursive universe")
                require(key not in third_rows_seen, "duplicate third-recursive identity")
                third_rows_seen.add(key)
                third_core_by_second[key] = ranked_labels(data["H1"])
                third_claim_by_second[key] = data["closed"] == "1"
            elif line.startswith("THIRD;"):
                data = fields(line)
                key = line_second_key(data)
                require(key in third_recursive_second, "third pivot escaped recursive universe")
                pivot_key = (
                    *key,
                    int(data["z"]),
                    ints(data["zearlier"]),
                )
                require(
                    pivot_key not in third_pivots_by_second[key],
                    "duplicate third-pivot identity",
                )
                third_pivots_by_second[key].add(pivot_key)
                derived_closed = F(data["margin"]) > 0
                require(
                    (data["closed"] == "1") == derived_closed,
                    f"{pivot_key}: THIRD closed bit disagrees with exact margin",
                )
                third_pivot_closed_by_second[key].append(derived_closed)
        require(third_rows_seen == third_recursive_second, "third recursion is incomplete")
        require(
            sum(map(len, third_pivots_by_second.values())) == 117,
            "third-pivot count changed",
        )
        for second, core in third_core_by_second.items():
            expected_thirds = {
                (*second, center, tuple(core[:index]))
                for index, center in enumerate(core)
            }
            require(
                third_pivots_by_second[second] == expected_thirds,
                f"{second}: THIRD pivots do not exactly realize the ordered H1 core",
            )
            require(
                third_claim_by_second[second]
                == all(third_pivot_closed_by_second[second]),
                f"{second}: third-recursive closed bit disagrees with THIRD pivots",
            )
        third_closed_second = {
            second
            for second in third_recursive_second
            if all(third_pivot_closed_by_second[second])
        }
        third_open_second = third_recursive_second - third_closed_second
        require(
            len(third_closed_second) == 48 and len(third_open_second) == 1,
            "derived third-recursive partition changed",
        )
        upgraded_second |= third_closed_second
        remaining_second -= upgraded_second
        deep_closed_first = {
            first
            for first, seconds in failed_seconds_by_first.items()
            if not (seconds & remaining_second)
        }
        require(
            deep_closed_first <= recursive_open,
            "deep closure escaped first-level recursive residual",
        )
    all_closed = pair_closed | recursive_closed | deep_closed_first

    completed_failed_parents = {
        key
        for key, group in child_groups.items()
        if any(not closed for _, closed in group)
        and {
            child for child, closed in group if not closed
        }
        <= all_closed
    }
    T2 = T0 | completed_failed_parents
    require(
        len(T2) == len(T0) + len(completed_failed_parents),
        "T0/new parent overlap",
    )

    P = branch_closed_keys
    E = {parent_key(row) for row in parents if row["pair_exception"]}
    H = X.finite_h1_keys(parents)
    require(
        (len(P), len(E), len(H)) == (279, 52, 3_090),
        "P/E/H route sizes changed",
    )

    q_bodies = set(literal_line(THM2913_OUT, "closed_roots="))
    require(len(q_bodies) == 38, "THM-2913 repaired root bank changed")
    Q: set[tuple[object, ...]] = set()
    for body in q_bodies:
        candidates = [
            parent_key(row)
            for row in parents
            if row["body"] == body
            and not row["pair_exception"]
            and not row["branch_closed"]
        ]
        require(len(candidates) == 1, f"{body}: Q parent is not unique")
        Q.add(candidates[0])
    require(len(Q) == 38, "Q key count changed")

    if hasattr(X, "proved_unions"):
        _, _, baseline314, source_q_bodies, baseline351 = X.proved_unions()
        require(source_q_bodies == q_bodies, "THM-2913 source/result Q bank changed")
    else:
        _, _, baseline314 = X.proved_union_through_2912()
        baseline351 = baseline314 | q_bodies
    require(
        len(baseline314) == 314
        and len(q_bodies & baseline314) == 1
        and len(baseline351) == 351,
        "canonical union through THM-2913 changed",
    )

    route_sets = {
        "PT2": P | T2,
        "EPT2": E | P | T2,
        "HPT2": H | P | T2,
        "EHPT2": E | H | P | T2,
        "QEHPT2": Q | E | H | P | T2,
    }
    root_sets = {
        name: root_projection(rows_by_body, routes)
        for name, routes in route_sets.items()
    }
    route_summaries: dict[str, tuple[int, ...]] = {}
    for name, roots in root_sets.items():
        union = baseline351 | roots
        route_summaries[name] = (
            len(route_sets[name]),
            len(roots),
            len(roots & baseline351),
            len(roots - baseline351),
            len(union),
            3_432 - len(union),
        )

    branch_atoms = Counter(
        "".join(
            label
            for label, route in (
                ("E", E),
                ("H", H),
                ("P", P),
                ("Q", Q),
                ("T", T2),
            )
            if parent_key(row) in route
        )
        or "-"
        for row in parents
    )
    newly_completed_roots = {
        body
        for body, failures in failure_by_body.items()
        if failures <= all_closed
    }
    remaining_first = all_failures - all_closed
    remaining_children_by_parent = Counter(
        sum(
            child in remaining_first
            for child, closed in group
            if not closed
        )
        for key, group in child_groups.items()
        if any(child in remaining_first for child, _ in group)
    )
    remaining_children_by_root = Counter(
        len(failures & remaining_first)
        for failures in failure_by_body.values()
        if failures & remaining_first
    )

    if args.grandchild_ledger is not None:
        endpoint_body = (2, 4, 6, 8, 10, 12, 14)
        endpoint_parent = (endpoint_body, 1, 22, (22,))
        endpoint_first = (*endpoint_parent, 18, ())
        endpoint_second = (*endpoint_first, 26, ())
        require(
            remaining_first == {endpoint_first}
            and remaining_second == {endpoint_second},
            "deep residual is not the unique doubled-AP endpoint branch",
        )
        require(endpoint_parent in E - T2, "endpoint parent is not E-only")
        require(
            route_summaries
            == {
                "PT2": (11_841, 3_410, 330, 3_080, 3_431, 1),
                "EPT2": (11_842, 3_411, 330, 3_081, 3_432, 0),
                "HPT2": (11_841, 3_410, 330, 3_080, 3_431, 1),
                "EHPT2": (11_842, 3_411, 330, 3_081, 3_432, 0),
                "QEHPT2": (11_842, 3_411, 330, 3_081, 3_432, 0),
            },
            "deep route composition changed",
        )
        require(
            remaining_children_by_parent == Counter({1: 1})
            and remaining_children_by_root == Counter({1: 1}),
            "deep residual histogram changed",
        )
        require(
            branch_atoms
            == Counter(
                {
                    "E": 1,
                    "ET": 51,
                    "HP": 215,
                    "HQT": 1,
                    "HT": 2_874,
                    "P": 64,
                    "QT": 37,
                    "T": 8_599,
                }
            ),
            "deep branch atoms changed",
        )

    root_lines: list[str] = []
    for name, roots in root_sets.items():
        for body in sorted(roots):
            root_lines.append(f"{name}_ROOT={','.join(map(str, body))}")
    root_semantic = hashlib.sha256(
        (
            "LRC14/THM2915-failure-recursion/route-roots/v1\n"
            + "\n".join(root_lines)
            + "\n"
        ).encode()
    ).hexdigest()
    if args.grandchild_ledger is not None:
        require(
            root_semantic == EXPECTED_ROOT_SEMANTIC_SHA256,
            "deep route-root semantic digest changed",
        )
    if args.root_ledger is not None:
        args.root_ledger.write_text(
            "LRC14 THM2915 failure-recursion route roots\n"
            + "\n".join(root_lines)
            + "\n"
            + f"semantic_sha256={root_semantic}\n"
            + "scope=complete seven-body/six-external route composition through "
            + "THM2913;not unrestricted LRC14\n",
            encoding="utf-8",
            newline="\n",
        )

    summary: dict[str, object] = {
        "children": 51_222,
        "top4_closed_children": 46_356,
        "top4_failed_children": len(all_failures),
        "pair_or_hunter_closed": len(pair_closed),
        "recursive_closed": len(recursive_closed),
        "recursive_open": len(recursive_open),
        "deep_closed_first": len(deep_closed_first),
        "deep_remaining_second": len(remaining_second),
        "T0_rows": len(T0),
        "new_T2_rows": len(completed_failed_parents),
        "T2_rows": len(T2),
        "remaining_open_rows": 11_563 - len(T2),
        "newly_completed_roots": len(newly_completed_roots),
        "newly_completed_overlap351": len(newly_completed_roots & baseline351),
        "newly_completed_additive351": len(newly_completed_roots - baseline351),
        "remaining_failed_roots": len(failure_by_body) - len(newly_completed_roots),
        "Q_keys": len(Q),
        "Q_intersection_T0": len(Q & T0),
        "Q_intersection_T2": len(Q & T2),
        "Q_not_T2": len(Q - T2),
        "E_intersection_T2": len(E & T2),
        "H_intersection_T2": len(H & T2),
        "P_intersection_T2": len(P & T2),
        "routes": route_summaries,
        "root_digests": {
            name: digest_bodies(name, roots)
            for name, roots in root_sets.items()
        },
        "branch_atoms": tuple(sorted(branch_atoms.items())),
        "dependency_graph": (
            "THM2904-parent-hostile-centres",
            "THM2915-child-top4",
            "exact-B2-or-G4",
            "ordered-second-centre-top3",
            "exact-B2-plus-q1-or-G3",
            "ordered-third-centre-B2",
            "THM2907-endpoint-E",
            "THM2920-positive-control",
            "THM2913-baseline-composition",
        ),
        "remaining_children_by_parent": tuple(
            sorted(remaining_children_by_parent.items())
        ),
        "remaining_children_by_root": tuple(
            sorted(remaining_children_by_root.items())
        ),
        "root_semantic_sha256": root_semantic,
    }
    rendered = "\n".join(f"{key}={value}" for key, value in summary.items()) + "\n"
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
