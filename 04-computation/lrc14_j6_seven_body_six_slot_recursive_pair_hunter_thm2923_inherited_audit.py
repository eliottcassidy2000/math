#!/usr/bin/env python3
"""Exact inherited-slice audit against THM-2913 and THM-2920.

The two canonical computations are independent implementations of the same
marked 4->3 Hunter recursion on narrow strata.  This checker joins their full
branch identities to the all-centre ledgers, compares exact values after
normalizing harmless tied-rank labels, and verifies that THM-2920's 296 roots
are contained in the all-centre deep route.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2913_LEDGER = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.ledger.out"
)
THM2920_LEDGER = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.ledger.out"
)
THM2920_OUTPUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.out"
)
EXPECTED_SHA256 = {
    THM2913_LEDGER: "5419b87511bbf51c43a2bed9647e82cb5178ad99ac8f667d0e66318caa049632",
    THM2920_LEDGER: "9c96a24c90c07c69d96b86f29355f9a86599c9d9174e0260dd952aa357f7d7f1",
    THM2920_OUTPUT: "1a38fd441dfd77a4f5d30d45d3160febc33d2d4eeb6247b223f10a1e31a8aefb",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def fields(line: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in line.rstrip().split(";")[1:])


def ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def first_key(data: dict[str, str]) -> tuple[object, ...]:
    return (
        ints(data["E"]),
        int(data["rank"]),
        int(data["a"]),
        ints(data["P"]),
        int(data["x"]),
        ints(data["earlier"]),
    )


def all_second_key(data: dict[str, str]) -> tuple[object, ...]:
    return (
        *first_key(data),
        int(data["y"]),
        ints(data["yearlier"]),
    )


def ranked_values(text: str) -> tuple[F, ...]:
    return tuple(sorted((F(item.split(":", 1)[1]) for item in text.split(",")), reverse=True))


def canonical_children(path: Path) -> dict[tuple[object, ...], dict[str, str]]:
    rows: dict[tuple[object, ...], dict[str, str]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("CHILD;"):
            continue
        data = fields(line)
        key = first_key(data)
        require(key not in rows, f"{path.name}: duplicate CHILD identity")
        rows[key] = data
    return rows


def all_children(path: Path) -> dict[tuple[object, ...], dict[str, str]]:
    rows: dict[tuple[object, ...], dict[str, str]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        data = fields(line)
        key = first_key(data)
        require(key not in rows, "duplicate all-centre PAIR identity")
        rows[key] = data
    require(len(rows) == 4_866, "all-centre PAIR universe changed")
    return rows


def canonical_second_rows(
    path: Path,
    children: dict[tuple[object, ...], dict[str, str]],
) -> tuple[
    dict[tuple[object, ...], dict[str, str]],
    dict[tuple[object, ...], tuple[object, ...]],
]:
    simple_first: dict[tuple[object, ...], tuple[object, ...]] = {}
    for key in children:
        simple = (key[0], key[1], key[2], key[4])
        require(simple not in simple_first, f"{path.name}: ambiguous simple CHILD key")
        simple_first[simple] = key
    rows: dict[tuple[object, ...], dict[str, str]] = {}
    simple_second: dict[tuple[object, ...], tuple[object, ...]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("SECOND;"):
            continue
        data = fields(line)
        simple = (ints(data["E"]), int(data["rank"]), int(data["a"]), int(data["x"]))
        require(simple in simple_first, f"{path.name}: SECOND parent missing")
        first = simple_first[simple]
        key = (*first, int(data["y"]), ints(data["earlier"]))
        require(key not in rows, f"{path.name}: duplicate SECOND identity")
        rows[key] = data
        short = (*simple, int(data["y"]))
        require(short not in simple_second, f"{path.name}: ambiguous simple SECOND key")
        simple_second[short] = key
    return rows, simple_second


def all_second_rows(path: Path) -> dict[tuple[object, ...], dict[str, str]]:
    rows: dict[tuple[object, ...], dict[str, str]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("SECOND;"):
            continue
        data = fields(line)
        key = all_second_key(data)
        require(key not in rows, "duplicate all-centre SECOND identity")
        rows[key] = data
    require(len(rows) == 6_172, "all-centre SECOND universe changed")
    return rows


def audit_slice(
    name: str,
    ledger: Path,
    expected_children: int,
    expected_seconds: int,
    all_first: dict[tuple[object, ...], dict[str, str]],
    all_second: dict[tuple[object, ...], dict[str, str]],
) -> tuple[
    dict[tuple[object, ...], dict[str, str]],
    dict[tuple[object, ...], tuple[object, ...]],
    Counter[str],
]:
    children = canonical_children(ledger)
    seconds, simple_seconds = canonical_second_rows(ledger, children)
    require(
        len(children) == expected_children and len(seconds) == expected_seconds,
        f"{name}: inherited universe changed",
    )
    route_map = {
        "pair": "pair",
        "hunter": "hunter",
        "toothpick": "recursive",
        "deep": "recursive",
    }
    routes: Counter[str] = Counter()
    for key, canonical in children.items():
        require(key in all_first, f"{name}: CHILD identity absent from all-centre ledger")
        broad = all_first[key]
        route = route_map[canonical["route"]]
        routes[route] += 1
        require(broad["route"] == route, f"{name}: CHILD route mismatch")
        require(
            (
                F(canonical["h"]),
                int(canonical["r"]),
                ranked_values(canonical["top4"]),
                F(canonical["B2"]),
                F(canonical["pairmargin"]),
                F(canonical["G4"]),
                F(canonical["G4margin"]),
                canonical["lambda"],
            )
            == (
                F(broad["h"]),
                int(broad["r"]),
                ranked_values(broad["top4"]),
                F(broad["B2"]),
                F(broad["pairmargin"]),
                F(broad["G4"]),
                F(broad["G4margin"]),
                broad["lambda"],
            ),
            f"{name}: CHILD exact data mismatch",
        )
    for key, canonical in seconds.items():
        require(key in all_second, f"{name}: SECOND identity absent from all-centre ledger")
        broad = all_second[key]
        require(
            (
                F(canonical["h"]),
                int(canonical["r"]),
                ranked_values(canonical["top3"]),
                F(canonical["margin"]),
                canonical["closed"],
            )
            == (
                F(broad["h"]),
                int(broad["r"]),
                ranked_values(broad["top3"]),
                F(broad["margin"]),
                broad["closed"],
            ),
            f"{name}: SECOND exact data mismatch",
        )
    return children, simple_seconds, routes


def literal_line(path: Path, prefix: str) -> object:
    values = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(values) == 1, f"{path.name}: expected one {prefix}")
    return ast.literal_eval(values[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-ledger", type=Path, required=True)
    parser.add_argument("--recursive-ledger", type=Path, required=True)
    parser.add_argument("--grandchild-ledger", type=Path, required=True)
    parser.add_argument("--root-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    for path, expected in EXPECTED_SHA256.items():
        require(file_sha256(path) == expected, f"{path.name}: canonical input changed")

    broad_children = all_children(args.pair_ledger)
    broad_seconds = all_second_rows(args.recursive_ledger)
    children2913, _, routes2913 = audit_slice(
        "THM2913",
        THM2913_LEDGER,
        42,
        29,
        broad_children,
        broad_seconds,
    )
    children2920, simple_second2920, routes2920 = audit_slice(
        "THM2920",
        THM2920_LEDGER,
        367,
        283,
        broad_children,
        broad_seconds,
    )
    require(
        routes2920 == Counter({"pair": 149, "hunter": 111, "recursive": 107}),
        "THM2920 first-stage classification changed",
    )

    broad_grandchildren: dict[tuple[object, ...], dict[str, str]] = {}
    for line in args.grandchild_ledger.read_text().splitlines():
        if not line.startswith("GRANDCHILD;"):
            continue
        data = fields(line)
        key = all_second_key(data)
        require(key not in broad_grandchildren, "duplicate all-centre GRANDCHILD identity")
        broad_grandchildren[key] = data
    require(len(broad_grandchildren) == 228, "all-centre GRANDCHILD universe changed")

    deep_rows: list[tuple[tuple[object, ...], dict[str, str]]] = []
    for line in THM2920_LEDGER.read_text().splitlines():
        if not line.startswith("DEEP;"):
            continue
        data = fields(line)
        short = (
            ints(data["E"]),
            int(data["rank"]),
            int(data["a"]),
            int(data["x"]),
            int(data["y"]),
        )
        require(short in simple_second2920, "THM2920 DEEP parent missing")
        deep_rows.append((simple_second2920[short], data))
    require(len(deep_rows) == 5, "THM2920 deep-row count changed")
    deep_routes: Counter[str] = Counter()
    for key, canonical in deep_rows:
        require(key in broad_grandchildren, "THM2920 DEEP row absent from broad ledger")
        broad = broad_grandchildren[key]
        expected_route = {
            "pair3": "pair_single",
            "hunter3": "hunter3",
        }[canonical["route"]]
        deep_routes[expected_route] += 1
        require(
            broad["route3"] == expected_route
            and F(canonical["h"]) == F(broad["h"])
            and int(canonical["r"]) == int(broad["r"])
            and F(canonical["B2"]) == F(broad["B2"])
            and F(canonical["pair3margin"]) == F(broad["partitionmargin"])
            and F(canonical["G3"]) == F(broad["G3"])
            and F(canonical["G3margin"]) == F(broad["G3margin"]),
            "THM2920 deep exact data mismatch",
        )
    require(
        deep_routes == Counter({"pair_single": 3, "hunter3": 2}),
        "THM2920 deep classification changed",
    )

    roots2920 = set(literal_line(THM2920_OUTPUT, "closed_roots="))
    pt2_roots = {
        ints(line.removeprefix("PT2_ROOT="))
        for line in args.root_ledger.read_text().splitlines()
        if line.startswith("PT2_ROOT=")
    }
    require(
        len(roots2920) == 296
        and len(pt2_roots) == 3_410
        and roots2920 <= pt2_roots,
        "THM2920 roots are not contained in the all-centre deep route",
    )

    lines = [
        f"thm2913=(children={len(children2913)},seconds=29,routes={tuple(sorted(routes2913.items()))})",
        f"thm2920=(children={len(children2920)},seconds=283,routes={tuple(sorted(routes2920.items()))})",
        f"thm2920_deep={tuple(sorted(deep_routes.items()))}",
        f"thm2920_roots=(slice=296,pt2=3410,contained={len(roots2920 & pt2_roots)})",
        "normalized_tie_policy=coverage-values-and-margins",
    ]
    digest = hashlib.sha256(
        (
            "LRC14/THM2915-failure-recursion/inherited-slices/v1\n"
            + "\n".join(lines)
            + "\n"
        ).encode()
    ).hexdigest()
    rendered = (
        "LRC14 THM2915-failure recursion inherited-slice audit\n"
        + "\n".join(lines)
        + "\n"
        + f"semantic_sha256={digest}\n"
        + "mode=LOCKED\n"
        + "all_exact_controls=PASS\n"
    )
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
