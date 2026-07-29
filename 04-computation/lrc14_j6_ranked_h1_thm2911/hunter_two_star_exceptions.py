#!/usr/bin/env python3
"""Exact relative-Hunter audit of the two K>=20 star exceptions.

The ordered H1 pivot closes every five-set except sets assigned to one
explicit first pivot on each of two branches.  For those sets, enumerate
only the combinations whose exact pivot-star invoice is at least the
carrier mass, then replace the star by a maximum-weight spanning tree.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCOUT = HERE / "scout.py"
SCOUT_SHA = "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
TARGETS = (
    ((1, 2, 5, 7, 8, 9, 11), 5, 25),
    ((1, 2, 5, 7, 9, 11, 12), 3, 32),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_scout():
    require(
        hashlib.sha256(SCOUT.read_bytes()).hexdigest() == SCOUT_SHA,
        "scout source changed",
    )
    spec = importlib.util.spec_from_file_location("h1_scout", SCOUT)
    require(spec is not None and spec.loader is not None, "cannot load scout")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_scout()


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def maximum_spanning_tree(
    vertices: tuple[int, ...],
    pair_weight,
) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex: int) -> int:
        root = vertex
        while parent[root] != root:
            root = parent[root]
        while parent[vertex] != vertex:
            nxt = parent[vertex]
            parent[vertex] = root
            vertex = nxt
        return root

    edges = sorted(
        (
            (pair_weight(x, y), min(x, y), max(x, y))
            for index, x in enumerate(vertices)
            for y in vertices[index + 1 :]
        ),
        key=lambda row: (-row[0], row[1], row[2]),
    )
    tree = []
    total = F(0)
    for weight, x, y in edges:
        root_x = find(x)
        root_y = find(y)
        if root_x == root_y:
            continue
        parent[root_y] = root_x
        total += weight
        tree.append((x, y, weight))
        if len(tree) == len(vertices) - 1:
            break
    require(len(tree) == len(vertices) - 1, "MST construction failed")
    return total, tuple(tree)


def audit_target(
    body: tuple[int, ...],
    rank: int,
    expected_pivot: int,
) -> dict[str, object]:
    root = S.profile_root_task((body, 15_000, 0, True))
    row = next(record for record in root["rows"] if record["rank"] == rank)
    require(row["cutoff"] <= 15_000, "target left finite window")
    require(not row["depth1_closed"], "target is no longer a star exception")
    unresolved = row["depth1_unresolved"]
    require(
        len(unresolved) == 1 and unresolved[0][1] == expected_pivot,
        "canonical unresolved pivot changed",
    )

    carrier = S.S.R.subtract_local(
        S.S.G.T.CORE.good_norm(body)[0],
        row["apex"],
    )
    require(carrier == S.S.R.subtract_local(
        S.S.G.T.CORE.good_norm(body)[0],
        row["apex"],
    ), "literal carrier reconstruction changed")
    h = S.residual_mass(carrier)
    rows = tuple(row["core_rows"])
    coverage = {speed: value for value, speed in rows}
    pivot_index = next(
        index for index, (_, speed) in enumerate(rows)
        if speed == expected_pivot
    )
    require(pivot_index == 0, "exceptional pivot is no longer first")
    pivot = expected_pivot
    pivot_residual = S.S.R.subtract_local(carrier, pivot)
    pivot_mass = S.residual_mass(pivot_residual)

    pair_cache: dict[tuple[int, int], F] = {}
    residual_cache: dict[int, tuple[list[tuple[F, F]], F]] = {}

    def pair_weight(x: int, y: int) -> F:
        key = (min(x, y), max(x, y))
        if key in pair_cache:
            return pair_cache[key]
        if x not in residual_cache:
            residual_x = S.S.R.subtract_local(carrier, x)
            residual_cache[x] = (
                residual_x,
                S.residual_mass(residual_x),
            )
        residual_x, mass_x = residual_cache[x]
        residual_xy = S.S.R.subtract_local(residual_x, y)
        child_gain = mass_x - S.residual_mass(residual_xy)
        intersection = coverage[y] - child_gain
        require(
            F(0) <= intersection <= min(coverage[x], coverage[y]),
            f"bad pair intersection {x},{y}",
        )
        pair_cache[key] = intersection
        return intersection

    suffix = []
    for _, speed in rows[pivot_index + 1 :]:
        child = S.S.R.subtract_local(pivot_residual, speed)
        child_coverage = pivot_mass - S.residual_mass(child)
        require(
            child_coverage == coverage[speed] - pair_weight(pivot, speed),
            "star inclusion-exclusion identity failed",
        )
        suffix.append((child_coverage, speed))
    suffix.sort(key=lambda item: (-item[0], item[1]))

    counts = {
        "nodes": 0,
        "star_bound_prunes": 0,
        "star_safe_leaves": 0,
        "star_hostile_sets": 0,
        "hunter_repairs": 0,
        "hunter_hard_sets": 0,
    }
    digest = hashlib.sha256(b"LRC14/j6/H1-relative-Hunter/v1\n")
    max_hunter: F | None = None
    max_hunter_row = None
    min_extra_credit: F | None = None

    def visit(
        start: int,
        needed: int,
        chosen: tuple[int, ...],
        child_sum: F,
    ) -> None:
        nonlocal max_hunter, max_hunter_row, min_extra_credit
        counts["nodes"] += 1
        if needed == 0:
            vertices = tuple(sorted((pivot, *chosen)))
            singleton_sum = sum((coverage[x] for x in vertices), F(0))
            star_credit = sum(
                (pair_weight(pivot, y) for y in chosen),
                F(0),
            )
            star_invoice = singleton_sum - star_credit
            require(
                star_invoice == coverage[pivot] + child_sum,
                "star invoice mismatch",
            )
            if star_invoice < h:
                counts["star_safe_leaves"] += 1
                return
            counts["star_hostile_sets"] += 1
            tree_credit, tree = maximum_spanning_tree(
                vertices,
                pair_weight,
            )
            require(tree_credit >= star_credit, "MST is below pivot star")
            hunter_invoice = singleton_sum - tree_credit
            extra_credit = tree_credit - star_credit
            if min_extra_credit is None or extra_credit < min_extra_credit:
                min_extra_credit = extra_credit
            if hunter_invoice < h:
                counts["hunter_repairs"] += 1
            else:
                counts["hunter_hard_sets"] += 1
            candidate = (
                hunter_invoice,
                vertices,
                star_invoice,
                tree_credit,
                tree,
            )
            if (
                max_hunter_row is None
                or (hunter_invoice, vertices)
                > (max_hunter_row[0], max_hunter_row[1])
            ):
                max_hunter = hunter_invoice
                max_hunter_row = candidate
            digest.update(
                (
                    f"Q={vertices};star={ftext(star_invoice)};"
                    f"tree={ftext(tree_credit)};"
                    f"psi={ftext(hunter_invoice)}\n"
                ).encode()
            )
            return
        if len(suffix) - start < needed:
            return
        upper = child_sum + sum(
            (value for value, _ in suffix[start : start + needed]),
            F(0),
        )
        if coverage[pivot] + upper < h:
            counts["star_bound_prunes"] += 1
            return
        last = len(suffix) - needed
        for index in range(start, last + 1):
            value, speed = suffix[index]
            visit(
                index + 1,
                needed - 1,
                (*chosen, speed),
                child_sum + value,
            )

    visit(0, 4, (), F(0))
    require(counts["star_hostile_sets"] > 0, "exception has no hostile set")
    require(
        counts["hunter_repairs"] + counts["hunter_hard_sets"]
        == counts["star_hostile_sets"],
        "Hunter classification lost a set",
    )
    require(max_hunter is not None and max_hunter_row is not None, "no maximum")
    return {
        "body": body,
        "rank": rank,
        "apex": row["apex"],
        "prefix": row["prefix"],
        "h": h,
        "H": len(rows),
        "pivot": pivot,
        **counts,
        "pair_evaluations": len(pair_cache),
        "max_hunter": max_hunter,
        "margin": h - max_hunter,
        "max_row": max_hunter_row,
        "min_extra_credit": min_extra_credit,
        "ledger": digest.hexdigest(),
    }


def main() -> None:
    results = [audit_target(*target) for target in TARGETS]
    print("K>=20 ordered-H1 relative-Hunter exception audit")
    for row in results:
        print(
            f"E={row['body']};rank={row['rank']};a={row['apex']};"
            f"P={row['prefix']};h={ftext(row['h'])};H={row['H']};"
            f"pivot={row['pivot']};nodes={row['nodes']};"
            f"star_prunes={row['star_bound_prunes']};"
            f"star_safe_leaves={row['star_safe_leaves']};"
            f"star_hostile={row['star_hostile_sets']};"
            f"hunter_repairs={row['hunter_repairs']};"
            f"hunter_hard={row['hunter_hard_sets']};"
            f"pairs={row['pair_evaluations']};"
            f"maxPsi={ftext(row['max_hunter'])};"
            f"margin={ftext(row['margin'])};"
            f"min_extra_credit={ftext(row['min_extra_credit'])};"
            f"maxrow={row['max_row']};ledger={row['ledger']}"
        )
    require(
        all(row["hunter_hard_sets"] == 0 for row in results),
        "a star exception survives the relative-Hunter cap",
    )
    aggregate = hashlib.sha256(b"LRC14/j6/H1-relative-Hunter/aggregate/v1\n")
    for row in results:
        aggregate.update((row["ledger"] + "\n").encode())
    print(f"aggregate_sha256={aggregate.hexdigest()}")
    print("scope=two K>=20 star exceptions;finite exact;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
