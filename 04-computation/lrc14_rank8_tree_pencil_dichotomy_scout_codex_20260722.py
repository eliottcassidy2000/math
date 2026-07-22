#!/usr/bin/env python3
"""Exact bounded scout for THM-2103's tree-or-affine-pencil dichotomy.

Two hostile 3^8 coefficient cubes are exhausted with cached THM-2099 edge
weights.  Every row has guard g=(1,0), eight transverse columns (x_i,y_i),
and distinct y_i, so a sufficiently large direction (1,D) gives an odd guard
and distinct positive terminal specializations.

This is bounded evidence only.  It does not replace the missing relation-
height tail in THM-2103.

Tournament Analysis: the exact pair weight is the observable.  Sorting edges
by weight supplies Kruskal's graphic-matroid maximum; orienting tied edges by
labels would add no invariant.  We report weight ties, low-edge degrees, the
MST, and the signed-affine gauge for exceptional rows.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from itertools import product
from pathlib import Path


HERE = Path(__file__).resolve().parent
EDGE_SCRIPT = HERE / "lrc14_mixed_torus_edge_spectrum_scout_codex_20260722.py"
SPEC = importlib.util.spec_from_file_location("mixed_edge", EDGE_SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load THM-2099 edge provider")
EDGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EDGE)

GUARD = (1, 0)
BUDGET = F(5, 49)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def det(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[1] - a[1] * b[0]


def circle_distance(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def small_safe_point(
    vectors: tuple[tuple[int, int], ...], max_denominator: int = 14
) -> tuple[F, F, tuple[F, ...]] | None:
    """Find a proof-readable rational mixed-safe point in a small grid."""
    for denominator in range(2, max_denominator + 1):
        for x_num in range(1, denominator):
            x = F(x_num, denominator)
            if circle_distance(x) <= F(1, 7):
                continue
            for y_num in range(denominator):
                y = F(y_num, denominator)
                distances = tuple(circle_distance(a * x + b * y) for a, b in vectors)
                if all(distance > F(1, 14) for distance in distances):
                    return x, y, distances
    return None


def signed_affine_pencil(
    vectors: tuple[tuple[int, int], ...],
) -> tuple[int, ...] | None:
    """Return a sign gauge for a nondegenerate affine pencil, if one exists."""
    n = len(vectors)
    for bits in range(1 << (n - 1)):
        signs = (1,) + tuple(-1 if bits & (1 << (i - 1)) else 1 for i in range(1, n))
        points = tuple((sign * v[0], sign * v[1]) for sign, v in zip(signs, vectors))
        direction = None
        for point in points[1:]:
            delta = (point[0] - points[0][0], point[1] - points[0][1])
            if delta != (0, 0):
                direction = delta
                break
        if direction is None:
            continue
        if any(
            det(direction, (point[0] - points[0][0], point[1] - points[0][1])) != 0
            for point in points[1:]
        ):
            continue
        if det(direction, GUARD) == 0:
            continue
        if det(direction, points[0]) == 0:
            continue
        return signs
    return None


def maximum_tree_from_edges(
    n: int, edges: list[tuple[F, int, int]]
) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    parent = list(range(n))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    chosen: list[tuple[int, int, F]] = []
    for weight, i, j in sorted(edges, reverse=True):
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
            chosen.append((i, j, weight))
            if len(chosen) == n - 1:
                break
    require(len(chosen) == n - 1, "complete graph failed to span")
    return sum((edge[2] for edge in chosen), F(0)), tuple(chosen)


def cached_cube(
    name: str,
    y_values: tuple[int, ...],
    x_values: tuple[int, ...],
) -> dict[str, object]:
    require(len(y_values) == 8 and len(set(y_values)) == 8, "need eight distinct y-values")
    cache: dict[tuple[int, int, int, int], F] = {}
    for i in range(8):
        for j in range(i + 1, 8):
            for x_i in x_values:
                for x_j in x_values:
                    cache[(i, x_i, j, x_j)] = EDGE.edge_data(
                        (x_i, y_values[i]), (x_j, y_values[j])
                    )[0]

    minimum_nonpencil: tuple[F, tuple[tuple[int, int], ...] | None, object] = (
        F(10),
        None,
        None,
    )
    at_or_below = 0
    pencil_at_or_below = 0
    nonpencil_at_or_below = 0
    exceptions: list[tuple[F, tuple[tuple[int, int], ...], object]] = []
    minimum_all = F(10)
    for assignment in product(x_values, repeat=8):
        vectors = tuple((assignment[i], y_values[i]) for i in range(8))
        edges = [
            (cache[(i, assignment[i], j, assignment[j])], i, j)
            for i in range(8)
            for j in range(i + 1, 8)
        ]
        tau, tree = maximum_tree_from_edges(8, edges)
        minimum_all = min(minimum_all, tau)
        gauge = signed_affine_pencil(vectors)
        if tau <= BUDGET:
            at_or_below += 1
            if gauge is None:
                nonpencil_at_or_below += 1
                exceptions.append((tau, vectors, tree))
            else:
                pencil_at_or_below += 1
        if gauge is None and tau < minimum_nonpencil[0]:
            minimum_nonpencil = (tau, vectors, tree)

    require(at_or_below == pencil_at_or_below + nonpencil_at_or_below, "count split failed")
    print(f"{name}: rows={len(x_values)**8}")
    print(
        f"  minimum-all={EDGE.fmt(minimum_all)}; tau<=5/49={at_or_below}; "
        f"pencil={pencil_at_or_below}; nonpencil={nonpencil_at_or_below}"
    )
    tau, vectors, tree = minimum_nonpencil
    require(vectors is not None, "no nonpencil row found")
    print(
        f"  minimum-nonpencil={EDGE.fmt(tau)}; "
        f"surplus={EDGE.fmt(tau-BUDGET)}; vectors={vectors}"
    )
    print("  tree=" + str(tuple((i, j, EDGE.fmt(w)) for i, j, w in tree)))
    for index, (bad_tau, bad_vectors, bad_tree) in enumerate(
        sorted(exceptions, key=lambda row: (row[0], row[1])), start=1
    ):
        witness = small_safe_point(bad_vectors)
        require(witness is not None, f"exception {index} lacks a small safe point")
        x, y, distances = witness
        print(
            f"  exception-{index}: tau={EDGE.fmt(bad_tau)}; "
            f"vectors={bad_vectors}; safe=({EDGE.fmt(x)},{EDGE.fmt(y)}); "
            f"terminal-distances={tuple(EDGE.fmt(value) for value in distances)}"
        )
        print(
            "    tree="
            + str(tuple((i, j, EDGE.fmt(weight)) for i, j, weight in bad_tree))
        )
    return {
        "at_or_below": at_or_below,
        "pencil": pencil_at_or_below,
        "nonpencil": nonpencil_at_or_below,
        "minimum_nonpencil": tau,
        "exceptions": tuple(exceptions),
    }


def main() -> None:
    print("THM-2103 RANK-EIGHT TREE-OR-PENCIL EXACT BOUNDED SCOUT")
    consecutive = cached_cube("consecutive-y cube", tuple(range(1, 9)), (-1, 0, 1))
    dyadic = cached_cube("dyadic-y cube", tuple(2**k for k in range(8)), (0, 1, 2))
    require(consecutive["nonpencil"] == 0, "consecutive cube refuted target")
    require(dyadic["nonpencil"] == 6, "dyadic counterexample count changed")
    require(consecutive["minimum_nonpencil"] > BUDGET, "consecutive strict gap failed")
    require(dyadic["minimum_nonpencil"] < BUDGET, "dyadic no-go disappeared")
    print("PASS -- 2*3^8 exact rows refute (T); all six no-go rows have small safe phases")


if __name__ == "__main__":
    main()
