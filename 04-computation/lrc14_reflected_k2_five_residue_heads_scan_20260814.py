#!/usr/bin/env python3
"""Exact scan of the five residual reflected-residue projected-k2 heads."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import sys
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path


HERE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else HERE_ROOT
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
PAIR = ROOT / "04-computation/lrc_general_reflected_pair_mass_thm3352.py"
FLOOR_THEOREM = ROOT / "01-canon/theorems/THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails.md"
EXPECTED = {
    MEDIAN: "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276",
    PAIR: "afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b",
    FLOOR_THEOREM: "bd5237632ff70293f4701674ce72fd95e3255a0aee05a62a1986ff62df773633",
}
SHAPES = (
    (1, 2, 3, 4, 6),
    (1, 2, 3, 4, 12),
    (1, 2, 3, 6, 12),
    (1, 2, 4, 6, 8),
    (1, 2, 4, 6, 12),
)
EDGES = tuple(combinations(range(5), 2))
TARGET = F(25, 91)
EXPECTED_SEMANTIC = "9a3a9ab7b779c990d41f908469159f52c02a5e09053330b45269fdf03b5f7685"
_MEDIAN_MODULE = None
_PAIR_MODULE = None


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mst(weights):
    parent = list(range(5))

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F()
    chosen = []
    for weight, i, j in sorted(weights, reverse=True):
        a, b = root(i), root(j)
        if a == b:
            continue
        parent[a] = b
        total += weight
        chosen.append((weight, i, j))
        if len(chosen) == 4:
            return total, tuple(chosen)
    raise RuntimeError("K5 disconnected")


def worker(item):
    global _MEDIAN_MODULE, _PAIR_MODULE
    body_index, body, ruler = item
    if _MEDIAN_MODULE is None:
        _MEDIAN_MODULE = load("k2_head_worker_median", MEDIAN)
    if _PAIR_MODULE is None:
        _PAIR_MODULE = load("k2_head_worker_pair", PAIR)
    median = _MEDIAN_MODULE
    pair = _PAIR_MODULE
    actual, ranges = median.R.safe_cell_ranges(body)
    require(actual == ruler, (body, actual, ruler))
    cells = tuple(cell for left, right in ranges for cell in range(left, right))
    cell = cells[len(cells) // 2]
    require(median.R.body_cell_is_safe(ruler, body, cell), (body, cell))
    failures = 0
    weakest = None
    shape_weakest = {}
    checks = 0
    for shape in SHAPES:
        cache = {}
        for slots in permutations(range(6), 5):
            debt = sum(
                (F(body[slots[i]], 7 * (shape[i] * ruler - body[slots[i]])) for i in range(5)),
                F(),
            )
            weights = []
            for i, j in EDGES:
                key = (slots[i], shape[i], slots[j], shape[j])
                if key not in cache:
                    cache[key] = pair.mass(
                        ruler, cell,
                        body[slots[i]], shape[i],
                        body[slots[j]], shape[j],
                    )
                weights.append((cache[key], i, j))
            credit, tree = mst(weights)
            margin = F(2, 7) - debt + credit - TARGET
            record = (margin, body_index, body, ruler, cell, shape, slots, debt, credit, tree)
            if margin <= 0:
                failures += 1
            if weakest is None or record < weakest:
                weakest = record
            prior = shape_weakest.get(shape)
            if prior is None or record < prior:
                shape_weakest[shape] = record
            checks += 1
        require(len(cache) == 300, (body, shape, len(cache)))
    require(checks == 3600, (body, checks))
    return body_index, checks, failures, weakest, tuple(shape_weakest[shape] for shape in SHAPES)


def interval_intersection(first, second):
    i = j = 0
    total = F()
    while i < len(first) and j < len(second):
        total += max(F(), min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", default=str(ROOT))
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()
    require(Path(args.root).resolve() == ROOT.resolve(), (args.root, ROOT))
    for path, expected in EXPECTED.items():
        require(sha(path) == expected, (path, sha(path)))
    median = load("k2_heads_median", MEDIAN)
    bodies = median.body_universe()
    require(len(bodies) == 649, len(bodies))
    tasks = tuple((index, body, ruler) for index, (body, ruler) in enumerate(bodies))
    if args.processes == 1:
        rows = tuple(worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(args.processes) as pool:
            rows = tuple(pool.imap(worker, tasks, chunksize=1))
    require(tuple(row[0] for row in rows) == tuple(range(649)), "body order")
    checks = sum(row[1] for row in rows)
    failures = sum(row[2] for row in rows)
    require(checks == 649 * 5 * 720 == 2_336_400, checks)
    require(failures == 0, failures)
    weakest = min(row[3] for row in rows)
    shape_weakest = tuple(min(row[4][i] for row in rows) for i in range(5))

    # Separate literal-interval control on every edge of the weakest row
    # and on the resulting five-clause union.
    pair = load("k2_heads_pair_control", PAIR)
    (_margin, _bi, body, ruler, cell, shape, slots, debt, credit, tree) = weakest
    clauses = tuple(
        median.R.reflected_level_arcs(ruler, body[slots[i]], shape[i], cell)
        for i in range(5)
    )
    edge_controls = 0
    for i, j in EDGES:
        literal = interval_intersection(clauses[i], clauses[j])
        fast = pair.mass(ruler, cell, body[slots[i]], shape[i], body[slots[j]], shape[j])
        require(literal == fast, (i, j, literal, fast))
        edge_controls += 1
    union = median.R.merge_intervals(tuple(arc for clause in clauses for arc in clause))
    actual_safe = 1 - median.R.interval_mass(union)
    hunter_safe = F(2, 7) - debt + credit
    require(actual_safe >= hunter_safe > TARGET, (actual_safe, hunter_safe))
    require(hunter_safe - TARGET == weakest[0], (hunter_safe, weakest[0]))

    packet = (tuple((path.name, expected) for path, expected in EXPECTED.items()), SHAPES,
              checks, failures, weakest, shape_weakest, edge_controls, actual_safe)
    semantic = hashlib.sha256(repr(packet).encode("ascii")).hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    require(
        weakest[:7]
        == (F(1107040893783422917, 37818341551113357325), 556,
            (2, 6, 7, 8, 9, 12), 7056, 3780, (1, 2, 3, 4, 6),
            (1, 4, 0, 3, 2)),
        weakest,
    )
    compact_shape_weakest = tuple(
        (record[5], record[0], record[2], record[4], record[6])
        for record in shape_weakest
    )
    print("LRC14 REFLECTED-RESIDUE PROJECTED-K2 FIVE-HEAD EXACT SCAN")
    print("scope=649_upper_median_bodies;five_reflected_residues;five_primitive_scale1_heads")
    print("universe", len(bodies), len(SHAPES), 720, "checks", checks, "failures", failures)
    print("weakest", weakest)
    print("shape_minima", compact_shape_weakest)
    print("literal_control", edge_controls, "actual_safe", actual_safe, "hunter_safe", hunter_safe)
    print("semantic_sha256", semantic)
    print("conclusion=all_five_scale1_heads_close_by_exact_common_cell_K5_Hunter_tree")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
