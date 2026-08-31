#!/usr/bin/env python3
"""Exact branch-and-bound set-cover audit for THM-4296 response atlases.

This promotes the scratch generator that produced the recorded endpoint-636
and r632 outputs.  The only interface extension is a maintained-packet path
fallback, so the recorded atlases can be replayed without copying or renaming.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix


MAINTAINED_INPUTS = {
    "LOCAL_100": pathlib.Path("endpoint636_atlas/LOCAL_100.tsv"),
    "LOCAL_256": pathlib.Path("endpoint636_atlas/LOCAL_256.tsv"),
    "COMMON_ACTIVE": pathlib.Path("endpoint636_atlas/COMMON_ACTIVE.tsv"),
    "CARRIER_UNION": pathlib.Path("endpoint636_atlas/CARRIER_UNION.tsv"),
    "r632_mixed89_atlas": pathlib.Path("r632_atlas/mixed89.tsv"),
}


def read_classes(path: pathlib.Path):
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pattern = int(row["a_hex"], 16) | (int(row["b_hex"], 16) << 64)
            rows.append((pattern, row["least_mask"], int(row["count"])))
    return rows


def class_path(root: pathlib.Path, name: str) -> pathlib.Path:
    legacy = root / f"{name}.tsv"
    if legacy.is_file() or name not in MAINTAINED_INPUTS:
        return legacy
    return root / MAINTAINED_INPUTS[name]


def maximal_rows(rows):
    ordered = sorted(rows, key=lambda row: (-row[0].bit_count(), row[1]))
    maximal = []
    for row in ordered:
        if any(row[0] & ~prior[0] == 0 for prior in maximal):
            continue
        maximal.append(row)
    return maximal


def rational_dual(rows, full):
    bits = [bit for bit in range(101) if (full >> bit) & 1]
    matrix = np.array(
        [[(pattern >> bit) & 1 for bit in bits] for pattern, _, _ in rows],
        dtype=float,
    )
    result = linprog(
        -np.ones(len(bits)),
        A_ub=csc_matrix(matrix),
        b_ub=np.ones(len(rows)),
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = {
        bit: Fraction(float(value)).limit_denominator(1_000_000)
        for bit, value in zip(bits, result.x)
        if value > 1e-10
    }
    for pattern, _, _ in rows:
        load = sum(value for bit, value in weights.items() if (pattern >> bit) & 1)
        if load > 1:
            raise RuntimeError(f"rationalized dual violation {load}")
    return weights


class Search:
    def __init__(self, rows, full, weights):
        self.rows = rows
        self.patterns = [row[0] for row in rows]
        self.full = full
        self.weights = weights
        self.by_bit = [[] for _ in range(101)]
        for index, pattern in enumerate(self.patterns):
            for bit in range(101):
                if (pattern >> bit) & 1:
                    self.by_bit[bit].append(index)
        self.memo = set()
        self.nodes = 0
        self.prune_gain = 0
        self.prune_dual = 0
        self.prune_memo = 0

    def lower_bound(self, remaining):
        if not remaining:
            return 0
        maximum = max((pattern & remaining).bit_count() for pattern in self.patterns)
        gain_bound = math.ceil(remaining.bit_count() / maximum)
        weight = sum(value for bit, value in self.weights.items() if (remaining >> bit) & 1)
        dual_bound = math.ceil(weight)
        return max(gain_bound, dual_bound)

    def dfs(self, covered, depth):
        self.nodes += 1
        if covered == self.full:
            return []
        key = (covered, depth)
        if key in self.memo:
            self.prune_memo += 1
            return None
        remaining = self.full & ~covered
        if depth == 0:
            self.memo.add(key)
            return None
        maximum = max((pattern & remaining).bit_count() for pattern in self.patterns)
        if remaining.bit_count() > depth * maximum:
            self.prune_gain += 1
            self.memo.add(key)
            return None
        weight = sum(value for bit, value in self.weights.items() if (remaining >> bit) & 1)
        if weight > depth:
            self.prune_dual += 1
            self.memo.add(key)
            return None

        uncovered_bits = [bit for bit in range(101) if (remaining >> bit) & 1]
        chosen_bit = min(
            uncovered_bits,
            key=lambda bit: (
                sum(1 for index in self.by_bit[bit] if self.patterns[index] & remaining),
                -max((self.patterns[index] & remaining).bit_count() for index in self.by_bit[bit]),
                bit,
            ),
        )
        choices = sorted(
            self.by_bit[chosen_bit],
            key=lambda index: (
                -(self.patterns[index] & remaining).bit_count(),
                self.rows[index][1],
            ),
        )
        seen_next = set()
        for index in choices:
            joined = covered | self.patterns[index]
            if joined in seen_next:
                continue
            seen_next.add(joined)
            answer = self.dfs(joined, depth - 1)
            if answer is not None:
                return [index] + answer
        self.memo.add(key)
        return None


def audit(root: pathlib.Path, name: str, expected: int) -> None:
    all_rows = read_classes(class_path(root, name))
    full = 0
    for pattern, _, _ in all_rows:
        full |= pattern
    rows = maximal_rows(all_rows)
    weights = rational_dual(all_rows, full)
    exact_sum = sum(weights.values(), Fraction())
    print(
        "MODE",
        name,
        "UNIVERSE",
        full.bit_count(),
        "CLASSES",
        len(all_rows),
        "MAXIMAL",
        len(rows),
        "DUAL",
        exact_sum,
    )
    for depth in (expected - 1, expected):
        search = Search(rows, full, weights)
        answer = search.dfs(0, depth)
        print(
            " DEPTH",
            depth,
            "FOUND",
            int(answer is not None),
            "NODES",
            search.nodes,
            "MEMO",
            len(search.memo),
            "PRUNE_GAIN",
            search.prune_gain,
            "PRUNE_DUAL",
            search.prune_dual,
            "PRUNE_MEMO",
            search.prune_memo,
        )
        if answer is not None:
            joined = 0
            for index in answer:
                pattern, mask, count = rows[index]
                joined |= pattern
                print(
                    "  WITNESS",
                    mask,
                    "COVER",
                    pattern.bit_count(),
                    "CLASS_COUNT",
                    count,
                    "A",
                    f"{pattern & ((1 << 64) - 1):016x}",
                    "B",
                    f"{pattern >> 64:016x}",
                )
            if joined != full:
                raise RuntimeError("reported witness does not cover universe")
        if depth == expected - 1 and answer is not None:
            raise RuntimeError("expected lower boundary failed")
        if depth == expected and answer is None:
            raise RuntimeError("expected optimum witness absent")


def main() -> None:
    # Freeze transcript newlines across Windows and POSIX replays.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    root = pathlib.Path(sys.argv[1])
    expected = {
        "LOCAL_100": 8,
        "LOCAL_256": 6,
        "COMMON_ACTIVE": 15,
        "CARRIER_UNION": 14,
        "R632_LOCAL_100": 3,
        "R632_LOCAL_256": 10,
        "r632_mixed89_atlas": 3,
    }
    names = sys.argv[2:] or [
        name for name in expected if class_path(root, name).is_file()
    ]
    if not names:
        raise RuntimeError("no recognized response atlases found")
    for name in names:
        audit(root, name, expected[name])


if __name__ == "__main__":
    main()
