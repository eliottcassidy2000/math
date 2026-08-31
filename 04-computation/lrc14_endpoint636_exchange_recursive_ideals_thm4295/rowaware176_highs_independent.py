#!/usr/bin/env python3
"""Independent HiGHS audit of the 14-mask downstream optimum."""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


EXPECTED_SHA256 = "cfe73a7d08795dde301af82629b3a6003b5afd6c889ca8da017826fd6afd4e21"


def load(path: Path) -> list[tuple[int, int]]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SHA256:
        raise ValueError("176-response atlas identity changed")
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pattern = (
                int(row["w0"], 16)
                | (int(row["w1"], 16) << 64)
                | (int(row["w2"], 16) << 128)
            )
            rows.append((pattern, int(row["least_mask_hex"], 16)))
    if len(rows) != 2620:
        raise ValueError("176-response class count changed")
    return rows


def maximal(rows: list[tuple[int, int]]) -> list[tuple[int, int]]:
    keep = []
    for i, row in enumerate(rows):
        pattern = row[0]
        if not any(
            i != j and pattern != other[0] and pattern & ~other[0] == 0
            for j, other in enumerate(rows)
        ):
            keep.append(row)
    if len(keep) != 379:
        raise ValueError("maximal response-class count changed")
    return keep


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: rowaware176_highs_independent.py ATLAS.tsv")
    rows = maximal(load(Path(sys.argv[1])))
    masks = len(rows)
    misses = 75
    variables = masks + misses
    objective = np.zeros(variables)
    objective[masks:] = 1.0
    matrix = []
    lower = []
    upper = []

    for bit in range(101):
        row = np.zeros(variables)
        for index, (pattern, _) in enumerate(rows):
            row[index] = (pattern >> bit) & 1
        matrix.append(row)
        lower.append(1.0)
        upper.append(np.inf)

    cardinality = np.zeros(variables)
    cardinality[:masks] = 1.0
    matrix.append(cardinality)
    lower.append(-np.inf)
    upper.append(14.0)

    for local, bit in enumerate(range(101, 176)):
        row = np.zeros(variables)
        for index, (pattern, _) in enumerate(rows):
            row[index] = (pattern >> bit) & 1
        row[masks + local] = 1.0
        matrix.append(row)
        lower.append(1.0)
        upper.append(np.inf)

    result = milp(
        c=objective,
        integrality=np.ones(variables),
        bounds=Bounds(np.zeros(variables), np.ones(variables)),
        constraints=LinearConstraint(np.asarray(matrix), lower, upper),
        options={"time_limit": 300.0, "mip_rel_gap": 0.0, "presolve": True},
    )
    if not result.success or result.fun is None or result.x is None:
        raise RuntimeError(result.message)
    chosen = [i for i in range(masks) if result.x[i] > 0.5]
    missed = [101 + i for i in range(misses) if result.x[masks + i] > 0.5]
    union = 0
    for index in chosen:
        union |= rows[index][0]
    if (
        len(chosen) != 14
        or (union & ((1 << 101) - 1)) != (1 << 101) - 1
        or len(missed) != 37
        or round(result.fun) != 37
    ):
        raise ValueError("independent optimum/witness changed")

    print("THM4295_ROWAWARE176_HIGHS_INDEPENDENT_V1")
    print(
        f"ATLAS_SHA256 {EXPECTED_SHA256} MAXIMAL {masks} "
        f"VARIABLES {variables}"
    )
    print(
        f"STATUS {result.status} SUCCESS {result.success} "
        f"OBJECTIVE {result.fun:.0f} MIP_GAP {getattr(result, 'mip_gap', None)} "
        f"NODE_COUNT {getattr(result, 'mip_node_count', None)}"
    )
    print("CHOSEN " + " ".join(f"{rows[i][1]:08x}" for i in chosen))
    print("MISSED_BITS " + " ".join(map(str, missed)))
    print("VERDICT PASS INDEPENDENT_HIGHS_EXACT14_DOWNSTREAM_OPTIMUM")


if __name__ == "__main__":
    main()
