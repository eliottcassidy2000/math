#!/usr/bin/env python3
"""Import-free exact set-cover audit for the 28 labelled (100,629) bodies."""

from __future__ import annotations

import csv
from fractions import Fraction
import pathlib
import sys

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix


def read_rows(path: pathlib.Path):
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pattern = int(row["a_hex"], 16)
            rows.append(
                {
                    "pattern": pattern,
                    "mask": row["least_mask"],
                    "rank": int(row["least_rank"]),
                    "count8": int(row["count8"]),
                    "count9": int(row["count9"]),
                    "count": int(row["count"]),
                    "mass": int(row["mass"]),
                    "ticks": int(row["margin_ticks63"]),
                }
            )
    return rows


def maximal_rows(rows):
    ordered = sorted(
        rows,
        key=lambda row: (
            -row["pattern"].bit_count(),
            -int(row["rank"] == 8),
            row["mask"],
        ),
    )
    maximal = []
    for row in ordered:
        if any(row["pattern"] & ~prior["pattern"] == 0 for prior in maximal):
            continue
        maximal.append(row)
    return maximal


def rational_dual(rows):
    matrix = np.array(
        [[(row["pattern"] >> bit) & 1 for bit in range(28)] for row in rows],
        dtype=float,
    )
    result = linprog(
        -np.ones(28),
        A_ub=csc_matrix(matrix),
        b_ub=np.ones(len(rows)),
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = {
        bit: Fraction(float(value)).limit_denominator(1_000_000)
        for bit, value in enumerate(result.x)
        if value > 1e-10
    }
    maximum_load = Fraction()
    for row in rows:
        load = sum(
            (weight for bit, weight in weights.items()
             if (row["pattern"] >> bit) & 1),
            Fraction(),
        )
        maximum_load = max(maximum_load, load)
        if load > 1:
            raise RuntimeError(f"rational dual violation {load}")
    return weights, maximum_load


def main():
    path = pathlib.Path(sys.argv[1])
    rows = read_rows(path)
    full = 0
    for row in rows:
        full |= row["pattern"]
        if row["ticks"] < 0:
            raise RuntimeError("atlas contains inactive representative")
    if full != (1 << 28) - 1:
        raise RuntimeError("response quotient misses an obligation")
    maximum = max(row["pattern"].bit_count() for row in rows)
    maximal = maximal_rows(rows)
    if maximum != 13:
        raise RuntimeError("maximum response size changed")

    weights, maximum_load = rational_dual(rows)
    dual_total = sum(weights.values(), Fraction())
    if dual_total != Fraction(7, 2) or maximum_load > 1:
        raise RuntimeError("exact rational dual changed")

    # This witness was selected from the exact response quotient with the
    # rank-eight representative preferred whenever its response class has one.
    witness_masks = ["002ac4c0", "3882a082", "0041c325", "08c28e40"]
    by_mask = {row["mask"]: row for row in rows}
    answer = [by_mask[mask] for mask in witness_masks]

    joined = 0
    multiplicity = [0] * 28
    for row in answer:
        joined |= row["pattern"]
        for bit in range(28):
            multiplicity[bit] += (row["pattern"] >> bit) & 1
    if joined != full:
        raise RuntimeError("witness does not cover all obligations")

    print("R629_EXACT_MIXED_SET_COVER_V1")
    print(
        "UNIVERSE", 28,
        "CLASSES", len(rows),
        "MAXIMAL_CLASSES", len(maximal),
        "MAXIMUM_CLASS_COVER", maximum,
    )
    print("DUAL_TOTAL", dual_total, "SUPPORT", len(weights),
          "MAXIMUM_CLASS_LOAD", maximum_load)
    for bit, weight in weights.items():
        print("DUAL_WEIGHT", bit, weight)
    print("LOWER_BOUND", 4, "REASON", "exact-rational-dual-total-7/2")
    for ordinal, row in enumerate(answer):
        print(
            "WITNESS", ordinal,
            row["mask"],
            "RANK", row["rank"],
            "COVER", row["pattern"].bit_count(),
            "RESPONSE", f"{row['pattern']:07x}",
            "COUNT8", row["count8"],
            "COUNT9", row["count9"],
            "MASS", row["mass"],
            "MARGIN_TICKS63", row["ticks"],
        )
    print(
        "MULTIPLICITY_ONE", sum(value == 1 for value in multiplicity),
        "MULTIPLICITY_TWO", sum(value == 2 for value in multiplicity),
        "MULTIPLICITY_THREE", sum(value == 3 for value in multiplicity),
    )
    print("MINIMUM_ADDITIONS 4")
    print("VERDICT PASS EXACT_MIXED_MINIMUM_FOUR")


if __name__ == "__main__":
    main()
