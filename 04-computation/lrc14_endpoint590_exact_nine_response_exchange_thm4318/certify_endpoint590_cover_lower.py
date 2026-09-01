#!/usr/bin/env python3
"""Seek and exactly validate an integer-dual lower bound for the cover size."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

from solve_endpoint590_response_cover import build_matrix, read_signatures, require


def exact_loads(
    records: list[tuple[int, int, int]], weights: list[int]
) -> list[int]:
    loads: list[int] = []
    for reply, _, _ in records:
        load = 0
        work = reply
        while work:
            least = work & -work
            load += weights[least.bit_length() - 1]
            work ^= least
        loads.append(load)
    return loads


def find_certificate(records, matrix):
    for denominator in range(1, 17):
        result = milp(
            c=-np.ones(100),
            integrality=np.ones(100),
            bounds=Bounds(np.zeros(100), np.full(100, denominator)),
            constraints=LinearConstraint(
                matrix.T,
                np.full(len(records), -np.inf),
                np.full(len(records), denominator),
            ),
            options={"mip_rel_gap": 0.0, "time_limit": 60, "presolve": True},
        )
        if result.x is None:
            continue
        weights = [int(round(value)) for value in result.x]
        loads = exact_loads(records, weights)
        total = sum(weights)
        if min(weights) >= 0 and max(loads) <= denominator and total > 7 * denominator:
            return denominator, weights, loads, "INTEGER_MILP_EXACTLY_RECHECKED"

    lp = linprog(
        -np.ones(100),
        A_ub=matrix.T,
        b_ub=np.ones(len(records)),
        bounds=(0, None),
        method="highs",
    )
    require(lp.success, f"dual LP failed: {lp.message}")
    for denominator in range(17, 10001):
        scaled = lp.x * denominator
        rounded = np.rint(scaled)
        if np.max(np.abs(scaled - rounded)) > 1e-7:
            continue
        weights = [int(value) for value in rounded]
        loads = exact_loads(records, weights)
        total = sum(weights)
        if min(weights) >= 0 and max(loads) <= denominator and total > 7 * denominator:
            return denominator, weights, loads, "LP_RATIONALIZED_EXACTLY_RECHECKED"
    raise RuntimeError("no integer-dual lower-eight certificate found")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signatures", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = read_signatures(args.signatures)
    matrix = build_matrix(records)
    denominator, weights, loads, method = find_certificate(records, matrix)
    total = sum(weights)
    maximum_load = max(loads)
    lower = math.ceil(total / denominator)
    require(lower == 8 and total > 7 * denominator,
            "certificate does not prove cover lower bound eight")
    require(maximum_load <= denominator,
            "a realized response violates the exact dual constraint")

    args.weights.parent.mkdir(parents=True, exist_ok=True)
    with args.weights.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["failure_ordinal", "weight"])
        for index, weight in enumerate(weights):
            if weight:
                writer.writerow([index, weight])
    args.output.write_text(
        "ENDPOINT590_RESPONSE_COVER_INTEGER_DUAL_V1\n"
        f"SIGNATURES {len(records)} DENOMINATOR {denominator}\n"
        f"WEIGHT_TOTAL {total} MAX_RESPONSE_LOAD {maximum_load}\n"
        f"LOWER_BOUND {lower} METHOD {method}\n"
        "SCOPE EXACT_OVER_COMPLETE_REALIZED_RANK8_RANK9_RESPONSE_ATLAS\n"
        "VERDICT PASS\n",
        encoding="ascii",
    )
    print(args.output.read_text(encoding="ascii"), end="")


if __name__ == "__main__":
    main()
