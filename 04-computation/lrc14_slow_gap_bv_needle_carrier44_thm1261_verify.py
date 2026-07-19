#!/usr/bin/env python3
"""Independent cumulative-integral verifier for THM-1261."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path


SCHEMA = "lrc14-slow-gap-bv-needle-carrier44-thm1261-v1"
CARRIER = 44
BINS = 200
ALLOWED_MULTIPLIERS = {12, 16, 20, 28, 44, 60, 76, 108}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"independent THM-1261 verification failed: {message}")


def cumulative(y: Q) -> Q:
    integer = y.numerator // y.denominator
    return Q(integer, 7) + min(y - integer, Q(1, 7))


def comb_load(k: int, speed: int, weights: list[int], total: int) -> Q:
    left = Q(14 * k + 1, 14 * CARRIER)
    width = Q(6, 7 * CARRIER * BINS)
    values = [
        cumulative(speed * (left + j * width) + Q(1, 14))
        for j in range(BINS + 1)
    ]
    weighted = weights[-1] * values[-1] - weights[0] * values[0]
    weighted += sum(
        (weights[j - 1] - weights[j]) * values[j]
        for j in range(1, BINS)
    )
    return weighted / (speed * width * total)


def check_row(row: dict) -> tuple[Q, Q, Q, int]:
    k, multiplier, weights = row.get("k"), row.get("multiplier"), row.get("weights", [])
    require(row.get("a") == CARRIER and k in range(22), "row address")
    require(row.get("bins") == BINS and len(weights) == BINS, "row shape")
    require(multiplier in ALLOWED_MULTIPLIERS, "multiplier")
    require(all(w >= 0 for w in weights), "weights")
    total = sum(weights)
    require(total > 0 and total == row.get("total"), "mass")
    cutoff = multiplier * CARRIER
    low, argmax_speed = max(
        (comb_load(k, speed, weights, total), speed)
        for speed in range(CARRIER + 1, cutoff + 1)
    )
    variation_mass = (
        weights[0]
        + sum(abs(y - x) for x, y in zip(weights, weights[1:]))
        + weights[-1]
    )
    width = Q(6, 7 * CARRIER * BINS)
    variation = Q(variation_mass, total) / width
    tail = Q(1, 7) + Q(3, 49 * (cutoff + 1)) * variation
    stored = (
        Q(row["low_num"], row["low_den"]),
        Q(row["tail_num"], row["tail_den"]),
        Q(row["variation_num"], row["variation_den"]),
    )
    require((low, tail, variation) == stored, "stored invariants")
    require(low < Q(1, 6) and tail < Q(1, 6), "strict loads")
    return low, tail, variation, argmax_speed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    document = json.loads(args.certificate.read_text())
    require(document.get("schema") == SCHEMA, "schema")
    require(document.get("carrier") == CARRIER, "carrier")
    require(document.get("reflection") == "k maps to 43-k modulo 44", "reflection")
    rows = document.get("rows", [])
    keys = [(row.get("a"), row.get("k")) for row in rows]
    require(
        set(keys) == {(CARRIER, k) for k in range(22)}
        and len(keys) == len(set(keys)),
        "phase rows",
    )
    checked = [
        (row, *check_row(row)) for row in sorted(rows, key=lambda item: item["k"])
    ]
    count = Counter(row["multiplier"] for row, _, _, _, _ in checked)
    low_row, max_low, _, _, low_speed = max(checked, key=lambda item: item[1])
    tail_row, _, max_tail, tail_variation, _ = max(checked, key=lambda item: item[2])
    print("THM-1261 CARRIER-44 INDEPENDENT CUMULATIVE REPLAY")
    print(
        "VERIFIED_INDEPENDENT_EXACT", f"rows={len(rows)}", "carrier=44",
        "multipliers=" + ",".join(f"{m}:{count[m]}" for m in sorted(count)),
    )
    print(
        f"max_low={max_low} gap_to_1/6={Q(1,6)-max_low} "
        f"at_k={low_row['k']} M={low_row['multiplier']} "
        f"argmax_d={low_speed} quotient_residue={divmod(low_speed,CARRIER)}"
    )
    print(
        f"max_tail={max_tail} gap_to_1/6={Q(1,6)-max_tail} "
        f"at_k={tail_row['k']} M={tail_row['multiplier']}"
    )
    print("hardest_tail_unit_gap_TV=" + str(Q(6, 7*CARRIER) * tail_variation))
    print("RESULT=no six-comb cover at carrier 44")


if __name__ == "__main__":
    main()
