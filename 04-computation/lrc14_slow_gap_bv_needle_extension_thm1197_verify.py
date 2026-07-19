#!/usr/bin/env python3
"""Independent exact verifier for the separate THM-1197 extension.

This pure-standard-library replay does not import the LP generator or its
clipped-tooth integrator.  It integrates each comb using the cumulative
periodic danger-length function instead.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path


SCHEMA = "lrc14-slow-gap-bv-needle-extension-thm1197-v1"
MIN_A = 31
MAX_A = 40
ALLOWED_MULTIPLIERS = {12, 16, 20, 28, 44, 60, 76, 108}


def cumulative(y: Q) -> Q:
    integer = y.numerator // y.denominator
    fraction = y - integer
    return Q(integer, 7) + min(fraction, Q(1, 7))


def comb_load(a: int, k: int, b: int, weights: list[int], total: int) -> Q:
    bins = len(weights)
    left = Q(14 * k + 1, 14 * a)
    width = Q(6, 7 * a * bins)
    values = [
        cumulative(b * (left + j * width) + Q(1, 14))
        for j in range(bins + 1)
    ]
    weighted = weights[-1] * values[-1] - weights[0] * values[0]
    weighted += sum(
        (weights[j - 1] - weights[j]) * values[j]
        for j in range(1, bins)
    )
    return weighted / (b * width * total)


def check_row(row: dict) -> tuple[Q, Q, Q]:
    a, k = row["a"], row["k"]
    multiplier = row["multiplier"]
    weights = row["weights"]
    if (
        row["bins"] != 200
        or len(weights) != 200
        or multiplier not in ALLOWED_MULTIPLIERS
        or any(weight < 0 for weight in weights)
    ):
        raise RuntimeError(("shape", a, k))
    total = sum(weights)
    if total <= 0 or total != row["total"]:
        raise RuntimeError(("mass", a, k))
    cutoff = multiplier * a
    low = max(
        comb_load(a, k, b, weights, total)
        for b in range(a + 1, cutoff + 1)
    )
    variation_mass = (
        weights[0]
        + sum(abs(y - x) for x, y in zip(weights, weights[1:]))
        + weights[-1]
    )
    width = Q(6, 7 * a * 200)
    variation = Q(variation_mass, total) / width
    tail = Q(1, 7) + Q(3, 49 * (cutoff + 1)) * variation
    stored = (
        Q(row["low_num"], row["low_den"]),
        Q(row["tail_num"], row["tail_den"]),
        Q(row["variation_num"], row["variation_den"]),
    )
    if (low, tail, variation) != stored:
        raise RuntimeError(("stored invariant", a, k))
    if not (low < Q(1, 6) and tail < Q(1, 6)):
        raise RuntimeError(("strict bound", a, k, low, tail))
    return low, tail, variation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    document = json.loads(args.certificate.read_text())
    if (
        document.get("schema") != SCHEMA
        or document.get("min_a") != MIN_A
        or document.get("max_a") != MAX_A
    ):
        raise RuntimeError("schema/range mismatch")
    expected = {
        (a, k)
        for a in range(MIN_A, MAX_A + 1)
        for k in range((a + 1) // 2)
    }
    keys = [(row["a"], row["k"]) for row in document["rows"]]
    if set(keys) != expected or len(keys) != len(set(keys)):
        raise RuntimeError("phase representatives")
    lows: list[Q] = []
    tails: list[Q] = []
    for index, row in enumerate(document["rows"], 1):
        low, tail, _ = check_row(row)
        lows.append(low)
        tails.append(tail)
        if index % 30 == 0:
            print(f"checked={index}", flush=True)
    count = Counter(row["multiplier"] for row in document["rows"])
    print(
        "VERIFIED_INDEPENDENT_EXACT",
        f"rows={len(document['rows'])}",
        f"carriers={MIN_A}..{MAX_A}",
        "multipliers=" + ",".join(f"{key}:{count[key]}" for key in sorted(count)),
    )
    print(f"max_low={max(lows)} gap_to_1/6={Q(1, 6) - max(lows)}")
    print(f"max_tail={max(tails)} gap_to_1/6={Q(1, 6) - max(tails)}")


if __name__ == "__main__":
    main()
