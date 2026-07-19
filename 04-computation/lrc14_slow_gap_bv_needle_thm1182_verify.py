#!/usr/bin/env python3
"""Pure-standard-library exact verifier for the THM-1182 certificate.

This verifier is independent of the generator's clipped-tooth integration.
It uses the exact cumulative danger-length function ``F`` instead.
"""

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


SCHEMA = "lrc14-slow-gap-bv-needle-thm1182-v1"


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
    # Telescope sum_j w_j (F_(j+1)-F_j), then divide only once.
    weighted = weights[-1] * values[-1] - weights[0] * values[0]
    weighted += sum(
        (weights[j - 1] - weights[j]) * values[j]
        for j in range(1, bins)
    )
    return weighted / (b * width * total)


def check_row(row: dict) -> None:
    a = row["a"]
    k = row["k"]
    multiplier = row["multiplier"]
    bins = row["bins"]
    weights = row["weights"]
    cutoff = multiplier * a
    if (
        bins != 200
        or len(weights) != bins
        or multiplier not in (12, 16, 20, 28)
        or any(weight < 0 for weight in weights)
    ):
        raise RuntimeError(("shape", a, k))
    total = sum(weights)
    if total != row["total"] or total <= 0:
        raise RuntimeError(("mass", a, k))

    low = max(
        (comb_load(a, k, b, weights, total) for b in range(a + 1, cutoff + 1)),
        default=Q(),
    )
    variation_mass = (
        weights[0]
        + sum(abs(y - x) for x, y in zip(weights, weights[1:]))
        + weights[-1]
    )
    width = Q(6, 7 * a * bins)
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
        raise RuntimeError(("bound", a, k, low, tail))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    document = json.loads(args.certificate.read_text())
    if document.get("schema") != SCHEMA:
        raise RuntimeError("certificate schema mismatch")

    expected = {
        (a, k)
        for a in range(1, document["max_a"] + 1)
        for k in range((a + 1) // 2)
    }
    actual = {(row["a"], row["k"]) for row in document["rows"]}
    if actual != expected or len(actual) != len(document["rows"]):
        raise RuntimeError("phase representatives")

    # k -> a-1-k has exactly one fixed point (a-1)/2 iff a is odd.
    for a in range(1, document["max_a"] + 1):
        fixed = [k for k in range(a) if (a - 1 - k) % a == k]
        expected_fixed = [(a - 1) // 2] if a % 2 else []
        if fixed != expected_fixed:
            raise RuntimeError(("reflection", a, fixed))

    for index, row in enumerate(document["rows"], 1):
        check_row(row)
        if index % 40 == 0:
            print(f"checked={index}", flush=True)
    print(
        f"VERIFIED_EXACT rows={len(document['rows'])} "
        f"max_a={document['max_a']} constant=3/49 threshold=1/6"
    )


if __name__ == "__main__":
    main()
