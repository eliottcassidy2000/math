#!/usr/bin/env python3
"""Independent exact verifier for the one-row THM-1194 certificate.

The discovery code used clipped teeth.  This verifier instead uses the exact
cumulative danger-length function

    F(y) = floor(y)/7 + min(frac(y), 1/7).

It also checks the [1,6]-core interval and the embedded 78-slow gap.
"""

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


SCHEMA = "lrc1194-top7-m6-bv-v1"


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


def interval_min_norm(left: Q, right: Q) -> Q:
    start = left.numerator // left.denominator - 1
    stop = right.numerator // right.denominator + 2
    if any(left <= integer <= right for integer in range(start, stop + 1)):
        return Q()
    return min(
        [abs(left - integer) for integer in range(start, stop + 1)]
        + [abs(right - integer) for integer in range(start, stop + 1)]
    )


def verify(document: dict) -> dict:
    if document.get("schema") != SCHEMA or len(document.get("rows", [])) != 1:
        raise RuntimeError("certificate schema or row count")
    row = document["rows"][0]
    a, k = row["a"], row["k"]
    multiplier, bins = row["multiplier"], row["bins"]
    weights = row["weights"]
    if (a, k, multiplier, bins) != (78, 11, 140, 200):
        raise RuntimeError("unexpected certified horn")
    if len(weights) != bins or any(weight < 0 for weight in weights):
        raise RuntimeError("invalid weights")
    total = sum(weights)
    if total != row["total"] or total <= 0:
        raise RuntimeError("invalid total mass")

    cutoff = multiplier * a
    low, low_speed = max(
        (comb_load(a, k, b, weights, total), b)
        for b in range(a + 1, cutoff + 1)
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
        raise RuntimeError("stored exact invariant mismatch")
    if not (low < Q(1, 6) and tail < Q(1, 6)):
        raise RuntimeError("one-comb mass cap failed")

    core_interval = Q(11, 84), Q(13, 84)
    slow_gap = Q(14 * k + 1, 14 * a), Q(14 * k + 13, 14 * a)
    if not (
        core_interval[0] <= slow_gap[0]
        <= slow_gap[1] <= core_interval[1]
    ):
        raise RuntimeError("slow gap is not inside the core interval")
    core_minima = [
        interval_min_norm(speed * core_interval[0], speed * core_interval[1])
        for speed in range(1, 7)
    ]
    if min(core_minima) < Q(1, 14):
        raise RuntimeError("[1,6]-core is not safe")

    return {
        "low": low,
        "low_speed": low_speed,
        "low_margin": Q(1, 6) - low,
        "tail": tail,
        "tail_margin": Q(1, 6) - tail,
        "variation": variation,
        "core_minima": core_minima,
        "core_interval": core_interval,
        "slow_gap": slow_gap,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    result = verify(json.loads(args.certificate.read_text()))
    print("VERIFIED_EXACT THM-1194 rows=1 a=78 k=11 multiplier=140 bins=200")
    print(f"max_low={result['low']} at_b={result['low_speed']}")
    print(f"low_margin={result['low_margin']}")
    print(f"tail={result['tail']}")
    print(f"tail_margin={result['tail_margin']}")
    print(f"variation={result['variation']}")
    print("core_minima=" + ",".join(map(str, result["core_minima"])))
    print(f"core_interval={result['core_interval'][0]},{result['core_interval'][1]}")
    print(f"slow_gap={result['slow_gap'][0]},{result['slow_gap'][1]}")
    print("a79_search_failure_is_telemetry_only=true")


if __name__ == "__main__":
    main()
