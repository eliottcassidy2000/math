#!/usr/bin/env python3
"""Exact divisor-universal support audit for six projected z1=229 prefixes.

This companion is intentionally independent of the THM-3113 screen engine.
It reconstructs six explicit complete-cell carriers directly from their weak
endpoint inequalities and checks every divisor of each ruler by two separate
projection aggregations (a residue set and a byte histogram).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction
from math import ceil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z229_six_prefix_divisor_universal_support.py"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z229_six_prefix_divisor_universal_support.out"
)

CASES = (
    ((1, 4, 5, 9, 11, 14), 194040, (234, 243), 34442),
    ((1, 5, 6, 8, 9, 14), 35280, (234, 246), 6812),
    ((1, 5, 9, 11, 12, 14), 194040, (234, 243), 35738),
    ((1, 8, 10, 11, 12, 14), 129360, (260, 312), 25358),
    ((2, 5, 8, 9, 11, 14), 388080, (234, 246), 73130),
    ((2, 8, 10, 11, 12, 14), 129360, (260, 364), 25798),
)

EXPECTED_CELL_SHA256 = (
    "4a0831c8736b30d822f646e570c2cd5a11c97d4ad7fbf1de14c60b2fdbf4a38e",
    "53cecffe9a47e74ebded39d293e5de79e803481473bc0b04c47e448924186ca5",
    "ce4e7355f63f2576bfaeafcd9d84f32a1157a69cd889f7a0fc0bad3708797e6d",
    "1fceea75eb8835bc70d8aef1a0c0cd9ee5ec4066fd5dce83f34a0c109b64fadc",
    "8b96033a89cc4ea040a41786855d155b509ddb6b0e8ca3295d0833e30a47b996",
    "a2fbb9979e653e4857a94e5ac9657118d871cc5f0c9c494c9abeddcbe48b58cb",
)
EXPECTED_SEMANTIC_SHA256 = "9f58b66c5ed2930023e1023efeb9fb9d873a856afe895ff2768a31889bd9500a"
EXPECTED_TOTALS = (792, 744, 48, 0)
EXPECTED_EXACT_MODULI = (2, 3, 4, 5, 8, 9, 10, 11, 15, 22)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def divisors(n):
    low = []
    high = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            low.append(d)
            if d * d != n:
                high.append(n // d)
        d += 1
    answer = tuple(low + high[::-1])
    require(answer[0] == 1 and answer[-1] == n, (n, answer))
    require(all(n % d == 0 for d in answer), ("nondivisor", n))
    return answer


def whole_cell_safe(j, label, ruler):
    """Weak closed-cell safety for a strict-open 1/14 danger comb."""
    residue = (j * label) % ruler
    return 14 * residue >= ruler and 14 * (residue + label) <= 13 * ruler


def build_carrier(body, ruler, lows):
    labels = body + (229,) + lows
    return tuple(
        j
        for j in range(ruler)
        if all(whole_cell_safe(j, label, ruler) for label in labels)
    )


def set_projection(cells, modulus):
    return frozenset(j % modulus for j in cells)


def histogram_projection(cells, modulus):
    flags = bytearray(modulus)
    for j in cells:
        flags[j % modulus] = 1
    return frozenset(index for index, flag in enumerate(flags) if flag)


def run_audit():
    minimum_density = min(
        Fraction(expected_size, ruler)
        for _, ruler, _, expected_size in CASES
    )
    require(minimum_density == Fraction(17221, 97020), minimum_density)
    bulk_margin = (
        (minimum_density - Fraction(1, 7)) * 25 - Fraction(6, 7)
    )
    require(bulk_margin == Fraction(173, 19404) and bulk_margin > 0, bulk_margin)

    carrier_rows = []
    divisor_semantic = []
    total_divisors = 0
    total_density = 0
    total_exact = 0
    failures = 0
    exact_moduli = set()
    global_minimum = None

    for index, (body, ruler, lows, expected_size) in enumerate(CASES):
        cells = build_carrier(body, ruler, lows)
        require(len(cells) == expected_size, (body, len(cells), expected_size))
        cell_sha = hashlib.sha256(repr(cells).encode()).hexdigest()
        if EXPECTED_CELL_SHA256[index] != "PENDING":
            require(cell_sha == EXPECTED_CELL_SHA256[index], (body, cell_sha))

        rows = []
        exact_needed = []
        full_supports = 0
        density_closed = 0
        for modulus in divisors(ruler):
            if modulus == 1:
                continue
            support_a = set_projection(cells, modulus)
            support_b = histogram_projection(cells, modulus)
            require(support_a == support_b, (body, modulus, "projection mismatch"))
            support = len(support_a)
            capacity = ceil(modulus / 7)
            fibre_size = ruler // modulus
            coarse = ceil(len(cells) / fibre_size)
            gap = support - capacity
            require(gap > 0, (body, modulus, support, capacity))
            if support == modulus:
                full_supports += 1
            if coarse > capacity:
                density_closed += 1
            else:
                exact_needed.append(
                    (modulus, support, capacity, coarse, support == modulus)
                )
                exact_moduli.add(modulus)
                require(support == modulus, (body, modulus, "tie not full"))
            row = (modulus, support, capacity, coarse, support == modulus)
            rows.append(row)
            divisor_semantic.append((body, ruler, lows, row))
            failures += gap <= 0
            minimum_row = (gap, modulus, support, capacity, coarse)
            if global_minimum is None or minimum_row < global_minimum[0]:
                global_minimum = (minimum_row, body)

        total_divisors += len(rows)
        total_density += density_closed
        total_exact += len(exact_needed)
        carrier_rows.append(
            (
                body,
                ruler,
                lows,
                len(cells),
                cell_sha,
                len(rows),
                full_supports,
                density_closed,
                tuple(exact_needed),
                min(rows, key=lambda row: (row[1] - row[2], row[0])),
            )
        )

    totals = (total_divisors, total_density, total_exact, failures)
    require(totals == EXPECTED_TOTALS, totals)
    require(tuple(sorted(exact_moduli)) == EXPECTED_EXACT_MODULI, exact_moduli)
    require(global_minimum is not None and global_minimum[0][0] == 1, global_minimum)

    # Hostile boundary: total carrier cardinality alone cannot settle d=2.
    body, ruler, lows, expected_size = CASES[0]
    hostile = tuple(2 * index for index in range(expected_size))
    require(hostile[-1] < ruler, ("hostile exceeds ruler", hostile[-1], ruler))
    hostile_support = set_projection(hostile, 2)
    actual_support = set_projection(build_carrier(body, ruler, lows), 2)
    require(len(hostile_support) == 1, hostile_support)
    require(len(actual_support) == 2, actual_support)
    require(ceil(2 / 7) == 1, "d=2 capacity")

    # MISTAKE-334 regression: centered beta(28)=3 is not translated kappa(28)=4.
    beta_28 = 2 * ((28 - 1) // 14) + 1
    kappa_28 = ceil(28 / 7)
    require((beta_28, kappa_28) == (3, 4), (beta_28, kappa_28))

    semantic = (
        tuple(carrier_rows),
        tuple(divisor_semantic),
        totals,
        minimum_density,
        bulk_margin,
        tuple(sorted(exact_moduli)),
        global_minimum,
        (len(hostile_support), len(actual_support), beta_28, kappa_28),
    )
    semantic_sha = hashlib.sha256(repr(semantic).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "PENDING":
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)
    return {
        "carrier_rows": tuple(carrier_rows),
        "totals": totals,
        "minimum_density": minimum_density,
        "bulk_margin": bulk_margin,
        "exact_moduli": tuple(sorted(exact_moduli)),
        "global_minimum": global_minimum,
        "hostile": (len(hostile_support), len(actual_support)),
        "beta_kappa_28": (beta_28, kappa_28),
        "semantic_sha": semantic_sha,
    }


def render(result):
    lines = [
        "LRC14 z229 six-prefix divisor-universal complete-cell support",
        f"source_sha256={lf_sha256(SOURCE)}",
        "scope=the six explicit prefix carriers are reconstructed directly;"
        "THM-3113 screen exhaustion and the cap-227 consequence are not assumed",
    ]
    for row in result["carrier_rows"]:
        (
            body,
            ruler,
            lows,
            cells,
            cell_sha,
            divisor_count,
            full_supports,
            density_closed,
            exact_needed,
            minimum,
        ) = row
        lines.append(
            "CARRIER;E=%s;L=%d;lows=%s;cells=%d;cells_sha256=%s;"
            "divisors=%d;full=%d;density_closed=%d;exact_needed=%d;"
            "minimum=(d:%d,support:%d,kappa:%d,coarse:%d)"
            % (
                body,
                ruler,
                lows,
                cells,
                cell_sha,
                divisor_count,
                full_supports,
                density_closed,
                len(exact_needed),
                minimum[0],
                minimum[1],
                minimum[2],
                minimum[3],
            )
        )
        lines.append(f"EXACT;E={body};rows={exact_needed}")
    lines.extend(
        [
            "bulk=minimum_density:%s;uniform_threshold:d>=25;"
            "threshold_margin:%s"
            % (result["minimum_density"], result["bulk_margin"]),
            "summary=divisor_pairs:%d;density_closed:%d;exact_needed:%d;failures:%d;"
            "exact_moduli:%s"
            % (*result["totals"], result["exact_moduli"]),
            "hostile=d2_same_size_even_support:%d;actual_support:%d;kappa:1"
            % result["hostile"],
            "translated_control=d28_beta:%d;kappa:%d"
            % result["beta_kappa_28"],
            "direction=for every d>1 dividing L, support exceeds ceil(d/7);"
            "THM-2984 therefore makes the projected residual full for every fourth nonaligned tail",
            "aligned_boundary=d=1 is the four-aligned/three-drift branch closed by THM-2928",
            "candidate_boundary=THM-3113 must still prove screen exhaustion, identify these six prefixes,"
            "and close z1=228;this artifact alone does not lower the proved cap 229",
            f"semantic_sha256={result['semantic_sha']}",
            "all_exact_controls=PASS",
        ]
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = run_audit()
    args.output.write_text(render(result))


if __name__ == "__main__":
    main()
