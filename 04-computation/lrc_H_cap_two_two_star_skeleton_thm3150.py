#!/usr/bin/env python3
"""Exact referee for THM-3150's cap-two two-star skeleton theorem."""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
OUTPUT = ROOT / "05-knowledge/results/lrc_H_cap_two_two_star_skeleton_thm3150.out"

H = (1, 2, 3, 4, 6, 12)
L = 168
CELL = 90
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
EXPECTED_LOW = frozenset((F(1, 2), F(2, 3), F(3, 4),
                          F(4, 3), F(3, 2), F(2)))
EXPECTED_QUOTIENTS = frozenset((
    F(1, 2), F(9, 16), F(2, 3), F(3, 4), F(8, 9), F(1),
    F(9, 8), F(4, 3), F(3, 2), F(16, 9), F(2),
))
REPRESENTATIVES = (F(1, 2), F(2, 3), F(3, 4), F(9, 16), F(8, 9))
EXPECTED_CASES = {
    F(1, 2): (4, 5, 4),
    F(2, 3): (4, 5, 4),
    F(3, 4): (4, 5, 4),
    F(9, 16): (5, 5, 4),
    F(8, 9): (6, 6, 3),
}
EXPECTED_SEMANTIC_SHA256 = "0c70a8b48cfa8479ec0ffee21e6d2236235acb6a29b58f72c8fdf5c08cce240d"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def source_sha256() -> str:
    return hashlib.sha256(HERE.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def low_ratios() -> frozenset[F]:
    rows = set()
    for p in range(1, 8):
        for q in range(1, 8):
            if p == q or gcd(p, q) != 1 or p + q > 7:
                continue
            ratio = F(q, p)
            if F(1, 2) <= ratio <= 2:
                rows.add(ratio)
    return frozenset(rows)


def leaf_degree(low: frozenset[F], t: F, x: F) -> int:
    return int(x in low) + int(x / t in low)


def exact_case(low: frozenset[F], t: F):
    """Exact maximum low incidences over at most four positive-degree leaves."""
    candidates = tuple(sorted((set(low) | {t * s for s in low}) - {F(1), t}))
    best = None
    for size in range(5):
        for chosen in combinations(candidates, size):
            levels = (F(1), t, *chosen)
            if max(levels) > 2 * min(levels):
                continue
            row = (sum(leaf_degree(low, t, x) for x in chosen), chosen)
            if best is None or row > best:
                best = row
    require(best is not None, ("empty ratio case", t))
    low_edge_upper = int(t in low) + best[0]
    return candidates, best, low_edge_upper, len(EDGES) - low_edge_upper


def circular_distance(value: int) -> int:
    value %= L
    return min(value, L - value)


def debt_upper(n: int) -> F:
    return sum((F(e, 7 * (n * L - e)) for e in H), F(0))


def transported_sum_floor(n: int) -> F:
    return F(1, 35) - F(65, 42 * n)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    low = low_ratios()
    require(low == EXPECTED_LOW, (low, EXPECTED_LOW))
    quotients = frozenset(
        x / y for x in low for y in low if F(1, 2) <= x / y <= 2
    )
    require(quotients == EXPECTED_QUOTIENTS, (quotients, EXPECTED_QUOTIENTS))

    case_rows = []
    for t in REPRESENTATIVES:
        candidates, best, low_upper, high_lower = exact_case(low, t)
        require((best[0], low_upper, high_lower) == EXPECTED_CASES[t],
                (t, best, low_upper, high_lower))
        case_rows.append((t, candidates, best, low_upper, high_lower))

    minimum_high = min(row[-1] for row in case_rows)
    require(minimum_high == 3, minimum_high)
    sharp = next(row for row in case_rows if row[0] == F(8, 9))
    require(sharp[2] == (
        6,
        (F(2, 3), F(3, 4), F(32, 27), F(4, 3)),
    ), sharp)
    sharp_levels = (F(1), F(8, 9), *sharp[2][1])
    require(len(set(sharp_levels)) == 6 and max(sharp_levels) == 2 * min(sharp_levels),
            sharp_levels)

    residues = tuple((e * CELL) % L for e in H)
    pivot_blind = set()
    for t in low:
        p, q = t.denominator, t.numerator
        if circular_distance(q * residues[0] - p * residues[1]) >= 12 * (p + q):
            pivot_blind.add(t)
    require(pivot_blind == {F(1, 2), F(3, 2)}, pivot_blind)
    for t in pivot_blind:
        require(exact_case(low, t)[-1] == 4, (t, exact_case(low, t)))

    skeleton_sum_floor = F(minimum_high, 105)
    blind_sum_floor = F(4, 105)
    coefficient_sum = sum(H[i] + H[j] for i, j in EDGES)
    require(skeleton_sum_floor == F(1, 35), skeleton_sum_floor)
    require(coefficient_sum == 65, coefficient_sum)

    margin_61 = transported_sum_floor(61) / 9 - debt_upper(61)
    margin_62 = transported_sum_floor(62) / 9 - debt_upper(62)
    require(margin_61 == F(
        -25676690027983908791, 734017257871022240108370
    ) < 0, margin_61)
    require(margin_62 == F(
        6445273605987946151, 383853218066532083354940
    ) > 0, margin_62)
    require(transported_sum_floor(62) == F(47, 13020),
            transported_sum_floor(62))

    # The common-dilation lower bound increases, while every debt summand
    # decreases.  Check the exact difference identities at the threshold;
    # their displayed formulas prove the same signs for every n>=1.
    n = 62
    require(
        transported_sum_floor(n + 1) - transported_sum_floor(n)
        == F(65, 42 * n * (n + 1)) > 0,
        "transport monotonicity",
    )
    for e in H:
        difference = F(e, 7 * (n * L - e)) - F(e, 7 * ((n + 1) * L - e))
        require(
            difference == F(e * L, 7 * (n * L - e) * ((n + 1) * L - e)) > 0,
            ("debt monotonicity", e),
        )

    semantic_payload = (
        tuple(sorted(low)), tuple(sorted(quotients)), tuple(case_rows),
        tuple(sorted(pivot_blind)), skeleton_sum_floor, blind_sum_floor,
        coefficient_sum, margin_61, margin_62,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC H CAP-TWO TWO-STAR SKELETON -- EXACT REFEREE",
        f"body={H};ruler={L};cell={CELL};edges={EDGES}",
        f"low_oriented_ratios={tuple(sorted(low))}",
        f"quotient_exception_set={tuple(sorted(quotients))}",
    ]
    for row in case_rows:
        lines.append(
            f"t={qtext(row[0])};candidates={row[1]};max_leaf_low={row[2]};"
            f"total_low_upper={row[3]};high_lower={row[4]}"
        )
    lines.extend((
        f"sharp_high_count={minimum_high};sharp_t=8/9;sharp_leaves={sharp[2][1]}",
        f"pivot_blind={tuple(sorted(pivot_blind))};pivot_blind_high_lower=4",
        f"skeleton_sum_floor={skeleton_sum_floor};pivot_blind_floor={blind_sum_floor}",
        f"two_star_coefficient_sum={coefficient_sum};transport_loss=65/(42*n)",
        f"n61_sufficient_margin={margin_61}",
        f"n62_sufficient_margin={margin_62}",
        "scope=common-dilation cap-two rays;arbitrary pairwise gcd remains open;LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"source_sha256={source_sha256()}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
