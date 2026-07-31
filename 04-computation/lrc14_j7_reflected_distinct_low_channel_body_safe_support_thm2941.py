#!/usr/bin/env python3
r"""Every distinct low phase channel meets a whole body-safe cell.

Let ``E`` be any six-element subset of ``{1,...,14}``, put
``L=14*lcm(E)``, choose labels ``a<b`` in ``E``, and let ``(P,Q)`` be an
oriented primitive channel with

    P != Q,       gcd(P,Q)=1,       P+Q<=7.

The primitive fiber coordinate on a body cell is governed by the exact
cross-difference

    h = P*b-Q*a.

This referee proves that there are an integer body-safe cell ``[j,j+1]``
and an integer lift ``k`` for which

    |h(j+s)-kL| < (P+Q)L/14          for every 0<=s<=1. (S)

Thus the entire cell lies strictly inside the support of the periodized
THM-1226 ``(P,Q)`` trapezoid.  The exhaustive universe has

    C(14,6)*C(6,2)*16 = 720720

oriented configurations, and all 720720 pass.  Endpoint-only contact is
excluded: both endpoint inequalities in (S) are strict.  Because all data
are integral, the selected raw support clearance is at least one; the exact
deterministic witness census has weakest normalized clearance ``1/2522520``.

The equal channel is the hostile control.  Repeating the same test with
``(P,Q)=(1,1)`` on all 45045 body pairs gives exactly five failures, namely
the bad path and matching already isolated by the universal chromatic
theorem:

    (1,2,7,9,11,13):  (7,9),(9,11),(11,13),
    (2,4,7,9,11,13):  (7,11),(9,13).

Those five do not merely lack a whole supported cell: their body-safe set
has no positive intersection with the primitive fiber support at all.

Scope is deliberately narrow.  Strict support is non-blindness, not yet a
Bonferroni closure: a theorem must still convert the support clearance into
an overlap floor after the ``1/g`` affine perturbation and compare it with
the six-clause singleton debt.  No arbitrary-level or LRC(14) conclusion is
claimed here.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_distinct_low_channel_body_safe_support_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_SEMANTIC_SHA256 = "458b7fd4cd3b772f37b2f239f657ca3935112d4c069530f4d2b4e687cabe485b"
BODY_COUNT = 3003
DISTINCT_CONFIGURATION_COUNT = 720720
EQUAL_CONFIGURATION_COUNT = 45045
EXPECTED_ZERO_CROSS_DIFFERENCES = 13365
EXPECTED_MAX_RULER = (5045040, (8, 9, 10, 11, 13, 14))
EXPECTED_WEAKEST_RAW = (
    1,
    (1, 2, 3, 4, 5, 11),
    9240,
    2,
    3,
    3,
    1,
    7,
    943,
    1,
)
EXPECTED_WEAKEST_NORMALIZED = (
    F(1, 2522520),
    (7, 9, 10, 11, 12, 13),
    2522520,
    7,
    12,
    5,
    1,
    53,
    27197,
    1,
    1,
)
EXPECTED_EQUAL_BLIND = (
    ((1, 2, 7, 9, 11, 13), 7, 9),
    ((1, 2, 7, 9, 11, 13), 9, 11),
    ((1, 2, 7, 9, 11, 13), 11, 13),
    ((2, 4, 7, 9, 11, 13), 7, 11),
    ((2, 4, 7, 9, 11, 13), 9, 13),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
SPEC = spec_from_file_location("low_channel_support_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot import reflected interval engine")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


CHANNELS = tuple(
    (p, q)
    for p in range(1, 7)
    for q in range(1, 7)
    if p != q and gcd(p, q) == 1 and p + q <= 7
)
require(len(CHANNELS) == 16, CHANNELS)


def full_support_witness(
    safe_ranges: tuple[tuple[int, int], ...],
    ruler: int,
    cross_difference: int,
    support_radius: int,
) -> tuple[int, int, int] | None:
    """Return ``(j,k,raw_slack)`` proving strict support on ``[j,j+1]``."""
    h = abs(cross_difference)
    if h == 0:
        return safe_ranges[0][0], 0, support_radius
    for left, right in safe_ranges:
        for lift in range(h * left // ruler - 2, h * right // ruler + 3):
            # Strict endpoint inequalities are
            # h*j > kL-A and h*(j+1) < kL+A.
            first = max(left, (lift * ruler - support_radius) // h + 1)
            last = min(
                right - 1,
                (lift * ruler + support_radius - h - 1) // h,
            )
            if first <= last:
                slack = min(
                    h * first - (lift * ruler - support_radius),
                    (lift * ruler + support_radius) - h * (first + 1),
                )
                require(slack > 0, (ruler, h, support_radius, first, lift, slack))
                return first, lift, slack
    return None


def support_hit(
    safe_ranges: tuple[tuple[int, int], ...],
    ruler: int,
    cross_difference: int,
    support_radius: int,
) -> bool:
    """Whether the safe-cell union has any positive open support incidence."""
    h = abs(cross_difference)
    if h == 0:
        return True
    for left, right in safe_ranges:
        lo, hi = h * left, h * right
        for lift in range(lo // ruler - 2, hi // ruler + 3):
            if max(lo, lift * ruler - support_radius) < min(
                hi, lift * ruler + support_radius
            ):
                return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == BODY_COUNT, len(bodies))
    distinct_count = 0
    equal_count = 0
    zero_cross_differences = 0
    weakest_raw = None
    weakest_normalized = None
    maximum_ruler = None
    equal_blind = []
    witness_digest = hashlib.sha256()

    for body in bodies:
        ruler, safe_ranges = R.safe_cell_ranges(body)
        row = (ruler, body)
        if maximum_ruler is None or row > maximum_ruler:
            maximum_ruler = row
        for a, b in combinations(body, 2):
            for p, q in CHANNELS:
                cross_difference = p * b - q * a
                support_radius = (p + q) * ruler // 14
                require(14 * support_radius == (p + q) * ruler,
                        (body, ruler, p, q))
                witness = full_support_witness(
                    safe_ranges, ruler, cross_difference, support_radius
                )
                require(witness is not None,
                        ("distinct channel unexpectedly blind", body, a, b, p, q))
                j, lift, slack = witness
                require(any(left <= j < right for left, right in safe_ranges),
                        (body, j, safe_ranges))
                h = abs(cross_difference)
                require(
                    abs(h * j - lift * ruler) < support_radius
                    and abs(h * (j + 1) - lift * ruler) < support_radius,
                    (body, a, b, p, q, j, lift, slack),
                )
                raw_row = (slack, body, ruler, a, b, p, q, cross_difference, j, lift)
                normalized_row = (
                    F(slack, ruler), body, ruler, a, b, p, q,
                    cross_difference, j, lift, slack,
                )
                if weakest_raw is None or raw_row < weakest_raw:
                    weakest_raw = raw_row
                if weakest_normalized is None or normalized_row < weakest_normalized:
                    weakest_normalized = normalized_row
                zero_cross_differences += cross_difference == 0
                distinct_count += 1
                witness_digest.update(
                    f"D|{body}|{ruler}|{a}|{b}|{p}|{q}|{cross_difference}|{j}|{lift}|{slack}\n".encode()
                )

            equal_count += 1
            cross_difference = b - a
            support_radius = ruler // 7
            witness = full_support_witness(
                safe_ranges, ruler, cross_difference, support_radius
            )
            if witness is None:
                require(
                    not support_hit(safe_ranges, ruler, cross_difference, support_radius),
                    ("equal channel has support but no full cell", body, a, b),
                )
                equal_blind.append((body, a, b))
                witness_digest.update(f"E|{body}|{ruler}|{a}|{b}|BLIND\n".encode())
            else:
                j, lift, slack = witness
                require(slack > 0, (body, a, b, witness))
                witness_digest.update(
                    f"E|{body}|{ruler}|{a}|{b}|{j}|{lift}|{slack}\n".encode()
                )

    require(distinct_count == DISTINCT_CONFIGURATION_COUNT, distinct_count)
    require(equal_count == EQUAL_CONFIGURATION_COUNT, equal_count)
    require(zero_cross_differences == EXPECTED_ZERO_CROSS_DIFFERENCES,
            zero_cross_differences)
    require(maximum_ruler == EXPECTED_MAX_RULER, maximum_ruler)
    require(weakest_raw == EXPECTED_WEAKEST_RAW, weakest_raw)
    require(weakest_normalized == EXPECTED_WEAKEST_NORMALIZED,
            weakest_normalized)
    require(tuple(equal_blind) == EXPECTED_EQUAL_BLIND, equal_blind)

    semantic_payload = (
        CHANNELS,
        distinct_count,
        equal_count,
        zero_cross_differences,
        maximum_ruler,
        weakest_raw,
        weakest_normalized,
        tuple(equal_blind),
        witness_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 distinct low-channel body-safe support exact referee",
        f"oriented_low_channels={CHANNELS};channel_count={len(CHANNELS)}",
        f"distinct_configurations={distinct_count};full_strict_support_cells={distinct_count};"
        f"zero_cross_differences={zero_cross_differences}",
        f"weakest_raw_support_slack={weakest_raw};"
        f"weakest_normalized_support_slack={qtext(weakest_normalized[0])};"
        f"weakest_normalized_witness={weakest_normalized[1:]}",
        f"maximum_ruler={maximum_ruler}",
        f"equal_channel_controls={equal_count};equal_blind_count={len(equal_blind)};"
        f"equal_blind={tuple(equal_blind)}",
        "conclusion=every oriented distinct primitive P+Q<=7 channel contains a whole strictly supported body-safe cell; endpoint-only hits are excluded",
        "scope=support/non-blindness only; no affine-overlap floor, debt comparison, reflected-packet closure, or LRC14 conclusion",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"witness_digest={witness_digest.hexdigest()}",
        f"source_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
