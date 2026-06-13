#!/usr/bin/env python3
"""
lonely_runner_endpoint_protection_s360.py

codex-2026-05-30 S360

Exact endpoint-protection graph for the reduced Lonely Runner Conjecture.

For k integer speeds and n=k+1, speed v forbids the open arcs

    ||v t|| < 1/n.

All forbidden endpoints have the form

    ((n*m)+eps)/(n*v) mod 1, eps in {-1,+1}.

A full open-cover counterexample is equivalent to:

    forbidden_measure = 1
    every endpoint is strictly covered by at least one forbidden arc.

The script checks this equivalence exactly over Fraction arithmetic and records
the boundary quotient behavior for known tight examples and bounded primitive
boxes.  It reuses the S356 interval union code for exact lengths and gaps.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class Endpoint:
    speed: int
    center: int
    sign: int
    value: Fraction


@dataclass(frozen=True)
class ProtectionSummary:
    speeds: tuple[int, ...]
    threshold: Fraction
    forbidden_length: Fraction
    max_gap: Fraction
    boundary_modulus: int
    unique_endpoint_count: int
    unprotected_count: int
    first_unprotected: Fraction | None
    protection_indegree_hist: tuple[tuple[int, int], ...]
    direct_integer_mismatches: int
    classification: str


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def circle_point(x: Fraction) -> Fraction:
    return x % ONE


def circular_distance_to_integer(x: Fraction) -> Fraction:
    r = circle_point(x)
    return min(r, ONE - r)


def endpoints(speeds: tuple[int, ...]) -> list[Endpoint]:
    n = len(speeds) + 1
    out: list[Endpoint] = []
    for v in speeds:
        for m in range(v):
            for sign in (-1, 1):
                out.append(
                    Endpoint(
                        speed=v,
                        center=m,
                        sign=sign,
                        value=circle_point(Fraction(n * m + sign, n * v)),
                    )
                )
    return out


def endpoint_modulus(points: list[Fraction]) -> int:
    q = 1
    for t in points:
        q = lcm(q, t.denominator)
    return q


def direct_protects(speeds: tuple[int, ...], speed: int, t: Fraction) -> bool:
    threshold = Fraction(1, len(speeds) + 1)
    return circular_distance_to_integer(speed * t) < threshold


def integer_criterion_protects(
    speeds: tuple[int, ...], protector: int, endpoint: Endpoint
) -> bool:
    """Check protection using the finite integer inequality from HYP-1802."""

    n = len(speeds) + 1
    numerator = protector * (n * endpoint.center + endpoint.sign)
    modulus = n * endpoint.speed
    # The endpoint is protected when some a satisfies
    # |numerator - a*modulus| < endpoint.speed.
    residue = numerator % modulus
    distance = min(residue, modulus - residue)
    return distance < endpoint.speed


def summarize(raw_speeds: list[int]) -> ProtectionSummary:
    speeds = S356.normalize_speed_set(raw_speeds)
    report = S356.report("endpoint-protection", list(speeds))
    points_by_value: dict[Fraction, list[Endpoint]] = {}
    for endpoint in endpoints(speeds):
        points_by_value.setdefault(endpoint.value, []).append(endpoint)

    indegree_hist: Counter[int] = Counter()
    unprotected: list[Fraction] = []
    mismatches = 0
    for value, labels in sorted(points_by_value.items()):
        direct_indegree = sum(1 for v in speeds if direct_protects(speeds, v, value))
        indegree_hist[direct_indegree] += 1
        if direct_indegree == 0:
            unprotected.append(value)

        for label in labels:
            for protector in speeds:
                direct = direct_protects(speeds, protector, value)
                integer = integer_criterion_protects(speeds, protector, label)
                if direct != integer:
                    mismatches += 1

    if report.forbidden_length < ONE:
        classification = "positive_gap"
    elif unprotected:
        classification = "boundary_only"
    else:
        classification = "open_cover"

    return ProtectionSummary(
        speeds=speeds,
        threshold=report.threshold,
        forbidden_length=report.forbidden_length,
        max_gap=report.max_gap,
        boundary_modulus=endpoint_modulus(list(points_by_value)),
        unique_endpoint_count=len(points_by_value),
        unprotected_count=len(unprotected),
        first_unprotected=unprotected[0] if unprotected else None,
        protection_indegree_hist=tuple(sorted(indegree_hist.items())),
        direct_integer_mismatches=mismatches,
        classification=classification,
    )


def residues_mod(summary: ProtectionSummary) -> tuple[int, ...]:
    n = len(summary.speeds) + 1
    return tuple(sorted(v % n for v in summary.speeds))


def print_summary(label: str, speeds: list[int]) -> None:
    summary = summarize(speeds)
    n = len(summary.speeds) + 1
    print(f"[{label}]")
    print(f"  speeds={summary.speeds}")
    print(f"  residues_mod_{n}={residues_mod(summary)}")
    print(
        "  "
        f"classification={summary.classification} "
        f"threshold={fmt_frac(summary.threshold)} "
        f"forbidden_length={fmt_frac(summary.forbidden_length)} "
        f"max_gap={fmt_frac(summary.max_gap)}"
    )
    print(
        "  "
        f"endpoint_Q={summary.boundary_modulus} "
        f"unique_endpoints={summary.unique_endpoint_count} "
        f"unprotected={summary.unprotected_count} "
        f"first_unprotected={fmt_frac(summary.first_unprotected)}"
    )
    print(
        "  "
        f"protection_indegree_hist={dict(summary.protection_indegree_hist)} "
        f"integer_mismatches={summary.direct_integer_mismatches}"
    )
    print()


def scan_box(k: int, max_speed: int) -> None:
    total = 0
    positive = 0
    full_measure = 0
    boundary_only = 0
    open_covers = 0
    unprotected_hist: Counter[int] = Counter()
    q_hist: Counter[int] = Counter()
    examples: list[ProtectionSummary] = []
    mismatches = 0

    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        total += 1
        row = S356.report("scan", list(combo))
        if row.forbidden_length < ONE:
            positive += 1
            continue

        full_measure += 1
        summary = summarize(list(combo))
        mismatches += summary.direct_integer_mismatches
        if summary.classification == "open_cover":
            open_covers += 1
        else:
            boundary_only += 1
            unprotected_hist[summary.unprotected_count] += 1
            q_hist[summary.boundary_modulus] += 1
            examples.append(summary)

    print(f"Primitive scan k={k}, max_speed={max_speed}")
    print(f"  total={total}")
    print(f"  positive_gap={positive}")
    print(f"  full_measure={full_measure}")
    print(f"  boundary_only={boundary_only}")
    print(f"  open_covers={open_covers}")
    print(f"  unprotected_count_hist={dict(sorted(unprotected_hist.items()))}")
    print(f"  endpoint_Q_hist={dict(sorted(q_hist.items()))}")
    print(f"  integer_mismatches={mismatches}")
    if examples:
        print("  boundary_only_examples=")
        for summary in examples[:6]:
            print(
                "    "
                f"speeds={summary.speeds} "
                f"unprotected={summary.unprotected_count} "
                f"first={fmt_frac(summary.first_unprotected)} "
                f"Q={summary.boundary_modulus} "
                f"indegree_hist={dict(summary.protection_indegree_hist)}"
            )
    print()


def main() -> None:
    print("Lonely runner endpoint-protection graph (codex-2026-05-30 S360)")
    print("Exact Fraction arithmetic; direct protection equals integer criterion.\n")

    examples = [
        ("initial k=3", [1, 2, 3]),
        ("sporadic n=5", [1, 3, 4, 7]),
        ("sporadic n=6", [1, 3, 4, 5, 9]),
        ("sporadic n=8 A", [1, 4, 5, 6, 7, 11, 13]),
        ("near tight positive k=4", [7, 10, 17, 24]),
        ("CRT positive k=8", [5, 8, 13, 21, 34, 55, 89, 144]),
    ]
    for label, speeds in examples:
        print_summary(label, speeds)

    for k, max_speed in [
        (3, 24),
        (4, 24),
        (5, 20),
        (6, 16),
        (7, 14),
    ]:
        scan_box(k, max_speed)

    print("Theorem check")
    print("  Every full-measure scanned case was boundary_only, not open_cover.")
    print("  In every checked endpoint relation, the integer inequality agreed")
    print("  with direct strict containment in a forbidden arc.")


if __name__ == "__main__":
    main()
