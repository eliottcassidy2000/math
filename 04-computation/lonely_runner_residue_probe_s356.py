#!/usr/bin/env python3
"""
lonely_runner_residue_probe_s356.py

codex-2026-05-30 S356

Exact rational interval probe for reading the lonely runner conjecture as a
quotient-gap problem.  For k nonzero integer speeds, the reduced lonely runner
threshold is 1/(k+1).  Runner v forbids the interval union

    dist(v t, Z) < 1/(k+1).

After pulling this back to t in R/Z, every boundary is a rational point with
denominator dividing (k+1) * v.  The conjectural witness is a gap in the union
of these pulled-back forbidden fibers.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm
import random


ONE = Fraction(1, 1)


@dataclass(frozen=True)
class GapReport:
    label: str
    speeds: tuple[int, ...]
    threshold: Fraction
    components: int
    forbidden_length: Fraction
    gap_count: int
    max_gap: Fraction
    witness: Fraction | None
    boundary_witness_count: int
    boundary_witness: Fraction | None
    boundary_modulus: int
    max_gap_boundary_residues: tuple[int, int] | None


def normalize_speed_set(speeds: list[int]) -> tuple[int, ...]:
    cleaned = sorted({abs(v) for v in speeds if v})
    if not cleaned:
        raise ValueError("speed set must contain at least one nonzero speed")
    g = cleaned[0]
    for v in cleaned[1:]:
        g = gcd(g, v)
    return tuple(v // g for v in cleaned)


def split_circle_interval(lo: Fraction, hi: Fraction) -> list[tuple[Fraction, Fraction]]:
    """Return pieces of [lo, hi] modulo 1 as intervals in [0,1]."""
    if hi - lo >= ONE:
        return [(Fraction(0), ONE)]

    # Move lo into [0,1).  Python // is floor division for Fractions.
    shift = lo // ONE
    lo -= shift
    hi -= shift

    if lo < 0:
        lo += ONE
        hi += ONE

    if hi <= ONE:
        return [(lo, hi)]
    return [(lo, ONE), (Fraction(0), hi - ONE)]


def forbidden_intervals(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    k = len(speeds)
    threshold = Fraction(1, k + 1)
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        radius = threshold / v
        for a in range(v):
            center = Fraction(a, v)
            intervals.extend(split_circle_interval(center - radius, center + radius))
    return intervals


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged: list[tuple[Fraction, Fraction]] = []
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        # The forbidden intervals are open.  Touching endpoints are valid
        # lonely boundary candidates, so only strictly overlapping intervals
        # should be merged.
        if lo < cur_hi:
            if hi > cur_hi:
                cur_hi = hi
        else:
            merged.append((cur_lo, cur_hi))
            cur_lo, cur_hi = lo, hi
    merged.append((cur_lo, cur_hi))

    return merged


def circular_gaps(components: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not components:
        return [(Fraction(0), ONE)]
    if len(components) == 1 and components[0] == (Fraction(0), ONE):
        return []

    comps = sorted(components)
    gaps: list[tuple[Fraction, Fraction]] = []
    for (_, hi), (next_lo, _) in zip(comps, comps[1:]):
        if next_lo > hi:
            gaps.append((hi, next_lo))
    last_hi = comps[-1][1]
    first_lo = comps[0][0]
    wrap_gap = first_lo + ONE - last_hi
    if wrap_gap > 0 and wrap_gap < ONE:
        gaps.append((last_hi, first_lo + ONE))
    return gaps


def midpoint_on_circle(gap: tuple[Fraction, Fraction]) -> Fraction:
    mid = (gap[0] + gap[1]) / 2
    if mid >= ONE:
        mid -= ONE
    return mid


def boundary_modulus(components: list[tuple[Fraction, Fraction]]) -> int:
    q = 1
    for lo, hi in components:
        q = lcm(q, lo.denominator, hi.denominator)
    return q


def residue_pair(q: int, gap: tuple[Fraction, Fraction]) -> tuple[int, int]:
    lo = gap[0]
    hi = gap[1]
    if hi >= ONE:
        hi -= ONE
    return (int((lo * q) % q), int((hi * q) % q))


def circle_point(x: Fraction) -> Fraction:
    return x % ONE


def is_lonely_witness(speeds: tuple[int, ...], t: Fraction) -> bool:
    threshold = Fraction(1, len(speeds) + 1)
    for v in speeds:
        r = circle_point(v * t)
        if min(r, ONE - r) < threshold:
            return False
    return True


def report(label: str, raw_speeds: list[int]) -> GapReport:
    speeds = normalize_speed_set(raw_speeds)
    raw_intervals = forbidden_intervals(speeds)
    components = merge_intervals(raw_intervals)
    gaps = circular_gaps(components)
    max_gap_tuple = max(gaps, key=lambda g: g[1] - g[0]) if gaps else None
    q = boundary_modulus(components)
    boundary_candidates = sorted(
        {
            circle_point(point)
            for interval in raw_intervals
            for point in interval
        }
    )
    boundary_witnesses = [t for t in boundary_candidates if is_lonely_witness(speeds, t)]
    return GapReport(
        label=label,
        speeds=speeds,
        threshold=Fraction(1, len(speeds) + 1),
        components=len(components),
        forbidden_length=sum(hi - lo for lo, hi in components),
        gap_count=len(gaps),
        max_gap=(max_gap_tuple[1] - max_gap_tuple[0]) if max_gap_tuple else Fraction(0),
        witness=midpoint_on_circle(max_gap_tuple) if max_gap_tuple else None,
        boundary_witness_count=len(boundary_witnesses),
        boundary_witness=boundary_witnesses[0] if boundary_witnesses else None,
        boundary_modulus=q,
        max_gap_boundary_residues=residue_pair(q, max_gap_tuple) if max_gap_tuple else None,
    )


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def deterministic_samples() -> list[tuple[str, list[int]]]:
    rng = random.Random(356)
    samples: list[tuple[str, list[int]]] = []
    for k in range(2, 10):
        samples.append((f"initial segment k={k}", list(range(1, k + 1))))
    samples.extend(
        [
            ("powers of 2 k=6", [1, 2, 4, 8, 16, 32]),
            ("odd arithmetic k=6", [1, 3, 5, 7, 9, 11]),
            ("small primes k=8", [2, 3, 5, 7, 11, 13, 17, 19]),
            ("mixed CRT k=8", [5, 8, 13, 21, 34, 55, 89, 144]),
        ]
    )
    for k, limit in [(4, 24), (5, 32), (7, 50), (9, 70)]:
        speeds = rng.sample(range(1, limit + 1), k)
        samples.append((f"random356 k={k}", speeds))
    return samples


def main() -> None:
    print("Lonely runner quotient-gap probe (codex-2026-05-30 S356)")
    print("Reduced form: k moving runners, threshold = 1/(k+1).")
    print("All interval arithmetic is exact over Fraction.\n")

    rows = [report(label, speeds) for label, speeds in deterministic_samples()]
    for row in rows:
        ratio = row.max_gap / row.threshold if row.threshold else Fraction(0)
        print(f"[{row.label}]")
        print(f"  speeds={row.speeds}")
        print(
            "  "
            f"threshold={fmt_frac(row.threshold)} "
            f"components={row.components} "
            f"forbidden_length={fmt_frac(row.forbidden_length)} "
            f"gap_count={row.gap_count}"
        )
        print(
            "  "
            f"max_gap={fmt_frac(row.max_gap)} ({fmt_float(row.max_gap)}) "
            f"max_gap/threshold={fmt_float(ratio)} "
            f"witness={fmt_frac(row.witness)}"
        )
        print(
            "  "
            f"boundary_witness_count={row.boundary_witness_count} "
            f"boundary_witness={fmt_frac(row.boundary_witness)}"
        )
        print(
            "  "
            f"boundary_modulus={row.boundary_modulus} "
            f"max_gap_boundary_residues={row.max_gap_boundary_residues}"
        )
        print()

    no_positive_gap = [row.label for row in rows if row.max_gap == 0]
    no_witness = [
        row.label
        for row in rows
        if row.max_gap == 0 and row.boundary_witness_count == 0
    ]
    print("Summary")
    print(f"  sample_count={len(rows)}")
    print(f"  positive_gap_count={len(rows) - len(no_positive_gap)}")
    print(f"  zero_positive_gap_labels={no_positive_gap}")
    print(f"  positive_or_boundary_witness_count={len(rows) - len(no_witness)}")
    print(f"  no_witness_labels={no_witness}")
    best = max(rows, key=lambda r: r.max_gap / r.threshold)
    tight = min((r for r in rows if r.max_gap > 0), key=lambda r: r.max_gap / r.threshold)
    print(
        "  widest_relative_gap="
        f"{best.label} ratio={fmt_float(best.max_gap / best.threshold)} "
        f"witness={fmt_frac(best.witness)}"
    )
    print(
        "  tightest_positive_gap="
        f"{tight.label} ratio={fmt_float(tight.max_gap / tight.threshold)} "
        f"boundary_residues={tight.max_gap_boundary_residues}"
    )


if __name__ == "__main__":
    main()
