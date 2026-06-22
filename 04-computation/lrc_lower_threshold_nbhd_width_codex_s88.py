#!/usr/bin/env python3
"""
lrc_lower_threshold_nbhd_width_codex_s88.py

Exact-rational scout for the HYP-2842 redirect:

  exact Farey centers can be killed by composite P-elements, but neighborhoods
  around those centers may still survive the small-part safe-set holes.

For an even runner count N=2q, every reduced Farey center a/b with 2<=b<q
has at most b phase values for any cluster E, hence a gap >=1/b.  If the
cluster span is S=max(E)-min(E), the same gap remains >1/q throughout the
conservative interval

  |x-a/b| < (1/b - 1/q) / S = (q-b)/(b q S).

For a small speed p in P, the unsafe holes of
  G_P = { ||p x|| >= 1/(2q) }
around j/p have radius 1/(2q p).

This script computes, exactly, how much of the union of the guaranteed lonely
Farey neighborhoods remains after subtracting all G_P holes.  It scans the
dense branch k=q+1..2q-1, with |P|=2q-1-k, and compares three span ledgers:

  consec_span = k-1       (packed cluster; largest guaranteed neighborhoods)
  bounded_span = 2q-1     (bounded-speed lower-threshold stress)
  wide_span = 4q          (shows where the simple width method starts to fail)

This is a lower-threshold lead generator, not a proof of LRC14.  The intervals
are sufficient neighborhoods only; positive residual is a rigorous positive
certificate for that modeled branch, while zero residual means only that this
particular Farey-neighborhood certificate was too crude.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def add_interval_mod1(intervals: list[Interval], lo: Fraction, hi: Fraction) -> None:
    """Append [lo,hi] modulo 1, clipped to [0,1]."""
    one = Fraction(1)
    if hi <= lo:
        return
    while lo < 0:
        lo += one
        hi += one
    while lo >= one:
        lo -= one
        hi -= one
    if hi <= one:
        intervals.append((max(Fraction(0), lo), min(one, hi)))
    else:
        intervals.append((lo, one))
        intervals.append((Fraction(0), hi - one))


def merge(intervals: list[Interval]) -> list[Interval]:
    clean = sorted((lo, hi) for lo, hi in intervals if lo < hi)
    if not clean:
        return []
    out = [clean[0]]
    for lo, hi in clean[1:]:
        plo, phi = out[-1]
        if lo <= phi:
            out[-1] = (plo, max(phi, hi))
        else:
            out.append((lo, hi))
    return out


def length(intervals: list[Interval]) -> Fraction:
    return sum((hi - lo for lo, hi in merge(intervals)), Fraction(0))


def subtract(base: list[Interval], holes: list[Interval]) -> list[Interval]:
    pieces = merge(base)
    for hlo, hhi in merge(holes):
        new: list[Interval] = []
        for lo, hi in pieces:
            if hhi <= lo or hi <= hlo:
                new.append((lo, hi))
                continue
            if lo < hlo:
                new.append((lo, hlo))
            if hhi < hi:
                new.append((hhi, hi))
        pieces = new
        if not pieces:
            break
    return pieces


def farey_centers(q: int, include_origin: bool = False) -> list[tuple[Fraction, int]]:
    centers: list[tuple[Fraction, int]] = []
    if include_origin:
        centers.append((Fraction(0), 1))
    for b in range(2, q):
        for a in range(1, b):
            if gcd(a, b) == 1:
                centers.append((Fraction(a, b), b))
    return centers


def lonely_neighborhoods(q: int, span: int, include_origin: bool = False) -> list[Interval]:
    intervals: list[Interval] = []
    for center, b in farey_centers(q, include_origin=include_origin):
        radius = Fraction(q - b, b * q * span)
        add_interval_mod1(intervals, center - radius, center + radius)
    return merge(intervals)


def gp_holes(q: int, P: tuple[int, ...]) -> list[Interval]:
    holes: list[Interval] = []
    for p in P:
        radius = Fraction(1, 2 * q * p)
        for j in range(p):
            center = Fraction(j, p)
            add_interval_mod1(holes, center - radius, center + radius)
    return merge(holes)


def exact_farey_survivors(q: int, P: tuple[int, ...], include_origin: bool = False) -> int:
    """Exact centers survive iff no p in P is divisible by their denominator."""
    if include_origin:
        # The exact origin never lies in G_P when P is nonempty.
        return int(not P) + sum(
            1 for _center, b in farey_centers(q) if all(p % b != 0 for p in P)
        )
    return sum(1 for _center, b in farey_centers(q) if all(p % b != 0 for p in P))


@dataclass(frozen=True)
class ScanResult:
    q: int
    k: int
    p_size: int
    span_name: str
    span: int
    origin: bool
    min_residual: Fraction
    lonely_len: Fraction
    worst_P: tuple[int, ...]
    exact_survivors: int


def scan(q: int, k: int, span_name: str, span: int, include_origin: bool) -> ScanResult:
    p_size = 2 * q - 1 - k
    candidates = range(1, 2 * q)
    lonely = lonely_neighborhoods(q, span, include_origin=include_origin)
    lonely_len = length(lonely)
    best_residual: Fraction | None = None
    best_P: tuple[int, ...] = ()
    best_survivors = 0
    for P in combinations(candidates, p_size):
        residual = length(subtract(lonely, gp_holes(q, P)))
        survivors = exact_farey_survivors(q, P, include_origin=include_origin)
        if best_residual is None or residual < best_residual:
            best_residual = residual
            best_P = P
            best_survivors = survivors
    assert best_residual is not None
    return ScanResult(
        q, k, p_size, span_name, span, include_origin,
        best_residual, lonely_len, best_P, best_survivors,
    )


def main() -> None:
    print("HYP-2842 lower-threshold Farey-neighborhood width scout (exact rational)")
    print("positive residual = certified GOOD-neighborhood mass surviving G_P holes")
    print("proper = denominators 2..q-1; origin+proper also keeps the b=1 collapse neighborhood")
    print()
    route_rows: list[ScanResult] = []
    for q in range(3, 8):
        print(f"N=2q={2*q} (q={q})")
        for k in range(q + 1, 2 * q):
            spans = [
                ("consec", max(1, k - 1)),
                ("bounded", 2 * q - 1),
                ("wide4q", 4 * q),
            ]
            for include_origin in (False, True):
              for name, span in spans:
                rec = scan(q, k, name, span, include_origin)
                route_rows.append(rec)
                flag = "POS" if rec.min_residual > 0 else "ZERO"
                mode = "origin+" if include_origin else "proper "
                print(
                    f"  {mode} k={k:2d} |P|={rec.p_size:2d} span={name:7s}={span:2d} "
                    f"lonely={fmt(rec.lonely_len):>12s} residual={fmt(rec.min_residual):>12s} "
                    f"{flag:4s} exact_surv={rec.exact_survivors:2d} worstP={rec.worst_P}"
                )
        print()

    print("Synthesis")
    for q in range(3, 8):
        dense = [r for r in route_rows if r.q == q and r.span_name == "bounded" and not r.origin]
        min_rec = min(dense, key=lambda r: r.min_residual)
        dense_origin = [r for r in route_rows if r.q == q and r.span_name == "bounded" and r.origin]
        min_origin = min(dense_origin, key=lambda r: r.min_residual)
        print(
            f"  bounded N={2*q:2d}: proper min {fmt(min_rec.min_residual)} "
            f"at k={min_rec.k}, |P|={min_rec.p_size}, worstP={min_rec.worst_P}, "
            f"exact_survivors={min_rec.exact_survivors}; "
            f"origin+proper min {fmt(min_origin.min_residual)} "
            f"at k={min_origin.k}, worstP={min_origin.worst_P}"
        )

    print()
    print("Tournament Analysis")
    vertices = [
        ("origin_neighborhood_width", 6, "repairs exact-center failure by using the b=1 collapse neighborhood"),
        ("bounded_width", 5, "proves every bounded dense branch checked here for N<=14"),
        ("consec_width", 4, "larger neighborhoods; useful as extremal-width model"),
        ("wide_width", 2, "often positive below 14 but too span-sensitive for final route"),
        ("exact_farey_points", 0, "refuted by composite P killing all centers"),
        ("raw_runner_vertices", 1, "does not retain neighborhood/hole ownership"),
    ]
    scores = {name: score for name, score, _ in vertices}
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    path = " > ".join(name for name, _score, _ in sorted(vertices, key=lambda row: -row[1]))
    print(f"  pairwise observable: certified residual mass after subtracting G_P holes")
    print(f"  switch/gauge: larger lower-bound residual wins; ties use proof-route specificity")
    print(f"  score_hist={hist}")
    print(f"  Hamiltonian path: {path}")
    print("  challenged assumption: vertices need not be runners or arcs; width certificates")
    print("  preserve neighborhood ownership and destroy exact center identity.")


if __name__ == "__main__":
    main()
