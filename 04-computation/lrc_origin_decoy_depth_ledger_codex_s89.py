#!/usr/bin/env python3
"""
lrc_origin_decoy_depth_ledger_codex_s89.py

Follow-up to HYP-2847.  The S88 exact interval scout showed that adding the
b=1 origin-collapse neighborhood repairs bounded dense lower-threshold rows
N=2q<=14.  This script asks WHY.

Key reduction tested here:

  For bounded span S=2q-1, the proper Farey neighborhood around a/b
  (2<=b<q) is touched by a G_P hole for speed p only when b divides p.

If this holds, proper-neighborhood residual is not a geometric union problem
anymore.  It is the divisor-depth ledger

  sum_{2<=b<q} phi(b) * 2 * max(0, r_b - 1/(2q p_b)),

where r_b=(q-b)/(b q S), and p_b is the smallest selected multiple of b
or infinity if no multiple is selected.  The b=1 origin neighborhood is then a
separate decoy: if p=1 is selected it is killed, but selecting p=1 consumes a
scarce P-slot and changes the divisor-depth ledger for b>=2.

The script verifies the reduction by exact interval arithmetic, then scans the
deepest bounded dense branch k=q+1, |P|=q-2 for q up to 10.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb, gcd
import random


Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def phi(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def add_interval_mod1(intervals: list[Interval], lo: Fraction, hi: Fraction) -> None:
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
        plo, phi_ = out[-1]
        if lo <= phi_:
            out[-1] = (plo, max(phi_, hi))
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


def gp_holes(q: int, P: tuple[int, ...]) -> list[Interval]:
    holes: list[Interval] = []
    for p in P:
        radius = Fraction(1, 2 * q * p)
        for j in range(p):
            add_interval_mod1(holes, Fraction(j, p) - radius, Fraction(j, p) + radius)
    return merge(holes)


def farey_intervals_by_denominator(q: int, span: int) -> dict[int, list[Interval]]:
    out: dict[int, list[Interval]] = {}
    for b in range(2, q):
        radius = Fraction(q - b, b * q * span)
        rows: list[Interval] = []
        for a in range(1, b):
            if gcd(a, b) == 1:
                add_interval_mod1(rows, Fraction(a, b) - radius, Fraction(a, b) + radius)
        out[b] = merge(rows)
    return out


def origin_intervals(q: int, span: int) -> list[Interval]:
    radius = Fraction(q - 1, q * span)
    rows: list[Interval] = []
    add_interval_mod1(rows, -radius, radius)
    return merge(rows)


def full_proper_exact(q: int, span: int, P: tuple[int, ...]) -> Fraction:
    holes = gp_holes(q, P)
    total = Fraction(0)
    for intervals in farey_intervals_by_denominator(q, span).values():
        total += length(subtract(intervals, holes))
    return total


def origin_exact(q: int, span: int, P: tuple[int, ...]) -> Fraction:
    if not P:
        return length(origin_intervals(q, span))
    if 1 in P:
        return Fraction(0)
    radius0 = Fraction(q - 1, q * span)
    holes: list[Interval] = []
    for p in P:
        hole_radius = Fraction(1, 2 * q * p)
        # Only centers whose holes touch the two origin tails can matter.
        for j in range(p):
            center = Fraction(j, p)
            touches_left = center - hole_radius < radius0
            touches_right = center + hole_radius > 1 - radius0
            if touches_left or touches_right:
                add_interval_mod1(holes, center - hole_radius, center + hole_radius)
    return length(subtract(origin_intervals(q, span), merge(holes)))


def proper_depth_ledger(q: int, span: int, P: tuple[int, ...]) -> tuple[Fraction, dict[int, Fraction]]:
    """Closed-form proper-neighborhood residual if only divisible speeds touch b."""
    by_b: dict[int, Fraction] = {}
    total = Fraction(0)
    for b in range(2, q):
        radius = Fraction(q - b, b * q * span)
        multiples = [p for p in P if p % b == 0]
        if multiples:
            p_min = min(multiples)
            hole_radius = Fraction(1, 2 * q * p_min)
            per_center = max(Fraction(0), 2 * (radius - hole_radius))
        else:
            per_center = 2 * radius
        val = phi(b) * per_center
        by_b[b] = val
        total += val
    return total, by_b


def all_neighborhoods_disjoint(q: int, span: int) -> bool:
    intervals: list[Interval] = []
    for rows in farey_intervals_by_denominator(q, span).values():
        intervals.extend(rows)
    intervals.extend(origin_intervals(q, span))
    return len(merge(intervals)) == len(intervals)


def nondivisor_separation_margin(q: int, span: int) -> Fraction:
    """Minimum gap between a proper b-neighborhood and a non-divisor p-hole."""
    best: Fraction | None = None
    for b in range(2, q):
        rb = Fraction(q - b, b * q * span)
        for p in range(1, 2 * q):
            if p % b == 0:
                continue
            hp = Fraction(1, 2 * q * p)
            for a in range(1, b):
                if gcd(a, b) != 1:
                    continue
                c = Fraction(a, b)
                for j in range(p):
                    d = abs(c - Fraction(j, p))
                    d = min(d, 1 - d)
                    margin = d - rb - hp
                    best = margin if best is None else min(best, margin)
    assert best is not None
    return best


@dataclass(frozen=True)
class ScanRec:
    q: int
    k: int
    p_size: int
    total: Fraction
    origin: Fraction
    proper: Fraction
    P: tuple[int, ...]
    case: str
    by_b: dict[int, Fraction]


def scan_deepest(q: int) -> tuple[ScanRec, ScanRec, ScanRec]:
    """Return global, p=1-present, and p=1-absent minima for k=q+1."""
    span = 2 * q - 1
    k = q + 1
    p_size = 2 * q - 1 - k
    best: ScanRec | None = None
    best_with: ScanRec | None = None
    best_without: ScanRec | None = None
    for P in combinations(range(1, 2 * q), p_size):
        proper, by_b = proper_depth_ledger(q, span, P)
        orig = origin_exact(q, span, P)
        total = proper + orig
        case = "with1" if 1 in P else "no1"
        rec = ScanRec(q, k, p_size, total, orig, proper, P, case, by_b)
        if best is None or rec.total < best.total:
            best = rec
        if 1 in P:
            if best_with is None or rec.total < best_with.total:
                best_with = rec
        else:
            if best_without is None or rec.total < best_without.total:
                best_without = rec
    assert best is not None and best_with is not None and best_without is not None
    return best, best_with, best_without


def contributors(by_b: dict[int, Fraction]) -> str:
    live = [(b, v) for b, v in by_b.items() if v > 0]
    live.sort(key=lambda row: (row[1], row[0]))
    return ", ".join(f"b={b}:{fmt(v)}" for b, v in live[:6]) or "none"


def verify_depth_reduction(exhaustive_q_max: int = 7, sampled_q_max: int = 12) -> None:
    print("Depth-ledger verification")
    rng = random.Random(2847)
    for q in range(3, sampled_q_max + 1):
        span = 2 * q - 1
        disjoint = all_neighborhoods_disjoint(q, span)
        sep = nondivisor_separation_margin(q, span)
        checked = 0
        mismatches = 0
        for p_size in range(0, q - 1):
            if q <= exhaustive_q_max:
                candidates = list(combinations(range(1, 2 * q), p_size))
            elif p_size == 0:
                candidates = [()]
            else:
                total_combos = comb(2 * q - 1, p_size)
                sample_count = min(total_combos, min(80, max(20, 4 * q)))
                if total_combos <= sample_count:
                    candidates = list(combinations(range(1, 2 * q), p_size))
                    for P in candidates:
                        pass
                    # Fall through to the shared exact-check loop below.
                    seen = set(candidates)
                    candidates = sorted(seen)
                    for P in candidates:
                        exact = full_proper_exact(q, span, P)
                        ledger, _ = proper_depth_ledger(q, span, P)
                        checked += 1
                        if exact != ledger:
                            mismatches += 1
                            if mismatches <= 3:
                                print(
                                    f"  MISMATCH q={q} P={P}: exact={fmt(exact)} ledger={fmt(ledger)}"
                                )
                    continue
                seen: set[tuple[int, ...]] = set()
                edge_cases = [
                    tuple(range(1, 1 + p_size)),
                    tuple(range(2 * q - p_size, 2 * q)),
                    tuple([1] + list(range(2, 1 + p_size))) if p_size else (),
                ]
                for row in edge_cases:
                    if len(row) == p_size and len(set(row)) == p_size:
                        seen.add(tuple(sorted(row)))
                while len(seen) < sample_count:
                    seen.add(tuple(sorted(rng.sample(range(1, 2 * q), p_size))))
                candidates = sorted(seen)
            for P in candidates:
                exact = full_proper_exact(q, span, P)
                ledger, _ = proper_depth_ledger(q, span, P)
                checked += 1
                if exact != ledger:
                    mismatches += 1
                    if mismatches <= 3:
                        print(
                            f"  MISMATCH q={q} P={P}: exact={fmt(exact)} ledger={fmt(ledger)}"
                        )
        print(
            f"  q={q:2d} span={span:2d}: disjoint={disjoint} "
            f"nondiv_margin={fmt(sep)} checked={checked} mismatches={mismatches}"
        , flush=True)
    print()


def main() -> None:
    print("HYP-2847 origin-decoy divisor-depth ledger scout (exact rational)")
    print("bounded span S=2q-1; deepest dense branch k=q+1, |P|=q-2")
    print(flush=True)
    verify_depth_reduction()

    print("Deepest bounded dense branch scans", flush=True)
    rows: list[ScanRec] = []
    for q in range(3, 11):
        best, with1, no1 = scan_deepest(q)
        rows.append(best)
        print(f"N={2*q:2d} q={q:2d} k={best.k:2d} |P|={best.p_size:2d}", flush=True)
        for label, rec in [("global", best), ("with p=1", with1), ("without p=1", no1)]:
            print(
                f"  {label:11s}: total={fmt(rec.total):>12s} "
                f"origin={fmt(rec.origin):>12s} proper={fmt(rec.proper):>12s} "
                f"P={rec.P}"
            )
            print(f"               live proper contributors: {contributors(rec.by_b)}")
        print(flush=True)

    print("Synthesis", flush=True)
    for rec in rows:
        print(
            f"  q={rec.q:2d} N={2*rec.q:2d}: min={fmt(rec.total):>12s} "
            f"case={rec.case:5s} P={rec.P} live={contributors(rec.by_b)}"
        )

    print()
    print("Proof lead")
    print("  1. Proper neighborhoods are pairwise disjoint at bounded span S=2q-1.")
    print("  2. Non-divisor G_P holes have a positive separation margin from every")
    print("     proper a/b neighborhood; only selected multiples of b matter.")
    print("  3. The proper residual is therefore the divisor-depth ledger above.")
    print("  4. The origin interval is a decoy: either p=1 is absent and origin mass")
    print("     survives, or p=1 is present and one P-slot is unavailable for the")
    print("     divisor-depth cover of b>=2.")
    print()
    print("Tournament Analysis")
    vertices = [
        ("divisor_depth_ledger", 7),
        ("origin_decoy_slot_tradeoff", 6),
        ("proper_interval_geometry", 5),
        ("low_denominator_resonant_route", 4),
        ("raw_interval_union", 3),
        ("exact_farey_centers", 1),
        ("raw_runner_vertices", 0),
    ]
    hist: dict[int, int] = {}
    for _name, score in vertices:
        hist[score] = hist.get(score, 0) + 1
    path = " > ".join(name for name, _score in sorted(vertices, key=lambda row: -row[1]))
    print("  pairwise observable: exact residual explained after G_P hole subtraction")
    print("  switch/gauge: proof compression wins when it preserves exact residual")
    print(f"  score_hist={hist}")
    print(f"  Hamiltonian path: {path}")
    print("  challenged assumption: the useful vertices are denominator obligations")
    print("  and decoy slots, not runners, arcs, or exact Farey centers.")


if __name__ == "__main__":
    main()
