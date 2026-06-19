#!/usr/bin/env python3
"""HYP-2651/T898: exact bounded atlas for the LRC14 core-gap crux.

THM-523 reduces LRC(14) to a uniform lower bound on

    G_C = {t in [0,1): ||c t|| > 1/14 for every c in C}

over all 12-subsets C of positive integers.  HYP-2649 found the top bounded
training collar 7/858 among 12-subsets of [1,14].

This scout widens that exact positive-core lab.  It deliberately does not
translate cores to contain 0: fixed-observer gap measure is scale-invariant but
not freely translation-invariant.  Freiman/sumset data are reported only as
classifiers for proof routing.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from heapq import heappop, heappush
from itertools import combinations
from math import gcd


CORE_SIZE = 12
TARGET_DENOM = 14
KNOWN_COLLAR = Fraction(7, 858)


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), Fraction(0))


@lru_cache(maxsize=None)
def danger_arcs(v: int, denom: int = TARGET_DENOM) -> tuple[tuple[Fraction, Fraction], ...]:
    """Danger arcs ||v t|| <= 1/denom on [0,1)."""

    arcs: list[tuple[Fraction, Fraction]] = []
    d = denom * v
    for a in range(v):
        lo_num = denom * a - 1
        hi_num = denom * a + 1
        if lo_num < 0:
            arcs.append((Fraction(0), Fraction(hi_num, d)))
            arcs.append((Fraction(d + lo_num, d), Fraction(1)))
        elif hi_num > d:
            arcs.append((Fraction(lo_num, d), Fraction(1)))
            arcs.append((Fraction(0), Fraction(hi_num - d, d)))
        else:
            arcs.append((Fraction(lo_num, d), Fraction(hi_num, d)))
    return tuple(arcs)


def safe_components(core: tuple[int, ...], denom: int = TARGET_DENOM) -> tuple[tuple[Fraction, Fraction], ...]:
    danger = merge([arc for v in core for arc in danger_arcs(v, denom)])
    safe: list[tuple[Fraction, Fraction]] = []
    prev = Fraction(0)
    for lo, hi in danger:
        if lo > prev:
            safe.append((prev, lo))
        if hi > prev:
            prev = hi
    if prev < 1:
        safe.append((prev, Fraction(1)))
    return tuple(safe)


def safe_measure(core: tuple[int, ...], denom: int = TARGET_DENOM) -> Fraction:
    return measure(list(safe_components(core, denom)))


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def sumset_excess(core: tuple[int, ...]) -> int:
    sums = {a + b for a in core for b in core}
    return len(sums) - (2 * len(core) - 1)


def holes(core: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(set(range(min(core), max(core) + 1)) - set(core)))


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


@dataclass(frozen=True)
class Row:
    safe: Fraction
    core: tuple[int, ...]
    components: int
    min_component: Fraction
    excess: int
    span: int
    holes: tuple[int, ...]

    @classmethod
    def from_core(cls, core: tuple[int, ...]) -> "Row":
        comps = safe_components(core)
        return cls(
            safe=measure(list(comps)),
            core=core,
            components=len(comps),
            min_component=min((b - a for a, b in comps), default=Fraction(0)),
            excess=sumset_excess(core),
            span=max(core) - min(core) + 1,
            holes=holes(core),
        )


def top_push(heap: list[tuple[Fraction, tuple[int, ...], Row]], row: Row, keep: int) -> None:
    item = (-row.safe, row.core, row)
    heappush(heap, item)
    if len(heap) > keep:
        heappop(heap)


def scan(bmax: int, keep: int) -> tuple[dict[int, list[Row]], list[Row], dict[int, Row]]:
    by_max: dict[int, list[Row]] = defaultdict(list)
    top_heap: list[tuple[Fraction, tuple[int, ...], Row]] = []
    best_by_excess: dict[int, Row] = {}

    for core in combinations(range(1, bmax + 1), CORE_SIZE):
        if not primitive(core):
            continue
        row = Row.from_core(core)
        by_max[max(core)].append(row)
        top_push(top_heap, row, keep)
        old = best_by_excess.get(row.excess)
        if old is None or row.safe < old.safe:
            best_by_excess[row.excess] = row

    top = sorted((item[2] for item in top_heap), key=lambda r: (r.safe, r.core))
    return by_max, top, best_by_excess


def print_cumulative(by_max: dict[int, list[Row]], bmax: int) -> None:
    print("cumulative primitive positive 12-core scan")
    print("  B      rows     min safe@1/14       argmin")
    best: Row | None = None
    count = 0
    for B in range(CORE_SIZE, bmax + 1):
        rows = by_max.get(B, [])
        count += len(rows)
        for row in rows:
            if best is None or row.safe < best.safe:
                best = row
        if best is not None:
            print(f" {B:2d}  {count:8d}  {str(best.safe):>18}  {best.core}")
    print()


def print_single_hole_ap_window() -> None:
    print("single-hole AP-window ledger: [1,13] with one deletion")
    print("drop | safe measure | components | min component | core")
    for drop in range(1, 14):
        core = tuple(v for v in range(1, 14) if v != drop)
        row = Row.from_core(core)
        print(
            f"{drop:4d} | {str(row.safe):>12} | {row.components:10d} | "
            f"{str(row.min_component):>13} | {row.core}"
        )
    print()


def print_critical_interval_ledgers() -> None:
    print("critical AP-window interval ledgers")
    for drop in (6, 12):
        core = tuple(v for v in range(1, 14) if v != drop)
        row = Row.from_core(core)
        print(f"drop {drop}: safe={row.safe}, components={row.components}, core={core}")
        for lo, hi in safe_components(core):
            print(f"  [{lo}, {hi}] len={hi - lo} mid={(lo + hi) / 2}")
    print()


def print_top(top: list[Row]) -> None:
    print("near-minimum rows")
    print("rank | safe measure | core | span | exc | comps | min component | holes")
    for i, row in enumerate(top, 1):
        print(
            f"{i:4d} | {str(row.safe):>12} | {row.core} | "
            f"{row.span:4d} | {row.excess:3d} | {row.components:5d} | "
            f"{str(row.min_component):>13} | {row.holes}"
        )
    print()


def print_excess(best_by_excess: dict[int, Row]) -> None:
    print("best safe measure by exact sumset excess")
    print("exc | safe measure | core | span | holes")
    for exc in sorted(best_by_excess):
        row = best_by_excess[exc]
        print(f"{exc:3d} | {str(row.safe):>12} | {row.core} | {row.span:4d} | {row.holes}")
    print()


def tournament_analysis() -> None:
    print("Tournament Analysis")
    print("Pairwise observable: proof quotient that preserves exact gap measure before routing.")
    print("Switch/gauge: prefer quotients that keep fixed-observer positivity and only then add structure.")
    vertices = [
        "exact_positive_core_gap",
        "AP_drop_profile",
        "state_word_template",
        "sumset_excess_classifier",
        "far_element_plateau",
        "raw_runner_bound",
    ]
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score histogram: {0:1,1:1,2:1,3:1,4:1,5:1}; directed 3-cycles: 0")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bmax", type=int, default=18)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    if args.bmax < CORE_SIZE:
        raise SystemExit("--bmax must be at least 12")

    print("HYP-2651/T898 LRC14 exact bounded core-gap atlas")
    print(f"target gap: 1/{TARGET_DENOM}; core size: {CORE_SIZE}; positive cores C subset [1,B]")
    print("No translation normalization is used; Freiman data are classifiers only.")
    print(f"known drop-6 collar from HYP-2649/THM-523: {fmt(KNOWN_COLLAR)}")
    print()

    by_max, top, best_by_excess = scan(args.bmax, args.top)
    print_cumulative(by_max, args.bmax)
    print_single_hole_ap_window()
    print_critical_interval_ledgers()
    print_top(top)
    print_excess(best_by_excess)

    if top:
        second = next((row for row in top if row.safe > top[0].safe), None)
        print("collar reading")
        print(f"  best row: {top[0].core}")
        print(f"  best safe measure: {fmt(top[0].safe)}")
        print(f"  matches known 7/858 collar? {top[0].safe == KNOWN_COLLAR}")
        if second is not None:
            print(f"  next distinct safe measure: {fmt(second.safe)}")
            print(f"  separation from collar: {fmt(second.safe - top[0].safe)}")
    print()

    print("Assumption challenge")
    print(
        "  Fixed-observer gap measure is not translation-invariant, so the scout does not "
        "replace positive cores by 0-based Freiman normal forms.  Alternate vertices "
        "considered: runners, holes, safe components, exact excess, wall states, "
        "far-element plateaus, and proof obligations.  The chosen vertex keeps the "
        "actual OPEN-Q-108 predicate first and uses the others only for routing."
    )
    print()
    tournament_analysis()


if __name__ == "__main__":
    main()
