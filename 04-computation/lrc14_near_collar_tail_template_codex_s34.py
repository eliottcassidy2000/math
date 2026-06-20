#!/usr/bin/env python3
"""HYP-2654/T901 scout: near-collar AP-tail templates after THM-541.

THM-541 proves that the AP-window one-hole collar is uniquely minimized by
deleting 6.  The next obligation is subtler than "only the exact drop-6 row":
AP-tail rows can sit below the one-hole second value 426/35035 while still
retaining the drop-6 mouth geometry.

This script scans structured AP-tail pockets

    ({1,...,13} \\ holes) union replacements

where len(replacements)=len(holes)-1, so the core size remains 12.  It records
whether low-measure rows preserve the original drop-6 mouths or pay for any
damage by creating new mouth mass elsewhere.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd


TARGET_DENOM = 14
DROP6_COLLAR = Fraction(7, 858)
AP_SECOND = Fraction(426, 35035)
DROP6_CORE = tuple(v for v in range(1, 14) if v != 6)


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for lo, hi in sorted(intervals):
        if lo >= hi:
            continue
        if out and lo <= out[-1][1]:
            if hi > out[-1][1]:
                out[-1] = (out[-1][0], hi)
        else:
            out.append((lo, hi))
    return out


@lru_cache(maxsize=None)
def danger_arcs(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    arcs: list[tuple[Fraction, Fraction]] = []
    denom = TARGET_DENOM * speed
    for tooth in range(speed):
        left = TARGET_DENOM * tooth - 1
        right = TARGET_DENOM * tooth + 1
        if left < 0:
            arcs.append((Fraction(0), Fraction(right, denom)))
            arcs.append((Fraction(denom + left, denom), Fraction(1)))
        else:
            arcs.append((Fraction(left, denom), Fraction(right, denom)))
    return tuple(arcs)


def safe_components(core: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    danger: list[tuple[Fraction, Fraction]] = []
    for speed in core:
        danger.extend(danger_arcs(speed))
    safe: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merge(danger):
        if lo > cursor:
            safe.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        safe.append((cursor, Fraction(1)))
    return tuple(safe)


def measure(intervals: tuple[tuple[Fraction, Fraction], ...] | list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((hi - lo for lo, hi in intervals), Fraction(0))


def intersect_measure(
    xs: tuple[tuple[Fraction, Fraction], ...], ys: tuple[tuple[Fraction, Fraction], ...]
) -> Fraction:
    total = Fraction(0)
    i = j = 0
    while i < len(xs) and j < len(ys):
        lo = max(xs[i][0], ys[j][0])
        hi = min(xs[i][1], ys[j][1])
        if lo < hi:
            total += hi - lo
        if xs[i][1] < ys[j][1]:
            i += 1
        else:
            j += 1
    return total


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def sumset_excess(core: tuple[int, ...]) -> int:
    return len({a + b for a in core for b in core}) - (2 * len(core) - 1)


@dataclass(frozen=True)
class Row:
    safe: Fraction
    core: tuple[int, ...]
    holes: tuple[int, ...]
    replacements: tuple[int, ...]
    components: int
    old_survivor: Fraction
    old_damage: Fraction
    new_mass: Fraction
    excess: int


def row_for(holes: tuple[int, ...], replacements: tuple[int, ...], drop6_safe) -> Row:
    base = tuple(v for v in range(1, 14) if v not in holes)
    core = tuple(sorted(base + replacements))
    comps = safe_components(core)
    old_survivor = intersect_measure(comps, drop6_safe)
    return Row(
        safe=measure(comps),
        core=core,
        holes=holes,
        replacements=replacements,
        components=len(comps),
        old_survivor=old_survivor,
        old_damage=DROP6_COLLAR - old_survivor,
        new_mass=measure(comps) - old_survivor,
        excess=sumset_excess(core),
    )


def scan(bmax: int, max_holes: int, keep: int) -> tuple[list[Row], list[Row], int]:
    drop6_safe = safe_components(DROP6_CORE)
    rows_seen = 0
    top: list[Row] = []
    below_second: list[Row] = []
    for hcount in range(1, max_holes + 1):
        for holes in combinations(range(1, 14), hcount):
            for replacements in combinations(range(14, bmax + 1), hcount - 1):
                core = tuple(sorted(tuple(v for v in range(1, 14) if v not in holes) + replacements))
                if not primitive(core):
                    continue
                rows_seen += 1
                row = row_for(holes, replacements, drop6_safe)
                if row.safe < AP_SECOND:
                    below_second.append(row)
                top.append(row)
                top.sort(key=lambda r: (r.safe, r.core))
                if len(top) > keep:
                    top.pop()
    below_second.sort(key=lambda r: (r.safe, r.core))
    return top, below_second, rows_seen


def print_rows(rows: list[Row], title: str, limit: int | None = None) -> None:
    print(title)
    print("rank | safe | holes | repl | comps | old_surv | old_damage | new_mass | exc | core")
    selected = rows if limit is None else rows[:limit]
    for i, row in enumerate(selected, 1):
        print(
            f"{i:4d} | {str(row.safe):>12} | {row.holes!s:>12} | {row.replacements!s:>10} | "
            f"{row.components:5d} | {str(row.old_survivor):>9} | {str(row.old_damage):>10} | "
            f"{str(row.new_mass):>8} | {row.excess:3d} | {row.core}"
        )
    if limit is not None and len(rows) > limit:
        print(f"  ... {len(rows) - limit} more rows omitted")
    print()


def summarize_below_second(rows: list[Row]) -> None:
    print("below the AP one-hole second value")
    print(f"  threshold 426/35035: {fmt(AP_SECOND)}")
    print(f"  rows below threshold: {len(rows)}")
    missing6 = [r for r in rows if 6 in r.holes]
    print(f"  rows with 6 among AP holes: {len(missing6)}")
    print(f"  rows without 6 among AP holes: {len(rows) - len(missing6)}")
    if rows:
        min_without6 = min((r.safe for r in rows if 6 not in r.holes), default=None)
        if min_without6 is None:
            print("  no below-threshold row omits the drop-6 mouth")
        else:
            print(f"  best below-threshold row without 6: {fmt(min_without6)}")
    print()


def tournament_analysis(rows: list[Row]) -> None:
    vertices = [
        "drop6_mouth_retained",
        "old_mouth_damage",
        "new_mouth_mass",
        "hole_count",
        "sumset_excess",
        "raw_tail_size",
    ]
    print("Tournament Analysis")
    print("  vertices: proof-obligation observables for AP-tail near-collar rows")
    print("  pairwise observable: which observable separates rows below 426/35035 from the rest")
    print("  switch/gauge: classify by retained drop-6 mouth mass before scalar safe measure")
    print("  Hamiltonian path:")
    print("    " + " > ".join(vertices))
    print("  fingerprint: transitive priority tournament; empirical separator is presence of hole 6")
    if rows:
        fully_retained = sum(1 for r in rows if r.old_damage == 0)
        print(f"  below-threshold rows with undamaged old mouths: {fully_retained}/{len(rows)}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bmax", type=int, default=35)
    parser.add_argument("--max-holes", type=int, default=3)
    parser.add_argument("--keep", type=int, default=30)
    parser.add_argument("--below-limit", type=int, default=80)
    args = parser.parse_args()

    print("HYP-2654/T901 LRC14 near-collar AP-tail template scout")
    print(f"AP base: [1,13], target 1/{TARGET_DENOM}, bmax={args.bmax}, max_holes={args.max_holes}")
    print(f"drop-6 collar: {fmt(DROP6_COLLAR)}")
    print(f"AP one-hole second value: {fmt(AP_SECOND)}")
    print()

    top, below_second, rows_seen = scan(args.bmax, args.max_holes, args.keep)
    print(f"primitive AP-tail rows scanned: {rows_seen}")
    print()
    print_rows(top, "top AP-tail rows", limit=args.keep)
    summarize_below_second(below_second)
    print_rows(below_second, "rows below 426/35035", limit=args.below_limit)
    tournament_analysis(below_second)

    if any(6 not in row.holes for row in below_second):
        print("FAIL: found a below-threshold AP-tail row without the drop-6 hole.")
    else:
        print("PASS: in this AP-tail pocket, every row below 426/35035 keeps the drop-6 hole.")


if __name__ == "__main__":
    main()
