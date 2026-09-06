#!/usr/bin/env python3
"""Exact strict-versus-a.e. component scout for odd q=2 failure sets."""

from __future__ import annotations

from fractions import Fraction as Q
from functools import cache
from itertools import combinations
from math import gcd
from pathlib import Path
import importlib.util
import sys


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "lrc14_dyadic_physical_union_sharp_thm4449.py"
SPEC = importlib.util.spec_from_file_location("all_odd_q2", SOURCE)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)
DELTA = Q(1, 14)


def gap(z: Q) -> Q:
    z %= 1
    return min(z, 1 - z)


def failure(tails: tuple[int, int, int], x: Q) -> bool:
    return (
        any(gap(Q(t) * x) < DELTA for t in tails)
        and any(gap(Q(t) * (x + Q(1, 2))) < DELTA for t in tails)
    )


def merge_strict(intervals: list[M.Interval]) -> list[M.Interval]:
    """Merge positive overlaps but retain a deleted common endpoint."""
    out: list[list[Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not out or left >= out[-1][1]:
            out.append([left, right])
        elif right > out[-1][1]:
            out[-1][1] = right
    return [(left, right) for left, right in out]


def intersect_strict(a: list[M.Interval], b: list[M.Interval]) -> list[M.Interval]:
    i = j = 0
    out: list[M.Interval] = []
    while i < len(a) and j < len(b):
        left = max(a[i][0], b[j][0])
        right = min(a[i][1], b[j][1])
        if left < right:
            out.append((left, right))
        if a[i][1] < b[j][1]:
            i += 1
        elif b[j][1] < a[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge_strict(out)


@cache
def oriented_cross(i: int, j: int) -> tuple[M.Interval, ...]:
    return tuple(intersect_strict(M.danger(i), M.shift_half(M.danger(j))))


@cache
def strict_stats(tails: tuple[int, int, int]) -> tuple[Q, int, tuple[Q, ...]]:
    """Actual open-set component widths, retaining isolated wall holes."""
    pieces = [
        arc
        for i in tails
        for j in tails
        if i != j
        for arc in oriented_cross(i, j)
    ]
    f_strict = merge_strict(pieces)
    widths = [right - left for left, right in f_strict]
    # For odd tails F never contains zero: sheet one is exactly antipodal for
    # every tail there.  Thus no first/last circular join can occur.
    assert not failure(tails, Q(0))
    widths.sort(reverse=True)
    mass = sum(widths, Q(0))
    assert mass == M.literal_physical(tails), (tails, mass, M.literal_physical(tails))
    return widths[0] if widths else Q(0), len(widths), tuple(widths)


@cache
def ae_stats(tails: tuple[int, int, int]) -> tuple[Q, int]:
    a = M.merge([arc for t in tails for arc in M.danger(t)])
    f = M.intersect(a, M.shift_half(a))
    return max((right - left for left, right in f), default=Q(0)), len(f)


@cache
def literal_strict_stats(tails: tuple[int, int, int]) -> tuple[Q, int, tuple[Q, ...]]:
    """Independent wall-cell topology, including truth at shared walls."""
    families = [M.danger(t) for t in tails]
    shifted = [M.shift_half(family) for family in families]
    walls = sorted(
        {
            Q(0),
            Q(1),
            *(x for family in families + shifted for arc in family for x in arc),
        }
    )
    live = [failure(tails, (left + right) / 2) for left, right in zip(walls, walls[1:])]
    widths: list[Q] = []
    current = Q(0)
    previous_live = False
    for i, is_live in enumerate(live):
        cell_width = walls[i + 1] - walls[i]
        if is_live:
            connected = previous_live and failure(tails, walls[i])
            if current and not connected:
                widths.append(current)
                current = Q(0)
            current += cell_width
        elif current:
            widths.append(current)
            current = Q(0)
        previous_live = is_live
    if current:
        widths.append(current)
    # F(0) is false for all-odd tails, so there is no circular join.
    assert not failure(tails, Q(0))
    widths.sort(reverse=True)
    mass = sum(widths, Q(0))
    assert mass == M.literal_physical(tails)
    return widths[0] if widths else Q(0), len(widths), tuple(widths)


def main() -> None:
    height = int(sys.argv[1]) if len(sys.argv) > 1 else 99
    odds = list(range(1, height + 1, 2))
    primitive = [tails for tails in combinations(odds, 3) if gcd(*tails) == 1]
    all_rows = []
    unit_rows = []
    split_rows = []
    for tails in primitive:
        strict_longest, strict_count, _ = strict_stats(tails)
        ae_longest, ae_count = ae_stats(tails)
        row = (strict_longest, tails, strict_count, ae_longest, ae_count)
        all_rows.append(row)
        if all(t % 3 for t in tails):
            unit_rows.append(row)
        if strict_longest != ae_longest or strict_count != ae_count:
            split_rows.append(row)
    all_rows.sort(reverse=True)
    unit_rows.sort(reverse=True)
    split_rows.sort(reverse=True)
    print("Q2_ODD_COMPONENT_TOPOLOGY_PROBE")
    print("STATUS=FINITE_EXACT_DISCOVERY_ONLY")
    print(f"height={height} primitive_all_odd_rows={len(all_rows)}")
    print(f"primitive_odd_3unit_rows={len(unit_rows)} split_topology_rows={len(split_rows)}")
    print("ALL_ODD_TOP20_STRICT")
    for row in all_rows[:20]:
        print(row)
    print("ODD_3UNIT_TOP20_STRICT")
    for row in unit_rows[:20]:
        print(row)
    print("TOP_AE_VS_STRICT_SPLITS")
    for row in sorted(split_rows, key=lambda r: (r[3], r[1]), reverse=True)[:20]:
        print(row)
    controls = [(1, 9, 11), (1, 7, 13), (1, 9, 23)]
    print("CONTROLS")
    for tails in controls:
        print(f"T={tails} strict={strict_stats(tails)} ae={ae_stats(tails)}")
    print("PASS")


if __name__ == "__main__":
    main()
