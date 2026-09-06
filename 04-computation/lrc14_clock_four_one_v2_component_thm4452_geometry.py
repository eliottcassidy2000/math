#!/usr/bin/env python3
"""Exact finite scout for strict q=4 one-v2 tail failure components."""

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
SPEC = importlib.util.spec_from_file_location("dyadic_physical_union", SOURCE)
M = importlib.util.module_from_spec(SPEC)
if SPEC.loader is None:
    raise RuntimeError("missing dyadic geometry loader")
SPEC.loader.exec_module(M)

Interval = tuple[Q, Q]
DELTA = Q(1, 14)


def gap(z: Q) -> Q:
    z %= 1
    return min(z, 1 - z)


def need(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def merge_strict(intervals: list[Interval]) -> list[Interval]:
    out: list[list[Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not out or left >= out[-1][1]:
            out.append([left, right])
        elif right > out[-1][1]:
            out[-1][1] = right
    return [(left, right) for left, right in out]


def intersect_strict(a: list[Interval], b: list[Interval]) -> list[Interval]:
    i = j = 0
    out: list[Interval] = []
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


def shift(family: list[Interval], amount: Q) -> list[Interval]:
    out: list[Interval] = []
    amount %= 1
    for left, right in family:
        left += amount
        right += amount
        if left >= 1:
            out.append((left - 1, right - 1))
        elif right > 1:
            out.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            out.append((left, right))
    return merge_strict(out)


def circular_widths(family: list[Interval]) -> tuple[Q, ...]:
    if not family:
        return ()
    widths = [right - left for left, right in family]
    if len(family) > 1 and family[0][0] == 0 and family[-1][1] == 1:
        widths = [family[0][1] + 1 - family[-1][0], *widths[1:-1]]
    return tuple(sorted(widths, reverse=True))


@cache
def pair_failure(a: int, b: int) -> tuple[Interval, ...]:
    u = merge_strict([*M.danger(a), *M.danger(b)])
    return tuple(intersect_strict(u, shift(u, Q(1, 2))))


@cache
def q4_failure(r: int, a: int, b: int) -> tuple[Interval, ...]:
    """Physical x-set where tails (2r,a,b) kill all four quarter lifts."""
    p = list(pair_failure(a, b))
    even_owned = intersect_strict(M.danger(2 * r), shift(p, -Q(1, 4)))
    odd_owned = intersect_strict(shift(M.danger(2 * r), -Q(1, 4)), p)
    return tuple(merge_strict([*even_owned, *odd_owned]))


def direct_failure(r: int, a: int, b: int, x: Q) -> bool:
    tails = (2 * r, a, b)
    return all(any(gap(Q(t) * (x + Q(j, 4))) < DELTA for t in tails) for j in range(4))


def literal_stats(r: int, a: int, b: int) -> tuple[Q, int, tuple[Q, ...]]:
    tails = (2 * r, a, b)
    walls = {Q(0), Q(1)}
    for t in tails:
        for j in range(4):
            fam = shift(M.danger(t), -Q(j, 4))
            for left, right in fam:
                walls.add(left)
                walls.add(right)
    walls = sorted(walls)
    live = [direct_failure(r, a, b, (left + right) / 2) for left, right in zip(walls, walls[1:])]
    widths: list[Q] = []
    current = Q(0)
    previous_live = False
    for i, is_live in enumerate(live):
        cell_width = walls[i + 1] - walls[i]
        if is_live:
            connected = previous_live and direct_failure(r, a, b, walls[i])
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
    need(not direct_failure(r, a, b, Q(0)), ("live zero wall", r, a, b))
    widths.sort(reverse=True)
    return (widths[0] if widths else Q(0), len(widths), tuple(widths))


def main() -> None:
    height = int(sys.argv[1]) if len(sys.argv) > 1 else 51
    odds = list(range(1, height + 1, 2))
    rows = []
    rows3 = []
    pair_rows = []
    for a, b in combinations(odds, 2):
        pw = circular_widths(list(pair_failure(a, b)))
        pair_rows.append((pw[0] if pw else Q(0), (a, b), len(pw)))
        for r in odds:
            if gcd(r, a, b) != 1:
                continue
            f = list(q4_failure(r, a, b))
            widths = circular_widths(f)
            row = (widths[0] if widths else Q(0), (r, a, b), len(widths))
            rows.append(row)
            if r % 3 and a % 3 and b % 3:
                rows3.append(row)
    rows.sort(reverse=True)
    rows3.sort(reverse=True)
    pair_rows.sort(reverse=True)
    print("Q4_ONE_V2_COMPONENT_PROBE")
    print(f"height={height}; primitive_rows={len(rows)}; primitive_3unit_rows={len(rows3)}")
    print("ALL_ODD_TOP30")
    for row in rows[:30]:
        print(row)
    print("ODD_3UNIT_TOP30")
    for row in rows3[:30]:
        print(row)
    print("PAIR_PHYSICAL_TOP20")
    for row in pair_rows[:20]:
        print(row)
    for r, a, b in [(7, 1, 11), (1, 1, 9), (1, 9, 11)]:
        f = list(q4_failure(r, a, b))
        formula = (circular_widths(f)[0] if f else Q(0), len(circular_widths(f)), circular_widths(f))
        literal = literal_stats(r, a, b)
        need(formula == literal, ((r, a, b), formula, literal))
        print("CONTROL", (r, a, b), formula)
    print("PASS")


if __name__ == "__main__":
    main()
