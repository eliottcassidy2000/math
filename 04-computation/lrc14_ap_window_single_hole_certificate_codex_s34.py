#!/usr/bin/env python3
"""THM-541: exact certificate for the LRC14 AP-window single-hole collar.

For each one-hole 12-core

    C_e = {1,2,...,13} \\ {e},

this script computes the level-1/14 safe set

    G_C = {t in [0,1): ||c t|| > 1/14 for all c in C}

by exact rational interval arithmetic.  The certificate keeps the boundary
addresses: every safe component is a gap from a right danger wall

    R(v,a) = (14a+1)/(14v)

to a left danger wall

    L(w,b) = (14b-1)/(14w).

That signed wall-pair address is the proof object.  The scalar measure table is
only the final readout.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd


TARGET_DENOM = 14
KNOWN_COLLAR = Fraction(7, 858)
NEXT_COLLAR = Fraction(426, 35035)


@dataclass(frozen=True, order=True)
class Wall:
    """A danger-interval boundary at target 1/14."""

    point: Fraction
    speed: int
    tooth: int  # periodic tooth index; the wrapped left wall at 1 uses tooth=speed
    side: str  # "L" for (14a-1)/(14v), "R" for (14a+1)/(14v)

    def label(self) -> str:
        return f"{self.side}({self.speed},{self.tooth})"


@dataclass(frozen=True)
class Component:
    lo: Fraction
    hi: Fraction
    left_owners: tuple[Wall, ...]
    right_owners: tuple[Wall, ...]

    @property
    def length(self) -> Fraction:
        return self.hi - self.lo

    @property
    def midpoint(self) -> Fraction:
        return (self.lo + self.hi) / 2

    def primary_pair(self) -> tuple[Wall, Wall]:
        """Return a canonical R-to-L owner pair for the open safe gap."""

        left = min((w for w in self.left_owners if w.side == "R"), key=lambda w: (w.speed, w.tooth))
        right = min((w for w in self.right_owners if w.side == "L"), key=lambda w: (w.speed, w.tooth))
        return left, right

    def determinant_numerator(self) -> int:
        """Signed numerator before reducing the R-to-L wall-pair length.

        If lo=R(v,a) and hi=L(w,b), then

            hi-lo = (v(14b-1) - w(14a+1))/(14vw).
        """

        left, right = self.primary_pair()
        return left.speed * (TARGET_DENOM * right.tooth - 1) - right.speed * (
            TARGET_DENOM * left.tooth + 1
        )


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


def danger_arcs(speed: int) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    denom = TARGET_DENOM * speed
    for tooth in range(speed):
        left_num = TARGET_DENOM * tooth - 1
        right_num = TARGET_DENOM * tooth + 1
        if left_num < 0:
            arcs.append((Fraction(0), Fraction(right_num, denom)))
            arcs.append((Fraction(denom + left_num, denom), Fraction(1)))
        else:
            arcs.append((Fraction(left_num, denom), Fraction(right_num, denom)))
    return arcs


def walls(speed: int) -> list[Wall]:
    out: list[Wall] = []
    denom = TARGET_DENOM * speed
    for tooth in range(speed):
        left_num = TARGET_DENOM * tooth - 1
        right_num = TARGET_DENOM * tooth + 1
        if left_num < 0:
            out.append(Wall(Fraction(denom + left_num, denom), speed, speed, "L"))
            out.append(Wall(Fraction(right_num, denom), speed, tooth, "R"))
        else:
            out.append(Wall(Fraction(left_num, denom), speed, tooth, "L"))
            out.append(Wall(Fraction(right_num, denom), speed, tooth, "R"))
    return out


def safe_components(core: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    danger = merge([arc for speed in core for arc in danger_arcs(speed)])
    safe: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in danger:
        if lo > cursor:
            safe.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        safe.append((cursor, Fraction(1)))
    return tuple(safe)


def addressed_components(core: tuple[int, ...]) -> tuple[Component, ...]:
    wall_map: dict[Fraction, list[Wall]] = {}
    for speed in core:
        for wall in walls(speed):
            wall_map.setdefault(wall.point, []).append(wall)

    comps: list[Component] = []
    for lo, hi in safe_components(core):
        left_owners = tuple(sorted(wall_map.get(lo, [])))
        right_owners = tuple(sorted(wall_map.get(hi, [])))
        comp = Component(lo, hi, left_owners, right_owners)
        if not any(w.side == "R" for w in left_owners):
            raise AssertionError(f"left endpoint {lo} lacks a right-wall owner in {core}")
        if not any(w.side == "L" for w in right_owners):
            raise AssertionError(f"right endpoint {hi} lacks a left-wall owner in {core}")
        comps.append(comp)
    return tuple(comps)


def safe_measure(core: tuple[int, ...]) -> Fraction:
    return sum((hi - lo for lo, hi in safe_components(core)), Fraction(0))


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def drop_core(drop: int) -> tuple[int, ...]:
    return tuple(v for v in range(1, 14) if v != drop)


def inside_danger(point: Fraction, speed: int) -> bool:
    frac = (speed * point) % 1
    return min(frac, 1 - frac) <= Fraction(1, TARGET_DENOM)


def validate_core(core: tuple[int, ...]) -> None:
    comps = addressed_components(core)
    raw = safe_components(core)
    if tuple((c.lo, c.hi) for c in comps) != raw:
        raise AssertionError("addressed/raw component mismatch")
    for comp in comps:
        mid = comp.midpoint
        if any(inside_danger(mid, speed) for speed in core):
            raise AssertionError(f"unsafe midpoint {mid} in component {comp}")
        for endpoint in (comp.lo, comp.hi):
            if not any(endpoint == wall.point for speed in core for wall in walls(speed)):
                raise AssertionError(f"endpoint {endpoint} is not a wall")


def pair_label(comp: Component) -> str:
    left, right = comp.primary_pair()
    return f"{left.label()} -> {right.label()} det={comp.determinant_numerator()}"


def print_drop_table() -> None:
    rows = []
    for drop in range(1, 14):
        core = drop_core(drop)
        validate_core(core)
        comps = addressed_components(core)
        measure = safe_measure(core)
        rows.append(
            (
                measure,
                drop,
                len(comps),
                min((c.length for c in comps), default=Fraction(0)),
                max((c.length for c in comps), default=Fraction(0)),
                core,
            )
        )

    print("one-hole AP-window exact table")
    print("rank | drop | safe measure | delta from 7/858 | components | min arc | max arc | core")
    for rank, (measure, drop, count, min_arc, max_arc, core) in enumerate(sorted(rows), 1):
        print(
            f"{rank:4d} | {drop:4d} | {str(measure):>12} | "
            f"{str(measure - KNOWN_COLLAR):>16} | {count:10d} | "
            f"{str(min_arc):>8} | {str(max_arc):>8} | {core}"
        )
    print()

    best = min(rows)
    assert best[0] == KNOWN_COLLAR and best[1] == 6
    assert all(measure > KNOWN_COLLAR for measure, drop, *_ in rows if drop != 6)
    second = min(measure for measure, drop, *_ in rows if drop != 6)
    assert second == NEXT_COLLAR
    print("certificate readout")
    print(f"  unique minimum: drop 6 with {fmt(KNOWN_COLLAR)}")
    print(f"  next value:     {fmt(NEXT_COLLAR)}")
    print(f"  separation:     {fmt(NEXT_COLLAR - KNOWN_COLLAR)}")
    print()


def print_component_ledgers() -> None:
    print("addressed component ledgers")
    for drop in range(1, 14):
        core = drop_core(drop)
        comps = addressed_components(core)
        print(f"drop {drop}: core={core}, safe={safe_measure(core)}, components={len(comps)}")
        for comp in comps:
            left_owners = ",".join(w.label() for w in comp.left_owners)
            right_owners = ",".join(w.label() for w in comp.right_owners)
            print(
                f"  [{comp.lo}, {comp.hi}] len={comp.length} "
                f"owners=({left_owners})=>({right_owners}) {pair_label(comp)}"
            )
        print()


def print_drop6_signed_geometry() -> None:
    drop = 6
    core = drop_core(drop)
    comps = addressed_components(core)
    dets = [c.determinant_numerator() for c in comps]
    print("drop-6 signed-wall geometry")
    print("  The four surviving mouths are R-to-L gaps with determinant numerators:")
    print(f"  {dets}")
    for comp in comps:
        left, right = comp.primary_pair()
        print(
            f"  {left.label()} at {comp.lo} to {right.label()} at {comp.hi}: "
            f"det={comp.determinant_numerator()}, len={comp.length}"
        )
    print()


def print_tournament_analysis() -> None:
    rows = sorted((safe_measure(drop_core(drop)), drop) for drop in range(1, 14))
    score = {drop: 0 for drop in range(1, 14)}
    for a, b in combinations(range(1, 14), 2):
        ma = safe_measure(drop_core(a))
        mb = safe_measure(drop_core(b))
        if ma < mb:
            score[a] += 1
        elif mb < ma:
            score[b] += 1
    hist = Counter(score.values())
    print("Tournament Analysis")
    print("  vertices: deleted AP-window positions e=1..13")
    print("  pairwise observable: lower exact safe measure wins")
    print("  switch/gauge: retain R/L boundary-owner addresses before scalar measure")
    print("  Hamiltonian path:")
    print("    " + " > ".join(f"drop-{drop}" for _, drop in rows))
    print(f"  score histogram: {dict(sorted(hist.items()))}")
    print("  directed 3-cycles: 0 (strict total order by exact measure)")
    print()


def print_frontier_check() -> None:
    print("nearby two-delete/one-replacement sanity check")
    print("  This is not used in THM-541; it checks that the theorem is local.")
    best = None
    for missing in combinations(range(1, 14), 2):
        base = tuple(v for v in range(1, 14) if v not in missing)
        for replacement in range(14, 41):
            core = tuple(sorted(base + (replacement,)))
            if not primitive(core):
                continue
            val = safe_measure(core)
            if best is None or (val, core, missing, replacement) < best:
                best = (val, core, missing, replacement)
    if best is None:
        raise AssertionError("empty frontier check")
    val, core, missing, replacement = best
    print(f"  best primitive row with two AP holes and replacement <=40: {core}")
    print(f"  missing={missing}, replacement={replacement}, safe={fmt(val)}")
    print(f"  collar gap above 7/858: {fmt(val - KNOWN_COLLAR)}")
    print()


def main() -> None:
    print("THM-541 LRC14 AP-window single-hole collar certificate")
    print(f"target: 1/{TARGET_DENOM}")
    print("safe components are addressed as R(v,a)->L(w,b) signed wall gaps")
    print()
    print_drop_table()
    print_component_ledgers()
    print_drop6_signed_geometry()
    print_tournament_analysis()
    print_frontier_check()
    print("PASS: exact finite certificate proves the AP-window single-hole collar.")


if __name__ == "__main__":
    main()
