#!/usr/bin/env python3
"""Exact first fixed-point recurrence for the THM-2693 raw delayed word.

The promoted carrier uses the base map y -> {13y} and is nilpotent at depth
four.  This scout changes only the base multiplier and asks for the cheapest
possible recurrent control: an interior fixed point of y -> {b y}.

It reconstructs the union of THM-2623's two raw guard sectors, independently
checks the guard-free factorization used by THM-2693, and exhausts every fixed
point k/(b-1) for 2 <= b <= 34.  The first surviving multiplier is b=18,
with exactly y=4/17 and 13/17.  The first point is the old THM-789 heavy-phase
hostile.  This is a foreign-base design target, not a physical LRC handoff.
"""

from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_successor_halfcell_carry_no_go_thm2623 as prior


P = 13
THRESHOLD = Fraction(1, 14)
U0 = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge(intervals):
    out = []
    for left, right in sorted(intervals):
        require(left < right, "empty interval entered raw delayed union")
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def distance(value):
    value %= 1
    return min(value, 1 - value)


def locate(point, intervals, starts, grid):
    scaled = point * grid
    index = bisect_right(starts, scaled) - 1
    if index < 0:
        return None
    left, right = intervals[index]
    if not left <= scaled < right:
        return None
    return (
        Fraction(left, grid),
        Fraction(right, grid),
        Fraction(scaled - left, grid),
        Fraction(right - scaled, grid),
    )


def fixed_points(multiplier, intervals, starts, grid):
    points = []
    for numerator in range(multiplier - 1):
        point = Fraction(numerator, multiplier - 1)
        component = locate(point, intervals, starts, grid)
        if component is not None:
            points.append((point, component))
    return tuple(points)


def main():
    module, *_ = cross.build_carrier_data()
    grid = cross.T
    sectors = prior.sector_words(module)
    raw_union = merge(sectors[0] + sectors[1])

    guard_free = module.make_comb(module.C2, 182, -13, 13)
    for index in module.UNIT_IDX:
        guard_free = module.subtract_comb(
            guard_free, module.W[index], 182, -13, 13
        )
    guard_free = module.subtract_comb(
        guard_free, module.C3, 182, -13, 13
    )
    require(raw_union == guard_free,
            "raw sector union stopped equalling the guard-free word")
    require((len(sectors[0]), len(sectors[1]), len(raw_union))
            == (32_578, 14_906, 47_484),
            "raw delayed component census changed")

    starts = [left for left, _ in raw_union]
    atlas = {
        multiplier: fixed_points(multiplier, raw_union, starts, grid)
        for multiplier in range(2, 35)
    }
    first_multiplier = next(
        multiplier for multiplier, points in atlas.items() if points
    )
    require(first_multiplier == 18,
            "first recurrent foreign multiplier changed")
    require(tuple(point for point, _ in atlas[18])
            == (Fraction(4, 17), Fraction(13, 17)),
            "slope-18 fixed-point pair changed")
    require(all(not atlas[multiplier] for multiplier in range(2, 18)),
            "a smaller multiplier acquired a fixed delayed point")

    point = Fraction(4, 17)
    reflected = Fraction(13, 17)
    require((18 * point) % 1 == point
            and (18 * reflected) % 1 == reflected,
            "slope-18 points stopped being fixed")

    target_distance = distance(module.C2 * point)
    unit_distances = tuple(
        distance(module.W[index] * point)
        for index in module.UNIT_IDX
    )
    high_distance = distance(module.C3 * point)
    require(target_distance == Fraction(1, 17) < THRESHOLD,
            "target tooth lost the 1/17 danger residue")
    require(unit_distances
            == tuple(Fraction(value, 17) for value in (5, 6, 7, 8, 8)),
            "ordinary safety residue word changed")
    require(min(unit_distances) > THRESHOLD
            and high_distance == Fraction(2, 17) > THRESHOLD,
            "fixed point left the strict guard-free word")

    u0_clearances = tuple(distance(speed * point) for speed in U0)
    require(u0_clearances == tuple(
        Fraction(value, 17) for value in (4, 8, 5, 3, 6, 2, 2, 6, 7, 3)
    ), "THM-789 heavy-phase clearance word changed")

    component = locate(point, raw_union, starts, grid)
    reflected_component = locate(reflected, raw_union, starts, grid)
    require(component is not None and reflected_component is not None,
            "fixed points lost their strict components")
    require(component[2:] == tuple(reversed(reflected_component[2:])),
            "reflection stopped interchanging component margins")

    congruent_family = tuple(1 + 17 * j for j in range(1, 6))
    require(all((multiplier * point) % 1 == point
                for multiplier in congruent_family),
            "b=1 mod 17 fixed family changed")

    print("LRC14 raw delayed-word foreign-slope fixed-point scout")
    print("status=VERIFIED-EXACT SCOUT; alternate base only")
    print("raw_sector_components=safe:32578,danger:14906,union:47484")
    print("fixed_multiplier_search=2..34")
    print("first_multiplier_with_interior_fixed_point=18")
    print("smaller_multipliers_2..17=NONE")
    print(f"slope18_fixed_points={(point, reflected)}")
    print(f"point_component={component[:2]}; margins={component[2:]}")
    print(f"reflected_component={reflected_component[:2]}; margins={reflected_component[2:]}")
    print(f"delayed_factor_distances=target:{target_distance}; units:{unit_distances}; high:{high_distance}")
    print(f"thm789_U0_clearances={u0_clearances}; minimum={min(u0_clearances)}")
    print(f"fixed_multiplier_family_b_eq_1_mod17={congruent_family}")
    print("mechanism=13^3*4=-1 mod17; five ordinary units have distances 5,6,7,8,8; high speed has distance2")
    print("scope=no physical multiplier-18 handoff, endpoint transport, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
