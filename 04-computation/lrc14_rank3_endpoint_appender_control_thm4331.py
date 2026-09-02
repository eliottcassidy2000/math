#!/usr/bin/env python3
"""Exact global-wall control for the THM-4331 single-appender sidecar.

This primary path constructs the complete reduced wall arrangement, classifies
every open cell at its exact midpoint, and merges adjacent safe cells.  It is
deliberately independent of the successive interval-intersection audit.
"""

from __future__ import annotations

from fractions import Fraction as Q
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
K = (10, 80, 85, 95, 120, 143, 145, 168)
PAIR = (50, 70)
AUXILIARY = 193
TAILS = (1, 9)
EXPECTED_S = (1, 9, 20, 100, 140, 160, 170, 190, 240, 286, 290, 336)
EXPECTED_COMPONENT = (Q(227, 476), Q(937, 1960))

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(distance(speed * phase) for speed in speeds)


def walls(speeds: tuple[int, ...]) -> tuple[Q, ...]:
    result = {Q(0), Q(1)}
    for speed in speeds:
        for k in range(speed):
            result.add(Q(14 * k + 1, 14 * speed))
            result.add(Q(14 * k + 13, 14 * speed))
    return tuple(sorted(result))


def global_wall_components(
    speeds: tuple[int, ...],
) -> tuple[tuple[Q, ...], tuple[tuple[Q, Q], ...], int]:
    cuts = walls(speeds)
    safe_cells: list[tuple[Q, Q]] = []
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left + right) / 2
        if clearance(speeds, midpoint) < DELTA:
            continue
        require(clearance(speeds, left) >= DELTA, ("left wall", left, right))
        require(clearance(speeds, right) >= DELTA, ("right wall", left, right))
        safe_cells.append((left, right))

    components: list[tuple[Q, Q]] = []
    for left, right in safe_cells:
        if components and components[-1][1] == left:
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return cuts, tuple(components), len(safe_cells)


def nonstrict_integer_cutoff(value: Q) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def strict_integer_cutoff(value: Q) -> int:
    return value.numerator // value.denominator + 1


def main() -> None:
    body = tuple(sorted((*K, AUXILIARY)))
    body_mask = sum(1 << POOL.index(value) for value in body)
    require(body_mask == 0x011CD402, ("body mask", hex(body_mask)))

    base = tuple(sorted((*K, *PAIR)))
    speeds = tuple(sorted((*TAILS, *(2 * value for value in base))))
    require(speeds == EXPECTED_S, ("physical speeds", speeds))
    require(len(speeds) == 12 and len(set(speeds)) == 12, "speed typing")

    cuts, components, safe_cell_count = global_wall_components(speeds)
    mass = sum((right - left for left, right in components), Q(0))
    require(len(cuts) == 3882, ("wall count", len(cuts)))
    require(safe_cell_count == 362, ("safe cell count", safe_cell_count))
    require(len(components) == 362, ("component count", len(components)))
    require(mass == Q(64917367, 577007200), ("mass", mass))
    require(EXPECTED_COMPONENT in components, "distinguished component")

    left, right = EXPECTED_COMPONENT
    width = right - left
    left_owners = tuple(v for v in speeds if distance(v * left) == DELTA)
    right_owners = tuple(v for v in speeds if distance(v * right) == DELTA)
    require(width == Q(39, 33320), ("width", width))
    require((left.denominator, right.denominator) == (476, 1960), "denominators")
    require(left_owners == (170,), ("left owners", left_owners))
    require(right_owners == (140,), ("right owners", right_owners))

    simple_ratio = (Q(1, 7) - Q(1, min(left.denominator, right.denominator))) / width
    additive_ratio = (
        Q(1, 7) - Q(1, left.denominator) - Q(1, right.denominator)
    ) / width
    simple_cutoff = nonstrict_integer_cutoff(simple_ratio)
    additive_cutoff = strict_integer_cutoff(additive_ratio)
    coarse_cutoff = 27 * sum(speeds)
    reserve_cutoff = 27 * (sum(speeds) - 1)
    require(simple_ratio == Q(4690, 39), ("simple ratio", simple_ratio))
    require(additive_ratio == Q(4673, 39), ("additive ratio", additive_ratio))
    require(simple_cutoff == 121, ("simple cutoff", simple_cutoff))
    require(additive_cutoff == 120, ("additive cutoff", additive_cutoff))
    require(coarse_cutoff == 52434, ("coarse cutoff", coarse_cutoff))
    require(reserve_cutoff == 52407, ("reserve cutoff", reserve_cutoff))

    print("THM4331 RANK3 ENDPOINT APPENDER CONTROL PRIMARY")
    print(f"PROVENANCE pair={PAIR[0]},{PAIR[1]} body_mask={body_mask:08x} deleted={AUXILIARY}")
    print("S=" + ",".join(map(str, speeds)))
    print(f"walls={len(cuts)} safe_cells={safe_cell_count} components={len(components)}")
    print(f"mass={mass}")
    print(f"component={left},{right} width={width}")
    print(
        f"owners={left_owners[0]},{right_owners[0]} "
        f"denominators={left.denominator},{right.denominator}"
    )
    print(f"simple_ratio={simple_ratio} simple_cutoff={simple_cutoff}")
    print(
        f"additive_ratio={additive_ratio} "
        f"additive_strict_cutoff={additive_cutoff}"
    )
    print(f"coarse_component_bound={coarse_cutoff}")
    print(f"global_reserve_bound={reserve_cutoff}")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
