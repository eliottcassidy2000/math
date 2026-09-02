#!/usr/bin/env python3
"""Independent exact interval audit for the THM-4331 rank-three control.

Unlike the primary, this path never constructs a global wall arrangement or
classifies midpoints.  It intersects each speed's closed safe-band list with
the current exact interval list, one speed at a time.
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


def safe_bands(speed: int) -> tuple[tuple[Q, Q], ...]:
    return tuple(
        (Q(14 * k + 1, 14 * speed), Q(14 * k + 13, 14 * speed))
        for k in range(speed)
    )


def merge(intervals: list[tuple[Q, Q]]) -> tuple[tuple[Q, Q], ...]:
    result: list[tuple[Q, Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect_lists(
    first: tuple[tuple[Q, Q], ...],
    second: tuple[tuple[Q, Q], ...],
) -> tuple[tuple[Q, Q], ...]:
    intersections: list[tuple[Q, Q]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            intersections.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge(intersections)


def successive_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    current = ((Q(0), Q(1)),)
    for speed in speeds:
        current = intersect_lists(current, safe_bands(speed))
        require(
            all(current[i][1] < current[i + 1][0] for i in range(len(current) - 1)),
            ("disjoint stage", speed),
        )
    return current


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

    components = successive_components(speeds)
    mass = sum((right - left for left, right in components), Q(0))
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

    # Independently pin the active safe bands at the two endpoint owners.
    require(left == Q(14 * 81 + 1, 14 * 170), "left owner wall")
    require(right == Q(14 * 67 - 1, 14 * 140), "right owner wall")

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

    print("THM4331 RANK3 ENDPOINT APPENDER CONTROL INDEPENDENT")
    print(f"PROVENANCE pair={PAIR[0]},{PAIR[1]} body_mask={body_mask:08x} deleted={AUXILIARY}")
    print("S=" + ",".join(map(str, speeds)))
    print(f"stages={len(speeds)} components={len(components)}")
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
