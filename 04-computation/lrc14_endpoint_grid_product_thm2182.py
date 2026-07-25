#!/usr/bin/env python3
"""Independent exact checks for THM-2182."""

from fractions import Fraction
from math import lcm


RADIUS = Fraction(1, 14)
CORE = (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12)
TAIL_1 = (3, 12)
TAIL_2 = (4, 6)
ENDPOINT_MODULUS = 55440

EXPECTED_CORE = Fraction(883, 6930)
EXPECTED_TAIL_1 = Fraction(3, 4)
EXPECTED_TAIL_2 = Fraction(16, 21)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def is_safe(speed: int, time: Fraction) -> bool:
    phase = (speed * time) % 1
    return min(phase, 1 - phase) >= RADIUS


def danger_arcs(
    speeds: tuple[int, ...],
) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        half_width = Fraction(1, 14 * speed)
        for index in range(speed):
            left = (Fraction(index, speed) - half_width) % 1
            right = left + 2 * half_width
            if right <= 1:
                arcs.append((left, right))
            else:
                arcs.append((left, Fraction(1)))
                arcs.append((Fraction(0), right - 1))
    return arcs


def safe_measure_by_union(speeds: tuple[int, ...]) -> Fraction:
    arcs = sorted(danger_arcs(speeds))
    danger = Fraction(0)
    left: Fraction | None = None
    right: Fraction | None = None
    for next_left, next_right in arcs:
        if left is None:
            left, right = next_left, next_right
        elif next_left <= right:
            right = max(right, next_right)
        else:
            danger += right - left
            left, right = next_left, next_right
    if left is not None:
        danger += right - left
    return 1 - danger


def boundary_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    points = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for index in range(speed):
            points.add((Fraction(14 * index - 1, 14 * speed)) % 1)
            points.add((Fraction(14 * index + 1, 14 * speed)) % 1)
    return tuple(sorted(points))


def safe_measure_by_cells(speeds: tuple[int, ...]) -> Fraction:
    points = boundary_points(speeds)
    return sum(
        (
            right - left
            for left, right in zip(points, points[1:])
            if all(
                is_safe(speed, (left + right) / 2)
                for speed in speeds
            )
        ),
        Fraction(0),
    )


def main() -> None:
    required_modulus = lcm(*(14 * speed for speed in CORE))
    require(
        ENDPOINT_MODULUS == required_modulus,
        "incorrect endpoint modulus",
    )

    values: dict[str, Fraction] = {}
    for name, speeds, expected in (
        ("core", CORE, EXPECTED_CORE),
        ("tail_1", TAIL_1, EXPECTED_TAIL_1),
        ("tail_2", TAIL_2, EXPECTED_TAIL_2),
    ):
        union_value = safe_measure_by_union(speeds)
        cell_value = safe_measure_by_cells(speeds)
        require(union_value == expected, f"{name} union mismatch")
        require(cell_value == expected, f"{name} cell mismatch")
        require(union_value == cell_value, f"{name} evaluators disagree")
        values[name] = union_value

    row_1 = values["core"] * values["tail_1"]
    row_2 = values["core"] * values["tail_2"]
    gap = row_2 - row_1

    require(row_1 == Fraction(883, 9240), "first row mismatch")
    require(row_2 == Fraction(7064, 72765), "second row mismatch")
    require(gap == Fraction(883, 582120), "gap mismatch")
    require(
        sum(Fraction(1, ENDPOINT_MODULUS * c) for c in TAIL_1)
        == sum(Fraction(1, ENDPOINT_MODULUS * c) for c in TAIL_2)
        == Fraction(5, 12 * ENDPOINT_MODULUS),
        "reciprocal sidecar mismatch",
    )

    print("THM-2182 ENDPOINT-GRID PRODUCT -- exact audit")
    print(f"endpoint_modulus={ENDPOINT_MODULUS}")
    print(f"core_measure={values['core']}")
    print(f"tail_1={TAIL_1}, tail_1_measure={values['tail_1']}")
    print(f"tail_2={TAIL_2}, tail_2_measure={values['tail_2']}")
    print(f"row_1_measure={row_1}")
    print(f"row_2_measure={row_2}")
    print(f"gap={gap}")
    print("evaluators=interval_union,boundary_cells")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
