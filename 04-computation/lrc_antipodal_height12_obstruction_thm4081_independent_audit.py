#!/usr/bin/env python3
"""Independent original-theta wall/cell audit for THM-4081.

Unlike the primary script, this audit never doubles theta and never uses a
parity-dependent multiplier.  It constructs every equality wall of the two
literal constraints ||d theta|| >= 1/14 and ||d(theta+1/2)|| >= 1/14,
tests all walls and all intervening open cells, and classifies all subsets of
{1,...,12} by their exact wall/cell failure masks.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction


DELTA = Fraction(1, 14)
FULL = tuple(range(1, 13))
DSTAR = (1, 3, 4, 5, 7, 8, 9, 11, 12)
OPTIONAL_SHADOWS = (2, 6, 10)
ENDPOINT_HOSTILE = (1, 3, 5, 8, 11, 13, 23, 36)
DELETION_THETA = {
    1: Fraction(1, 28),
    3: Fraction(15, 98),
    4: Fraction(18, 77),
    5: Fraction(4, 21),
    7: Fraction(1, 14),
    8: Fraction(5, 42),
    9: Fraction(4, 35),
    11: Fraction(29, 126),
    12: Fraction(4, 49),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_values(speed: int, theta: Fraction) -> tuple[Fraction, Fraction]:
    return (
        circle_norm(speed * theta),
        circle_norm(speed * (theta + Fraction(1, 2))),
    )


def literal_safe(speed: int, theta: Fraction) -> bool:
    first, second = literal_values(speed, theta)
    return first >= DELTA and second >= DELTA


def equality_walls(body: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls = {Fraction(0), Fraction(1)}
    for speed in body:
        for half in (0, 1):
            for integer in range(speed):
                for sign in (-1, 1):
                    wall = ((Fraction(integer) + sign * DELTA) / speed - Fraction(half, 2)) % 1
                    walls.add(wall)
    return tuple(sorted(walls))


def wall_and_cell_points(body: tuple[int, ...]) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    walls = equality_walls(body)
    midpoints = tuple((left + right) / 2 for left, right in zip(walls, walls[1:]))
    return walls, midpoints


def failure_masks(
    body: tuple[int, ...],
    points: tuple[Fraction, ...],
    bit_positions: dict[int, int],
) -> tuple[tuple[int, ...], int]:
    masks: list[int] = []
    phase_gates = 0
    for theta in points:
        mask = 0
        for speed in body:
            values = literal_values(speed, theta)
            phase_gates += 2
            if values[0] < DELTA or values[1] < DELTA:
                mask |= 1 << bit_positions[speed]
        masks.append(mask)
    return tuple(masks), phase_gates


def full_mask(body: tuple[int, ...]) -> int:
    mask = 0
    for speed in body:
        mask |= 1 << (speed - 1)
    return mask


def main() -> None:
    dstar_walls, dstar_midpoints = wall_and_cell_points(DSTAR)
    require(len(dstar_walls) == 186, f"D* wall count changed: {len(dstar_walls)}")
    require(len(dstar_midpoints) == 185, "D* cell count changed")
    dstar_points = dstar_walls[:-1] + dstar_midpoints
    dstar_bits = {speed: index for index, speed in enumerate(DSTAR)}
    dstar_failures, dstar_phase_gates = failure_masks(DSTAR, dstar_points, dstar_bits)
    require(all(mask != 0 for mask in dstar_failures), "D* has a safe original-theta wall or cell")
    require(dstar_phase_gates == 6660, "D* phase-gate count changed")

    full_walls, full_midpoints = wall_and_cell_points(FULL)
    require(len(full_walls) == 222, f"full wall count changed: {len(full_walls)}")
    require(len(full_midpoints) == 221, "full cell count changed")
    full_points = full_walls[:-1] + full_midpoints
    full_bits = {speed: speed - 1 for speed in FULL}
    full_failures, full_phase_gates = failure_masks(FULL, full_points, full_bits)
    require(full_phase_gates == 10608, f"full phase-gate count changed: {full_phase_gates}")

    dstar_mask = full_mask(DSTAR)
    expected_masks = tuple(
        dstar_mask | full_mask(tuple(OPTIONAL_SHADOWS[index] for index in range(3) if optional_bits & (1 << index)))
        for optional_bits in range(8)
    )
    obstruction_masks: list[int] = []
    histogram = {size: 0 for size in range(13)}
    for subset_mask in range(1 << 12):
        is_obstruction = all(subset_mask & failed_mask for failed_mask in full_failures)
        if is_obstruction:
            obstruction_masks.append(subset_mask)
            histogram[subset_mask.bit_count()] += 1
    require(tuple(obstruction_masks) == tuple(sorted(expected_masks)), "direct [1,12] iff classification changed")
    require({size: count for size, count in histogram.items() if count} == {9: 1, 10: 3, 11: 3, 12: 1}, "direct histogram changed")

    deletion_survivor_phase_gates = 0
    deletion_strict_failures = 0
    for omitted, theta in DELETION_THETA.items():
        for speed in FULL:
            values = literal_values(speed, theta)
            if speed == omitted:
                require(min(values) < DELTA, f"omitted speed did not fail directly: {omitted}")
                deletion_strict_failures += 1
            else:
                require(values[0] >= DELTA and values[1] >= DELTA, f"direct deletion witness failed omit={omitted},d={speed}")
                deletion_survivor_phase_gates += 2
    require(deletion_survivor_phase_gates == 198, "direct deletion survivor count changed")
    require(deletion_strict_failures == 9, "direct deletion failure count changed")

    hostile_walls, hostile_midpoints = wall_and_cell_points(ENDPOINT_HOSTILE)
    require(len(hostile_walls) == 302, f"hostile wall count changed: {len(hostile_walls)}")
    require(len(hostile_midpoints) == 301, "hostile cell count changed")
    hostile_bits = {speed: index for index, speed in enumerate(ENDPOINT_HOSTILE)}
    hostile_circle_walls = hostile_walls[:-1]
    hostile_wall_failures, hostile_wall_phase_gates = failure_masks(ENDPOINT_HOSTILE, hostile_circle_walls, hostile_bits)
    hostile_midpoint_failures, hostile_midpoint_phase_gates = failure_masks(ENDPOINT_HOSTILE, hostile_midpoints, hostile_bits)
    hostile_safe_walls = tuple(
        theta for theta, mask in zip(hostile_circle_walls, hostile_wall_failures) if mask == 0
    )
    expected_hostile_theta = tuple(
        Fraction(index, 14) for index in tuple(range(1, 7)) + tuple(range(8, 14))
    )
    require(hostile_safe_walls == expected_hostile_theta, f"hostile singleton set changed: {hostile_safe_walls}")
    require(all(mask != 0 for mask in hostile_midpoint_failures), "hostile acquired a safe open cell")
    require(hostile_wall_phase_gates + hostile_midpoint_phase_gates == 9632, "hostile phase-gate count changed")

    odd_dilation_identity_gates = 0
    odd_dilation_survivor_phase_gates = 0
    odd_dilation_strict_failures = 0
    for dilation in (3, 5, 7):
        for speed in DSTAR:
            for half in (0, 1):
                offset_numerator = speed * half * (dilation - 1)
                require(offset_numerator % 2 == 0, "odd-dilation offset lost integrality")
                for theta in DELETION_THETA.values():
                    dilated_value = circle_norm(dilation * speed * (theta / dilation + Fraction(half, 2)))
                    original_value = circle_norm(speed * (theta + Fraction(half, 2)))
                    require(dilated_value == original_value, "literal odd-dilation identity failed")
                    odd_dilation_identity_gates += 1
        for omitted, theta in DELETION_THETA.items():
            lifted = theta / dilation
            for speed in DSTAR:
                values = literal_values(dilation * speed, lifted)
                if speed == omitted:
                    require(min(values) < DELTA, "odd-dilated omitted speed did not fail")
                    odd_dilation_strict_failures += 1
                else:
                    require(values[0] >= DELTA and values[1] >= DELTA, "odd-dilated direct witness failed")
                    odd_dilation_survivor_phase_gates += 2
    require(odd_dilation_identity_gates == 486, "odd-dilation identity count changed")
    require(odd_dilation_survivor_phase_gates == 432, "odd-dilation survivor count changed")
    require(odd_dilation_strict_failures == 27, "odd-dilation failure count changed")

    even_dilation_theta = Fraction(1, 26)
    even_dilation_phase_gates = 0
    for speed in DSTAR:
        values = literal_values(2 * speed, even_dilation_theta)
        require(values[0] >= DELTA and values[1] >= DELTA, f"direct even-dilation witness failed d={2 * speed}")
        even_dilation_phase_gates += 2
    require(even_dilation_phase_gates == 18, "even-dilation gate count changed")

    semantic = {
        "Dstar_cells": len(dstar_midpoints),
        "Dstar_walls": len(dstar_walls),
        "deletion_theta": {str(speed): str(theta) for speed, theta in DELETION_THETA.items()},
        "full_cells": len(full_midpoints),
        "full_walls": len(full_walls),
        "hostile_theta": [str(theta) for theta in hostile_safe_walls],
        "obstruction_masks": obstruction_masks,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4081 independent literal original-theta wall/cell audit")
    print(f"Dstar_wall_endpoints={len(dstar_walls)} Dstar_circle_walls={len(dstar_walls) - 1} Dstar_open_cells={len(dstar_midpoints)} Dstar_test_points={len(dstar_points)} Dstar_phase_gates={dstar_phase_gates}")
    print(f"full_wall_endpoints={len(full_walls)} full_circle_walls={len(full_walls) - 1} full_open_cells={len(full_midpoints)} full_test_points={len(full_points)} full_phase_gates={full_phase_gates}")
    print(f"subsets_classified={1 << 12} obstruction_count={len(obstruction_masks)} cardinality_histogram=9:1,10:3,11:3,12:1")
    print(f"deletion_survivor_phase_gates={deletion_survivor_phase_gates} deletion_strict_failures={deletion_strict_failures}")
    print(f"endpoint_hostile_wall_endpoints={len(hostile_walls)} endpoint_hostile_circle_walls={len(hostile_walls) - 1} endpoint_hostile_open_cells={len(hostile_midpoints)} safe_walls={len(hostile_safe_walls)} safe_open_cells=0 hostile_phase_gates={hostile_wall_phase_gates + hostile_midpoint_phase_gates}")
    print(f"odd_dilation_identity_gates={odd_dilation_identity_gates} odd_dilation_survivor_phase_gates={odd_dilation_survivor_phase_gates} odd_dilation_strict_failures={odd_dilation_strict_failures}")
    print(f"even_dilation_phase_gates={even_dilation_phase_gates} even_dilation_theta={even_dilation_theta}")
    print(f"semantic_sha256={digest}")
    print("PASS: independent literal-theta obstruction, iff subset classification, endpoints, and dilation boundary")


if __name__ == "__main__":
    main()
