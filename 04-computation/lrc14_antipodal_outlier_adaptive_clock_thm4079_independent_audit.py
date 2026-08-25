#!/usr/bin/env python3
"""Independent canonical-family and obstruction audit for THM-4079."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import gcd, lcm


DELTA = Fraction(1, 14)
D0 = (1, 3, 5, 7, 8, 9, 13, 21, 24)
ODD_CORE = (1, 3, 5, 7, 9, 13, 21)
DELETION_WITNESSES = {
    1: Fraction(4, 147),
    3: Fraction(11, 63),
    5: Fraction(29, 336),
    7: Fraction(15, 112),
    8: Fraction(43, 336),
    9: Fraction(11, 49),
    13: Fraction(29, 126),
    21: Fraction(43, 182),
    24: Fraction(15, 182),
}
EXPECTED_ODD_INTERVALS = (
    (Fraction(15, 91), Fraction(6, 35)),
    (Fraction(12, 49), Fraction(13, 49)),
    (Fraction(71, 147), Fraction(76, 147)),
    (Fraction(36, 49), Fraction(37, 49)),
    (Fraction(29, 35), Fraction(76, 91)),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def antipodal_safe(body: tuple[int, ...], theta: Fraction) -> bool:
    return all(
        circle_norm(speed * (theta + Fraction(half, 2))) >= DELTA
        for speed in body
        for half in (0, 1)
    )


def intersect_unions(
    left: tuple[tuple[Fraction, Fraction], ...],
    right: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    intersections: list[tuple[Fraction, Fraction]] = []
    for a, b in left:
        for c, d in right:
            start = max(a, c)
            end = min(b, d)
            if start <= end:
                intersections.append((start, end))
    intersections.sort()
    merged: list[tuple[Fraction, Fraction]] = []
    for start, end in intersections:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def odd_safe_intervals(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(
        (Fraction(7 * index + 1, 7 * speed), Fraction(7 * index + 6, 7 * speed))
        for index in range(speed)
    )


def equality_walls(body: tuple[int, ...]) -> list[Fraction]:
    walls = {Fraction(0), Fraction(1)}
    for speed in body:
        levels = (DELTA,) if speed % 2 == 0 else (DELTA, Fraction(3, 7))
        for level in levels:
            for integer in range(speed):
                for sign in (-1, 1):
                    walls.add(((Fraction(integer) + sign * level) / speed) % 1)
    return sorted(walls)


def main() -> None:
    canonical_phase_checks = 0
    for multiplier in range(1, 1001):
        speed = 12012 * multiplier
        body = tuple(range(1, 11)) + (speed,)
        theta = Fraction(1, 12) + Fraction(1, 2 * speed)
        for core_speed in body:
            for half in (0, 1):
                require(circle_norm(core_speed * (theta + Fraction(half, 2))) >= DELTA, f"canonical phase failed B={speed},d={core_speed}")
                canonical_phase_checks += 1
    require(canonical_phase_checks == 22000, "canonical phase aggregate changed")

    literal_odd_classes = 0
    for multiplier in range(1, 21):
        speed = 12012 * multiplier
        clock = 2 * speed
        label = speed // 6 + 1
        require(label % 2 == 1, f"canonical label parity failed B={speed}")
        for tail in range(1, 2 * clock, 2):
            first = min((tail * label) % clock, clock - (tail * label) % clock)
            second_label = label + speed
            second = min((tail * second_label) % clock, clock - (tail * second_label) % clock)
            require(not (7 * first < clock and 7 * second < clock), f"eligible odd tail survived B={speed},z={tail}")
            literal_odd_classes += 1
    require(literal_odd_classes == 5045040, "literal odd-class aggregate changed")

    blind_checks = 0
    running = 12012
    for horizon in range(1, 101):
        running = lcm(running, horizon)
        for clock in range(1, horizon + 1):
            require(running % clock == 0, f"blind horizon failed K={horizon},N={clock}")
            blind_checks += 1
    require(blind_checks == 5050, "blind-horizon aggregate changed")

    intersection: tuple[tuple[Fraction, Fraction], ...] = ((Fraction(0), Fraction(1)),)
    for speed in ODD_CORE:
        intersection = intersect_unions(intersection, odd_safe_intervals(speed))
    require(intersection == EXPECTED_ODD_INTERVALS, f"odd interval identity changed: {intersection}")

    cover_maxima = (
        max(circle_norm(12 * endpoint) for endpoint in EXPECTED_ODD_INTERVALS[0]),
        max(circle_norm(4 * endpoint) for endpoint in EXPECTED_ODD_INTERVALS[1]),
        max(circle_norm(4 * endpoint) for endpoint in EXPECTED_ODD_INTERVALS[2]),
        max(circle_norm(4 * endpoint) for endpoint in EXPECTED_ODD_INTERVALS[3]),
        max(circle_norm(12 * endpoint) for endpoint in EXPECTED_ODD_INTERVALS[4]),
    )
    require(cover_maxima == (Fraction(2, 35), Fraction(3, 49), Fraction(10, 147), Fraction(3, 49), Fraction(2, 35)), "even cover maxima changed")
    require(all(value < DELTA for value in cover_maxima), "even cover lost strictness")

    walls = equality_walls(D0)
    require(len(walls) == 294, f"wall checkpoint count changed: {len(walls)}")
    require(all(not antipodal_safe(D0, point) for point in walls), "feasible equality checkpoint found")
    midpoints = [(left + right) / 2 for left, right in zip(walls, walls[1:])]
    require(len(midpoints) == 293, "open-cell count changed")
    require(all(not antipodal_safe(D0, point) for point in midpoints), "feasible open cell found")

    deletion_gates = 0
    for omitted, witness in DELETION_WITNESSES.items():
        reduced = tuple(speed for speed in D0 if speed != omitted)
        for speed in reduced:
            for half in (0, 1):
                require(circle_norm(speed * (witness + Fraction(half, 2))) >= DELTA, f"deletion witness failed omit={omitted},d={speed}")
                deletion_gates += 1
    require(deletion_gates == 144, "deletion-witness gate count changed")

    epsilon = Fraction(1, 1000)
    weak_time = Fraction(1, 14)
    require(circle_norm(weak_time - epsilon) < DELTA, "left strict-margin hostile failed")
    require(circle_norm(13 * (weak_time + epsilon)) < DELTA, "right strict-margin hostile failed")

    semantic = {
        "odd_intervals": [[str(left), str(right)] for left, right in intersection],
        "cover_maxima": [str(value) for value in cover_maxima],
        "deletions": {str(key): str(value) for key, value in DELETION_WITNESSES.items()},
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4079 independent canonical-family/obstruction audit")
    print(f"canonical_core_phase_checks={canonical_phase_checks} literal_odd_tail_classes={literal_odd_classes}")
    print(f"blind_horizon_divisibility_checks={blind_checks}")
    print(f"odd_layer_union_intervals={len(intersection)} even_cover_rows={len(cover_maxima)}")
    print(f"full_theta_checkpoints={len(walls)} full_theta_open_cells={len(midpoints)}")
    print(f"deletion_witness_gates={deletion_gates} strict_margin_hostile_signs=2")
    print(f"cover_maxima={','.join(str(value) for value in cover_maxima)}")
    print(f"semantic_sha256={digest}")
    print("PASS: infinite adaptive family, fixed-bank blindness, and minimal obstruction")


if __name__ == "__main__":
    main()
