#!/usr/bin/env python3
"""General exact construction audit for THM-4079.

The script stress-tests continuous one-outlier absorption and both adaptive
odd-grid clock constructions over a deterministic universe of strict rational
antipodal bodies. All arithmetic is Fraction-exact.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import gcd


DELTA = Fraction(1, 14)
BODY_SPEED_MAX = 40
DENOMINATOR_MAX = 30
CASES_PER_BODY = 21


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional(value: Fraction) -> Fraction:
    return value % 1


def circle_norm(value: Fraction) -> Fraction:
    residue = fractional(value)
    return min(residue, 1 - residue)


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def in_band(x: Fraction, intervals: tuple[tuple[Fraction, Fraction], ...]) -> bool:
    residue = fractional(x)
    return any(left <= residue <= right for left, right in intervals)


def nearest_band_shift(x: Fraction, intervals: tuple[tuple[Fraction, Fraction], ...]) -> Fraction:
    residue = fractional(x)
    if in_band(residue, intervals):
        return Fraction(0)
    candidates: list[Fraction] = []
    for left, right in intervals:
        for endpoint in (left, right):
            for winding in (-1, 0, 1):
                candidates.append(endpoint + winding - residue)
    return min(candidates, key=lambda value: (abs(value), value))


def parity_band(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    if speed % 2 == 0:
        return ((DELTA, 1 - DELTA),)
    return ((DELTA, Fraction(3, 7)), (Fraction(4, 7), 1 - DELTA))


def nearest_odd(value: Fraction) -> int:
    floor_value = value.numerator // value.denominator
    candidates = [integer for integer in range(floor_value - 2, floor_value + 4) if integer % 2]
    return min(candidates, key=lambda integer: (abs(Fraction(integer) - value), integer))


def safe_pair(body: tuple[int, ...], time: Fraction) -> int:
    gates = 0
    for speed in body:
        for half in (0, 1):
            require(circle_norm(speed * (time + Fraction(half, 2))) >= DELTA, f"unsafe phase speed={speed},t={time},half={half}")
            gates += 1
    return gates


def strict_body(time: Fraction) -> tuple[tuple[int, ...], Fraction]:
    body = tuple(
        speed
        for speed in range(1, BODY_SPEED_MAX + 1)
        if min(
            circle_norm(speed * time),
            circle_norm(speed * (time + Fraction(1, 2))),
        ) > DELTA
    )
    require(body, f"empty strict body at t={time}")
    margin = min(
        circle_norm(speed * (time + Fraction(half, 2))) - DELTA
        for speed in body
        for half in (0, 1)
    )
    require(margin > 0, f"nonpositive strict margin at t={time}")
    return body, margin


def main() -> None:
    tagged_times = [
        Fraction(numerator, denominator)
        for denominator in range(3, DENOMINATOR_MAX + 1)
        for numerator in range(1, (denominator + 1) // 2)
        if gcd(numerator, denominator) == 1
    ]
    require(len(tagged_times) == 138, "strict rational body universe changed")

    continuous_cases = 0
    clock4_cases = 0
    clock2_cases = 0
    safe_phase_gates = 0
    body_rows: list[dict[str, int | str]] = []
    max_scaled_continuous_shift = Fraction(0)
    max_scaled_clock4_shift = Fraction(0)

    for time in tagged_times:
        body, margin = strict_body(time)
        radius = max(body)
        continuous_start = max(BODY_SPEED_MAX + 1, ceil_fraction(Fraction(radius, 14) / margin))
        clock4_start = max(BODY_SPEED_MAX + 1, ceil_fraction(Fraction(radius, 4) / margin))
        clock2_start = max(BODY_SPEED_MAX + 1, ceil_fraction(Fraction(radius, 2) / margin))
        if clock2_start % 2:
            clock2_start += 1

        for speed in range(continuous_start, continuous_start + CASES_PER_BODY):
            shift_phase = nearest_band_shift(speed * time, parity_band(speed))
            require(abs(shift_phase) <= DELTA, f"band-distance failed at t={time},B={speed}")
            shifted_time = time + shift_phase / speed
            safe_phase_gates += safe_pair(body + (speed,), shifted_time)
            max_scaled_continuous_shift = max(max_scaled_continuous_shift, abs(speed * (shifted_time - time)))
            continuous_cases += 1

        for speed in range(clock4_start, clock4_start + CASES_PER_BODY):
            clock = 4 * speed
            label = nearest_odd(clock * time)
            shifted_time = Fraction(label, clock)
            require(abs(label - clock * time) <= 1, f"odd quarter-grid approximation failed at t={time},B={speed}")
            safe_phase_gates += safe_pair(body + (speed,), shifted_time)
            require(label % clock != (label + 2 * speed) % clock, "quarter-grid labels collapsed")
            max_scaled_clock4_shift = max(max_scaled_clock4_shift, abs(speed * (shifted_time - time)))
            clock4_cases += 1

        for offset in range(CASES_PER_BODY):
            speed = clock2_start + 2 * offset
            clock = 2 * speed
            label = nearest_odd(clock * time)
            shifted_time = Fraction(label, clock)
            require(abs(label - clock * time) <= 1, f"odd half-grid approximation failed at t={time},B={speed}")
            safe_phase_gates += safe_pair(body + (speed,), shifted_time)
            require(speed % 2 == 0, "half-grid parity lost")
            clock2_cases += 1

        body_rows.append(
            {
                "t": str(time),
                "size": len(body),
                "R": radius,
                "eta": str(margin),
                "continuous_start": continuous_start,
                "clock4_start": clock4_start,
                "clock2_start": clock2_start,
            }
        )

    mesh_points = 0
    max_even_band_distance = Fraction(0)
    max_odd_band_distance = Fraction(0)
    for denominator in range(1, 101):
        for numerator in range(denominator):
            phase = Fraction(numerator, denominator)
            max_even_band_distance = max(max_even_band_distance, abs(nearest_band_shift(phase, parity_band(2))))
            max_odd_band_distance = max(max_odd_band_distance, abs(nearest_band_shift(phase, parity_band(1))))
            mesh_points += 1

    require(continuous_cases == 2898, "continuous case aggregate changed")
    require(clock4_cases == 2898, "4B case aggregate changed")
    require(clock2_cases == 2898, "2B case aggregate changed")
    require(safe_phase_gates == 568512, f"safe-phase aggregate changed: {safe_phase_gates}")
    require(mesh_points == 5050, "mesh aggregate changed")
    require(max_even_band_distance == DELTA, "even band covering radius changed")
    require(max_odd_band_distance == DELTA, "odd band covering radius changed")
    require(max_scaled_continuous_shift == DELTA, "continuous shift boundary changed")
    require(max_scaled_clock4_shift == Fraction(1, 4), "quarter-grid shift boundary changed")

    digest = hashlib.sha256(
        json.dumps(body_rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4079 primary general-construction audit")
    print(f"strict_rational_bodies={len(tagged_times)}")
    print(f"continuous_threshold_cases={continuous_cases} clock_4B_threshold_cases={clock4_cases} even_clock_2B_threshold_cases={clock2_cases}")
    print(f"all_safe_phase_inequalities={safe_phase_gates} phase_band_mesh_points={mesh_points}")
    print(f"max_band_distances=even:{max_even_band_distance},odd:{max_odd_band_distance}")
    print(f"max_scaled_shifts=continuous:{max_scaled_continuous_shift},clock4B:{max_scaled_clock4_shift}")
    print(f"semantic_sha256={digest}")
    print("PASS: continuous absorption and adaptive 4B/2B constructions")


if __name__ == "__main__":
    main()
