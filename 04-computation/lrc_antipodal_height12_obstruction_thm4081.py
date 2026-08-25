#!/usr/bin/env python3
"""Fraction-exact transformed-circle interval audit for THM-4081.

This is the primary audit.  It works only in the doubled coordinate t=2 theta:
an odd speed d asks ||d t|| >= 1/7, while an even speed d=2e asks
||e t|| >= 1/14.  Closed rational interval intersections retain isolated
endpoint solutions; this is essential for the eight-speed hostile example.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction


DELTA = Fraction(1, 14)
FULL = tuple(range(1, 13))
DSTAR = (1, 3, 4, 5, 7, 8, 9, 11, 12)
ODD_CORE = (1, 3, 5, 7, 9, 11)
OPTIONAL_SHADOWS = (2, 6, 10)
ENDPOINT_HOSTILE = (1, 3, 5, 8, 11, 13, 23, 36)
DELETION_T = {
    1: Fraction(1, 14),
    3: Fraction(15, 49),
    4: Fraction(36, 77),
    5: Fraction(8, 21),
    7: Fraction(1, 7),
    8: Fraction(5, 21),
    9: Fraction(8, 35),
    11: Fraction(29, 63),
    12: Fraction(8, 49),
}
EXPECTED_STAGES = (
    ((Fraction(1, 7), Fraction(6, 7)),),
    (
        (Fraction(1, 7), Fraction(2, 7)),
        (Fraction(8, 21), Fraction(13, 21)),
        (Fraction(5, 7), Fraction(6, 7)),
    ),
    (
        (Fraction(1, 7), Fraction(6, 35)),
        (Fraction(8, 35), Fraction(2, 7)),
        (Fraction(3, 7), Fraction(4, 7)),
        (Fraction(5, 7), Fraction(27, 35)),
        (Fraction(29, 35), Fraction(6, 7)),
    ),
    (
        (Fraction(8, 49), Fraction(6, 35)),
        (Fraction(8, 35), Fraction(13, 49)),
        (Fraction(22, 49), Fraction(27, 49)),
        (Fraction(36, 49), Fraction(27, 35)),
        (Fraction(29, 35), Fraction(41, 49)),
    ),
    (
        (Fraction(8, 49), Fraction(6, 35)),
        (Fraction(5, 21), Fraction(13, 49)),
        (Fraction(29, 63), Fraction(34, 63)),
        (Fraction(36, 49), Fraction(16, 21)),
        (Fraction(29, 35), Fraction(41, 49)),
    ),
    (
        (Fraction(8, 49), Fraction(13, 77)),
        (Fraction(5, 21), Fraction(20, 77)),
        (Fraction(36, 77), Fraction(41, 77)),
        (Fraction(57, 77), Fraction(16, 21)),
        (Fraction(64, 77), Fraction(41, 49)),
    ),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def transformed_multiplier(speed: int) -> int:
    return speed if speed % 2 else speed // 2


def transformed_threshold(speed: int) -> Fraction:
    return Fraction(1, 7) if speed % 2 else DELTA


def transformed_safe(speed: int, time: Fraction) -> bool:
    return circle_norm(transformed_multiplier(speed) * time) >= transformed_threshold(speed)


def safe_intervals(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    multiplier = transformed_multiplier(speed)
    denominator = 7 if speed % 2 else 14
    return tuple(
        (
            Fraction(denominator * index + 1, denominator * multiplier),
            Fraction(denominator * index + denominator - 1, denominator * multiplier),
        )
        for index in range(multiplier)
    )


def intersect_unions(
    left: tuple[tuple[Fraction, Fraction], ...],
    right: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    pieces: list[tuple[Fraction, Fraction]] = []
    for left_start, left_end in left:
        for right_start, right_end in right:
            start = max(left_start, right_start)
            end = min(left_end, right_end)
            if start <= end:
                pieces.append((start, end))
    pieces.sort()
    merged: list[tuple[Fraction, Fraction]] = []
    for start, end in pieces:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def joint_safe_set(body: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    safe: tuple[tuple[Fraction, Fraction], ...] = ((Fraction(0), Fraction(1)),)
    for speed in body:
        safe = intersect_unions(safe, safe_intervals(speed))
        if not safe:
            break
    return safe


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def max_norm_on_interval(multiplier: int, start: Fraction, end: Fraction) -> Fraction:
    candidates = [start, end]
    first = ceil_fraction(2 * multiplier * start)
    last = (2 * multiplier * end).numerator // (2 * multiplier * end).denominator
    for index in range(first, last + 1):
        candidates.append(Fraction(index, 2 * multiplier))
    return max(circle_norm(multiplier * point) for point in candidates)


def body_mask(body: tuple[int, ...]) -> int:
    mask = 0
    for speed in body:
        mask |= 1 << (speed - 1)
    return mask


def main() -> None:
    require(DSTAR == tuple(speed for speed in FULL if speed % 4 != 2), "canonical residue description changed")
    require(ODD_CORE == tuple(speed for speed in FULL if speed % 2), "odd core changed")
    require(OPTIONAL_SHADOWS == tuple(speed for speed in FULL if speed % 4 == 2), "shadow set changed")

    parity_shadow_gates = 0
    for even_speed in OPTIONAL_SHADOWS:
        odd_speed = even_speed // 2
        require(transformed_multiplier(even_speed) == transformed_multiplier(odd_speed), "shadow multiplier mismatch")
        require(transformed_threshold(even_speed) < transformed_threshold(odd_speed), "shadow threshold not weaker")
        require(
            all(
                any(outer_start <= inner_start and inner_end <= outer_end for outer_start, outer_end in safe_intervals(even_speed))
                for inner_start, inner_end in safe_intervals(odd_speed)
            ),
            "exact parity-shadow interval containment failed",
        )
        parity_shadow_gates += 3

    stages: list[tuple[tuple[Fraction, Fraction], ...]] = []
    running: tuple[tuple[Fraction, Fraction], ...] = ((Fraction(0), Fraction(1)),)
    generated_odd_intervals = 0
    for speed in ODD_CORE:
        layer = safe_intervals(speed)
        generated_odd_intervals += len(layer)
        running = intersect_unions(running, layer)
        stages.append(running)
    require(tuple(stages) == EXPECTED_STAGES, f"odd-core interval anatomy changed: {stages}")
    require(generated_odd_intervals == 36, "odd interval count changed")

    final_intervals = EXPECTED_STAGES[-1]
    cover_speeds = (12, 8, 4, 8, 12)
    cover_maxima = tuple(
        max_norm_on_interval(transformed_multiplier(speed), start, end)
        for speed, (start, end) in zip(cover_speeds, final_intervals)
    )
    expected_maxima = (
        Fraction(1, 49),
        Fraction(1, 21),
        Fraction(5, 77),
        Fraction(1, 21),
        Fraction(1, 49),
    )
    require(cover_maxima == expected_maxima, f"strict-cover maxima changed: {cover_maxima}")
    require(all(value < DELTA for value in cover_maxima), "an even cover lost strictness")
    require(not joint_safe_set(DSTAR), "D* ceased to obstruct")

    deletion_survivor_gates = 0
    omitted_failure_values: list[Fraction] = []
    for omitted, witness in DELETION_T.items():
        for speed in FULL:
            if speed != omitted:
                require(transformed_safe(speed, witness), f"strong deletion witness failed omit={omitted},d={speed}")
                deletion_survivor_gates += 1
        omitted_value = circle_norm(transformed_multiplier(omitted) * witness)
        require(omitted_value < transformed_threshold(omitted), f"omitted speed did not strictly fail: {omitted}")
        omitted_failure_values.append(omitted_value)
    require(deletion_survivor_gates == 99, "deletion survivor count changed")

    dstar_mask = body_mask(DSTAR)
    expected_obstruction_masks = tuple(
        dstar_mask | body_mask(tuple(OPTIONAL_SHADOWS[index] for index in range(3) if optional_bits & (1 << index)))
        for optional_bits in range(8)
    )
    obstruction_masks: list[int] = []
    histogram = {size: 0 for size in range(13)}
    for mask in range(1 << 12):
        body = tuple(speed for speed in FULL if mask & (1 << (speed - 1)))
        if not joint_safe_set(body):
            obstruction_masks.append(mask)
            histogram[len(body)] += 1
    require(tuple(obstruction_masks) == tuple(sorted(expected_obstruction_masks)), "[1,12] obstruction classification changed")
    require({size: count for size, count in histogram.items() if count} == {9: 1, 10: 3, 11: 3, 12: 1}, "obstruction histogram changed")

    expected_hostile_points = tuple((Fraction(index, 7), Fraction(index, 7)) for index in range(1, 7))
    hostile_safe_set = joint_safe_set(ENDPOINT_HOSTILE)
    require(hostile_safe_set == expected_hostile_points, f"endpoint hostile changed: {hostile_safe_set}")

    odd_dilation_identity_gates = 0
    odd_dilation_witness_gates = 0
    for dilation in (1, 3, 5, 7, 9):
        for speed in FULL:
            require((dilation * speed) % 2 == speed % 2, "odd dilation changed parity")
            require(transformed_multiplier(dilation * speed) == dilation * transformed_multiplier(speed), "odd dilation multiplier identity failed")
            require(transformed_threshold(dilation * speed) == transformed_threshold(speed), "odd dilation threshold identity failed")
            odd_dilation_identity_gates += 3
        for omitted, witness in DELETION_T.items():
            lifted = witness / dilation
            for speed in FULL:
                if speed != omitted:
                    require(transformed_safe(dilation * speed, lifted), "odd-dilated strong witness failed")
                    odd_dilation_witness_gates += 1

    even_dilation_time = Fraction(1, 13)
    even_dilation_safe_gates = 0
    for speed in DSTAR:
        doubled = 2 * speed
        require(transformed_safe(doubled, even_dilation_time), f"even-dilation witness failed d={doubled}")
        require(circle_norm(speed * even_dilation_time) >= Fraction(1, 13), "even-dilation margin changed")
        even_dilation_safe_gates += 1

    semantic = {
        "Dstar": list(DSTAR),
        "cover_maxima": [str(value) for value in cover_maxima],
        "deletion_t": {str(speed): str(time) for speed, time in DELETION_T.items()},
        "final_intervals": [[str(start), str(end)] for start, end in final_intervals],
        "hostile_t": [str(start) for start, end in hostile_safe_set if start == end],
        "obstruction_masks": obstruction_masks,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4081 primary transformed-circle interval audit")
    print(f"canonical_Dstar={','.join(str(speed) for speed in DSTAR)} parity_shadow_gates={parity_shadow_gates}")
    print(f"odd_layers={len(ODD_CORE)} generated_odd_intervals={generated_odd_intervals} stage_component_counts={','.join(str(len(stage)) for stage in stages)}")
    print("final_odd_intervals=" + ";".join(f"[{start},{end}]" for start, end in final_intervals))
    print("strict_even_cover=" + ";".join(f"d{speed}:{maximum}" for speed, maximum in zip(cover_speeds, cover_maxima)))
    print(f"deletion_survivor_gates={deletion_survivor_gates} omitted_strict_failures={len(omitted_failure_values)}")
    print(f"subsets_classified={1 << 12} obstruction_count={len(obstruction_masks)} cardinality_histogram=9:1,10:3,11:3,12:1")
    print(f"endpoint_hostile_speeds={len(ENDPOINT_HOSTILE)} exact_safe_t_singletons={len(hostile_safe_set)} positive_safe_intervals=0")
    print(f"odd_dilation_identity_gates={odd_dilation_identity_gates} odd_dilation_witness_gates={odd_dilation_witness_gates}")
    print(f"even_dilation_safe_gates={even_dilation_safe_gates} even_dilation_t={even_dilation_time}")
    print(f"semantic_sha256={digest}")
    print("PASS: D* obstruction, strong minimality, [1,12] iff classification, and dilation boundary")


if __name__ == "__main__":
    main()
