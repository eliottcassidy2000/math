#!/usr/bin/env python3
"""Exact rational circle geometry and three-lift masks for THM-4442.

The safe set is G_A={x in R/Z: ||a*x||>=1/14 for every a in A}.  Positive
components are reconstructed both as complements of merged open danger teeth
and as exact endpoint-arrangement cells, and the two representations must
agree.  The lift routines retain all three labelled sheets under y=3*x.
"""

from __future__ import annotations

from fractions import Fraction as Q


THRESHOLD = Q(1, 14)

def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def circle_distance(x: Q) -> Q:
    residue = x % 1
    return min(residue, 1 - residue)


def safe_at(speeds: tuple[int, ...], x: Q) -> bool:
    return all(circle_distance(v * x) >= THRESHOLD for v in speeds)


def danger_intervals(speed: int) -> tuple[tuple[Q, Q], ...]:
    """Positive-length clipped danger teeth on [0,1]."""
    answer: list[tuple[Q, Q]] = []
    radius = Q(1, 14 * speed)
    for k in range(speed + 1):
        left = max(Q(0), Q(k, speed) - radius)
        right = min(Q(1), Q(k, speed) + radius)
        if left < right:
            answer.append((left, right))
    return tuple(answer)


def merge_intervals(
    intervals: list[tuple[Q, Q]] | tuple[tuple[Q, Q], ...]
) -> tuple[tuple[Q, Q], ...]:
    """Merge for positive-measure geometry; isolated safe points are omitted."""
    merged: list[list[Q]] = []
    for left, right in sorted(intervals):
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return tuple((left, right) for left, right in merged)


def safe_components_by_union(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    danger = merge_intervals(
        [interval for speed in speeds for interval in danger_intervals(speed)]
    )
    safe: list[tuple[Q, Q]] = []
    cursor = Q(0)
    for left, right in danger:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        safe.append((cursor, Q(1)))
    return tuple(safe)


def speed_endpoints(speed: int) -> tuple[Q, ...]:
    endpoints = {Q(0), Q(1)}
    radius = Q(1, 14 * speed)
    for k in range(speed + 1):
        for sign in (-1, 1):
            point = Q(k, speed) + sign * radius
            if 0 <= point <= 1:
                endpoints.add(point)
    return tuple(sorted(endpoints))


def safe_components_by_arrangement(
    speeds: tuple[int, ...]
) -> tuple[tuple[Q, Q], ...]:
    endpoints = sorted(
        {point for speed in speeds for point in speed_endpoints(speed)}
    )
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        if safe_at(speeds, midpoint):
            if cells and cells[-1][1] == left:
                cells[-1] = (cells[-1][0], right)
            else:
                cells.append((left, right))
    return tuple(cells)


def interval_measure(intervals: tuple[tuple[Q, Q], ...]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def exact_safe_geometry(
    speeds: tuple[int, ...]
) -> tuple[tuple[tuple[Q, Q], ...], Q]:
    by_union = safe_components_by_union(speeds)
    by_arrangement = safe_components_by_arrangement(speeds)
    require(by_union == by_arrangement, f"geometry disagreement at {speeds}")
    return by_union, interval_measure(by_union)


def fmt_interval(interval: tuple[Q, Q]) -> str:
    left, right = interval
    return f"[{left},{right}]"


def endpoint_owners(speeds: tuple[int, ...], point: Q) -> tuple[str, ...]:
    owners: list[str] = []
    for speed in speeds:
        radius = Q(1, 14 * speed)
        for k in range(speed + 1):
            if point == Q(k, speed) - radius:
                owners.append(f"{speed}:{k}:L")
            if point == Q(k, speed) + radius:
                owners.append(f"{speed}:{k}:R")
    return tuple(owners)


def lift_safe_mask(y: Q, tails: tuple[int, int, int]) -> int:
    mask = 0
    for sheet in range(3):
        x = (y + sheet) / 3
        if safe_at(tails, x):
            mask |= 1 << sheet
    return mask


def pulled_tail_endpoints(tails: tuple[int, int, int]) -> set[Q]:
    """All y-walls induced by x_j=(y+j)/3 for the three tails."""
    endpoints = {Q(0), Q(1)}
    for tail in tails:
        for sheet in range(3):
            for k in range(tail + 1):
                for sign in (-1, 1):
                    y = Q(3, tail) * (Q(k) + sign * THRESHOLD) - sheet
                    if 0 <= y <= 1:
                        endpoints.add(y)
    return endpoints


def component_mask_pieces(
    components: tuple[tuple[Q, Q], ...], tails: tuple[int, int, int]
) -> tuple[tuple[int, Q, Q, int], ...]:
    walls = pulled_tail_endpoints(tails)
    answer: list[tuple[int, Q, Q, int]] = []
    for component_index, (component_left, component_right) in enumerate(components):
        local = sorted(
            {component_left, component_right}
            | {x for x in walls if component_left < x < component_right}
        )
        for left, right in zip(local, local[1:]):
            mask = lift_safe_mask((left + right) / 2, tails)
            if (
                answer
                and answer[-1][0] == component_index
                and answer[-1][2] == left
                and answer[-1][3] == mask
            ):
                old_index, old_left, _, old_mask = answer[-1]
                answer[-1] = (old_index, old_left, right, old_mask)
            else:
                answer.append((component_index, left, right, mask))
    return tuple(answer)


def first_safe_lift_cell(
    components: tuple[tuple[Q, Q], ...], tails: tuple[int, int, int]
) -> tuple[int, Q, Q, int] | None:
    """Return a positive quotient cell with at least one tail-safe lift."""
    walls = pulled_tail_endpoints(tails)
    for component_index, (component_left, component_right) in enumerate(components):
        local = sorted(
            {component_left, component_right}
            | {x for x in walls if component_left < x < component_right}
        )
        for left, right in zip(local, local[1:]):
            mask = lift_safe_mask((left + right) / 2, tails)
            if mask:
                return component_index, left, right, mask
    return None
