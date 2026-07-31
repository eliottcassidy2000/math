"""Small self-contained exact interval engine for the scratch j=7 audits."""

from __future__ import annotations

from fractions import Fraction as F


RULER = 14 * 360_360


def merge_intervals(
    rows: list[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    ordered = sorted((left, right) for left, right in rows if left < right)
    if not ordered:
        return ()
    merged: list[list[F]] = [[ordered[0][0], ordered[0][1]]]
    for left, right in ordered[1:]:
        if left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1][1] = right
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def carrier(body: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    danger: list[tuple[F, F]] = []
    for speed in body:
        radius = F(1, 14 * speed)
        danger.append((F(0), radius))
        danger.extend(
            (F(index, speed) - radius, F(index, speed) + radius)
            for index in range(1, speed)
        )
        danger.append((F(1) - radius, F(1)))
    union = merge_intervals(danger)
    safe: list[tuple[F, F]] = []
    cursor = F(0)
    for left, right in union:
        if cursor < left:
            safe.append((cursor, left))
        if right > cursor:
            cursor = right
    if cursor < 1:
        safe.append((cursor, F(1)))
    return tuple(safe)


def merge_integer_intervals(
    rows: list[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    ordered = sorted((left, right) for left, right in rows if left < right)
    if not ordered:
        return ()
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for left, right in ordered[1:]:
        if left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1][1] = right
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def integer_carrier(body: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    danger: list[tuple[int, int]] = []
    for speed in body:
        if RULER % (14 * speed):
            raise RuntimeError("inexact integer ruler")
        radius = RULER // (14 * speed)
        step = RULER // speed
        danger.append((0, radius))
        danger.extend(
            (index * step - radius, index * step + radius)
            for index in range(1, speed)
        )
        danger.append((RULER - radius, RULER))
    union = merge_integer_intervals(danger)
    safe: list[tuple[int, int]] = []
    cursor = 0
    for left, right in union:
        if cursor < left:
            safe.append((cursor, left))
        if right > cursor:
            cursor = right
    if cursor < RULER:
        safe.append((cursor, RULER))
    return tuple(safe)


def mass(rows: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in rows), F(0))


def unit_danger_primitive(position: F) -> F:
    """Integral from zero to ``position`` of the one-speed danger indicator."""
    cycles = position.numerator // position.denominator
    remainder = position - cycles
    return (
        F(cycles, 7)
        + min(remainder, F(1, 14))
        + max(F(0), remainder - F(13, 14))
    )


def singleton_coverage(
    rows: tuple[tuple[F, F], ...],
    label: int,
) -> F:
    return sum(
        (
            (
                unit_danger_primitive(label * right)
                - unit_danger_primitive(label * left)
            )
            / label
            for left, right in rows
        ),
        F(0),
    )


def integer_danger_primitive(numerator: int) -> int:
    """RULER times integral_0^(numerator/RULER) 1_{D_1}(s) ds."""
    cycles, remainder = divmod(numerator, RULER)
    return (
        cycles * (RULER // 7)
        + min(remainder, RULER // 14)
        + max(0, remainder - 13 * (RULER // 14))
    )


def integer_singleton_coverage(
    rows: tuple[tuple[int, int], ...],
    label: int,
) -> F:
    numerator = sum(
        integer_danger_primitive(label * right)
        - integer_danger_primitive(label * left)
        for left, right in rows
    )
    return F(numerator, RULER * label)
