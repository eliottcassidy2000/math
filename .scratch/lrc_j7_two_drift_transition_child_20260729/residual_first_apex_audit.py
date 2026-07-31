#!/usr/bin/env python3
"""Exact first-apex audit after deleting the bounded first drift.

Scratch theorem probe for the five-aligned/two-drift face of the literal
six-body LRC(14) carrier.  The universe and every comparison are exact:

* ``E`` runs through all six-subsets of ``{1,...,14}``;
* ``z1`` runs through all external nonaligned labels from 15 through the
  conservative inclusive rootwise multiplicity-excess/BV cap;
* ``delta_E(z1) >= (88/1365) h_E`` selects the rows not already given a
  finite ``z2`` cap by the scalar excess gap;
* the projected-safe wall gives
  ``z2 > (2275/18627) L_E``;
* THM-2893 on ``R= C_E \\ D_z1`` says one of the six remaining tails has
  label at most ``floor(36 r_R/(7 h_R))``.

The five aligned labels are multiples of ``L_E``.  Comparing those two
typed lower bounds with the residual first-apex cap is therefore a rigorous
row closure, not a heuristic search horizon.
"""

from __future__ import annotations

import hashlib
import math
from fractions import Fraction as F
from itertools import combinations


RULER = 14 * 360_360
BASE_LABEL = 15
ETA = F(88, 1365)
ALPHA = F(2275, 18627)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def merge_integer(
    rows: list[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    rows = sorted((left, right) for left, right in rows if left < right)
    if not rows:
        return ()
    out: list[list[int]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= out[-1][1]:
            out[-1][1] = max(out[-1][1], right)
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


def integer_carrier(body: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    danger: list[tuple[int, int]] = []
    for speed in body:
        require(RULER % (14 * speed) == 0, "integer ruler lost a body speed")
        radius = RULER // (14 * speed)
        step = RULER // speed
        danger.append((0, radius))
        danger.extend(
            (index * step - radius, index * step + radius)
            for index in range(1, speed)
        )
        danger.append((RULER - radius, RULER))
    occupied = merge_integer(danger)
    safe: list[tuple[int, int]] = []
    cursor = 0
    for left, right in occupied:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < RULER:
        safe.append((cursor, RULER))
    return tuple(safe)


def integer_primitive(numerator: int) -> int:
    cycles, residue = divmod(numerator, RULER)
    return (
        cycles * (RULER // 7)
        + min(residue, RULER // 14)
        + max(0, residue - 13 * (RULER // 14))
    )


def singleton_coverage(
    carrier: tuple[tuple[int, int], ...],
    speed: int,
) -> F:
    numerator = sum(
        integer_primitive(speed * right) - integer_primitive(speed * left)
        for left, right in carrier
    )
    return F(numerator, RULER * speed)


def merge_fraction(
    rows: list[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    rows = sorted((left, right) for left, right in rows if left < right)
    if not rows:
        return ()
    out: list[list[F]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= out[-1][1]:
            out[-1][1] = max(out[-1][1], right)
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


def danger(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, 14 * speed)
    rows: list[tuple[F, F]] = []
    for index in range(speed):
        center = F(index, speed)
        left, right = center - radius, center + radius
        if left < 0:
            rows.extend(((F(0), right), (left + 1, F(1))))
        elif right > 1:
            rows.extend(((left, F(1)), (F(0), right - 1)))
        else:
            rows.append((left, right))
    return merge_fraction(rows)


def subtract(
    carrier: tuple[tuple[F, F], ...],
    removed: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    first_removed = 0
    for left, right in carrier:
        cursor = left
        while (
            first_removed < len(removed)
            and removed[first_removed][1] <= left
        ):
            first_removed += 1
        index = first_removed
        while index < len(removed) and removed[index][0] < right:
            bad_left, bad_right = removed[index]
            if cursor < bad_left:
                out.append((cursor, min(right, bad_left)))
            cursor = max(cursor, bad_right)
            if cursor >= right:
                break
            index += 1
        if cursor < right:
            out.append((cursor, right))
    return tuple(out)


def fraction_carrier(
    integer_rows: tuple[tuple[int, int], ...],
) -> tuple[tuple[F, F], ...]:
    return tuple((F(left, RULER), F(right, RULER)) for left, right in integer_rows)


def mass(rows: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in rows), F(0))


def inclusive_integer_cap(bound: F) -> int:
    """Conservative integer cap for a non-strict analytic upper bound."""
    return bound.numerator // bound.denominator


def main() -> None:
    candidate_rows = 0
    high_excess: list[tuple[object, ...]] = []
    survivors: list[tuple[object, ...]] = []

    for body in combinations(range(1, 15), 6):
        carrier_i = integer_carrier(body)
        carrier_f = fraction_carrier(carrier_i)
        h = F(sum(right - left for left, right in carrier_i), RULER)
        components = len(carrier_i)
        canonical_l = 14 * math.lcm(*body)
        first_bound = F(12 * components, 49) / (ETA * h)
        first_cap = inclusive_integer_cap(first_bound)

        for first in range(BASE_LABEL, first_cap + 1):
            if first % canonical_l == 0:
                continue
            candidate_rows += 1
            delta = singleton_coverage(carrier_i, first) - h / 7
            if delta < ETA * h:
                continue

            residual = subtract(carrier_f, danger(first))
            residual_h = mass(residual)
            residual_r = len(residual)
            require(
                residual_h == h - singleton_coverage(carrier_i, first),
                f"{body},z1={first}: residual mass cross-check failed",
            )
            require(
                residual_h > 0 and residual_r > 0,
                f"{body},z1={first}: positive residual disappeared",
            )

            apex_ratio = F(36 * residual_r, 7) / residual_h
            apex_cap = apex_ratio.numerator // apex_ratio.denominator
            second_floor = max(
                first + 1,
                ALPHA.numerator * canonical_l // ALPHA.denominator + 1,
            )
            aligned_floor = canonical_l
            drift_first_possible = second_floor <= apex_cap
            aligned_first_possible = aligned_floor <= apex_cap
            row = (
                body,
                first,
                h,
                components,
                canonical_l,
                first_cap,
                delta,
                residual_h,
                residual_r,
                apex_cap,
                second_floor,
                aligned_floor,
                drift_first_possible,
                aligned_first_possible,
            )
            high_excess.append(row)
            if drift_first_possible or aligned_first_possible:
                survivors.append(row)

    require(
        candidate_rows == 626_787,
        f"analytic first-label universe changed: {candidate_rows}",
    )
    require(len(high_excess) == 4_084, "high-excess row count changed")
    require(
        len({row[0] for row in high_excess}) == 2_309,
        "high-excess root count changed",
    )
    require(max(row[1] for row in high_excess) == 66, "high-excess z1 max changed")
    require(len(survivors) == 257, "first-apex survivor count changed")
    require(
        len({row[0] for row in survivors}) == 183,
        "first-apex survivor root count changed",
    )
    require(max(row[1] for row in survivors) == 44, "surviving z1 max changed")

    typed_counts = {
        (drift, aligned): sum(
            row[12] == drift and row[13] == aligned for row in survivors
        )
        for drift in (False, True)
        for aligned in (False, True)
    }
    require(
        typed_counts
        == {
            (False, False): 0,
            (False, True): 0,
            (True, False): 257,
            (True, True): 0,
        },
        "typed first-apex split changed",
    )

    digest = hashlib.sha256(
        b"LRC14/five-aligned/two-drift/residual-first-apex/v1\n"
        + repr(tuple(survivors)).encode()
    ).hexdigest()
    require(
        digest == "ac39be269cbd8695fb1cdefdd5d2bb8ce5512f636fd93626795c1aba949c3362",
        f"survivor digest changed: {digest}",
    )

    minimum_margin = min(
        row[9] - row[10]
        for row in survivors
    )
    maximum_apex = max(row[9] for row in survivors)
    maximum_second_floor = max(row[10] for row in survivors)
    print("LRC14 five-aligned/two-drift residual first-apex audit")
    print(
        "universe=(E six-subset of 1..14,z1>=15,z1 mod L nonzero,"
        "z1 below rootwise inclusive multiplicity/BV cap,"
        "delta_E(z1)>=88h_E/1365)"
    )
    print(f"projected_safe_alpha={ALPHA}")
    print(
        f"candidate_first_rows={candidate_rows};"
        f"high_excess_rows={len(high_excess)};"
        f"high_excess_roots={len({row[0] for row in high_excess})}"
    )
    print(
        f"closed_by_typed_first_apex={len(high_excess)-len(survivors)};"
        f"survivors={len(survivors)};"
        f"survivor_roots={len({row[0] for row in survivors})};"
        f"survivor_max_z1={max(row[1] for row in survivors)}"
    )
    print(f"typed_survivor_counts={typed_counts}")
    print(
        f"maximum_residual_apex_cap={maximum_apex};"
        f"maximum_second_drift_floor={maximum_second_floor};"
        f"minimum_cap_minus_drift_floor={minimum_margin}"
    )
    print(f"survivor_digest={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
