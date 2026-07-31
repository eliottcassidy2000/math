#!/usr/bin/env python3
"""Exact projected-safe-set filter on the 257 first-apex k=5 rows.

Scratch only.  The row universe is reconstructed independently from
``residual_first_apex_audit.py`` using its exact integer carrier engine.
For every finite second-drift label left by that audit, this script uses

    P = phi_L(C_E \ (D_z1 union D_z2)),   phi_L(t)=Lt mod 1.

If five aligned combs complete the cover, then P is contained in their open
normalized union, whose measure is at most 1-u_5=887/1365.  Compact
containment is strict, so every row with mu(P)>=887/1365 is impossible.

The computation uses the equivalent exact fiber identity

    T \ P = intersection_(body-safe cells j)
            (E_z1(j) union E_z2(j)),

where E_z(j)={u in [0,1]: ||z(j+u)/L||<1/14}.  All endpoints and all
comparisons are rational.  The intersection stops as soon as its mass is at
most `u_5`.  Consequently `1` minus the current prefix intersection mass is
an exact rational **lower bound** for `mu(P)`, not necessarily its final
exact value.  That lower bound already certifies `mu(P)>=1-u_5`.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = (
    ROOT
    / ".scratch/lrc_j7_two_drift_transition_child_20260729"
    / "residual_first_apex_audit.py"
)
EXPECTED_SOURCE_SHA256 = (
    "8dfc811e44d6c32081b0f44e80377f8a8970ff2a17d10dffbd98b97bc5755825"
)
U5 = F(478, 1365)
ALIGNED_UNION_CAP = 1 - U5


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


spec = importlib.util.spec_from_file_location("first_apex_source", SOURCE)
require(spec is not None and spec.loader is not None, "cannot load first-apex source")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def merge(
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


def intersect(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def interval_mass(rows: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in rows), F(0))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def phase_danger(cell: int, speed: int, canonical_l: int) -> tuple[tuple[F, F], ...]:
    """Danger subset in normalized coordinate u on physical cell ``cell``."""

    start = F(speed * cell, canonical_l)
    finish = F(speed * (cell + 1), canonical_l)
    rows: list[tuple[F, F]] = []
    for tooth in range(floor_fraction(start) - 1, ceil_fraction(finish) + 2):
        left = max(
            F(0),
            F(canonical_l, speed) * (F(tooth) - F(1, 14)) - cell,
        )
        right = min(
            F(1),
            F(canonical_l, speed) * (F(tooth) + F(1, 14)) - cell,
        )
        if left < right:
            rows.append((left, right))
    return merge(rows)


def body_cells(
    carrier_i: tuple[tuple[int, int], ...],
    canonical_l: int,
) -> tuple[int, ...]:
    require(S.RULER % canonical_l == 0, "body grid does not divide ruler")
    scale = S.RULER // canonical_l
    cells: list[int] = []
    for left, right in carrier_i:
        require(left % scale == 0 and right % scale == 0, "non-grid carrier endpoint")
        cells.extend(range(left // scale, right // scale))
    require(cells, "body has no safe cells")
    return tuple(cells)


def projected_safe_lower_bound(
    cells: tuple[int, ...],
    canonical_l: int,
    first: int,
    second: int,
) -> tuple[F, int]:
    """Return a certified lower bound for mu(P), with its prefix length."""

    common_danger: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = merge(
            [
                *phase_danger(cell, first, canonical_l),
                *phase_danger(cell, second, canonical_l),
            ]
        )
        common_danger = intersect(common_danger, local_union)
        cells_used += 1
        if interval_mass(common_danger) <= U5:
            # The final common-danger intersection is contained in this
            # prefix intersection, so final mu(P) is at least the returned
            # value.  This is enough for the strict compact-in-open
            # contradiction.
            break
    return 1 - interval_mass(common_danger), cells_used


def enumerate_first_apex_survivors() -> list[tuple[object, ...]]:
    candidate_rows = 0
    high_excess = 0
    survivors: list[tuple[object, ...]] = []
    for body in combinations(range(1, 15), 6):
        carrier_i = S.integer_carrier(body)
        carrier_f = S.fraction_carrier(carrier_i)
        h = F(sum(right - left for left, right in carrier_i), S.RULER)
        components = len(carrier_i)
        canonical_l = 14 * math.lcm(*body)
        first_bound = F(12 * components, 49) / (S.ETA * h)
        # The THM-1094 component discrepancy is non-strict:
        # delta(w)<=6r/(49w).  Therefore an integral analytic boundary label
        # must be retained.  Five such equality rows occur in the universe;
        # all five fail the later exact high-excess test.
        first_cap = first_bound.numerator // first_bound.denominator
        for first in range(S.BASE_LABEL, first_cap + 1):
            if first % canonical_l == 0:
                continue
            candidate_rows += 1
            delta = S.singleton_coverage(carrier_i, first) - h / 7
            if delta < S.ETA * h:
                continue
            high_excess += 1
            residual = S.subtract(carrier_f, S.danger(first))
            residual_h = S.mass(residual)
            residual_r = len(residual)
            apex_ratio = F(36 * residual_r, 7) / residual_h
            apex_cap = apex_ratio.numerator // apex_ratio.denominator
            second_floor = max(
                first + 1,
                S.ALPHA.numerator * canonical_l // S.ALPHA.denominator + 1,
            )
            aligned_floor = canonical_l
            drift_possible = second_floor <= apex_cap
            aligned_possible = aligned_floor <= apex_cap
            if drift_possible or aligned_possible:
                survivors.append(
                    (
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
                        drift_possible,
                        aligned_possible,
                        carrier_i,
                    )
                )
    require(candidate_rows == 626_787, "first-apex candidate universe changed")
    require(high_excess == 4_084, "high-excess universe changed")
    require(len(survivors) == 257, "first-apex survivor universe changed")
    return survivors


def main() -> None:
    require(
        file_sha256(SOURCE) == EXPECTED_SOURCE_SHA256,
        f"first-apex source changed: {file_sha256(SOURCE)}",
    )
    first_rows = enumerate_first_apex_survivors()
    pair_rows = 0
    killed = 0
    surviving: list[tuple[object, ...]] = []
    maximum_cells_used = 0
    closest_kill: tuple[F, tuple[object, ...]] | None = None

    cell_cache: dict[tuple[int, ...], tuple[int, ...]] = {}
    for row in first_rows:
        (
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
            drift_possible,
            aligned_possible,
            carrier_i,
        ) = row
        cells = cell_cache.setdefault(
            body,
            body_cells(carrier_i, canonical_l),
        )
        for second in range(second_floor, apex_cap + 1):
            if second <= first or second % canonical_l == 0:
                continue
            pair_rows += 1
            projected_lower, cells_used = projected_safe_lower_bound(
                cells,
                canonical_l,
                first,
                second,
            )
            maximum_cells_used = max(maximum_cells_used, cells_used)
            pair = (
                body,
                first,
                second,
                canonical_l,
                h,
                components,
                residual_h,
                residual_r,
                apex_cap,
                second_floor,
                projected_lower,
                cells_used,
            )
            if projected_lower >= ALIGNED_UNION_CAP:
                killed += 1
                margin = projected_lower - ALIGNED_UNION_CAP
                if (
                    closest_kill is None
                    or margin < closest_kill[0]
                    or (margin == closest_kill[0] and pair < closest_kill[1])
                ):
                    closest_kill = (margin, pair)
            else:
                surviving.append(pair)

    digest = hashlib.sha256(
        b"LRC14/k5/projected-safe-set-pair-filter/v1\n"
        + repr(tuple(surviving)).encode()
    ).hexdigest()
    require(pair_rows == 42_912, "finite second-drift pair universe changed")
    require(killed == pair_rows, "projected-prefix filter no longer closes the bank")
    require(not surviving, "projected-prefix survivor bank is nonempty")
    require(maximum_cells_used == 22, "maximum prefix length changed")
    require(
        closest_kill
        == (
            F(1, 378_105),
            (
                (1, 2, 3, 5, 9, 10),
                24,
                277,
                1_260,
                F(136, 315),
                14,
                F(43, 126),
                26,
                391,
                154,
                F(180, 277),
                1,
            ),
        ),
        "minimum certificate margin row changed",
    )
    require(
        digest
        == "3f001616fc494fca31891589993823d9aab957b9a5b73b96159cc899361d27bd",
        "empty survivor digest changed",
    )
    print("LRC14 k5 exact projected-safe-set pair filter")
    print(f"first_apex_rows={len(first_rows)};candidate_pairs={pair_rows}")
    print(
        f"killed_by_projected_measure={killed};"
        f"surviving_pairs={len(surviving)};"
        f"surviving_roots={len({row[0] for row in surviving})}"
    )
    print(f"aligned_union_cap={ALIGNED_UNION_CAP};u5={U5}")
    print(f"maximum_cell_prefix_used_before_decision={maximum_cells_used}")
    print(f"minimum_certificate_margin={closest_kill}")
    print(f"survivor_digest={digest}")
    if surviving:
        print(f"survivors={tuple(surviving)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
