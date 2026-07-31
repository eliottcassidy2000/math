#!/usr/bin/env python3
"""Projected-safe-set filter for every exact subcritical k=5 drift pair.

For selected body cells ``j`` and normalized coordinate ``u`` define

    E_z(j) = {u : ||z(j+u)/L|| < 1/14}.

After the two drifts, the image under ``t -> Lt mod 1`` of their literal
residual is

    P = union_j (T \\ (E_z1(j) union E_z2(j))),

so

    T \\ P = intersection_j (E_z1(j) union E_z2(j)).

Five aligned tails can complete the cover only if the compact set ``P`` is
contained in their proper open normalized union.  The five-comb safe floor
therefore forces the strict inequality ``mu(P) < 887/1365``.  This script
tests that condition exactly on every exact excess-admissible pair from
``subcritical_exact_pair_bank.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import multiprocessing as mp
from fractions import Fraction as F
from itertools import combinations

import residual_first_apex_audit as A
import subcritical_exact_pair_bank as B


U5 = F(478, 1365)
ALIGNED_UNION_CAP = 1 - U5


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def mass(rows: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in rows), F(0))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def phase_danger(
    cell: int,
    speed: int,
    canonical_l: int,
) -> tuple[tuple[F, F], ...]:
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
    carrier: tuple[tuple[int, int], ...],
    canonical_l: int,
) -> tuple[int, ...]:
    require(A.RULER % canonical_l == 0, "body grid does not divide ruler")
    scale = A.RULER // canonical_l
    cells: list[int] = []
    for left, right in carrier:
        require(left % scale == 0 and right % scale == 0, "off-grid body endpoint")
        cells.extend(range(left // scale, right // scale))
    require(cells, "body carrier has no positive cells")
    return tuple(cells)


def filter_body(body: tuple[int, ...]) -> dict[str, object]:
    bank = B.profile(body)
    exact_pairs = bank["exact_pairs"]
    if not exact_pairs:
        return {
            "body": body,
            "pairs": 0,
            "killed": 0,
            "survivors": [],
            "maximum_cells": 0,
            "closest": None,
        }

    carrier = A.integer_carrier(body)
    canonical_l = exact_pairs[0][3]
    cells = body_cells(carrier, canonical_l)
    first_cache: dict[tuple[int, int], tuple[tuple[F, F], ...]] = {}
    killed = 0
    maximum_cells = 0
    closest: tuple[F, tuple[object, ...]] | None = None
    survivors: list[tuple[object, ...]] = []

    for pair in exact_pairs:
        first = pair[1]
        second = pair[2]
        common: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
        cells_used = 0
        for cell in cells:
            first_phase = first_cache.get((first, cell))
            if first_phase is None:
                first_phase = phase_danger(cell, first, canonical_l)
                first_cache[(first, cell)] = first_phase
            local_union = merge(
                [*first_phase, *phase_danger(cell, second, canonical_l)]
            )
            common = intersect(common, local_union)
            cells_used += 1
            if mass(common) <= U5:
                break
        projected_mass = 1 - mass(common)
        maximum_cells = max(maximum_cells, cells_used)
        row = (
            body,
            first,
            second,
            canonical_l,
            projected_mass,
            cells_used,
            len(cells),
        )
        if projected_mass >= ALIGNED_UNION_CAP:
            killed += 1
            margin = projected_mass - ALIGNED_UNION_CAP
            if (
                closest is None
                or margin < closest[0]
                or (margin == closest[0] and row < closest[1])
            ):
                closest = (margin, row)
        else:
            # No early stop occurred.  `common` is the exact full-cell
            # intersection, so its complement is the exact projected
            # residual, modulo measure-zero endpoint conventions.
            require(cells_used == len(cells), "uncertified early survivor")
            survivors.append((*row, common))

    return {
        "body": body,
        "pairs": len(exact_pairs),
        "killed": killed,
        "survivors": survivors,
        "maximum_cells": maximum_cells,
        "closest": closest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, mp.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        rows = [filter_body(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(filter_body, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])

    pair_count = sum(row["pairs"] for row in rows)
    killed = sum(row["killed"] for row in rows)
    survivors = sorted(
        survivor
        for row in rows
        for survivor in row["survivors"]
    )
    closest_rows = [row["closest"] for row in rows if row["closest"] is not None]
    closest = min(closest_rows) if closest_rows else None
    maximum_cells = max(row["maximum_cells"] for row in rows)
    digest = hashlib.sha256(
        b"LRC14/k5/subcritical-projected-pair-filter/v1\n"
        + repr(tuple(survivors)).encode()
    ).hexdigest()
    require(pair_count == 194_073, "subcritical exact-pair universe changed")
    require(killed == 194_073, "projected kill count changed")
    require(not survivors, "subcritical projected survivor appeared")
    require(maximum_cells == 871, "maximum projected prefix changed")
    require(
        closest
        == (
            F(1, 378_105),
            (
                (2, 3, 4, 5, 6, 9),
                22,
                554,
                2_520,
                F(180, 277),
                1,
                1_092,
            ),
        ),
        "closest projected kill changed",
    )
    require(
        digest == "18fedd7370fa6ac4a65cf748742345302de244d0c49b724f2317039a8e4f047c",
        f"projected survivor digest changed: {digest}",
    )

    print("LRC14 k5 subcritical exact projected-safe-set pair filter")
    print(
        f"exact_excess_pairs={pair_count};"
        f"killed_by_projected_measure={killed};"
        f"surviving_pairs={len(survivors)};"
        f"surviving_roots={len({row[0] for row in survivors})}"
    )
    print(
        f"aligned_union_cap={ALIGNED_UNION_CAP};"
        f"maximum_cells_used={maximum_cells};"
        f"closest_kill={closest}"
    )
    print(f"survivor_digest={digest}")
    if survivors:
        print(f"survivors={tuple(survivors)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
