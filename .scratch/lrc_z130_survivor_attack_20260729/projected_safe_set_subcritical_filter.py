#!/usr/bin/env python3
"""Projected-safe-set filter for the finite subcritical k=5 bank.

Scratch discovery companion to ``projected_safe_set_pair_filter.py``.
It reconstructs every first label with

    delta_E(z1) < (88/1365) h_E

that survives both the projected largest-drift floor and the rigorous
THM-1094 second-label cap.  It then applies the exact cell-prefix projected
safe-set certificate to every integer in the resulting (possibly relaxed)
second-label interval.  Scanning the whole analytic interval is stronger
as a universe audit than first filtering by the exact second excess.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from fractions import Fraction as F
from itertools import combinations

import projected_safe_set_pair_filter as P


S = P.S
ETA = F(88, 1365)
ALPHA = F(2275, 18627)


def projected_floor(canonical_l: int) -> int:
    wall = ALPHA * canonical_l
    return max(15, wall.numerator // wall.denominator + 1)


def enumerate_rows() -> list[tuple[object, ...]]:
    candidate_rows = 0
    high_excess_rows = 0
    subcritical_rows: list[tuple[object, ...]] = []
    for body in combinations(range(1, 15), 6):
        carrier_i = S.integer_carrier(body)
        h = F(sum(right - left for left, right in carrier_i), S.RULER)
        components = len(carrier_i)
        canonical_l = 14 * math.lcm(*body)
        target = ETA * h
        gamma = F(6 * components, 49)
        first_bound = F(2) * gamma / target
        first_cap = first_bound.numerator // first_bound.denominator
        for first in range(S.BASE_LABEL, first_cap + 1):
            if first % canonical_l == 0:
                continue
            candidate_rows += 1
            delta_first = S.singleton_coverage(carrier_i, first) - h / 7
            if delta_first >= target:
                high_excess_rows += 1
                continue
            gap = target - delta_first
            second_cap_q = gamma / gap
            second_cap = second_cap_q.numerator // second_cap_q.denominator
            second_floor = max(first + 1, projected_floor(canonical_l))
            if second_floor <= second_cap:
                subcritical_rows.append(
                    (
                        body,
                        first,
                        h,
                        components,
                        canonical_l,
                        first_cap,
                        delta_first,
                        gap,
                        second_floor,
                        second_cap,
                        carrier_i,
                    )
                )
    P.require(candidate_rows == 626_787, "candidate first-label universe changed")
    P.require(high_excess_rows == 4_084, "high-excess universe changed")
    P.require(
        len(subcritical_rows) == 2_290,
        f"subcritical row universe changed: {len(subcritical_rows)}",
    )
    return subcritical_rows


def interval_pair_count(rows: list[tuple[object, ...]]) -> int:
    count = 0
    for row in rows:
        canonical_l = row[4]
        second_floor = row[8]
        second_cap = row[9]
        count += sum(
            second % canonical_l != 0
            for second in range(second_floor, second_cap + 1)
        )
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-only", action="store_true")
    args = parser.parse_args()

    rows = enumerate_rows()
    pair_count = interval_pair_count(rows)
    print("LRC14 k5 projected subcritical finite-pair filter")
    print(
        f"subcritical_first_rows={len(rows)};"
        f"subcritical_roots={len({row[0] for row in rows})};"
        f"analytic_interval_pairs={pair_count};"
        f"maximum_z1={max(row[1] for row in rows)};"
        f"maximum_z2_cap={max(row[9] for row in rows)}"
    )
    if args.inventory_only:
        print("inventory_only=TRUE")
        return

    killed = 0
    surviving: list[tuple[object, ...]] = []
    maximum_prefix = 0
    minimum_certificate: tuple[F, tuple[object, ...]] | None = None
    cell_cache: dict[tuple[int, ...], tuple[int, ...]] = {}
    for row in rows:
        (
            body,
            first,
            h,
            components,
            canonical_l,
            first_cap,
            delta_first,
            gap,
            second_floor,
            second_cap,
            carrier_i,
        ) = row
        cells = cell_cache.setdefault(
            body,
            P.body_cells(carrier_i, canonical_l),
        )
        for second in range(second_floor, second_cap + 1):
            if second % canonical_l == 0:
                continue
            projected_lower, cells_used = P.projected_safe_lower_bound(
                cells,
                canonical_l,
                first,
                second,
            )
            maximum_prefix = max(maximum_prefix, cells_used)
            pair = (
                body,
                first,
                second,
                canonical_l,
                h,
                components,
                gap,
                second_floor,
                second_cap,
                projected_lower,
                cells_used,
            )
            if projected_lower >= P.ALIGNED_UNION_CAP:
                killed += 1
                margin = projected_lower - P.ALIGNED_UNION_CAP
                if (
                    minimum_certificate is None
                    or margin < minimum_certificate[0]
                    or (
                        margin == minimum_certificate[0]
                        and pair < minimum_certificate[1]
                    )
                ):
                    minimum_certificate = (margin, pair)
            else:
                surviving.append(pair)

    digest = hashlib.sha256(
        b"LRC14/k5/projected-safe-set-subcritical/v1\n"
        + repr(tuple(surviving)).encode()
    ).hexdigest()
    print(
        f"killed_by_projected_prefix={killed};"
        f"surviving_pairs={len(surviving)};"
        f"surviving_roots={len({row[0] for row in surviving})}"
    )
    print(f"maximum_cell_prefix={maximum_prefix}")
    print(f"minimum_certificate={minimum_certificate}")
    print(f"survivor_digest={digest}")
    if surviving:
        print(f"survivors={tuple(surviving)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
