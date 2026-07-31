#!/usr/bin/env python3
"""Exact pair bank for the projected-suffix subcritical k=5 rows.

The canonical projected suffix scan leaves 4,702 first-drift rows.  Of
these, 4,084 have ``delta(z1) >= eta*h`` and are handled by the companion
high-excess closure.  This script reconstructs the remaining 618 rows and
then replaces the suffix envelope by an exact scan of every finite ``z2``
allowed by

    z2 >= max(z1+1, floor((2275/18627)L)+1),
    z2 <= floor((6r/49)/(eta*h-delta(z1))).

The first inequality is the projected-safe largest-drift wall and the
second is the component discrepancy bound.  A pair is retained only when
its exact singleton excess pays the remaining safe-surplus gap.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import multiprocessing as mp
from fractions import Fraction as F
from itertools import combinations

import residual_first_apex_audit as A


HORIZON = 7_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def excess_numerator(
    carrier: tuple[tuple[int, int], ...],
    carrier_numerator: int,
    speed: int,
) -> int:
    """Numerator N with delta(speed)=N/(7*RULER*speed)."""
    covered = sum(
        A.integer_primitive(speed * right) - A.integer_primitive(speed * left)
        for left, right in carrier
    )
    return 7 * covered - carrier_numerator * speed


def excess_meets(
    numerator: int,
    speed: int,
    gap: F,
) -> bool:
    return (
        numerator * gap.denominator
        >= gap.numerator * 7 * A.RULER * speed
    )


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier = A.integer_carrier(body)
    carrier_numerator = sum(right - left for left, right in carrier)
    h = F(carrier_numerator, A.RULER)
    components = len(carrier)
    canonical_l = 14 * math.lcm(*body)
    gamma = F(6 * components, 49)
    first_bound = F(12 * components, 49) / (A.ETA * h)
    first_cap = A.inclusive_integer_cap(first_bound)
    wall_q = A.ALPHA * canonical_l
    wall = max(15, wall_q.numerator // wall_q.denominator + 1)

    analytic_rows = 0
    suffix_rows: list[tuple[object, ...]] = []
    exact_pairs: list[tuple[object, ...]] = []
    pair_candidates = 0

    for first in range(A.BASE_LABEL, first_cap + 1):
        if first % canonical_l == 0:
            continue
        first_numerator = excess_numerator(carrier, carrier_numerator, first)
        first_delta = F(first_numerator, 7 * A.RULER * first)
        lower = A.ETA * h
        if first_delta >= lower:
            continue
        gap = lower - first_delta
        second_floor = max(first + 1, wall)
        second_bound = gamma / gap
        second_cap = A.inclusive_integer_cap(second_bound)
        if second_floor > second_cap:
            continue
        analytic_rows += 1

        exact_suffix_hit = False
        row_pairs: list[tuple[object, ...]] = []
        for second in range(second_floor, second_cap + 1):
            if second % canonical_l == 0:
                continue
            pair_candidates += 1
            second_numerator = excess_numerator(
                carrier,
                carrier_numerator,
                second,
            )
            hit = excess_meets(second_numerator, second, gap)
            if second <= HORIZON and hit:
                exact_suffix_hit = True
            if hit:
                row_pairs.append(
                    (
                        body,
                        first,
                        second,
                        canonical_l,
                        h,
                        components,
                        first_delta,
                        gap,
                        second_floor,
                        second_cap,
                        F(second_numerator, 7 * A.RULER * second),
                    )
                )

        high_tail_start = max(HORIZON + 1, second_floor)
        high_tail = F(6 * components, 49 * high_tail_start)
        tail_suffix_hit = high_tail >= gap
        if exact_suffix_hit or tail_suffix_hit:
            suffix_rows.append(
                (
                    body,
                    first,
                    canonical_l,
                    h,
                    components,
                    first_delta,
                    gap,
                    second_floor,
                    second_cap,
                    exact_suffix_hit,
                    tail_suffix_hit,
                    high_tail_start,
                    high_tail,
                    len(row_pairs),
                )
            )
            exact_pairs.extend(row_pairs)
        else:
            require(
                not row_pairs,
                f"{body},z1={first}: exact pair escaped suffix predicate",
            )

    return {
        "body": body,
        "analytic_rows": analytic_rows,
        "suffix_rows": suffix_rows,
        "exact_pairs": exact_pairs,
        "pair_candidates": pair_candidates,
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
        profiles = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, roots, chunksize=1))
    profiles.sort(key=lambda row: row["body"])

    analytic_rows = sum(row["analytic_rows"] for row in profiles)
    pair_candidates = sum(row["pair_candidates"] for row in profiles)
    suffix_rows = sorted(
        row
        for profile_row in profiles
        for row in profile_row["suffix_rows"]
    )
    exact_pairs = sorted(
        row
        for profile_row in profiles
        for row in profile_row["exact_pairs"]
    )
    digest = hashlib.sha256(
        b"LRC14/k5/subcritical-exact-pair-bank/v1\n"
        + repr(tuple(suffix_rows)).encode()
        + b"\n"
        + repr(tuple(exact_pairs)).encode()
    ).hexdigest()
    require(analytic_rows == 2_290, "analytic subcritical row count changed")
    require(len(suffix_rows) == 618, "projected suffix row count changed")
    require(
        len({row[0] for row in suffix_rows}) == 435,
        "projected suffix root count changed",
    )
    require(pair_candidates == 7_218_110, "finite z2 universe changed")
    require(len(exact_pairs) == 194_073, "exact admissible pair count changed")
    require(
        len({row[0] for row in exact_pairs}) == 407,
        "exact admissible pair root count changed",
    )
    require(
        len({(row[0], row[1]) for row in exact_pairs}) == 590,
        "exact admissible (E,z1) count changed",
    )
    require(max(row[1] for row in suffix_rows) == 55, "subcritical z1 max changed")
    require(max(row[8] for row in suffix_rows) == 688_956, "z2 cap max changed")
    require(
        max(row[2] for row in exact_pairs) == 327_502,
        "admissible z2 max changed",
    )
    require(
        digest == "4827dabbc24f07bdf9e628e6b706fe66e7ed451c893cece993499a4cf13394c0",
        f"subcritical bank digest changed: {digest}",
    )
    suffix_mode_counts = {
        (exact_hit, tail_hit): sum(
            row[9] == exact_hit and row[10] == tail_hit
            for row in suffix_rows
        )
        for exact_hit in (False, True)
        for tail_hit in (False, True)
    }

    print("LRC14 k5 projected-suffix subcritical exact pair bank")
    print(
        "universe=(six-body roots,subcritical z1 under inclusive BV cap,"
        "projected-safe z2 floor,component-BV z2 cap)"
    )
    print(
        f"analytic_first_rows={analytic_rows};"
        f"projected_suffix_rows={len(suffix_rows)};"
        f"suffix_roots={len({row[0] for row in suffix_rows})}"
    )
    print(
        f"finite_z2_candidates={pair_candidates};"
        f"exact_excess_admissible_pairs={len(exact_pairs)};"
        f"pair_roots={len({row[0] for row in exact_pairs})};"
        f"pair_(E,z1)_rows={len({(row[0],row[1]) for row in exact_pairs})}"
    )
    print(f"suffix_mode_counts={suffix_mode_counts}")
    print(
        f"maximum_z1={max(row[1] for row in suffix_rows)};"
        f"maximum_z2_cap={max(row[8] for row in suffix_rows)};"
        f"maximum_admissible_z2={max((row[2] for row in exact_pairs),default=None)}"
    )
    print(f"bank_digest={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
