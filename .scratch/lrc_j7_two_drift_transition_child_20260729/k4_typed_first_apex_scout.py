#!/usr/bin/env python3
"""Typed first-apex scout for the four-aligned/three-drift k=4 frontier.

The canonical projected-suffix theorem proves that a hypothetical k=4 row
has first drift ``15 <= z1 <= 182``.  This script deliberately works on the
larger Cartesian superset of *all* such ``(E,z1)`` rows, rather than importing
the suffix survivor list.

After deleting ``z1``, six tails remain: four aligned labels, all at least
``L``, and two ordered drifts ``z2<z3``.  THM-2893 gives a residual first
apex at most

    A1 = floor(36 r1/(7 h1)).

Hence:

* if ``A1 < min(L,z1+1)``, the row is empty;
* if ``z1+1 <= A1 < L``, the witness must be ``z2``, so ``z2<=A1``;
* if ``L<=A1``, an aligned witness is also possible, but its multiplier is
  at most ``floor(A1/L)``.

All residual interval arithmetic is exact.  Counts are exploratory until
frozen below.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import residual_first_apex_audit as A


FIRST_MAX = 182
PROJECTED_ALPHA = F(2366, 21875)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def danger(speed: int) -> tuple[tuple[F, F], ...]:
    return A.danger(speed)


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier_i = A.integer_carrier(body)
    carrier = A.fraction_carrier(carrier_i)
    canonical_l = 14 * math.lcm(*body)
    projected_floor_q = PROJECTED_ALPHA * canonical_l
    projected_floor = max(
        15,
        projected_floor_q.numerator // projected_floor_q.denominator + 1,
    )
    rows: list[tuple[object, ...]] = []
    for first in range(A.BASE_LABEL, FIRST_MAX + 1):
        if first % canonical_l == 0:
            continue
        residual = A.subtract(carrier, danger(first))
        residual_h = A.mass(residual)
        residual_r = len(residual)
        require(residual_h > 0 and residual_r > 0, "z1 erased body carrier")
        ratio = F(36 * residual_r, 7) / residual_h
        apex_cap = ratio.numerator // ratio.denominator
        drift_floor = first + 1
        closed = apex_cap < min(canonical_l, drift_floor)
        drift_forced = drift_floor <= apex_cap < canonical_l
        aligned_possible = canonical_l <= apex_cap
        require(
            sum((closed, drift_forced, aligned_possible)) == 1,
            "typed first-apex trichotomy failed",
        )
        rows.append(
            (
                body,
                first,
                canonical_l,
                projected_floor,
                residual_h,
                residual_r,
                apex_cap,
                closed,
                drift_forced,
                aligned_possible,
                apex_cap // canonical_l,
            )
        )
    return {"body": body, "rows": rows}


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
    rows = sorted(row for item in profiles for row in item["rows"])
    closed = [row for row in rows if row[7]]
    forced = [row for row in rows if row[8]]
    aligned = [row for row in rows if row[9]]
    aligned_hist = tuple(sorted(Counter(row[10] for row in aligned).items()))
    digest = hashlib.sha256(
        b"LRC14/k4/z1<=182/typed-first-apex-superset/v1\n"
        + repr(tuple(rows)).encode()
    ).hexdigest()
    require(len(rows) == 504_503, "k4 z1 superset size changed")
    require(len(closed) == 0, "k4 first-gate closure count changed")
    require(len(forced) == 498_981, "k4 forced-z2 count changed")
    require(len(aligned) == 5_522, "k4 aligned-first count changed")
    require(
        len({row[0] for row in forced}) == 3_002,
        "k4 forced-z2 root count changed",
    )
    require(
        max(row[6] for row in forced) == 1_888,
        "k4 uniform forced z2 cap changed",
    )
    require(
        len({row[0] for row in aligned}) == 91,
        "k4 aligned-first root count changed",
    )
    require(
        aligned_hist
        == ((1, 4416), (2, 739), (3, 278), (4, 26), (5, 30), (6, 28), (7, 5)),
        "k4 aligned multiplier-cap histogram changed",
    )
    require(
        digest == "e46296a196084178db69742851d9f6c98397e6447453a769ad2ae0706ef10bf0",
        f"k4 typed first-apex digest changed: {digest}",
    )

    print("LRC14 k4 typed first-apex superset scout")
    print(
        "universe=(E six-subset of1..14,15<=z1<=182,z1 mod L nonzero);"
        "canonical projected-suffix survivors are a subset"
    )
    print(
        f"rows={len(rows)};closed={len(closed)};"
        f"z2_forced_bounded={len(forced)};"
        f"aligned_first_possible={len(aligned)}"
    )
    print(
        f"forced_roots={len({row[0] for row in forced})};"
        f"forced_max_z2_cap={max((row[6] for row in forced),default=None)};"
        f"aligned_roots={len({row[0] for row in aligned})};"
        f"aligned_multiplier_cap_histogram={aligned_hist}"
    )
    print(
        f"projected_z3_floor_min={min(row[3] for row in rows)};"
        f"projected_z3_floor_max={max(row[3] for row in rows)}"
    )
    print(f"row_digest={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
