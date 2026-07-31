#!/usr/bin/env python3
"""Recursive typed-apex cap for the k=4/three-drift frontier.

The canonical projected suffix theorem gives ``z1<=182``.  Work on the full
Cartesian superset of such rows.  After deleting ``z1``, repeatedly apply
THM-2893 to the literal residual.  At a state with ``s`` already selected
least aligned multipliers, the remaining tail count is ``p=6-s`` and

    A_s=floor(6 p r_s / (7 (7-p) h_s)).

Partition disjointly:

* ``z2<=A_s`` is a bounded drift branch;
* on the complement ``z2>A_s``, the low apex must be the least remaining
  aligned label.  Enumerate its multiplier in increasing order through
  ``floor(A_s/L)`` and recurse.

There are only four aligned labels, so every branch either closes or forces
an explicit ``z2`` cap after at most four aligned deletions.  Increasing
multiplier order is canonical and avoids permutation duplicates.
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def danger(speed: int) -> tuple[tuple[F, F], ...]:
    return A.danger(speed)


def apex_cap(
    residual: tuple[tuple[F, F], ...],
    remaining: int,
) -> int:
    h = A.mass(residual)
    components = len(residual)
    require(h > 0 and components > 0, "typed recursion reached empty residual")
    ratio = F(
        6 * remaining * components,
        7 * (7 - remaining),
    ) / h
    return ratio.numerator // ratio.denominator


class Ledger:
    def __init__(self) -> None:
        self.states: Counter[int] = Counter()
        self.aligned_children: Counter[int] = Counter()
        self.drift_branches: Counter[int] = Counter()
        self.closed: Counter[int] = Counter()
        self.maximum_z2_cap = 0
        self.maximum_multiplier = 0
        self.maximum_apex_by_p: dict[int, int] = {}
        self.rows: list[tuple[object, ...]] = []


def descend(
    ledger: Ledger,
    body: tuple[int, ...],
    first: int,
    canonical_l: int,
    residual: tuple[tuple[F, F], ...],
    remaining_aligned: int,
    previous_multiplier: int,
    drift_floor: int,
    prefix: tuple[int, ...],
) -> None:
    remaining = remaining_aligned + 2
    require(2 <= remaining <= 6, "wrong k4 remaining-tail count")
    cap = apex_cap(residual, remaining)
    ledger.states[remaining] += 1
    ledger.maximum_apex_by_p[remaining] = max(
        ledger.maximum_apex_by_p.get(remaining, 0),
        cap,
    )

    if drift_floor <= cap:
        ledger.drift_branches[remaining] += 1
        ledger.maximum_z2_cap = max(ledger.maximum_z2_cap, cap)
        ledger.rows.append(
            (
                body,
                first,
                canonical_l,
                remaining,
                prefix,
                drift_floor,
                cap,
                "Z2-BOUNDED",
            )
        )

    next_drift_floor = max(drift_floor, cap + 1)
    if remaining_aligned == 0:
        # The complementary branch z2>A_s has no aligned label left to pay
        # the first-apex gate, whether or not the bounded z2 branch exists.
        ledger.closed[remaining] += 1
        return

    first_multiplier = previous_multiplier + 1
    multiplier_cap = cap // canonical_l
    if first_multiplier > multiplier_cap:
        # The complement z2>A_s has neither a drift nor an aligned low apex.
        ledger.closed[remaining] += 1
        return

    for multiplier in range(first_multiplier, multiplier_cap + 1):
        ledger.aligned_children[remaining] += 1
        ledger.maximum_multiplier = max(ledger.maximum_multiplier, multiplier)
        child = A.subtract(
            residual,
            danger(canonical_l * multiplier),
        )
        require(
            A.mass(child) > 0 and len(child) > 0,
            "aligned prefix erased residual",
        )
        ledger.rows.append(
            (
                body,
                first,
                canonical_l,
                remaining,
                prefix + (multiplier,),
                next_drift_floor,
                cap,
                "ALIGNED-COMPLEMENT",
            )
        )
        descend(
            ledger,
            body,
            first,
            canonical_l,
            child,
            remaining_aligned - 1,
            multiplier,
            next_drift_floor,
            prefix + (multiplier,),
        )


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier_i = A.integer_carrier(body)
    carrier = A.fraction_carrier(carrier_i)
    canonical_l = 14 * math.lcm(*body)
    ledger = Ledger()
    for first in range(A.BASE_LABEL, FIRST_MAX + 1):
        if first % canonical_l == 0:
            continue
        residual = A.subtract(carrier, danger(first))
        descend(
            ledger,
            body,
            first,
            canonical_l,
            residual,
            remaining_aligned=4,
            previous_multiplier=0,
            drift_floor=first + 1,
            prefix=(),
        )
    return {"body": body, "ledger": ledger}


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

    states: Counter[int] = Counter()
    aligned_children: Counter[int] = Counter()
    drift_branches: Counter[int] = Counter()
    closed: Counter[int] = Counter()
    maximum_z2_cap = maximum_multiplier = 0
    maximum_apex_by_p: dict[int, int] = {}
    rows: list[tuple[object, ...]] = []
    for profile_row in profiles:
        ledger = profile_row["ledger"]
        states.update(ledger.states)
        aligned_children.update(ledger.aligned_children)
        drift_branches.update(ledger.drift_branches)
        closed.update(ledger.closed)
        maximum_z2_cap = max(maximum_z2_cap, ledger.maximum_z2_cap)
        maximum_multiplier = max(maximum_multiplier, ledger.maximum_multiplier)
        for remaining, cap in ledger.maximum_apex_by_p.items():
            maximum_apex_by_p[remaining] = max(
                maximum_apex_by_p.get(remaining, 0),
                cap,
            )
        rows.extend(ledger.rows)
    digest = hashlib.sha256(
        b"LRC14/k4/z1<=182/recursive-typed-apex/v1\n"
        + repr(tuple(rows)).encode()
    ).hexdigest()

    print("LRC14 k4 recursive typed-apex z2 cap scout")
    print(
        "universe=(all six-body roots,15<=z1<=182 nonaligned);"
        "canonical projected-suffix bank is a subset"
    )
    print("coefficients=(p6:36/7,p5:15/7,p4:8/7,p3:9/14,p2:12/35)")
    print(f"states_by_remaining={tuple(sorted(states.items(),reverse=True))}")
    print(
        "aligned_children_by_remaining="
        f"{tuple(sorted(aligned_children.items(),reverse=True))}"
    )
    print(
        "drift_branches_by_remaining="
        f"{tuple(sorted(drift_branches.items(),reverse=True))}"
    )
    print(f"closed_by_remaining={tuple(sorted(closed.items(),reverse=True))}")
    print(
        f"uniform_z2_cap={maximum_z2_cap};"
        f"maximum_aligned_multiplier={maximum_multiplier};"
        f"maximum_apex_by_remaining={tuple(sorted(maximum_apex_by_p.items(),reverse=True))}"
    )
    print(f"typed_recursion_digest={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
