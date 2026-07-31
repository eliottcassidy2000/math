#!/usr/bin/env python3
"""Recursive exact closure probe for the 257 high-excess two-drift rows.

This composes three proved necessary conditions:

1. the safe-surplus/BV cap and exact singleton excess select 4,084 possible
   escaping first-drift rows;
2. the projected-safe largest-drift wall and the six-tail first-apex bound
   leave 257 rows, on each of which the low apex is necessarily ``z2``;
3. after enumerating that now-bounded ``z2``, all remaining tails are
   aligned labels ``L*a``.  Reapply the cap-free THM-2893 first-apex gate at
   each literal residual, enumerate the least remaining multiplier, and
   recurse through the five aligned labels.

Every residual and every cutoff is rational and exact.  An empty terminal
would be a genuine a.e. cover candidate needing an endpoint audit; a
positive-mass terminal is rigorously not a pointwise cover.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import residual_first_apex_audit as A


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def danger(speed: int) -> tuple[tuple[F, F], ...]:
    return A.danger(speed)


def initial_survivors() -> list[tuple[object, ...]]:
    rows: list[tuple[object, ...]] = []
    for body in combinations(range(1, 15), 6):
        carrier_i = A.integer_carrier(body)
        carrier_f = A.fraction_carrier(carrier_i)
        h = F(sum(right - left for left, right in carrier_i), A.RULER)
        components = len(carrier_i)
        canonical_l = 14 * math.lcm(*body)
        first_bound = F(12 * components, 49) / (A.ETA * h)
        first_cap = A.inclusive_integer_cap(first_bound)
        for first in range(A.BASE_LABEL, first_cap + 1):
            if first % canonical_l == 0:
                continue
            delta = A.singleton_coverage(carrier_i, first) - h / 7
            if delta < A.ETA * h:
                continue
            residual = A.subtract(carrier_f, danger(first))
            residual_h = A.mass(residual)
            residual_r = len(residual)
            apex_ratio = F(36 * residual_r, 7) / residual_h
            apex_cap = apex_ratio.numerator // apex_ratio.denominator
            second_floor = max(
                first + 1,
                A.ALPHA.numerator * canonical_l // A.ALPHA.denominator + 1,
            )
            if second_floor <= apex_cap < canonical_l:
                rows.append(
                    (
                        body,
                        first,
                        canonical_l,
                        residual,
                        second_floor,
                        apex_cap,
                    )
                )
    return rows


class Ledger:
    def __init__(self) -> None:
        self.second_candidates = 0
        self.second_closed_at_aligned_gate = 0
        self.second_survivors = 0
        self.second_survivor_roots: set[tuple[int, ...]] = set()
        self.second_survivor_first_rows: set[tuple[tuple[int, ...], int]] = set()
        self.nodes_by_remaining: Counter[int] = Counter()
        self.children_by_remaining: Counter[int] = Counter()
        self.closed_without_child_by_remaining: Counter[int] = Counter()
        self.terminal_positive = 0
        self.terminal_zero_mass: list[tuple[object, ...]] = []
        self.maximum_multiplier = 0
        self.minimum_terminal_mass: tuple[F, tuple[object, ...]] | None = None
        self.digest_rows: list[tuple[object, ...]] = []


def first_apex_multiplier_cap(
    residual: tuple[tuple[F, F], ...],
    canonical_l: int,
    remaining: int,
) -> tuple[int, int]:
    """Return physical-label and multiplier caps for the next aligned apex."""
    h = A.mass(residual)
    components = len(residual)
    require(h > 0 and components > 0, "first-apex gate received empty residual")
    ratio = F(6 * remaining * components, 7 * (7 - remaining)) / h
    label_cap = ratio.numerator // ratio.denominator
    return label_cap, label_cap // canonical_l


def descend_aligned(
    ledger: Ledger,
    body: tuple[int, ...],
    first: int,
    second: int,
    canonical_l: int,
    residual: tuple[tuple[F, F], ...],
    remaining: int,
    previous_multiplier: int,
    prefix: tuple[int, ...],
) -> None:
    require(1 <= remaining <= 5, "invalid aligned recursion depth")
    ledger.nodes_by_remaining[remaining] += 1
    label_cap, multiplier_cap = first_apex_multiplier_cap(
        residual,
        canonical_l,
        remaining,
    )
    first_candidate = previous_multiplier + 1
    if first_candidate > multiplier_cap:
        ledger.closed_without_child_by_remaining[remaining] += 1
        return

    for multiplier in range(first_candidate, multiplier_cap + 1):
        ledger.children_by_remaining[remaining] += 1
        ledger.maximum_multiplier = max(ledger.maximum_multiplier, multiplier)
        child = A.subtract(residual, danger(canonical_l * multiplier))
        child_mass = A.mass(child)
        row = (
            body,
            first,
            second,
            canonical_l,
            remaining,
            prefix + (multiplier,),
            A.mass(residual),
            len(residual),
            label_cap,
            multiplier_cap,
            child_mass,
            len(child),
        )
        ledger.digest_rows.append(row)
        if remaining == 1:
            if child_mass == 0:
                ledger.terminal_zero_mass.append(row)
            else:
                ledger.terminal_positive += 1
                if (
                    ledger.minimum_terminal_mass is None
                    or child_mass < ledger.minimum_terminal_mass[0]
                    or (
                        child_mass == ledger.minimum_terminal_mass[0]
                        and row < ledger.minimum_terminal_mass[1]
                    )
                ):
                    ledger.minimum_terminal_mass = (child_mass, row)
            continue
        require(
            child_mass > 0 and len(child) > 0,
            "nonterminal aligned prefix unexpectedly covered the residual",
        )
        descend_aligned(
            ledger,
            body,
            first,
            second,
            canonical_l,
            child,
            remaining - 1,
            multiplier,
            prefix + (multiplier,),
        )


def main() -> None:
    initial = initial_survivors()
    require(len(initial) == 257, "initial typed survivor bank changed")
    ledger = Ledger()

    for body, first, canonical_l, residual_one, second_floor, apex_cap in initial:
        for second in range(second_floor, apex_cap + 1):
            require(second % canonical_l != 0, "bounded z2 became aligned")
            ledger.second_candidates += 1
            residual_two = A.subtract(residual_one, danger(second))
            require(
                A.mass(residual_two) > 0 and len(residual_two) > 0,
                "two drifts unexpectedly cover the body carrier",
            )
            _, multiplier_cap = first_apex_multiplier_cap(
                residual_two,
                canonical_l,
                5,
            )
            if multiplier_cap < 1:
                ledger.second_closed_at_aligned_gate += 1
                continue
            ledger.second_survivors += 1
            ledger.second_survivor_roots.add(body)
            ledger.second_survivor_first_rows.add((body, first))
            descend_aligned(
                ledger,
                body,
                first,
                second,
                canonical_l,
                residual_two,
                remaining=5,
                previous_multiplier=0,
                prefix=(),
            )

    require(ledger.second_candidates == 42_912, "bounded z2 universe changed")
    require(
        ledger.second_closed_at_aligned_gate == 39_913,
        "five-aligned first-apex closure count changed",
    )
    require(ledger.second_survivors == 2_999, "bounded z2 survivor count changed")
    require(
        len(ledger.second_survivor_roots) == 23,
        "bounded z2 survivor-root count changed",
    )
    require(
        len(ledger.second_survivor_first_rows) == 25,
        "bounded z2 survivor-(E,z1) count changed",
    )
    require(
        ledger.nodes_by_remaining == Counter({5: 2_999, 4: 3_138, 3: 907, 2: 79}),
        "aligned recursion node census changed",
    )
    require(
        ledger.children_by_remaining == Counter({5: 3_138, 4: 907, 3: 79}),
        "aligned recursion child census changed",
    )
    require(
        ledger.closed_without_child_by_remaining
        == Counter({4: 2_241, 3: 828, 2: 79}),
        "aligned recursion closure census changed",
    )
    require(
        ledger.terminal_positive == 0
        and not ledger.terminal_zero_mass
        and ledger.minimum_terminal_mass is None,
        "aligned recursion unexpectedly reached a terminal",
    )
    require(ledger.maximum_multiplier == 4, "aligned multiplier maximum changed")

    digest = hashlib.sha256(
        b"LRC14/five-aligned/two-drift/high-excess-recursion/v1\n"
        + repr(tuple(ledger.digest_rows)).encode()
    ).hexdigest()
    require(
        digest == "4ffaae21ccaf438b65ae64548fcf0091031e6ccaf29f2046d613fe90085272ca",
        f"aligned recursion digest changed: {digest}",
    )

    print("LRC14 high-excess two-drift recursive aligned closure")
    print(
        "universe=(257 typed residual rows,z2 in projected-safe-floor"
        "..six-tail-first-apex-cap,then five increasing aligned multipliers)"
    )
    print(
        f"z2_candidates={ledger.second_candidates};"
        f"closed_at_first_aligned_gate={ledger.second_closed_at_aligned_gate};"
        f"z2_survivors={ledger.second_survivors};"
        f"survivor_roots={len(ledger.second_survivor_roots)};"
        f"survivor_(E,z1)_rows={len(ledger.second_survivor_first_rows)}"
    )
    print(f"nodes_by_remaining={tuple(sorted(ledger.nodes_by_remaining.items()))}")
    print(
        f"children_by_remaining={tuple(sorted(ledger.children_by_remaining.items()))}"
    )
    print(
        "closed_without_child_by_remaining="
        f"{tuple(sorted(ledger.closed_without_child_by_remaining.items()))}"
    )
    print(
        f"terminal_positive_mass={ledger.terminal_positive};"
        f"terminal_zero_mass={len(ledger.terminal_zero_mass)};"
        f"maximum_multiplier={ledger.maximum_multiplier}"
    )
    print(f"minimum_terminal_mass={ledger.minimum_terminal_mass}")
    print(f"recursion_digest={digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
