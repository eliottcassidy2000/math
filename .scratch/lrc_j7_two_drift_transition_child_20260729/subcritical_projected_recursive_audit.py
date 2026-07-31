#!/usr/bin/env python3
"""Independent THM-2893 recursion on projected subcritical k=5 carriers.

For each exact excess-admissible drift pair, intersect the local two-drift
danger unions over body-safe cells until the complementary projected carrier
``P_s`` reaches mass ``887/1365``.  The full projected residual contains
``P_s``.  Therefore any five aligned completion would give five normalized
integer combs covering this exact prefix carrier.

This script does *not* use the sharp five-comb union cap to close that
carrier.  It invokes the cap-free THM-2893 first-apex gate recursively:

    next a <= floor(6 p r / (7 (7-p) h))

for ``p=5,4,...,1``, always retaining increasing multiplier order and
subtracting the literal rational danger comb.  It is an independent replay
of the projected-measure closure.
"""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import residual_first_apex_audit as A
import subcritical_exact_pair_bank as B
import subcritical_projected_pair_filter as P


TARGET_MASS = P.ALIGNED_UNION_CAP


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def normalized_danger(speed: int) -> tuple[tuple[F, F], ...]:
    return A.danger(speed)


def complement(
    rows: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    cursor = F(0)
    for left, right in rows:
        if cursor < left:
            out.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        out.append((cursor, F(1)))
    return tuple(out)


class Ledger:
    def __init__(self) -> None:
        self.nodes: Counter[int] = Counter()
        self.children: Counter[int] = Counter()
        self.closed_no_child: Counter[int] = Counter()
        self.terminal_positive = 0
        self.terminal_zero: list[tuple[object, ...]] = []
        self.max_multiplier = 0
        self.max_apex_cap = 0
        self.rows: list[tuple[object, ...]] = []


def descend(
    ledger: Ledger,
    carrier: tuple[tuple[F, F], ...],
    remaining: int,
    previous: int,
    prefix: tuple[int, ...],
    pair_key: tuple[object, ...],
) -> None:
    require(1 <= remaining <= 5, "bad projected recursion depth")
    h = P.mass(carrier)
    components = len(carrier)
    require(h > 0 and components > 0, "empty recursion carrier")
    ratio = F(
        6 * remaining * components,
        7 * (7 - remaining),
    ) / h
    apex_cap = ratio.numerator // ratio.denominator
    ledger.nodes[remaining] += 1
    ledger.max_apex_cap = max(ledger.max_apex_cap, apex_cap)
    if previous + 1 > apex_cap:
        ledger.closed_no_child[remaining] += 1
        return

    for multiplier in range(previous + 1, apex_cap + 1):
        ledger.children[remaining] += 1
        ledger.max_multiplier = max(ledger.max_multiplier, multiplier)
        child = A.subtract(carrier, normalized_danger(multiplier))
        child_mass = P.mass(child)
        row = (
            *pair_key,
            remaining,
            prefix + (multiplier,),
            h,
            components,
            apex_cap,
            child_mass,
            len(child),
        )
        ledger.rows.append(row)
        if child_mass == 0:
            ledger.terminal_zero.append(row)
            continue
        if remaining == 1:
            ledger.terminal_positive += 1
            continue
        descend(
            ledger,
            child,
            remaining - 1,
            multiplier,
            prefix + (multiplier,),
            pair_key,
        )


def audit_body(body: tuple[int, ...]) -> dict[str, object]:
    bank = B.profile(body)
    exact_pairs = bank["exact_pairs"]
    ledger = Ledger()
    if not exact_pairs:
        return {
            "body": body,
            "pairs": 0,
            "prefix_max": 0,
            "ledger": ledger,
        }

    carrier_i = A.integer_carrier(body)
    canonical_l = exact_pairs[0][3]
    cells = P.body_cells(carrier_i, canonical_l)
    first_cache: dict[tuple[int, int], tuple[tuple[F, F], ...]] = {}
    maximum_cells = 0

    for pair in exact_pairs:
        first, second = pair[1], pair[2]
        common: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
        cells_used = 0
        for cell in cells:
            first_phase = first_cache.get((first, cell))
            if first_phase is None:
                first_phase = P.phase_danger(cell, first, canonical_l)
                first_cache[(first, cell)] = first_phase
            local_union = P.merge(
                [*first_phase, *P.phase_danger(cell, second, canonical_l)]
            )
            common = P.intersect(common, local_union)
            cells_used += 1
            if 1 - P.mass(common) >= TARGET_MASS:
                break
        maximum_cells = max(maximum_cells, cells_used)
        prefix_carrier = complement(common)
        require(
            P.mass(prefix_carrier) >= TARGET_MASS,
            "projected prefix failed target mass",
        )
        descend(
            ledger,
            prefix_carrier,
            remaining=5,
            previous=0,
            prefix=(),
            pair_key=(body, first, second, canonical_l, cells_used),
        )

    return {
        "body": body,
        "pairs": len(exact_pairs),
        "prefix_max": maximum_cells,
        "ledger": ledger,
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
        audits = [audit_body(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            audits = list(pool.imap(audit_body, roots, chunksize=1))
    audits.sort(key=lambda row: row["body"])

    pairs = sum(row["pairs"] for row in audits)
    prefix_max = max(row["prefix_max"] for row in audits)
    nodes: Counter[int] = Counter()
    children: Counter[int] = Counter()
    closed: Counter[int] = Counter()
    terminal_positive = 0
    terminal_zero: list[tuple[object, ...]] = []
    maximum_multiplier = 0
    maximum_apex_cap = 0
    digest_rows: list[tuple[object, ...]] = []
    for audit in audits:
        ledger = audit["ledger"]
        nodes.update(ledger.nodes)
        children.update(ledger.children)
        closed.update(ledger.closed_no_child)
        terminal_positive += ledger.terminal_positive
        terminal_zero.extend(ledger.terminal_zero)
        maximum_multiplier = max(maximum_multiplier, ledger.max_multiplier)
        maximum_apex_cap = max(maximum_apex_cap, ledger.max_apex_cap)
        digest_rows.extend(ledger.rows)
    digest = hashlib.sha256(
        b"LRC14/k5/subcritical-projected-THM2893-recursion/v1\n"
        + repr(tuple(digest_rows)).encode()
    ).hexdigest()

    print("LRC14 k5 subcritical projected-prefix THM2893 recursion")
    print(
        f"exact_pairs={pairs};maximum_prefix_cells={prefix_max};"
        f"target_mass={TARGET_MASS}"
    )
    print(f"nodes_by_remaining={tuple(sorted(nodes.items()))}")
    print(f"children_by_remaining={tuple(sorted(children.items()))}")
    print(f"closed_no_child_by_remaining={tuple(sorted(closed.items()))}")
    print(
        f"terminal_positive={terminal_positive};"
        f"terminal_zero={len(terminal_zero)};"
        f"maximum_multiplier={maximum_multiplier};"
        f"maximum_apex_cap={maximum_apex_cap}"
    )
    print(f"recursion_digest={digest}")
    if terminal_zero:
        print(f"zero_terminals={tuple(terminal_zero)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
