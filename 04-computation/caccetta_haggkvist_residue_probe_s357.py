#!/usr/bin/env python3
"""
caccetta_haggkvist_residue_probe_s357.py

codex-2026-05-30 S357

Probe Caccetta-Haggkvist through cyclic Cayley digraphs.  For
Cay(Z/nZ, A), a directed cycle of length ell is the same as a zero-sum word
of length ell in A.  Hamidoune's Cayley proof can be read as a growth-residue
statement for B = A union {0}:

    if 0 is not in A, 2A, ..., kA, then |kB| >= k|B| - k + 1.

This script searches small cyclic groups for connection sets whose first
zero-sum length is as late as the CH bound ceil(n/|A|), and records whether
the Kemperman lower bound is tight one step before the first return.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import ceil, gcd


@dataclass(frozen=True)
class SetReport:
    n: int
    d: int
    A: tuple[int, ...]
    ch_bound: int
    first_zero: int
    growth: tuple[int, ...]
    kemperman_slack: tuple[int, ...]
    pre_return_residue: int
    is_unit_interval: bool


def sumset_mod(X: set[int], Y: set[int], n: int) -> set[int]:
    return {(x + y) % n for x in X for y in Y}


def first_zero_sum_length(n: int, A: tuple[int, ...], limit: int) -> tuple[int, tuple[int, ...]]:
    current = set(A)
    sizes = [len(current)]
    if 0 in current:
        return 1, tuple(sizes)
    for ell in range(2, limit + 1):
        current = sumset_mod(current, set(A), n)
        sizes.append(len(current))
        if 0 in current:
            return ell, tuple(sizes)
    return limit + 1, tuple(sizes)


def growth_profile_B(n: int, A: tuple[int, ...], limit: int) -> tuple[int, ...]:
    B = set(A) | {0}
    current = {0}
    sizes: list[int] = []
    for _ in range(1, limit + 1):
        current = sumset_mod(current, B, n)
        sizes.append(len(current))
    return tuple(sizes)


def unit_interval_test(n: int, A: tuple[int, ...]) -> bool:
    target = tuple(range(1, len(A) + 1))
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    aset = set(A)
    for u in units:
        image = tuple(sorted((u * a) % n for a in aset))
        if image == target:
            return True
    return False


def unit_canonical(n: int, A: tuple[int, ...]) -> tuple[int, ...]:
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    return min(tuple(sorted((u * a) % n for a in A)) for u in units)


def report_set(n: int, A: tuple[int, ...]) -> SetReport:
    d = len(A)
    ch_bound = ceil(n / d)
    first_zero, _ = first_zero_sum_length(n, A, ch_bound)
    growth = growth_profile_B(n, A, ch_bound)
    kemperman_slack = tuple(
        growth[j - 1] - (j * (d + 1) - j + 1)
        for j in range(1, ch_bound + 1)
    )
    pre = ch_bound - 1
    pre_return_residue = n - growth[pre - 1] if pre >= 1 else n - 1
    return SetReport(
        n=n,
        d=d,
        A=A,
        ch_bound=ch_bound,
        first_zero=first_zero,
        growth=growth,
        kemperman_slack=kemperman_slack,
        pre_return_residue=pre_return_residue,
        is_unit_interval=unit_interval_test(n, A),
    )


def cyclic_interval(n: int, d: int) -> tuple[int, ...]:
    return tuple(range(1, d + 1))


def enumerate_reports(max_n: int = 24) -> list[SetReport]:
    reports: list[SetReport] = []
    for n in range(4, max_n + 1):
        # The hard CH zone is ch_bound >= 3.  Above n/2, digons or triangles
        # dominate and the additive search is less informative.
        max_d = n // 3
        for d in range(1, max_d + 1):
            for A in combinations(range(1, n), d):
                reports.append(report_set(n, A))
    return reports


def summarize_by_bucket(reports: list[SetReport]) -> None:
    print("Caccetta-Haggkvist cyclic residue probe (codex-2026-05-30 S357)")
    print("Cay(Z/nZ,A): first directed cycle = first zero-sum word in A.")
    print("Enumerated all A subset {1,...,n-1} with |A| <= floor(n/3), n<=24.\n")

    buckets: dict[tuple[int, int], list[SetReport]] = {}
    for r in reports:
        buckets.setdefault((r.n, r.d), []).append(r)

    print("Bucket summary: max first-zero length against CH bound ceil(n/d)")
    for key in sorted(buckets):
        rows = buckets[key]
        n, d = key
        ch = rows[0].ch_bound
        max_first = max(r.first_zero for r in rows)
        tight = [r for r in rows if r.first_zero == ch]
        unit_tight = [r for r in tight if r.is_unit_interval]
        canon_tight = {unit_canonical(n, r.A) for r in tight}
        print(
            f"  n={n:2d} d={d:2d} sets={len(rows):6d} "
            f"bound={ch:2d} max_first={max_first:2d} "
            f"tight={len(tight):5d} canon_tight={len(canon_tight):4d} "
            f"unit_interval_tight={len(unit_tight):4d}"
        )


def print_extreme_examples(reports: list[SetReport]) -> None:
    tight = [r for r in reports if r.first_zero == r.ch_bound]
    non_interval_tight = [r for r in tight if not r.is_unit_interval]
    low_pre_residue = sorted(tight, key=lambda r: (r.pre_return_residue, r.n, r.d, r.A))[:12]
    highest_bound = sorted(tight, key=lambda r: (-r.ch_bound, r.n, r.d, r.A))[:12]

    print("\nTight examples with smallest pre-return growth residue")
    for r in low_pre_residue:
        print(
            f"  n={r.n:2d} d={r.d:2d} A={r.A} "
            f"bound={r.ch_bound} first={r.first_zero} "
            f"|jB|={r.growth} slack={r.kemperman_slack} "
            f"pre_residue={r.pre_return_residue} "
            f"unit_interval={r.is_unit_interval}"
        )

    print("\nHighest-bound tight examples")
    for r in highest_bound:
        print(
            f"  n={r.n:2d} d={r.d:2d} A={r.A} "
            f"bound={r.ch_bound} growth={r.growth} "
            f"pre_residue={r.pre_return_residue} "
            f"unit_interval={r.is_unit_interval}"
        )

    print("\nFirst non-interval tight examples")
    for r in non_interval_tight[:18]:
        print(
            f"  n={r.n:2d} d={r.d:2d} A={r.A} "
            f"bound={r.ch_bound} growth={r.growth} "
            f"slack={r.kemperman_slack} pre_residue={r.pre_return_residue}"
        )

    print("\nSummary")
    print(f"  total_sets={len(reports)}")
    print(f"  tight_sets={len(tight)}")
    print(f"  non_interval_tight_sets={len(non_interval_tight)}")
    print(
        "  tight_with_zero_pre_residue="
        f"{sum(1 for r in tight if r.pre_return_residue == 0)}"
    )
    print(
        "  tight_with_all_zero_kemperman_slack_until_return="
        f"{sum(1 for r in tight if all(s == 0 for s in r.kemperman_slack[: r.ch_bound - 1]))}"
    )
    violations = [r for r in reports if r.first_zero > r.ch_bound]
    print(f"  ch_violations={len(violations)}")


def main() -> None:
    reports = enumerate_reports(24)
    summarize_by_bucket(reports)
    print_extreme_examples(reports)

    print("\nCanonical CH-tight family A={1,...,d} in n=d(g-1)+1")
    for d, g in [(2, 4), (3, 5), (4, 6), (5, 6), (6, 5)]:
        n = d * (g - 1) + 1
        A = cyclic_interval(n, d)
        r = report_set(n, A)
        print(
            f"  d={d} g={g} n={n} A={A} "
            f"bound={r.ch_bound} first_zero={r.first_zero} "
            f"growth={r.growth} slack={r.kemperman_slack}"
        )


if __name__ == "__main__":
    main()
