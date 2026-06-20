#!/usr/bin/env python3
"""HYP-2666/T908: two-gate boundary currency scout for LRC14.

Recent work split into two complementary HYP-2664 notes:

* near AP-tail: shell-1 tower conservation removes most three-tail comb burden;
* far Delta: the intra-quadratic residual looks paid by p1 boundary mass.

This script tests whether those are the same proof order in two gauges.  It
widens the residual p1-tax bounded bank

    E' = {0} union A, |A|=7, A subset [1,B],

and stratifies by whether the shell-1 tower {1,2,4,8} is present.  The tested
certificate is

    p0(E') + (1/7+c) p1(E') <= cap_9,

with c in {1/4, 13/51, 1/3, 3/8}.  If dangerous p1-tax rows live in the shell-1-full
stratum while shell-1-damaged rows have slack, then HYP-2661 and the p1-tax
route are one two-gate boundary currency:

    shell-1 gate first, p1 boundary tax second.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm

from lrc14_residual_plateau_packet_codex_s39 import CAP9, fmt


SHELL1 = frozenset({1, 2, 4, 8})
INNER = tuple(range(1, 7))
COEFFS = (Fraction(1, 4), Fraction(13, 51), Fraction(1, 3), Fraction(3, 8))


@dataclass(frozen=True)
class BankRow:
    Ep: tuple[int, ...]
    p0: Fraction
    p1: Fraction
    values: tuple[Fraction, ...]

    @property
    def shell1_missing(self) -> tuple[int, ...]:
        return tuple(sorted(SHELL1 - set(self.Ep)))

    @property
    def shell1_full(self) -> bool:
        return not self.shell1_missing

    @property
    def min_allowed_c(self) -> Fraction | None:
        if self.p1 == 0:
            return None
        return (CAP9 - self.p0) / self.p1 - Fraction(1, 7)


def primitive(Ep: tuple[int, ...]) -> bool:
    return reduce(gcd, Ep, 0) == 1


def missed_distribution_fast(Ep: tuple[int, ...]) -> dict[int, Fraction]:
    """Exact missed-sector distribution using integer common-wall coordinates."""

    denom = 1
    for e in sorted(set(Ep)):
        if e:
            denom = lcm(denom, 7 * e)

    breakpoints = {0, denom}
    for e in sorted(set(Ep)):
        if not e:
            continue
        step = denom // (7 * e)
        for a in range(0, 7 * e + 1):
            breakpoints.add(a * step)

    ordered = sorted(breakpoints)
    totals = [0] * 7
    mid_denom = 2 * denom
    nonzero = tuple(e for e in sorted(set(Ep)) if e)

    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        mid_num = lo + hi
        hit = [False] * 7
        hit[0] = True
        for e in nonzero:
            rem = (e * mid_num) % mid_denom
            sector = (7 * rem) // mid_denom
            hit[sector] = True
        missed = sum(1 for s in INNER if not hit[s])
        totals[missed] += hi - lo

    return {i: Fraction(totals[i], denom) for i in range(7)}


def qdict(values: list[Fraction]) -> dict[str, Fraction]:
    values = sorted(values)
    if not values:
        return {}

    def q(p: float) -> Fraction:
        return values[min(len(values) - 1, int(round(p * (len(values) - 1))))]

    return {
        "min": values[0],
        "q25": q(0.25),
        "median": q(0.50),
        "q75": q(0.75),
        "q90": q(0.90),
        "q99": q(0.99),
        "max": values[-1],
    }


def build_rows(bound: int) -> list[BankRow]:
    rows: list[BankRow] = []
    for comb in combinations(range(1, bound + 1), 7):
        Ep = (0,) + comb
        if not primitive(Ep):
            continue
        dist = missed_distribution_fast(Ep)
        p0 = dist[0]
        p1 = dist[1]
        values = tuple(p0 + (Fraction(1, 7) + c) * p1 for c in COEFFS)
        rows.append(BankRow(Ep=Ep, p0=p0, p1=p1, values=values))
    return rows


def summarize(rows: list[BankRow], label: str) -> None:
    print(f"{label}")
    print(f"  rows={len(rows)}")
    if not rows:
        return
    for coeff, idx in zip(COEFFS, range(len(COEFFS))):
        max_row = max(rows, key=lambda r: r.values[idx])
        violations = sum(1 for r in rows if r.values[idx] > CAP9)
        slack = CAP9 - max_row.values[idx]
        print(
            f"  c={coeff}: violations={violations}, max={fmt(max_row.values[idx])}, "
            f"slack={fmt(slack)}, row={max_row.Ep}, shell1_missing={max_row.shell1_missing}"
        )
    allowed = [r.min_allowed_c for r in rows if r.min_allowed_c is not None]
    min_row = min((r for r in rows if r.min_allowed_c is not None), key=lambda r: r.min_allowed_c)
    print(
        f"  min allowed c={fmt(min(allowed))}, row={min_row.Ep}, "
        f"p0={fmt(min_row.p0)}, p1={fmt(min_row.p1)}, shell1_missing={min_row.shell1_missing}"
    )
    print(f"  p1 quantiles: { {k: str(v) for k, v in qdict([r.p1 for r in rows]).items()} }")
    print()


def missing_packet_summary(rows: list[BankRow]) -> None:
    packets: dict[tuple[int, ...], list[BankRow]] = defaultdict(list)
    for row in rows:
        packets[row.shell1_missing].append(row)

    print(f"shell-1 packet summary at c={COEFFS[-1]}")
    print("missing | rows | max_value | slack | row | p0 | p1")
    for missing, bucket in sorted(packets.items(), key=lambda item: (len(item[0]), item[0])):
        max_row = max(bucket, key=lambda r: r.values[-1])
        slack = CAP9 - max_row.values[-1]
        label = "full" if not missing else str(missing)
        print(
            f"{label:18s} | {len(bucket):5d} | {str(max_row.values[-1]):>12s} | "
            f"{str(slack):>12s} | {max_row.Ep!s:>30s} | {max_row.p0} | {max_row.p1}"
        )
    print()


def top_rows(rows: list[BankRow], keep: int) -> None:
    ranked = sorted(rows, key=lambda r: r.values[-1], reverse=True)[:keep]
    print(f"top {keep} rows by c={COEFFS[-1]} p1-tax value")
    print("rank | value | slack | shell1_missing | Ep | p0 | p1 | allowed_c")
    for i, row in enumerate(ranked, 1):
        value = row.values[-1]
        print(
            f"{i:3d} | {str(value):>13s} | {str(CAP9 - value):>13s} | "
            f"{str(row.shell1_missing):>14s} | {row.Ep!s:>34s} | "
            f"{row.p0} | {row.p1} | {row.min_allowed_c}"
        )
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=18)
    parser.add_argument("--keep", type=int, default=30)
    args = parser.parse_args()

    rows = build_rows(args.bound)
    full = [r for r in rows if r.shell1_full]
    damaged = [r for r in rows if not r.shell1_full]

    print("HYP-2666/T908 LRC14 two-gate boundary currency scout")
    print(f"bank: E'={{0}} union 7-subsets of [1,{args.bound}], primitive")
    print(f"cap9={fmt(CAP9)}")
    print(f"coefficients={', '.join(str(c) for c in COEFFS)}")
    print()

    summarize(rows, "all rows")
    summarize(full, "shell-1-full rows")
    summarize(damaged, "shell-1-damaged rows")
    missing_packet_summary(rows)
    top_rows(rows, args.keep)

    print("invariant reading")
    print("1. Shell-1 status is a gate, not a scalar. Damaged rows should route")
    print("   through HYP-2661 before far residual estimates are invoked.")
    print("2. The p1-tax inequality is the far-side boundary currency. In this")
    print(f"   bounded bank, coefficients up to c={COEFFS[-1]} are tested after splitting by shell-1 packet.")
    print("3. The proof order is two-gate: shell-1 conservation first, then p1")
    print("   boundary tax for the remaining shell-1-full or nonlocal rows.")
    print()
    print("Tournament Analysis")
    print("  vertices: shell1_gate, p1_boundary_tax, missing_packet, cap_slack, raw_value")
    print("  observable: preservation of p0+(1/7+c)p1 <= cap9 through c=3/8")
    print("  switch/gauge: split by shell-1 missing packet before scalarizing p1 tax")
    print("  Hamiltonian path: shell1_gate > p1_boundary_tax > missing_packet > cap_slack > raw_value")
    print("  challenged assumption: far Delta and AP-tail shell gates are separate; they are two gauges of boundary currency.")


if __name__ == "__main__":
    main()
