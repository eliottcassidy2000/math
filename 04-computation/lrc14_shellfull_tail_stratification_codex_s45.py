#!/usr/bin/env python3
"""HYP-2672 scout: shell-full tail stratification after HYP-2671.

HYP-2671 localizes the shell-full new-speed 1/3 constant at the dyadic block

    E'=(0,1,2,4,8,12,16,20), w=24.

This continuation asks what remains once that block is treated as an addressed
exception.  The script extends the exact shell-full quotient to B=36 and
separates:

* the finite high pocket max(E') <= 14;
* the new-speed/dyadic-block band max(E') <= 20;
* the intermediate rows 21 <= max(E') <= 24;
* the far tail max(E') > 24, where HYP-2670 suggested a p1/4 barrier.

Only exact Fraction arithmetic is used for the raw Delta^+/p1 ratios.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction

from lrc14_shellfull_packet_gap_codex_s44 import (
    LightRow,
    fold_profile,
    fold_recip_mass,
    odd_carry_profile,
    scan_shellfull_light,
    shell_extras,
)


BMAX = 36
W_EXTRA = 8


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def fold_mass(row: LightRow) -> Fraction:
    return fold_recip_mass(fold_profile(row.Ep))


def tax_gap(row: LightRow, c: Fraction) -> Fraction:
    return c * row.w * row.p1 - row.raw


def band_name(row: LightRow) -> str:
    m = max(row.Ep)
    if m <= 14:
        return "finite <=14"
    if m <= 20:
        return "dyadic/new 15..20"
    if m <= 24:
        return "intermediate 21..24"
    if m <= 30:
        return "tail 25..30"
    return "tail 31..36"


def print_band_table(rows: list[LightRow]) -> None:
    bands = (
        "finite <=14",
        "dyadic/new 15..20",
        "intermediate 21..24",
        "tail 25..30",
        "tail 31..36",
        "all >14",
        "all >20",
        "all >24",
        "all >30",
    )

    print("band table")
    print("band | rows | max ratio | gap 1/3 | gap 3/10 | gap 1/4 | >1/3 | >3/10 | >1/4 | max row")
    for band in bands:
        if band == "all >14":
            group = [r for r in rows if max(r.Ep) > 14]
        elif band == "all >20":
            group = [r for r in rows if max(r.Ep) > 20]
        elif band == "all >24":
            group = [r for r in rows if max(r.Ep) > 24]
        elif band == "all >30":
            group = [r for r in rows if max(r.Ep) > 30]
        else:
            group = [r for r in rows if band_name(r) == band]
        if not group:
            continue
        top = max(group, key=lambda r: (r.actual, r.raw))
        print(
            f"{band:20s} | {len(group):5d} | {fmt(top.actual):>20} | "
            f"{Fraction(1,3)-top.actual} | {Fraction(3,10)-top.actual} | "
            f"{Fraction(1,4)-top.actual} | "
            f"{sum(r.actual > Fraction(1,3) for r in group):4d} | "
            f"{sum(r.actual > Fraction(3,10) for r in group):5d} | "
            f"{sum(r.actual > Fraction(1,4) for r in group):4d} | "
            f"w={top.w}, E'={top.Ep}"
        )
    print()


def print_frontier(rows: list[LightRow], predicate, title: str, limit: int = 24) -> None:
    group = [r for r in rows if predicate(r)]
    print(title)
    print(f"rows={len(group)}")
    if not group:
        print()
        return
    print("rank | ratio | w | max(E') | band | E' | extras | fold_recip | odd-carry | gap 1/3 | gap 1/4")
    for rank, row in enumerate(sorted(group, key=lambda r: (r.actual, r.raw), reverse=True)[:limit], 1):
        print(
            f"{rank:4d} | {fmt(row.actual):>20} | {row.w:2d} | {max(row.Ep):2d} | "
            f"{band_name(row):20s} | {row.Ep} | {shell_extras(row.Ep)} | "
            f"{fold_mass(row)} | {odd_carry_profile(row.Ep)} | "
            f"{Fraction(1,3)-row.actual} | {Fraction(1,4)-row.actual}"
        )
    print()


def print_by_max(rows: list[LightRow]) -> None:
    by_max: dict[int, list[LightRow]] = defaultdict(list)
    for row in rows:
        by_max[max(row.Ep)].append(row)

    print("best row by exact max(E') layer")
    print("max(E') | rows | max ratio | gap 3/10 | gap 1/4 | w | E' | fold_recip")
    for m in sorted(by_max):
        group = by_max[m]
        top = max(group, key=lambda r: (r.actual, r.raw))
        if m <= 13 and top.actual <= Fraction(1, 4):
            continue
        print(
            f"{m:7d} | {len(group):4d} | {fmt(top.actual):>20} | "
            f"{Fraction(3,10)-top.actual} | {Fraction(1,4)-top.actual} | "
            f"{top.w:2d} | {top.Ep} | {fold_mass(top)}"
        )
    print()


def print_tail_features(rows: list[LightRow]) -> None:
    tail = [r for r in rows if max(r.Ep) > 24]
    top_tail = sorted(tail, key=lambda r: (r.actual, r.raw), reverse=True)[:18]
    fold_bins: Counter[int] = Counter((fold_mass(r).numerator * 20) // fold_mass(r).denominator for r in tail)
    odd_patterns: Counter[tuple[tuple[int, int], ...]] = Counter(odd_carry_profile(r.Ep) for r in top_tail)

    print("tail feature guardrails")
    print(f"tail rows={len(tail)}")
    print(f"tail rows above 1/4={sum(r.actual > Fraction(1,4) for r in tail)}")
    print(f"tail rows above 3/10={sum(r.actual > Fraction(3,10) for r in tail)}")
    print(f"tail fold twentieth-bin histogram={tuple(sorted(fold_bins.items()))}")
    print(f"top-tail odd-carry patterns={tuple(odd_patterns.items())}")
    print()


def main() -> None:
    print("HYP-2672 shell-full tail stratification scout")
    print("exact Fraction arithmetic; reuses HYP-2670/HYP-2671 shell-full quotient")
    print(f"family: E'={{0}}+{{1,2,4,8}}+3 extras from [1,{BMAX}], w=max(E')+1..max(E')+{W_EXTRA}")
    rows = scan_shellfull_light(BMAX, W_EXTRA)
    print(f"rows={len(rows)}")
    print()

    print_band_table(rows)
    print_frontier(rows, lambda r: r.actual > Fraction(3, 10), "rows above 3/10")
    print_frontier(
        rows,
        lambda r: max(r.Ep) > 20 and r.actual > Fraction(1, 4),
        "post-dyadic rows with max(E')>20 above 1/4",
    )
    print_frontier(
        rows,
        lambda r: max(r.Ep) > 24 and r.actual > Fraction(1, 4),
        "far-tail rows with max(E')>24 above 1/4",
    )
    print_by_max(rows)
    print_tail_features(rows)

    print("Interpretation")
    print("  HYP-2671 should own the m=4 dyadic new-speed 1/3 block.")
    print("  HYP-2672 isolates the remaining work: finite high pocket, finite")
    print("  intermediate >1/4 rows, and a corrected tail target.")
    print("  The naive HYP-2670 far-tail statement max(E')>24 => p1/4 is false")
    print("  at B36: exactly one low-fold doubled-odd row above 1/4 appears.")
    print("  The surviving broad tail target in this scan is instead <3/10 after")
    print("  the HYP-2671 dyadic block is removed.")
    print()
    print("Tournament Analysis")
    print("  vertices: doubled_odd_tail_exception > tail_3/10_decay > intermediate_21_24_ledger > dyadic_block_exception > finite_B13_pocket")
    print("  observable: exact Delta^+/p1 after shell-1-full quotienting")
    print("  switch/gauge: remove the HYP-2671 dyadic block, then sort by max(E') band")
    print("  Hamiltonian path: doubled_odd_tail_exception > tail_3/10_decay > intermediate_21_24_ledger > dyadic_block_exception > finite_B13_pocket")
    print("  challenged assumption: one constant explains all shell-full rows; after HYP-2671 the remaining tail is a different obligation.")
    print("PASS: shell-full tail stratification scout complete.")


if __name__ == "__main__":
    main()
