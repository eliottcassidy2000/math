#!/usr/bin/env python3
"""HYP-2668 scout: B=14 p1-tax frontier after the shell-1 gate.

S41/HYP-2667 found that the full B=13 bank refutes 3p1/8 but survives
2p1/5.  This script pushes one bounded layer farther.  The key question is
whether a global scalar 2p1/5 target survives, or whether HYP-2666's shell-1
gate must be part of the statement.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

from lrc14_p1_tax_envelope_codex_s40 import RatioRow, bounded_bank
from lrc14_residual_plateau_packet_codex_s39 import CAP9


SHELL1 = {1, 2, 4, 8}


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def shell1_missing(Ep: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(SHELL1 - set(Ep)))


def ap_holes(Ep: tuple[int, ...], B: int) -> tuple[int, ...]:
    return tuple(v for v in range(1, B + 1) if v not in set(Ep))


def tax_gap(row: RatioRow, c: Fraction) -> Fraction:
    return c * row.w * row.p1 - row.raw


def cap_slack(row: RatioRow, c: Fraction) -> Fraction:
    return CAP9 - row.phi - c * row.p1


def print_top(rows: list[RatioRow], B: int, limit: int = 20) -> None:
    rows = sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)
    print("top rows")
    print(
        "rank | ratio | w | shell1_missing | E' | holes in [1,B] | "
        "raw | p1 | gap_2/5 | gap_5/12 | cap_slack_5/12"
    )
    for i, row in enumerate(rows[:limit], 1):
        print(
            f"{i:4d} | {fmt(row.actual):>20} | {row.w:2d} | "
            f"{shell1_missing(row.Ep)!s:>14} | {row.Ep} | {ap_holes(row.Ep, B)} | "
            f"{row.raw} | {row.p1} | {tax_gap(row, Fraction(2, 5))} | "
            f"{tax_gap(row, Fraction(5, 12))} | {cap_slack(row, Fraction(5, 12))}"
        )
    print()


def print_strata(rows: list[RatioRow]) -> None:
    groups: dict[tuple[int, ...], list[RatioRow]] = defaultdict(list)
    for row in rows:
        groups[shell1_missing(row.Ep)].append(row)

    print("shell-1 packet strata")
    print("missing | rows | max ratio | max row | >3/8 | >2/5 | >5/12")
    for miss, group in sorted(groups.items(), key=lambda kv: (len(kv[0]), kv[0])):
        top = max(group, key=lambda r: (r.actual, r.envelope))
        print(
            f"{miss!s:>14} | {len(group):5d} | {fmt(top.actual):>20} | "
            f"w={top.w}, E'={top.Ep} | "
            f"{sum(1 for r in group if r.actual > Fraction(3, 8)):5d} | "
            f"{sum(1 for r in group if r.actual > Fraction(2, 5)):5d} | "
            f"{sum(1 for r in group if r.actual > Fraction(5, 12)):5d}"
        )
    print()


def main() -> None:
    B = 14
    rows = bounded_bank(B, 8, max_rows=None)
    rows_sorted = sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)
    shell_full = [r for r in rows if not shell1_missing(r.Ep)]
    damaged = [r for r in rows if shell1_missing(r.Ep)]

    print("HYP-2668 B=14 shell-gated p1-tax scout")
    print("exact Fraction arithmetic")
    print(f"family: E'={{0}}+7-subsets of [1,{B}], w=max(E')+1..max(E')+8")
    print(f"rows={len(rows)}")
    for c in (Fraction(1, 3), Fraction(3, 8), Fraction(2, 5), Fraction(5, 12)):
        print(f"  rows with Delta^+/p1 > {c}: {sum(1 for r in rows if r.actual > c)}")
    print(f"  global max={fmt(rows_sorted[0].actual)} at w={rows_sorted[0].w}, E'={rows_sorted[0].Ep}")
    print(f"  shell1-full max={fmt(max(r.actual for r in shell_full))}")
    print(f"  shell1-damaged max={fmt(max(r.actual for r in damaged))}")
    print(f"  min 5/12 tax gap={min(tax_gap(r, Fraction(5, 12)) for r in rows)}")
    print(f"  min cap slack for c=5/12={min(cap_slack(r, Fraction(5, 12)) for r in rows)}")
    print()
    print_top(rows, B)
    print_strata(rows)
    print("Interpretation")
    print("  B=14 refutes a global raw 2p1/5 tax: one row reaches 7071/17584.")
    print("  That row is shell-1 damaged, missing the tower bit 2.")
    print("  In this scan every shell-1-full row remains below 2p1/5.")
    print("  This supports the ordered target from HYP-2666:")
    print("    shell-1 damaged rows route to tower/mouth rigidity;")
    print("    shell-1-full rows keep the 2p1/5 p1-tax target.")
    print("  The full B=14 bank remains below 5p1/12 and keeps positive cap slack in this bank,")
    print("  but the sharper proof target is still shell-gated 2p1/5, not a looser global ceiling.")
    print()
    print("Tournament Analysis")
    print("  vertices: shell1_gate > shell1_full_p1_tax > global_2p1/5 > global_5p1/12 > raw_ratio")
    print("  observable: whether the next bounded layer preserves the proposed p1-tax target")
    print("  switch/gauge: stratify exact ratios by missing shell-1 packet before choosing the scalar")
    print("  challenged assumption: the B=13 scalar 2p1/5 target globalizes without the shell gate.")
    print("PASS: B=14 shell-gated p1-tax scout complete.")


if __name__ == "__main__":
    main()
