#!/usr/bin/env python3
"""
lrc_gap_debt_product_s398.py

codex-2026-05-31 S398

Audit the user's proposed product-law candidate:

    Archimedean gap * 2-adic endpoint debt ~= conserved.

For an LRC speed set V at denominator n=len(V)+1, use

    ArchGap(V) = max_gap(V) / (1/n)
    Debt(V)    = number of unprotected forbidden endpoints
    Product    = ArchGap(V) * Debt(V).

This is not meant to apply to the boundary-only initial segment, where
ArchGap=0 but visible endpoint debt is already the certificate.  It is meant
for the positive-gap debt-export branch: when a gate shrinks the visible
Archimedean gap, endpoint debt grows in a 2-adic quotient layer.  The raw
endpoint count is only the first debt proxy; the layer histogram records where
the 2-adic mass moved.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()


@dataclass(frozen=True)
class GapDebtRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    classification: str
    arch_gap: Fraction
    debt: int
    product: Fraction
    first_unprotected: Fraction | None
    layer_hist: tuple[tuple[int, int], ...]


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        raise ValueError(f"bad ladder n={n}, scale={scale}, skip={skip}")
    return speeds


def unprotected_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return tuple(
        point
        for point in sorted(endpoints)
        if not any(S360.direct_protects(speeds, speed, point) for speed in speeds)
    )


def debt_layer(n: int, point: Fraction) -> int:
    return point.denominator // gcd(point.denominator, n)


def layer_histogram(n: int, speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(debt_layer(n, point) for point in unprotected_points(speeds)).items()))


def summarize(label: str, speeds: tuple[int, ...]) -> GapDebtRow:
    summary = S360.summarize(list(speeds))
    n = len(summary.speeds) + 1
    arch_gap = summary.max_gap / summary.threshold
    debt = summary.unprotected_count
    return GapDebtRow(
        label=label,
        n=n,
        speeds=summary.speeds,
        classification=summary.classification,
        arch_gap=arch_gap,
        debt=debt,
        product=arch_gap * debt,
        first_unprotected=summary.first_unprotected,
        layer_hist=layer_histogram(n, summary.speeds),
    )


def print_rows(title: str, rows: tuple[GapDebtRow, ...]) -> None:
    print(title)
    print("-" * 78)
    print("  label                  class          gap/th      debt     product   first layer_hist")
    for row in rows:
        hist = ",".join(f"{layer}:{count}" for layer, count in row.layer_hist[:5])
        if len(row.layer_hist) > 5:
            hist += ",..."
        print(
            f"  {row.label:<22} {row.classification:<12} "
            f"{fmt_frac(row.arch_gap):>10} {row.debt:>7} "
            f"{fmt_frac(row.product):>11} {fmt_frac(row.first_unprotected):>7} {hist}"
        )
    print()


def print_doubling_checks(rows: tuple[GapDebtRow, ...]) -> None:
    print("DOUBLING CHECKS")
    print("-" * 78)
    print("  pair                 gap ratio  debt ratio  product ratio  note")
    for left, right in zip(rows, rows[1:]):
        if left.n != right.n:
            continue
        gap_ratio = right.arch_gap / left.arch_gap if left.arch_gap else None
        debt_ratio = Fraction(right.debt, left.debt) if left.debt else None
        product_ratio = right.product / left.product if left.product else None
        note = ""
        if product_ratio and product_ratio != 1:
            note = "layer tax"
        note_text = f"  {note}" if note else ""
        print(
            f"  {left.label:<10}->{right.label:<10} "
            f"{fmt_frac(gap_ratio):>9} {fmt_frac(debt_ratio):>11} "
            f"{fmt_frac(product_ratio):>14}{note_text}"
        )
    print()


def main() -> None:
    print("LRC GAP-DEBT PRODUCT AUDIT")
    print("=" * 78)
    print("ArchGap = max_gap / threshold.  Debt = unprotected endpoint count.")
    print("The layer histogram records the 2-adic quotient where the debt lives.")
    print()

    boundary_rows = tuple(
        summarize(f"n={n} initial", tuple(range(1, n))) for n in (14, 16, 18)
    )
    print_rows("BOUNDARY BRANCH: GAP=0 BUT DEBT IS A WITNESS", boundary_rows)

    n14_rows = (
        summarize("n14 d=7", ladder(14, 7, 6)),
        summarize("n14 d=14", ladder(14, 14, 6)),
    )
    n16_rows = (
        summarize("n16 d=2", ladder(16, 2, 14)),
        summarize("n16 d=4", ladder(16, 4, 14)),
        summarize("n16 d=8", ladder(16, 8, 14)),
        summarize("n16 d=16", ladder(16, 16, 14)),
    )
    n18_rows = (
        summarize("n18 d=9", ladder(18, 9, 8)),
        summarize("n18 d=18", ladder(18, 18, 8)),
    )

    print_rows("POSITIVE-GAP DEBT EXPORT ROWS", n14_rows + n16_rows + n18_rows)
    print_doubling_checks(n14_rows + n16_rows + n18_rows)

    print("EXACT PRODUCT PLATEAUS")
    print("-" * 78)
    print(f"  n=14: {fmt_frac(n14_rows[0].product)} = {fmt_frac(n14_rows[1].product)}")
    print(
        f"  n=16: {fmt_frac(n16_rows[0].product)} = {fmt_frac(n16_rows[1].product)}; "
        f"{fmt_frac(n16_rows[2].product)} = {fmt_frac(n16_rows[3].product)}"
    )
    print(f"  n=18: {fmt_frac(n18_rows[0].product)} = {fmt_frac(n18_rows[1].product)}")
    print("  n=16 phase step: (35/33)/(34/33) = 35/34")
    print()

    print("PROOF GRAMMAR")
    print("-" * 78)
    print("  boundary witness:      ArchGap=0, Debt>0")
    print("  positive-gap witness:  ArchGap>0, any Debt")
    print("  exported debt branch:  ArchGap*Debt >= c(n) > 0")
    print("  counterexample target: ArchGap=0 and Debt=0")
    print()
    print("A product lower bound belongs to the exported branch, not to the")
    print("initial boundary branch.  It says that if a repair shrinks the")
    print("Archimedean gap, 2-adic endpoint debt must grow quickly enough that")
    print("the pair cannot converge to (0,0).")


if __name__ == "__main__":
    main()
