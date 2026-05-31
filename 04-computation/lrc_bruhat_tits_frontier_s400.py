#!/usr/bin/env python3
"""
lrc_bruhat_tits_frontier_s400.py

codex-2026-05-31 S400

Model LRC endpoint debt as a finite Bruhat-Tits tree frontier.

For a rational endpoint x=a/b, the p-adic valuation sends x to a boundary
direction of the PGL_2(Q_p) tree.  If p divides the denominator past the base
denominator n, the endpoint lies deeper in the infinity cusp.  The relative
depth

    h_p(x;n) = max(0, v_p(b) - v_p(n))

is the denominator debt height beyond the base n-lattice.

The tree heuristic behind HYP-1866 is:

    real gap shrinks like p^(-h)
    raw endpoint debt grows like p^h * boundary_mass
    gap * raw_debt = (gap * p^h) * boundary_mass.

Thus a proof should either find a real gap, find exposed endpoint debt, or
show that the Bruhat-Tits frontier mass has a positive lower bound.
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
class Case:
    label: str
    speeds: tuple[int, ...]
    primes: tuple[int, ...]


@dataclass(frozen=True)
class CaseSummary:
    label: str
    n: int
    classification: str
    arch_gap: Fraction
    debt: int
    product: Fraction
    points: tuple[Fraction, ...]


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_hist(counter: Counter[int]) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


def vp(value: int, p: int) -> int:
    if value == 0:
        raise ValueError("vp(0) is not used in this finite endpoint audit")
    value = abs(value)
    out = 0
    while value % p == 0:
        out += 1
        value //= p
    return out


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


def summarize(case: Case) -> CaseSummary:
    summary = S360.summarize(list(case.speeds))
    n = len(summary.speeds) + 1
    arch_gap = summary.max_gap / summary.threshold
    debt = summary.unprotected_count
    return CaseSummary(
        label=case.label,
        n=n,
        classification=summary.classification,
        arch_gap=arch_gap,
        debt=debt,
        product=arch_gap * debt,
        points=unprotected_points(summary.speeds),
    )


def relative_depth(point: Fraction, n: int, p: int) -> int:
    return max(0, vp(point.denominator, p) - vp(n, p))


def infinity_depth(point: Fraction, p: int) -> int:
    return max(0, vp(point.denominator, p) - vp(point.numerator, p))


def unit_residue_in_infinity_chart(point: Fraction, p: int, modulus_power: int) -> int:
    """Return the unit residue of (1/x)/p^depth modulo p^modulus_power."""
    modulus = p**modulus_power
    if modulus == 1:
        return 0
    num = point.numerator
    den = point.denominator
    a = vp(num, p)
    b = vp(den, p)
    den_unit = den // (p**b)
    num_unit = num // (p**a)
    return (den_unit * pow(num_unit % modulus, -1, modulus)) % modulus


def bt_mass(points: tuple[Fraction, ...], n: int, p: int) -> Fraction:
    return sum(Fraction(1, p ** relative_depth(point, n, p)) for point in points)


def active_bt_mass(points: tuple[Fraction, ...], n: int, p: int) -> Fraction:
    return sum(
        Fraction(1, p ** relative_depth(point, n, p))
        for point in points
        if relative_depth(point, n, p) > 0
    )


def depth_hist(points: tuple[Fraction, ...], n: int, p: int) -> Counter[int]:
    return Counter(relative_depth(point, n, p) for point in points)


def unique_positive_depth(hist: Counter[int]) -> int | None:
    positive = [depth for depth, count in hist.items() if depth > 0 and count > 0]
    if len(positive) == 1:
        return positive[0]
    return None


def positive_depth_label(hist: Counter[int]) -> str:
    positive = [depth for depth, count in sorted(hist.items()) if depth > 0 and count > 0]
    if not positive:
        return "-"
    if len(positive) == 1:
        return str(positive[0])
    return "{" + ",".join(str(depth) for depth in positive) + "}"


def print_case_table(cases: tuple[Case, ...]) -> tuple[CaseSummary, ...]:
    rows = tuple(summarize(case) for case in cases)
    by_label = {row.label: row for row in rows}
    print("BRUHAT-TITS FRONTIER AUDIT FOR LRC ENDPOINT DEBT")
    print("=" * 88)
    print("h_p = max(0, v_p(den endpoint) - v_p(n)).")
    print("BT mass = sum p^(-h_p) over unprotected endpoints.")
    print()
    print("CASE SUMMARY")
    print("-" * 88)
    print("  label          class          gap/th      debt     product")
    for row in rows:
        print(
            f"  {row.label:<14} {row.classification:<12} "
            f"{fmt_frac(row.arch_gap):>10} {row.debt:>7} {fmt_frac(row.product):>11}"
        )
    print()

    print("P-ADIC FRONTIERS")
    print("-" * 88)
    print("  label          p  depth_hist        BT_mass  active_mass  gap*p^h")
    for case in cases:
        row = by_label[case.label]
        for p in case.primes:
            hist = depth_hist(row.points, row.n, p)
            h = unique_positive_depth(hist)
            scaled_gap = row.arch_gap * (p**h) if h is not None else None
            print(
                f"  {row.label:<14} {p:>2} {fmt_hist(hist):<16} "
                f"{fmt_frac(bt_mass(row.points, row.n, p)):>8} "
                f"{fmt_frac(active_bt_mass(row.points, row.n, p)):>12} "
                f"{fmt_frac(scaled_gap):>8}"
            )
    print()
    return rows


def print_translation_checks(rows: tuple[CaseSummary, ...]) -> None:
    by_label = {row.label: row for row in rows}
    pairs = (
        ("n14 d=7", "n14 d=14", 2),
        ("n16 d=2", "n16 d=4", 2),
        ("n16 d=4", "n16 d=8", 2),
        ("n16 d=8", "n16 d=16", 2),
        ("n18 d=9", "n18 d=18", 2),
    )
    print("TREE TRANSLATION CHECKS")
    print("-" * 88)
    print("  pair                 p  hL->hR  gap ratio  raw ratio  BT mass ratio")
    for left_label, right_label, p in pairs:
        left = by_label[left_label]
        right = by_label[right_label]
        left_hist = depth_hist(left.points, left.n, p)
        right_hist = depth_hist(right.points, right.n, p)
        gap_ratio = right.arch_gap / left.arch_gap
        raw_ratio = Fraction(right.debt, left.debt)
        mass_ratio = active_bt_mass(right.points, right.n, p) / active_bt_mass(
            left.points, left.n, p
        )
        print(
            f"  {left_label:<10}->{right_label:<10} {p:>2} "
            f"{positive_depth_label(left_hist):>3}->{positive_depth_label(right_hist):<5} "
            f"{fmt_frac(gap_ratio):>10} "
            f"{fmt_frac(raw_ratio):>10} {fmt_frac(mass_ratio):>14}"
        )
    print()


def print_n16_residue_frontier(rows: tuple[CaseSummary, ...]) -> None:
    by_label = {row.label: row for row in rows}
    print("N=16 DYADIC INFINITY-CUSP RESIDUES")
    print("-" * 88)
    print("Residues are units of (1/x)/2^depth modulo 16 for endpoints x.")
    for label in ("n16 d=2", "n16 d=4", "n16 d=8", "n16 d=16"):
        row = by_label[label]
        residue_counts: Counter[int] = Counter()
        depth_counts: Counter[int] = Counter()
        for point in row.points:
            depth = infinity_depth(point, 2)
            depth_counts[depth] += 1
            residue_counts[unit_residue_in_infinity_chart(point, 2, 4)] += 1
        residues = ",".join(
            f"{res}:{residue_counts[res]}" for res in sorted(residue_counts)
        )
        print(
            f"  {label:<10} abs_depths={fmt_hist(depth_counts):<8} "
            f"unit_residues_mod16={residues}"
        )
    print()


def main() -> None:
    cases = (
        Case("n14 initial", tuple(range(1, 14)), (2, 7)),
        Case("n14 d=7", ladder(14, 7, 6), (2, 7)),
        Case("n14 d=14", ladder(14, 14, 6), (2, 7)),
        Case("n16 initial", tuple(range(1, 16)), (2,)),
        Case("n16 d=2", ladder(16, 2, 14), (2,)),
        Case("n16 d=4", ladder(16, 4, 14), (2,)),
        Case("n16 d=8", ladder(16, 8, 14), (2,)),
        Case("n16 d=16", ladder(16, 16, 14), (2,)),
        Case("n18 initial", tuple(range(1, 18)), (2, 3)),
        Case("n18 d=9", ladder(18, 9, 8), (2, 3)),
        Case("n18 d=18", ladder(18, 18, 8), (2, 3)),
    )
    rows = print_case_table(cases)
    print_translation_checks(rows)
    print_n16_residue_frontier(rows)

    print("PROOF TRANSLATION")
    print("-" * 88)
    print("  LRC endpoint x=a/b     -> boundary point of P^1(Q_p)")
    print("  extra denominator p^h  -> depth h in the infinity cusp")
    print("  raw endpoint debt      -> horosphere population")
    print("  BT mass sum p^(-h)     -> normalized frontier measure")
    print("  real gap shrink        -> translation toward the cusp")
    print()
    print("A Bruhat-Tits proof would show that every repaired endpoint core")
    print("has positive normalized frontier mass in at least one p-tree.  Then")
    print("moving deeper can hide the real gap only by creating more endpoint")
    print("frontier, so the pair (Archimedean gap, p-adic debt) cannot both")
    print("vanish.")
    print()
    print("The n=16 rows are the clean single-tree lab: gap*2^h is constant")
    print("while BT mass is 17,17,35/2,35/2.  The 35/34 jump is therefore")
    print("a frontier-mass jump, not a failure of tree translation.  The n=14")
    print("and n=18 rows already live in a product of p-trees, so a final")
    print("proof should use an adelic product building rather than one tree.")


if __name__ == "__main__":
    main()
