#!/usr/bin/env python3
"""
lrc_adic_column_modes_s421.py

codex-2026-05-31 S421

Use the repo's 2-adic column-family picture for LRC.

Every natural number is n = 2^r * odd.  The odd roots form an x+2 chain,
and each root has a vertical x*2 chain.  Tournament history calls these:

    column/top-row step: odd -> odd+2
    row step:            n -> 2n

For LRC, the unit-boundary skeleton of {1,...,n-1} has size phi(n).
This script checks how that boundary debt and a designed largest-proper-divisor
quotient ladder behave along the same two modes.
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


MAX_N = 18


@dataclass(frozen=True)
class GridRow:
    n: int
    r: int
    odd: int
    column: int
    phi: int
    unit_density: Fraction
    lpd: int | None
    ladder_skip: int | None
    ladder_gap: Fraction | None
    ladder_debt: int | None
    ladder_product: Fraction | None
    bt2_mass: Fraction | None
    odd_bt_mass: Fraction | None


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def fmt_int(value: int | None) -> str:
    if value is None:
        return "-"
    return str(value)


def vp(value: int, p: int) -> int:
    value = abs(value)
    out = 0
    while value and value % p == 0:
        out += 1
        value //= p
    return out


def adic_row(n: int) -> tuple[int, int, int]:
    r = vp(n, 2)
    odd = n >> r
    column = (odd + 1) // 2
    return r, odd, column


def totient(n: int) -> int:
    return sum(1 for a in range(1, n) if gcd(a, n) == 1)


def proper_divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(2, n) if n % d == 0)


def largest_proper_divisor(n: int) -> int | None:
    divs = proper_divisors(n)
    return divs[-1] if divs else None


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def lpd_ladder(n: int, skip: int) -> tuple[int, ...] | None:
    d = largest_proper_divisor(n)
    if d is None:
        return None
    speeds = tuple(sorted({1} | {d * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        return None
    return speeds


def unprotected_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return tuple(
        point
        for point in sorted(endpoints)
        if not any(S360.direct_protects(speeds, speed, point) for speed in speeds)
    )


def relative_depth(point: Fraction, n: int, p: int) -> int:
    return max(0, vp(point.denominator, p) - vp(n, p))


def bt_mass(points: tuple[Fraction, ...], n: int, p: int) -> Fraction:
    return sum(Fraction(1, p ** relative_depth(point, n, p)) for point in points)


def designed_skip(n: int) -> int | None:
    _r, odd, _column = adic_row(n)
    if odd > 1:
        return odd - 1
    if n > 2:
        return n - 2
    return None


def designed_lpd_ladder_row(n: int) -> tuple[int, Fraction, int, Fraction, Fraction, Fraction] | None:
    d = largest_proper_divisor(n)
    if d is None:
        return None
    skip = designed_skip(n)
    if skip is None:
        return None
    speeds = lpd_ladder(n, skip)
    if speeds is None:
        return None
    summary = S360.summarize(list(speeds))
    gap = summary.max_gap / summary.threshold
    debt = summary.unprotected_count
    points = unprotected_points(speeds)
    odd = adic_row(n)[1]
    odd_p = smallest_prime_factor(odd)
    odd_mass = bt_mass(points, n, odd_p) if odd_p else Fraction(0)
    return skip, gap, debt, gap * debt, bt_mass(points, n, 2), odd_mass


def smallest_prime_factor(n: int) -> int | None:
    if n <= 1:
        return None
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1 if p == 2 else 2
    return n


def build_rows(max_n: int) -> dict[int, GridRow]:
    rows: dict[int, GridRow] = {}
    for n in range(2, max_n + 1):
        r, odd, column = adic_row(n)
        phi = totient(n)
        ladder = designed_lpd_ladder_row(n)
        rows[n] = GridRow(
            n=n,
            r=r,
            odd=odd,
            column=column,
            phi=phi,
            unit_density=Fraction(phi, n - 1),
            lpd=largest_proper_divisor(n),
            ladder_skip=ladder[0] if ladder else None,
            ladder_gap=ladder[1] if ladder else None,
            ladder_debt=ladder[2] if ladder else None,
            ladder_product=ladder[3] if ladder else None,
            bt2_mass=ladder[4] if ladder else None,
            odd_bt_mass=ladder[5] if ladder else None,
        )
    return rows


def chain_text(values: list[object]) -> str:
    return " -> ".join(str(value) for value in values)


def print_boundary_row_mode(rows: dict[int, GridRow]) -> None:
    print("BOUNDARY UNIT DEBT ALONG x*2 CHAINS")
    print("=" * 96)
    print("For initial segments, boundary debt = phi(n).  The first doubling is the seam.")
    print()
    print("  odd  n-chain              phi-chain            phi-ratios      phi/n-chain")
    for odd in range(1, 19, 2):
        ns: list[int] = []
        phis: list[int] = []
        densities: list[str] = []
        n = odd
        while n in rows:
            ns.append(n)
            phis.append(rows[n].phi)
            densities.append(fmt_frac(Fraction(rows[n].phi, n)))
            n *= 2
        if len(ns) < 2:
            continue
        ratios = [fmt_frac(Fraction(phis[i + 1], phis[i])) for i in range(len(phis) - 1)]
        print(
            f"  {odd:<3} {chain_text(ns):<20} {chain_text(phis):<20} "
            f"{chain_text(ratios):<15} {chain_text(densities)}"
        )
    print()


def print_first_even_column_mode(rows: dict[int, GridRow]) -> None:
    print("ODD ROOT x+2 CHAIN, FIRST EVEN CHILDREN 2*x")
    print("=" * 96)
    print("This is the LRC analogue of top-row column motion followed by first blowup.")
    print()
    print(
        "  odd  child n  lpd  skip  gap/th      debt   product   BT_2     BT_odd"
    )
    for odd in range(3, 19, 2):
        n = 2 * odd
        if n not in rows:
            continue
        row = rows[n]
        print(
            f"  {odd:<3} {n:>7} {fmt_int(row.lpd):>4} {fmt_int(row.ladder_skip):>5} "
            f"{fmt_frac(row.ladder_gap):>10} {fmt_int(row.ladder_debt):>7} "
            f"{fmt_frac(row.ladder_product):>9} {fmt_frac(row.bt2_mass):>8} "
            f"{fmt_frac(row.odd_bt_mass):>8}"
        )
    print()


def print_ladder_row_mode(rows: dict[int, GridRow]) -> None:
    print("DESIGNED QUOTIENT-LADDER ROW MODE n -> 2n")
    print("=" * 96)
    print("Compare designed largest-proper-divisor ladders inside each odd-root column.")
    print()
    print(
        "  odd  pair        gap ratio  debt ratio  product ratio  BT2 ratio  note"
    )
    for odd in range(1, 19, 2):
        n = odd
        while 2 * n in rows:
            left = rows.get(n)
            right = rows.get(2 * n)
            n *= 2
            if not left or not right or not left.ladder_gap or not right.ladder_gap:
                continue
            gap_ratio = right.ladder_gap / left.ladder_gap
            debt_ratio = Fraction(right.ladder_debt, left.ladder_debt)
            product_ratio = right.ladder_product / left.ladder_product
            bt_ratio = right.bt2_mass / left.bt2_mass if left.bt2_mass else None
            note = ""
            if left.r == 0:
                note = "first-blowup seam"
            elif product_ratio != 1:
                note = "frontier tax"
            print(
                f"  {odd:<3} {left.n:>2}->{right.n:<3} "
                f"{fmt_frac(gap_ratio):>10} {fmt_frac(debt_ratio):>11} "
                f"{fmt_frac(product_ratio):>14} {fmt_frac(bt_ratio):>10}  {note}"
            )
    print()


def print_pressure_ranking(rows: dict[int, GridRow]) -> None:
    print("SMALLEST DESIGNED FIRST-EVEN LADDERS BY VISIBLE GAP")
    print("=" * 96)
    first_even = [rows[2 * odd] for odd in range(3, 19, 2) if 2 * odd in rows]
    first_even = [row for row in first_even if row.ladder_gap is not None]
    first_even.sort(key=lambda row: (row.ladder_gap, -row.ladder_debt))
    print("  rank  n   odd  gap/th      debt   product   interpretation")
    for index, row in enumerate(first_even[:8], start=1):
        interp = "pure dyadic" if row.odd == 1 else f"first even child of {row.odd}"
        print(
            f"  {index:>4} {row.n:>3} {row.odd:>5} {fmt_frac(row.ladder_gap):>10} "
            f"{fmt_int(row.ladder_debt):>7} {fmt_frac(row.ladder_product):>9}  {interp}"
        )
    print()


def print_grid_snapshot(rows: dict[int, GridRow]) -> None:
    print("GRID SNAPSHOT")
    print("=" * 96)
    print("  n   (r,odd,k)  phi  phi/(n-1)  lpd  skip  ladder gap/debt/product")
    for n in range(2, MAX_N + 1):
        row = rows[n]
        print(
            f"  {n:>2}  ({row.r},{row.odd},{row.column})"
            f"{row.phi:>7} {fmt_float(row.unit_density):>11} "
            f"{fmt_int(row.lpd):>4} {fmt_int(row.ladder_skip):>5} "
            f"{fmt_frac(row.ladder_gap):>8}/{fmt_int(row.ladder_debt):<4}/"
            f"{fmt_frac(row.ladder_product)}"
        )
    print()


def main() -> None:
    rows = build_rows(MAX_N)
    print("LRC IN THE 2-ADIC COLUMN-FAMILY MODES")
    print("=" * 96)
    print("n = 2^r * odd.  Odd roots move by x+2; columns move by x*2.")
    print("Tournament blowup is the row step.  Mode-B/add-two is the top-row column step.")
    print()

    print_boundary_row_mode(rows)
    print_first_even_column_mode(rows)
    print_ladder_row_mode(rows)
    print_pressure_ranking(rows)
    print_grid_snapshot(rows)

    print("SYNTHESIS")
    print("=" * 96)
    print("1. Initial LRC boundary debt has the same seam shape as tournament pairs:")
    print("   odd -> first even keeps phi fixed, then every further row step doubles it.")
    print("2. The first-even children 2*odd are the natural LRC tests for x+2 motion")
    print("   among odd roots.  n=14 and n=18 are adjacent first-even children, while")
    print("   n=16 is the pure vertical column over odd root 1.")
    print("3. Row-mode proof objects should be p-adic: external n->2n and internal")
    print("   endpoint-depth translations both move down a Bruhat-Tits cusp.")
    print("4. Column-mode proof objects should be odd-payload transfer matrices:")
    print("   changing odd root by +2 changes which product-building trees carry debt.")


if __name__ == "__main__":
    main()
