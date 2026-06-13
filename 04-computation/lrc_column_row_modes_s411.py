#!/usr/bin/env python3
"""
lrc_column_row_modes_s411.py

codex-2026-05-31 S411

Use the old "natural numbers = odd +2 chains with even *2 chains" picture as
an LRC diagnostic.

Tournament history in this repo separates two size recursions:

* top-row column step: odd n -> n+2;
* doubling row step:  n -> 2n.

Every integer n = 2^r * b with b odd lives in that grid.  This script asks what
the same grid means for Lonely Runner threshold denominators.  The answer is
not just cosmetic:

* the initial-segment boundary witnesses are exactly the unit skeleton phi(n);
* along a fixed odd core b, the first doubling b -> 2b does not create new unit
  witnesses, while every later doubling doubles them;
* for even n, the largest-proper-divisor LRC quotient ladder is exactly the
  row parent n/2.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
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
S378 = SourceFileLoader(
    "natural_lrc_recursive_modes_s378",
    str(ROOT / "04-computation" / "natural_lrc_recursive_modes_s378.py"),
).load_module()


INITIAL_VERIFY_MAX = 40
LADDER_SCAN_MAX = 24
COLUMN_ODD_CORES = (3, 5, 7, 9, 15)


@dataclass(frozen=True)
class GridState:
    n: int
    row: int
    odd_core: int
    column_index: int


@dataclass(frozen=True)
class UnitRow:
    state: GridState
    phi: int
    unit_density: Fraction
    previous_row_phi: int | None
    row_ratio: Fraction | None
    seam_defect: int | None


@dataclass(frozen=True)
class LadderRow:
    n: int
    state: GridState
    lpd: int | None
    mode: str
    skip: int | None
    arch_gap: Fraction | None
    debt: int | None
    product: Fraction | None
    first_unprotected: Fraction | None
    first_layer: int | None
    layer_hist: tuple[tuple[int, int], ...]


def v2(value: int) -> int:
    out = 0
    while value and value % 2 == 0:
        out += 1
        value //= 2
    return out


def grid_state(n: int) -> GridState:
    row = v2(n)
    odd_core = n >> row
    return GridState(
        n=n,
        row=row,
        odd_core=odd_core,
        column_index=(odd_core + 1) // 2,
    )


def state_text(state: GridState) -> str:
    return f"2^{state.row}*{state.odd_core}"


def divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def proper_divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(2, n) if n % d == 0)


def largest_proper_divisor(n: int) -> int | None:
    divs = proper_divisors(n)
    return divs[-1] if divs else None


@lru_cache(maxsize=None)
def phi(n: int) -> int:
    return sum(1 for value in range(1, n + 1) if gcd(value, n) == 1)


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
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


def ladder_rank(speeds: tuple[int, ...]) -> tuple[Fraction, int, Fraction, tuple[int, ...]]:
    report = S356.report("ladder-rank", list(speeds))
    threshold = Fraction(1, len(speeds) + 1)
    return (
        report.max_gap / threshold,
        report.boundary_witness_count,
        1 - report.forbidden_length,
        report.speeds,
    )


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
    return tuple(
        sorted(Counter(debt_layer(n, point) for point in unprotected_points(speeds)).items())
    )


@lru_cache(maxsize=None)
def summarize_speeds(speeds: tuple[int, ...]):
    return S360.summarize(list(speeds))


@lru_cache(maxsize=None)
def initial_summary(n: int):
    return summarize_speeds(tuple(range(1, n)))


def unit_row(n: int) -> UnitRow:
    state = grid_state(n)
    previous_phi: int | None = None
    row_ratio: Fraction | None = None
    seam_defect: int | None = None
    if state.row > 0:
        previous_phi = phi(n // 2)
        row_ratio = Fraction(phi(n), previous_phi)
        seam_defect = 2 * previous_phi - phi(n)
    return UnitRow(
        state=state,
        phi=phi(n),
        unit_density=Fraction(phi(n), n - 1),
        previous_row_phi=previous_phi,
        row_ratio=row_ratio,
        seam_defect=seam_defect,
    )


def best_lpd_ladder(n: int) -> LadderRow:
    state = grid_state(n)
    lpd = largest_proper_divisor(n)
    if lpd is None:
        return LadderRow(
            n=n,
            state=state,
            lpd=None,
            mode="prime/top",
            skip=None,
            arch_gap=None,
            debt=None,
            product=None,
            first_unprotected=None,
            first_layer=None,
            layer_hist=(),
        )

    candidates = []
    for skip in range(1, n):
        speeds = ladder(n, lpd, skip)
        candidates.append((ladder_rank(speeds), skip, speeds))

    _rank, skip, speeds = min(candidates, key=lambda row: row[0])
    summary = summarize_speeds(speeds)
    arch_gap = summary.max_gap / summary.threshold
    hist = layer_histogram(n, speeds)
    first_layer = debt_layer(n, summary.first_unprotected) if summary.first_unprotected else None
    row_parent = n // 2 if state.row > 0 else None
    if row_parent == lpd:
        mode = "row-parent"
    elif state.row == 0:
        mode = "odd-divisor"
    else:
        mode = "mixed-divisor"
    return LadderRow(
        n=n,
        state=state,
        lpd=lpd,
        mode=mode,
        skip=skip,
        arch_gap=arch_gap,
        debt=summary.unprotected_count,
        product=arch_gap * summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        first_layer=first_layer,
        layer_hist=hist,
    )


def print_unit_column_law() -> None:
    print("1. UNIT SKELETON ON ODD +2 CHAINS AND EVEN *2 CHAINS")
    print("=" * 86)
    print("Initial LRC segments leave unit endpoints.  Along n=2^r*b, b odd:")
    print("  r=0 -> r=1: phi(2b)=phi(b), so the first doubling has seam defect phi(b)")
    print("  r>=1:       phi(2^(r+1)b)=2*phi(2^r b)")
    print()
    print("  odd_core row n      phi(n)  phi/(n-1)  row_ratio  seam_defect")
    for odd_core in COLUMN_ODD_CORES:
        for row in range(0, 5):
            n = odd_core * (1 << row)
            item = unit_row(n)
            print(
                f"  {odd_core:>8} {row:>3} {n:>5} {item.phi:>7} "
                f"{fmt_float(item.unit_density):>10} "
                f"{fmt_frac(item.row_ratio):>9} "
                f"{'-' if item.seam_defect is None else item.seam_defect:>11}"
            )
        print()


def print_initial_segment_verification() -> None:
    print("2. INITIAL SEGMENT LRC CHECK")
    print("=" * 86)
    print(
        f"Verified for 3<=n<={INITIAL_VERIFY_MAX}: initial-segment unprotected "
        "endpoint count equals phi(n)."
    )
    mismatches = []
    for n in range(3, INITIAL_VERIFY_MAX + 1):
        summary = initial_summary(n)
        if summary.unprotected_count != phi(n):
            mismatches.append((n, summary.unprotected_count, phi(n)))
    print(f"  mismatches: {mismatches or 'none'}")
    print()
    print("  n   state      class          boundary debt  phi(n)  first witness")
    for n in (7, 9, 11, 13, 14, 15, 16, 18, 20, 24, 28, 32):
        summary = initial_summary(n)
        print(
            f"  {n:>2}  {state_text(grid_state(n)):<9} {summary.classification:<14} "
            f"{summary.unprotected_count:>13} {phi(n):>7} "
            f"{fmt_frac(summary.first_unprotected):>14}"
        )
    print()


def print_mode_crossing_table() -> None:
    print("3. TOURNAMENT MODE B VS LRC ROW-PARENT DESCENT")
    print("=" * 86)
    print("For odd n, n->n-2 stays on the odd +2 chain.  For even n, n-2")
    print("usually changes odd core; the clean LRC multiplicative descent is n/2.")
    print()
    print("  n   state      n-2 state  top-row step?  n/2 state  lpd  lpd mode")
    for n in range(9, 25):
        state = grid_state(n)
        prev = grid_state(n - 2)
        half = grid_state(n // 2) if n % 2 == 0 else None
        lpd = largest_proper_divisor(n)
        top_row_step = state.row == 0 and prev.row == 0
        lpd_mode = "-"
        if lpd is not None:
            lpd_mode = "row" if n % 2 == 0 and lpd == n // 2 else "div"
        print(
            f"  {n:>2}  {state_text(state):<9} {state_text(prev):<9} "
            f"{'yes' if top_row_step else 'no ':>13} "
            f"{state_text(half) if half else '-':<9} "
            f"{'-' if lpd is None else lpd:>4}  {lpd_mode}"
        )
    print()


def print_lpd_ladder_table() -> None:
    print("4. LARGEST-PROPER-DIVISOR LADDER PRESSURE")
    print("=" * 86)
    print(
        "For even n, lpd(n)=n/2: the first quotient ladder is exactly the "
        "doubling-row parent.  For odd composites it is a divisor channel "
        "inside the odd top row."
    )
    print()
    print("  n   state      lpd mode        skip  gap/th      debt  product  first layer")
    for n in range(9, LADDER_SCAN_MAX + 1):
        row = best_lpd_ladder(n)
        if row.lpd is None:
            if n in (11, 13, 17, 19, 23, 29, 31):
                print(
                    f"  {n:>2}  {state_text(row.state):<9} {'-':>3} "
                    f"{row.mode:<11} {'-':>4} {'-':>8} {'-':>7} {'-':>8} {'-':>6}"
                )
            continue
        if n <= 24 or n in (27, 28, 32, 36):
            print(
                f"  {n:>2}  {state_text(row.state):<9} {row.lpd:>3} "
                f"{row.mode:<11} {row.skip:>4} {fmt_frac(row.arch_gap):>8} "
                f"{row.debt:>7} {fmt_frac(row.product):>8} "
                f"{fmt_frac(row.first_unprotected):>8} {row.first_layer:>5}"
            )
    print()


def print_column_ladder_families() -> None:
    print("5. ROW FAMILIES: SAME ODD CORE, DOUBLING THE DENOMINATOR")
    print("=" * 86)
    print("This follows the user's x*2 chains directly.  The lpd ladder on")
    print("the first even row is the odd core; later rows keep doubling that parent.")
    print()
    print("  odd_core  n-chain              gap/th sequence                 debt sequence")
    for odd_core in (3, 5, 7, 9):
        ns = [odd_core * (1 << row) for row in range(1, 1 + 4)]
        ns = [n for n in ns if n <= LADDER_SCAN_MAX and largest_proper_divisor(n)]
        rows = [best_lpd_ladder(n) for n in ns]
        gap_seq = ", ".join(fmt_frac(row.arch_gap) for row in rows)
        debt_seq = ", ".join(str(row.debt) for row in rows)
        n_seq = ", ".join(str(n) for n in ns)
        print(f"  {odd_core:>8}  {n_seq:<20} {gap_seq:<30} {debt_seq}")
    print()


def print_product_sum_overlay() -> None:
    print("6. PRODUCT-SUM OVERLAY NEAR THE LRC FRONTIER")
    print("=" * 86)
    minima = S378.enumerate_product_sum_minima(40)
    print("The +2/additive chain supplies slack; multiplication supplies sparse")
    print("factor cores.  At LRC denominator n, use arity k=n-1.")
    print()
    print("  n   state      k   min product  seed")
    for n in range(12, 25):
        best = minima[n - 1]
        print(
            f"  {n:>2}  {state_text(grid_state(n)):<9} {n-1:>2} "
            f"{best.product:>11}  {best.compact()}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 86)
    print("1. The LRC unit boundary skeleton is the same object as the odd top row:")
    print("   initial segments leave exactly phi(n) unit endpoints.")
    print("2. The first row doubling b->2b is a seam: it creates quotient/nonunit")
    print("   room but no new unit witnesses.  Later doublings double the unit")
    print("   skeleton.  This mirrors the tournament pair anomaly.")
    print("3. For even n, lpd(n)=n/2.  Thus the most visible multiplicative LRC")
    print("   descent is literally the row parent in the 2-adic grid.")
    print("4. The proof strategy suggested by the grid is two-dimensional:")
    print("   prove top-row +2 movement with unit/product-sum controls, and prove")
    print("   row movement with endpoint-debt transfer from n/2 to n.")
    print("5. In this language, n=14 is the first prime-core seam 7->14; n=16 is")
    print("   the pure dyadic row lab; n=18 is the first square-odd-core seam 9->18.")


def main() -> None:
    print("LRC column/row mode atlas (codex-2026-05-31 S411)")
    print()
    print_unit_column_law()
    print_initial_segment_verification()
    print_mode_crossing_table()
    print_lpd_ladder_table()
    print_column_ladder_families()
    print_product_sum_overlay()
    print_synthesis()


if __name__ == "__main__":
    main()
