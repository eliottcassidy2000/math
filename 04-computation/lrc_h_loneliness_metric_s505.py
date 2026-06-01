#!/usr/bin/env python3
"""
lrc_h_loneliness_metric_s505.py

codex-2026-06-01 S505

Explore H(T), the Hamiltonian-path count of the half-turn phase tournament,
as a global loneliness/spread metric for Lonely Runner rows.

The point is to separate two recursions:

* H sees the coarse half-turn phase shape: bunched/transitive cells have tiny H;
  near-regular/even-spread cells have large H.
* LRC endpoint debt sees the anchored 1/n boundary clock.  Along hard row
  ladders the endpoint debt doubles while the phase tournament can freeze or
  move only through a tiny corridor.

Tournament Analysis declaration:

* pairwise observable: clockwise half-turn phase difference between runners;
* switch/gauge: orient i -> j when j lies in the clockwise open half-circle
  from i, with collision/antipodal ties broken by the fixed speed-order
  Hamiltonian path;
* tie Hamiltonian path: numerical speed order, inherited from S24/S502.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from math import factorial, gcd, log2
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
S502 = SourceFileLoader(
    "lrc_tournament_clock_overlay_s502",
    str(ROOT / "04-computation" / "lrc_tournament_clock_overlay_s502.py"),
).load_module()


ONE = Fraction(1, 1)


@dataclass(frozen=True)
class HSnapshot:
    t: Fraction
    hamiltonian_paths: int
    random_ratio: Fraction
    log2_h: float
    max_circle_gap: Fraction
    score_width: int
    cyclic_triples: int


@dataclass(frozen=True)
class LadderRow:
    label: str
    n: int
    scale: int | None
    depth: int | None
    skip: int | None
    classification: str
    gap_ratio: Fraction
    endpoint_debt: int
    product: Fraction
    selected_kind: str
    selected_t: Fraction
    origin_margin: Fraction
    phase_walls_in_gap: int
    h_mid: HSnapshot
    h_corridor_cell_count: int
    h_corridor_values: tuple[int, ...]
    h_corridor_ratio_min: Fraction
    h_corridor_ratio_max: Fraction


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return S356.fmt_frac(value)


def fmt_dec(value: Fraction, places: int = 4) -> str:
    return f"{float(value):.{places}f}"


def primitive(values: tuple[int, ...]) -> bool:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g == 1


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        raise ValueError((n, scale, skip, speeds))
    return speeds


def row_masks(adj: list[list[bool]]) -> tuple[int, ...]:
    masks: list[int] = []
    for row in adj:
        mask = 0
        for j, value in enumerate(row):
            if value:
                mask |= 1 << j
        masks.append(mask)
    return tuple(masks)


@lru_cache(maxsize=None)
def hamiltonian_paths_from_masks(n: int, masks: tuple[int, ...]) -> int:
    full = (1 << n) - 1
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[((1 << v) * n) + v] = 1

    for mask in range(1 << n):
        base = mask * n
        for last in range(n):
            value = dp[base + last]
            if not value:
                continue
            allowed = masks[last] & ~mask
            while allowed:
                bit = allowed & -allowed
                nxt = bit.bit_length() - 1
                dp[((mask | bit) * n) + nxt] += value
                allowed -= bit
    return sum(dp[full * n : (full + 1) * n])


def h_snapshot(speeds: tuple[int, ...], t: Fraction) -> HSnapshot:
    speeds0 = (0,) + speeds
    adj = S502.phase_tournament(speeds0, t)
    h_value = hamiltonian_paths_from_masks(len(speeds0), row_masks(adj))
    random_ratio = Fraction(h_value * (2 ** (len(speeds0) - 1)), factorial(len(speeds0)))
    pos = S502.positions(speeds0, t)
    return HSnapshot(
        t=t,
        hamiltonian_paths=h_value,
        random_ratio=random_ratio,
        log2_h=log2(h_value),
        max_circle_gap=S502.max_circle_gap(pos),
        score_width=S502.score_width(adj),
        cyclic_triples=S502.directed_triangles(adj),
    )


def walls_inside_arc(
    walls: tuple[Fraction, ...], lo: Fraction, hi: Fraction
) -> tuple[Fraction, ...]:
    if hi <= lo:
        return tuple()
    out: list[Fraction] = []
    for wall in walls:
        shifted = wall
        while shifted <= lo:
            shifted += ONE
        if shifted < hi:
            out.append(shifted)
    return tuple(sorted(out))


def corridor_sample_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    kind, t, gap = S502.selected_time_and_gap(speeds)
    if kind != "gap-mid" or gap is None:
        return (t,)
    walls = S502.phase_walls((0,) + speeds)
    inside = walls_inside_arc(walls, gap[0], gap[1])
    cuts = (gap[0],) + inside + (gap[1],)
    return tuple(S502.circle((a + b) / 2) for a, b in zip(cuts, cuts[1:]))


def summarize_row(
    label: str,
    n: int,
    speeds: tuple[int, ...],
    scale: int | None,
    depth: int | None,
    skip: int | None,
) -> LadderRow:
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    kind, t, gap = S502.selected_time_and_gap(speeds)
    walls_in_gap = 0
    if gap is not None:
        walls_in_gap = len(walls_inside_arc(S502.phase_walls((0,) + speeds), gap[0], gap[1]))
    h_mid = h_snapshot(speeds, t)
    corridor = tuple(h_snapshot(speeds, sample) for sample in corridor_sample_times(speeds))
    ratios = tuple(row.random_ratio for row in corridor)
    return LadderRow(
        label=label,
        n=n,
        scale=scale,
        depth=depth,
        skip=skip,
        classification=summary.classification,
        gap_ratio=report.max_gap / report.threshold,
        endpoint_debt=summary.unprotected_count,
        product=(report.max_gap / report.threshold) * summary.unprotected_count,
        selected_kind=kind,
        selected_t=t,
        origin_margin=S502.origin_margin(speeds, t),
        phase_walls_in_gap=walls_in_gap,
        h_mid=h_mid,
        h_corridor_cell_count=len(corridor),
        h_corridor_values=tuple(sorted({row.hamiltonian_paths for row in corridor})),
        h_corridor_ratio_min=min(ratios),
        h_corridor_ratio_max=max(ratios),
    )


def build_rows(max_depth: int = 3) -> tuple[LadderRow, ...]:
    hard = {
        14: (7, 6),
        18: (9, 8),
    }
    rows: list[LadderRow] = []
    for n, (base_scale, skip) in hard.items():
        initial = tuple(range(1, n))
        rows.append(summarize_row(f"n{n} initial", n, initial, None, None, None))
        for depth in range(max_depth + 1):
            scale = base_scale * (2**depth)
            speeds = ladder(n, scale, skip)
            rows.append(
                summarize_row(
                    f"n{n} scale {scale}",
                    n,
                    speeds,
                    scale,
                    depth,
                    skip,
                )
            )
    return tuple(rows)


def print_method() -> None:
    print("LRC H-as-loneliness metric exploration (codex-2026-06-01 S505)")
    print("=" * 118)
    print("H(T) here is the Hamiltonian-path count of the half-turn phase tournament.")
    print("Normalize by the random-tournament baseline n!/2^(n-1):")
    print()
    print("  H_ratio = H(T) * 2^(n-1) / n!")
    print()
    print("Large H_ratio means many compatible circular orderings: global spread/no")
    print("clean hierarchy.  Small H_ratio means bunched/transitive-side structure.")
    print("This is a global loneliness entropy, not the anchored LRC endpoint witness.")
    print()


def print_ladder_table(rows: tuple[LadderRow, ...]) -> None:
    print("H METRIC ON LRC RECURSIVE ROWS")
    print("=" * 118)
    print(
        f"{'row':<14} {'scale':>5} {'d':>2} {'class':>12} {'gap/th':>8} "
        f"{'debt':>5} {'prod':>7} {'orig/th':>8} {'walls':>5} "
        f"{'scoreW':>6} {'gapMax':>8} {'H_mid':>18} {'H_ratio':>8}"
    )
    print("-" * 118)
    for row in rows:
        snap = row.h_mid
        print(
            f"{row.label:<14} {fmt(row.scale):>5} {fmt(row.depth):>2} "
            f"{row.classification:>12} {fmt(row.gap_ratio):>8} "
            f"{row.endpoint_debt:>5} {fmt(row.product):>7} "
            f"{fmt(row.origin_margin):>8} {row.phase_walls_in_gap:>5} "
            f"{snap.score_width:>6} {fmt(snap.max_circle_gap):>8} "
            f"{snap.hamiltonian_paths:>18} {fmt_dec(snap.random_ratio):>8}"
        )
    print()


def print_corridor_table(rows: tuple[LadderRow, ...]) -> None:
    print("H PROFILE ACROSS THE LONELY CORRIDOR")
    print("=" * 118)
    print(
        f"{'row':<14} {'cells':>5} {'H values in corridor':<44} "
        f"{'ratio range':<22} {'read'}"
    )
    print("-" * 118)
    for row in rows:
        values = ",".join(str(value) for value in row.h_corridor_values)
        if len(values) > 42:
            values = values[:39] + "..."
        if len(row.h_corridor_values) == 1:
            read = "phase shape frozen"
        elif row.h_corridor_ratio_max - row.h_corridor_ratio_min < Fraction(1, 20):
            read = "small H wiggle"
        else:
            read = "phase corridor visible"
        print(
            f"{row.label:<14} {row.h_corridor_cell_count:>5} {values:<44} "
            f"{fmt_dec(row.h_corridor_ratio_min)}..{fmt_dec(row.h_corridor_ratio_max):<12} {read}"
        )
    print()


def print_recursion(rows: tuple[LadderRow, ...]) -> None:
    print("RECURSIVE READ")
    print("=" * 118)
    for n in (14, 18):
        ladder_rows = [row for row in rows if row.n == n and row.scale is not None]
        h_sequence = " -> ".join(str(row.h_mid.hamiltonian_paths) for row in ladder_rows)
        ratio_sequence = " -> ".join(fmt_dec(row.h_mid.random_ratio) for row in ladder_rows)
        product_sequence = " -> ".join(fmt(row.product) for row in ladder_rows)
        wall_sequence = " -> ".join(str(row.phase_walls_in_gap) for row in ladder_rows)
        print(f"n={n}:")
        print(f"  H_mid:      {h_sequence}")
        print(f"  H_ratio:    {ratio_sequence}")
        print(f"  gap*debt:   {product_sequence}")
        print(f"  phase walls:{wall_sequence}")
    print()
    print("Interpretation:")
    print("  The LRC endpoint recursion and H recursion are coupled but not the same.")
    print("  gap*debt is the anchored endpoint ledger; H_ratio is the half-turn phase")
    print("  entropy.  When H freezes while debt doubles, the LRC recursion has moved")
    print("  in endpoint/denominator depth without changing the coarse phase shape.")


def main() -> None:
    print_method()
    rows = build_rows(max_depth=3)
    print_ladder_table(rows)
    print_corridor_table(tuple(row for row in rows if row.scale is not None))
    print_recursion(rows)


if __name__ == "__main__":
    main()
