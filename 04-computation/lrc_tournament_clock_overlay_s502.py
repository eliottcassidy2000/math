#!/usr/bin/env python3
"""Overlay LRC witness intervals with the half-turn tournament clock.

The S24 tournament clock records pairwise half-turn wall crossings

    t = m / (2 |s_i - s_j|)

for the runner positions, including the observer as speed 0.  LRC, by contrast,
uses anchored forbidden-boundary times

    ||v t|| = 1/n.

This script overlays the two exact wall systems for the n=14 and n=18 hard
rows.  The goal is to see whether a lonely interval is visible to the half-turn
clock, or whether it lies inside one clock cell and therefore requires a finer
threshold/endpoint gauge.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
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

ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class Row:
    label: str
    n: int
    speeds: tuple[int, ...]
    mode: str
    scale: int | None = None
    skip: int | None = None


@dataclass(frozen=True)
class Overlay:
    label: str
    n: int
    mode: str
    scale: int | None
    skip: int | None
    classification: str
    gap_ratio: Fraction
    unprotected: int
    selected_kind: str
    t: Fraction
    origin_margin: Fraction
    clock_on_wall: bool
    clock_cell_width: Fraction
    lonely_gap_width: Fraction
    phase_walls_in_lonely_gap: int
    lonely_boundaries: int
    lonely_boundaries_on_phase_wall: int
    max_circle_gap: Fraction
    score_width: int
    cyclic_triples: int
    largest_scc: int
    hamiltonian_paths: int | None


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return S356.fmt_frac(value)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def circular_distance_to_zero(value: Fraction) -> Fraction:
    value = circle(value)
    return min(value, ONE - value)


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = circle(b - a)
    return min(delta, ONE - delta)


def positions(speeds_with_observer: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in speeds_with_observer)


def phase_walls(speeds_with_observer: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls: set[Fraction] = set()
    for a, b in combinations(speeds_with_observer, 2):
        d = abs(a - b)
        if d == 0:
            continue
        for m in range(2 * d):
            walls.add(Fraction(m, 2 * d))
    return tuple(sorted(walls))


def count_walls_in_arc(walls: tuple[Fraction, ...], lo: Fraction, hi: Fraction) -> int:
    """Count clock walls in the open circular arc (lo, hi).  hi may exceed 1."""
    if hi <= lo:
        return 0
    if hi - lo >= ONE:
        return len(walls)
    lo0 = circle(lo)
    hi0 = lo0 + (hi - lo)
    if hi0 <= ONE:
        return bisect_left(walls, hi0) - bisect_right(walls, lo0)
    return (
        len(walls)
        - bisect_right(walls, lo0)
        + bisect_left(walls, hi0 - ONE)
    )


def clock_cell(walls: tuple[Fraction, ...], t: Fraction) -> tuple[bool, Fraction, Fraction]:
    """Return (on_wall, left_width, right_width) around t."""
    t = circle(t)
    idx = bisect_left(walls, t)
    on_wall = idx < len(walls) and walls[idx] == t
    if on_wall:
        prev_wall = walls[idx - 1] if idx > 0 else walls[-1] - ONE
        next_wall = walls[(idx + 1) % len(walls)]
        if next_wall <= t:
            next_wall += ONE
        return True, t - prev_wall, next_wall - t

    prev_wall = walls[idx - 1] if idx > 0 else walls[-1] - ONE
    next_wall = walls[idx] if idx < len(walls) else walls[0] + ONE
    return False, t - prev_wall, next_wall - t


def lonely_gaps(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    raw = S356.forbidden_intervals(speeds)
    return tuple(S356.circular_gaps(S356.merge_intervals(raw)))


def selected_time_and_gap(speeds: tuple[int, ...]) -> tuple[str, Fraction, tuple[Fraction, Fraction] | None]:
    gaps = lonely_gaps(speeds)
    if gaps:
        gap = max(gaps, key=lambda item: item[1] - item[0])
        return "gap-mid", S356.midpoint_on_circle(gap), gap
    summary = S360.summarize(list(speeds))
    if summary.first_unprotected is not None:
        return "boundary", summary.first_unprotected, None
    return "unit", Fraction(1, len(speeds) + 1), None


def lonely_boundary_witnesses(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    raw = S356.forbidden_intervals(speeds)
    candidates = sorted({S356.circle_point(point) for interval in raw for point in interval})
    return tuple(t for t in candidates if S356.is_lonely_witness(speeds, t))


def phase_tournament(speeds_with_observer: tuple[int, ...], t: Fraction) -> list[list[bool]]:
    pos = positions(speeds_with_observer, t)
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = circle(pos[j] - pos[i])
        if delta == 0 or delta == HALF:
            winner = i
        elif delta < HALF:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def score_width(adj: list[list[bool]]) -> int:
    scores = [sum(row) for row in adj]
    return max(scores) - min(scores)


def score_hist(adj: list[list[bool]]) -> str:
    counts = Counter(sum(row) for row in adj)
    return " ".join(f"{score}:{counts[score]}" for score in sorted(counts))


def directed_triangles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> tuple[int, ...]:
    n = len(adj)
    graph = [[] for _ in range(n)]
    reverse = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                graph[i].append(j)
                reverse[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in graph[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    sizes: list[int] = []
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


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
def hamiltonian_paths_cached(n: int, masks: tuple[int, ...]) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            allowed = masks[last] & ~mask
            while allowed:
                bit = allowed & -allowed
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += value
                allowed -= bit
    return sum(dp[-1])


def hamiltonian_paths(adj: list[list[bool]]) -> int | None:
    n = len(adj)
    if n > 14:
        return None
    return hamiltonian_paths_cached(n, row_masks(adj))


def max_circle_gap(pos: tuple[Fraction, ...]) -> Fraction:
    ordered = sorted(pos)
    gaps = [ordered[i + 1] - ordered[i] for i in range(len(ordered) - 1)]
    gaps.append(ordered[0] + ONE - ordered[-1])
    return max(gaps)


def origin_margin(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    threshold = Fraction(1, len(speeds) + 1)
    return min(circular_distance_to_zero(speed * t) for speed in speeds) / threshold


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


def rows() -> tuple[Row, ...]:
    # Established by the S481/S490/S493 hard-row audits.
    hard_skip = {14: 6, 18: 8}
    out: list[Row] = []
    for n in (14, 18):
        out.append(Row(f"n{n} initial", n, tuple(range(1, n)), "initial"))
        for mode, scale in (
            ("row-parent", n // 2),
            ("gate", n),
            ("double-gate", 2 * n),
        ):
            skip = hard_skip[n]
            out.append(Row(f"n{n} {mode}", n, ladder(n, scale, skip), mode, scale, skip))
    return tuple(out)


def overlay(row: Row) -> Overlay:
    report = S356.report(row.label, list(row.speeds))
    summary = S360.summarize(list(row.speeds))
    kind, t, lonely_gap = selected_time_and_gap(row.speeds)
    speeds0 = (0,) + row.speeds
    walls = phase_walls(speeds0)
    on_wall, left, right = clock_cell(walls, t)
    boundaries = lonely_boundary_witnesses(row.speeds)
    boundary_wall_count = sum(1 for b in boundaries if b in walls)
    gap_width = Fraction(0)
    walls_inside_gap = 0
    if lonely_gap is not None:
        gap_width = lonely_gap[1] - lonely_gap[0]
        walls_inside_gap = count_walls_in_arc(walls, lonely_gap[0], lonely_gap[1])

    adj = phase_tournament(speeds0, t)
    pos = positions(speeds0, t)
    return Overlay(
        label=row.label,
        n=row.n,
        mode=row.mode,
        scale=row.scale,
        skip=row.skip,
        classification=summary.classification,
        gap_ratio=report.max_gap / report.threshold,
        unprotected=summary.unprotected_count,
        selected_kind=kind,
        t=t,
        origin_margin=origin_margin(row.speeds, t),
        clock_on_wall=on_wall,
        clock_cell_width=min(left, right) if on_wall else left + right,
        lonely_gap_width=gap_width,
        phase_walls_in_lonely_gap=walls_inside_gap,
        lonely_boundaries=len(boundaries),
        lonely_boundaries_on_phase_wall=boundary_wall_count,
        max_circle_gap=max_circle_gap(pos),
        score_width=score_width(adj),
        cyclic_triples=directed_triangles(adj),
        largest_scc=scc_sizes(adj)[0],
        hamiltonian_paths=hamiltonian_paths(adj),
    )


def print_overlay_table(items: tuple[Overlay, ...]) -> None:
    print("LRC WITNESS TIMES OVERLAID ON THE HALF-TURN TOURNAMENT CLOCK")
    print("=" * 132)
    print(
        f"{'row':<18} {'mode':<11} {'scale':>5} {'skip':>4} {'class':>12} "
        f"{'gap/th':>9} {'debt':>6} {'t-kind':<8} {'t':>13} {'orig/th':>8} "
        f"{'clock':>6} {'cell/wall':>10} {'LRCgap':>10} {'walls':>5} "
        f"{'bdy wall':>9} {'H':>12}"
    )
    print("-" * 132)
    for item in items:
        clock_state = "wall" if item.clock_on_wall else "cell"
        bdy = f"{item.lonely_boundaries_on_phase_wall}/{item.lonely_boundaries}"
        print(
            f"{item.label:<18} {item.mode:<11} {fmt(item.scale):>5} {fmt(item.skip):>4} "
            f"{item.classification:>12} {fmt(item.gap_ratio):>9} "
            f"{item.unprotected:>6} {item.selected_kind:<8} {fmt(item.t):>13} "
            f"{fmt(item.origin_margin):>8} {clock_state:>6} "
            f"{fmt(item.clock_cell_width):>10} {fmt(item.lonely_gap_width):>10} "
            f"{item.phase_walls_in_lonely_gap:>5} {bdy:>9} "
            f"{fmt(item.hamiltonian_paths):>12}"
        )
    print()


def print_tournament_table(items: tuple[Overlay, ...]) -> None:
    print("HALF-TURN TOURNAMENT SHAPE AT THE LRC SELECTED TIME")
    print("=" * 104)
    print(
        f"{'row':<18} {'t-kind':<8} {'max circle gap':>14} {'score width':>12} "
        f"{'cyc3':>6} {'largest SCC':>11} {'read'}"
    )
    print("-" * 104)
    for item in items:
        if item.max_circle_gap > HALF:
            read = "bunched/transitive-side"
        elif item.score_width <= 1:
            read = "near-regular clock top"
        elif item.clock_on_wall:
            read = "clock-wall boundary"
        else:
            read = "interior circular cell"
        print(
            f"{item.label:<18} {item.selected_kind:<8} {fmt(item.max_circle_gap):>14} "
            f"{item.score_width:>12} {item.cyclic_triples:>6} "
            f"{item.largest_scc:>11} {read}"
        )
    print()


def print_synthesis(items: tuple[Overlay, ...]) -> None:
    positive = [item for item in items if item.lonely_gap_width > 0]
    crossing = [item for item in positive if item.phase_walls_in_lonely_gap > 0]
    initial = [item for item in items if item.mode == "initial"]
    print("SYNTHESIS")
    print("=" * 104)
    print(
        f"Positive-gap rows whose LRC gap crosses at least one half-turn wall: "
        f"{len(crossing)}/{len(positive)}"
    )
    for item in crossing:
        print(
            f"  {item.label}: LRC gap {fmt(item.lonely_gap_width)} crosses "
            f"{item.phase_walls_in_lonely_gap} wall(s), so it is a corridor through "
            f"{item.phase_walls_in_lonely_gap + 1} adjacent clock cells; midpoint "
            f"t={fmt(item.t)} is a {'wall' if item.clock_on_wall else 'cell interior'}."
        )
    print()
    print("Initial rows:")
    for item in initial:
        print(
            f"  {item.label}: {item.lonely_boundaries_on_phase_wall}/"
            f"{item.lonely_boundaries} lonely boundary witnesses are half-turn "
            f"clock walls; selected t={fmt(item.t)} has score width {item.score_width}."
        )
    print()
    print("Boundary alignment ratios on the hard ladders:")
    for item in positive:
        ratio = Fraction(item.lonely_boundaries_on_phase_wall, item.lonely_boundaries)
        print(
            f"  {item.label}: {item.lonely_boundaries_on_phase_wall}/"
            f"{item.lonely_boundaries} = {fmt(ratio)} LRC boundary witnesses "
            f"are also half-turn clock walls."
        )
    print()
    print(
        "Interpretation: the half-turn clock sees global circular spread.  The "
        "initial-segment unit witnesses are exactly clock-wall events.  The hard "
        "quotient-ladder lonely intervals are different: they are tiny corridors "
        "through a few adjacent circular-tournament cells, with fixed alignment "
        "ratios along the row-parent/gate/double-gate ladder.  LRC proof data "
        "should therefore keep two clocks: the coarse half-turn walk for "
        "spread/topology, and an anchored 1/n endpoint clock for the actual "
        "lonely boundary certificate."
    )


def main() -> None:
    print("LRC tournament-clock overlay (codex-2026-06-01 S502b)")
    print()
    items = tuple(overlay(row) for row in rows())
    print_overlay_table(items)
    print_tournament_table(items)
    print_synthesis(items)


if __name__ == "__main__":
    main()
