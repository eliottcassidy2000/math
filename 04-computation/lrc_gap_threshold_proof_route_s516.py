#!/usr/bin/env python3
"""Audit the gap-threshold proof route suggested by THM-382 and THM-383.

codex-2026-06-01 S516

THM-382 says that raw A000568/phase classes are too coarse, while threshold
decorated gap fibers certify all bounded n=3,4 systems in the S512 audit.
THM-383 says that equality walls must be compactified rather than ignored.

This script separates the part that is now essentially theorem-level from the
part that is still the Lonely Runner proof:

1. Local source-gap criterion:
   for any speed system and time, the stationary observer is lonely iff the two
   circular gaps adjacent to the observer are both at least 1/n.  This checks
   the criterion over exact wall/midpoint samples and records boundary-only
   witnesses.

2. Global forced-visit problem:
   after the local criterion, a counterexample is a compactified arithmetic
   walk avoiding all source-gap fibers.  The remaining proof must rule those
   walks out using sieve completeness and endpoint-pressure structure.

Tournament Analysis declaration:

Proof-route tournament
    vertices: candidate proof languages for the LRC/A000568 route.
    pairwise observable: seven route dimensions: exactness, boundary handling,
        compression, non-tautological content, proof payoff, compatibility with
        THM-369/THM-380, and projection risk.
    switch/gauge: route A points to route B when A wins a majority of the
        dimensions; projection risk is reversed because lower risk is better.
    tie Hamiltonian path: the listed route order.

Stored output:
    05-knowledge/results/lrc_gap_threshold_proof_route_s516.out
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
S514 = SourceFileLoader(
    "lrc_three_layer_stack_audit_s514",
    str(ROOT / "04-computation" / "lrc_three_layer_stack_audit_s514.py"),
).load_module()


Adj = list[list[int]]


def fmt_frac(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: float | Fraction | None, places: int = 3) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def frac(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def circular_distance_zero(x: Fraction) -> Fraction:
    y = frac(x)
    return min(y, 1 - y)


def positions(speeds: tuple[int, ...], time: Fraction) -> tuple[Fraction, ...]:
    return tuple(frac(Fraction(speed) * time) for speed in speeds)


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        g = 0
        for speed in moving:
            g = gcd(g, speed)
        if g == 1:
            out.append((0,) + moving)
    return out


def lrc_good(speeds: tuple[int, ...], time: Fraction) -> bool:
    threshold = Fraction(1, len(speeds))
    pos = positions(speeds, time)
    return all(circular_distance_zero(pos[idx]) >= threshold for idx in range(1, len(speeds)))


def observer_adjacent_gaps(speeds: tuple[int, ...], time: Fraction) -> tuple[Fraction, Fraction]:
    pos = positions(speeds, time)
    # Observer has index 0 and is ordered first among exact ties.  If any moving
    # runner is at the observer, one adjacent gap is therefore exactly zero.
    order = sorted(range(len(speeds)), key=lambda idx: (pos[idx], idx))
    k = order.index(0)
    prev_idx = order[(k - 1) % len(order)]
    next_idx = order[(k + 1) % len(order)]
    left_gap = frac(pos[0] - pos[prev_idx])
    right_gap = frac(pos[next_idx] - pos[0])
    return left_gap, right_gap


def source_gap_good(speeds: tuple[int, ...], time: Fraction) -> bool:
    threshold = Fraction(1, len(speeds))
    left_gap, right_gap = observer_adjacent_gaps(speeds, time)
    return left_gap >= threshold and right_gap >= threshold


def endpoint_walls(speeds: tuple[int, ...]) -> list[Fraction]:
    total_n = len(speeds)
    out: set[Fraction] = {Fraction(0), Fraction(1)}
    for speed in speeds[1:]:
        for m in range(speed):
            out.add(Fraction(m * total_n + 1, total_n * speed))
            out.add(Fraction(m * total_n + total_n - 1, total_n * speed))
    return sorted(time for time in out if 0 <= time <= 1)


def collision_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = set()
    for i, j in combinations(range(len(speeds)), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for m in range(d + 1):
            out.add(Fraction(m, d))
    return out


def pair_threshold_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    total_n = len(speeds)
    out: set[Fraction] = set()
    for i, j in combinations(range(len(speeds)), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for m in range(d):
            out.add(Fraction(m * total_n + 1, total_n * d))
            out.add(Fraction(m * total_n + total_n - 1, total_n * d))
    return {time for time in out if 0 <= time <= 1}


def compactified_gap_walls(speeds: tuple[int, ...]) -> list[Fraction]:
    out = set(endpoint_walls(speeds))
    out.update(collision_walls(speeds))
    out.update(pair_threshold_walls(speeds))
    return sorted(time for time in out if 0 <= time <= 1)


def sampled_times(walls: list[Fraction]) -> list[Fraction]:
    out: set[Fraction] = set(walls)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            out.add((left + right) / 2)
    return sorted(out)


def witness_profile(speeds: tuple[int, ...], walls: list[Fraction]) -> tuple[str, int, int, Fraction]:
    wall_set = set(walls)
    open_hits = 0
    wall_hits = 0
    best_margin: Fraction | None = None
    threshold = Fraction(1, len(speeds))
    for time in sampled_times(walls):
        margin = min(circular_distance_zero(pos) for pos in positions(speeds, time)[1:]) - threshold
        if best_margin is None or margin > best_margin:
            best_margin = margin
        if lrc_good(speeds, time):
            if time in wall_set:
                wall_hits += 1
            else:
                open_hits += 1
    assert best_margin is not None
    if open_hits and wall_hits:
        kind = "open+wall"
    elif open_hits:
        kind = "open"
    elif wall_hits:
        kind = "wall-only"
    else:
        kind = "none"
    return kind, open_hits, wall_hits, best_margin


def open_corridors(speeds: tuple[int, ...], walls: list[Fraction]) -> list[tuple[Fraction, Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction, Fraction]] = []
    for left, right in zip(walls, walls[1:]):
        if left >= right:
            continue
        mid = (left + right) / 2
        if lrc_good(speeds, mid):
            out.append((left, right, right - left))
    return out


def sieve_missing_denominators(speeds: tuple[int, ...]) -> tuple[int, ...]:
    total_n = len(speeds)
    moving = speeds[1:]
    return tuple(q for q in range(2, total_n + 1) if all(speed % q != 0 for speed in moving))


@dataclass(frozen=True)
class BoundedAudit:
    total_n: int
    max_speed: int
    systems: int
    sample_states: int
    criterion_mismatches: int
    fiber_status: Counter[tuple[int, int]]
    profiles: Counter[str]
    sieve_explained: int
    all_have_source_gap: int


def bounded_audit(total_n: int, max_speed: int) -> BoundedAudit:
    systems = primitive_speed_sets(total_n, max_speed)
    fiber_status: Counter[tuple[int, int]] = Counter()
    profiles: Counter[str] = Counter()
    criterion_mismatches = 0
    sample_states = 0
    sieve_explained = 0
    all_have_source_gap = 0
    for speeds in systems:
        walls = compactified_gap_walls(speeds)
        has_good = False
        for time in sampled_times(walls):
            sample_states += 1
            direct = lrc_good(speeds, time)
            source_gap = source_gap_good(speeds, time)
            if direct != source_gap:
                criterion_mismatches += 1
            left_gap, right_gap = observer_adjacent_gaps(speeds, time)
            threshold = Fraction(1, total_n)
            fiber = (int(left_gap >= threshold), int(right_gap >= threshold))
            fiber_status[(fiber, int(direct))] += 1
            has_good = has_good or direct
        kind, _, _, _ = witness_profile(speeds, walls)
        profiles[kind] += 1
        if sieve_missing_denominators(speeds):
            sieve_explained += 1
        if has_good:
            all_have_source_gap += 1
    return BoundedAudit(
        total_n=total_n,
        max_speed=max_speed,
        systems=len(systems),
        sample_states=sample_states,
        criterion_mismatches=criterion_mismatches,
        fiber_status=fiber_status,
        profiles=profiles,
        sieve_explained=sieve_explained,
        all_have_source_gap=all_have_source_gap,
    )


def print_bounded_audits() -> None:
    print("Part A. Exact compactified source-gap criterion audit")
    print("-" * 88)
    print(
        "total n  max speed  systems  states  mismatches  source-hit  "
        "sieve-explained  witness profiles"
    )
    for total_n, max_speed in ((3, 16), (4, 10), (5, 9), (6, 8)):
        audit = bounded_audit(total_n, max_speed)
        profile = ", ".join(f"{key}:{audit.profiles[key]}" for key in sorted(audit.profiles))
        print(
            f"{total_n:>7}  {max_speed:>9}  {audit.systems:>7}  "
            f"{audit.sample_states:>6}  {audit.criterion_mismatches:>10}  "
            f"{audit.all_have_source_gap:>10}/{audit.systems:<5}  "
            f"{audit.sieve_explained:>15}  {profile}"
        )
        good_fibers = sorted(
            fiber for (fiber, status), count in audit.fiber_status.items() if status and count
        )
        bad_fibers = sorted(
            fiber for (fiber, status), count in audit.fiber_status.items() if not status and count
        )
        print(f"         good fibers={good_fibers} bad fibers={bad_fibers}")
    print()


def hard_row_profiles() -> None:
    print("Part B. Selected n14/n16/n18 rows through the source-gap lens")
    print("-" * 88)
    print(
        "row                                  N speeds max  open corridors  "
        "wall hits  best margin  widest corridor"
    )
    for label, total_n, moving in S514.selected_rows():
        speeds = (0,) + moving
        walls = endpoint_walls(speeds)
        corridors = open_corridors(speeds, walls)
        kind, open_hits, wall_hits, best_margin = witness_profile(speeds, walls)
        widest = max((width for _, _, width in corridors), default=Fraction(0))
        print(
            f"{label:<36} {total_n:>2} {max(moving):>10}  "
            f"{len(corridors):>14}  {wall_hits:>9}  "
            f"{fmt_frac(best_margin):>11}  {fmt_frac(widest):>15}  {kind}"
        )
    print()


@dataclass(frozen=True)
class Route:
    name: str
    exactness: int
    boundary: int
    compression: int
    non_taut: int
    payoff: int
    sieve_pressure: int
    projection_risk: int


ROUTES = (
    Route("raw_A000568_phase", 2, 1, 5, 5, 1, 1, 5),
    Route("observer_source_marked", 5, 3, 4, 3, 3, 2, 2),
    Route("threshold_gap_source", 5, 5, 5, 2, 3, 3, 1),
    Route("compactified_gap_walk", 4, 5, 4, 4, 5, 4, 2),
    Route("endpoint_pressure_core", 5, 4, 3, 5, 4, 5, 2),
    Route("add_mult_row_stack", 3, 3, 3, 5, 3, 4, 3),
)


DIMENSIONS: tuple[tuple[str, Callable[[Route], int], bool], ...] = (
    ("exactness", lambda route: route.exactness, True),
    ("boundary", lambda route: route.boundary, True),
    ("compression", lambda route: route.compression, True),
    ("non_taut", lambda route: route.non_taut, True),
    ("payoff", lambda route: route.payoff, True),
    ("sieve_pressure", lambda route: route.sieve_pressure, True),
    ("projection_risk", lambda route: route.projection_risk, False),
)


def route_tournament() -> Adj:
    adj = [[0] * len(ROUTES) for _ in ROUTES]
    for i, j in combinations(range(len(ROUTES)), 2):
        wins_i = 0
        wins_j = 0
        for _, getter, high_is_good in DIMENSIONS:
            vi = getter(ROUTES[i])
            vj = getter(ROUTES[j])
            if vi == vj:
                continue
            if (vi > vj) == high_is_good:
                wins_i += 1
            else:
                wins_j += 1
        winner = i if wins_i >= wins_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def hamiltonian_paths(adj: Adj) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for start in range(n):
        dp[1 << start][start] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def triangle_count(adj: Adj) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: Adj) -> tuple[int, ...]:
    n = len(adj)
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w, bit in enumerate(adj[v]):
            if bit and not seen[w]:
                dfs(w)
        order.append(v)

    for vertex in range(n):
        if not seen[vertex]:
            dfs(vertex)

    reverse = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            reverse[j][i] = adj[i][j]

    seen = [False] * n
    sizes: list[int] = []
    for start in reversed(order):
        if seen[start]:
            continue
        queue = deque([start])
        seen[start] = True
        size = 0
        while queue:
            v = queue.popleft()
            size += 1
            for w, bit in enumerate(reverse[v]):
                if bit and not seen[w]:
                    seen[w] = True
                    queue.append(w)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def print_route_tournament() -> None:
    print("Part C. Proof-route Tournament Analysis")
    print("-" * 88)
    adj = route_tournament()
    scores = [sum(row) for row in adj]
    score_hist = Counter(scores)
    print(
        f"routes={len(ROUTES)} H={hamiltonian_paths(adj)} c3={triangle_count(adj)} "
        f"SCCs={scc_sizes(adj)} score_hist={dict(sorted(score_hist.items()))}"
    )
    ranked = sorted(((scores[idx], -idx, ROUTES[idx]) for idx in range(len(ROUTES))), reverse=True)
    print("ranked routes:")
    for score, _, route in ranked:
        print(
            f"  {score:>2}  {route.name:<26} "
            f"exact={route.exactness} boundary={route.boundary} "
            f"compress={route.compression} non_taut={route.non_taut} "
            f"payoff={route.payoff} sieve/pressure={route.sieve_pressure} "
            f"risk={route.projection_risk}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 88)
    print(
        "1. The local THM-382 fiber phenomenon has collapsed to a proof-level "
        "criterion: LRC-good status is exactly the closed threshold color of "
        "the two gaps adjacent to the observer."
    )
    print(
        "2. THM-383 remains essential because equality witnesses are common; "
        "the proof object is a compactified walk with wall points, not merely "
        "an open-cell circular menu."
    )
    print(
        "3. This does not prove LRC.  It moves the global burden to a sharper "
        "statement: every admissible compactified arithmetic walk must visit "
        "the source-gap fiber (1,1)."
    )
    print(
        "4. By THM-369 and THM-380, a source-gap-avoiding counterexample must "
        "also be sieve-complete and carry a nonempty owner-compatible pressure "
        "core with a directed cycle.  The next proof route is to show that the "
        "compactified source-gap walk forbids that conjunction."
    )


def main() -> None:
    print("LRC gap-threshold proof-route audit (S516)")
    print("=" * 88)
    print_bounded_audits()
    hard_row_profiles()
    print_route_tournament()
    print_synthesis()


if __name__ == "__main__":
    main()
