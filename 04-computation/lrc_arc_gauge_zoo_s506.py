#!/usr/bin/env python3
"""Try many LRC tournament arc gauges and score their loneliness signal.

codex-2026-06-01 S506c

Tournament Analysis declaration:

    observable data -> pairwise switch/arc criterion -> tournament -> invariant

For a Lonely Runner time slice, this script tries several ways to orient every
pair of runners (usually including the stationary observer).  The goal is not
to find one magic tournament immediately.  It is to separate useful loneliness
gauges from inert rankers:

* scalar rank gauges usually have H=1 and are readable but too compressed;
* half-turn phase H is a global spread/class meter;
* close-pair and blocker-relief switches keep edge-local crowding data;
* endpoint-protection gauges keep the finite LRC boundary object.

Stored output:
    05-knowledge/results/lrc_arc_gauge_zoo_s506.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import factorial, gcd, sqrt
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
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class SpeedRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    mode: str


@dataclass(frozen=True)
class GaugeMetrics:
    gauge: str
    vertices: int
    hamiltonian_paths: int | None
    h_ratio: Fraction | None
    score_width: int
    score_hist: str
    cyclic_triples: int
    largest_scc: int
    source_count: int
    sink_count: int
    observer_out: int | None


@dataclass(frozen=True)
class RowGaugeRecord:
    row: SpeedRow
    t_kind: str
    t: Fraction
    origin_margin: Fraction
    safe_vertices: int
    gap_ratio: Fraction
    endpoint_debt: int
    gauge: GaugeMetrics


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return S356.fmt_frac(value)


def fmt_dec(value: Fraction | float | None, places: int = 3) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def circle(value: Fraction) -> Fraction:
    return value % ONE


def dist0(value: Fraction) -> Fraction:
    value = circle(value)
    return min(value, ONE - value)


def clockwise(a: Fraction, b: Fraction) -> Fraction:
    return circle(b - a)


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    return min(clockwise(a, b), clockwise(b, a))


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


def rows() -> tuple[SpeedRow, ...]:
    return (
        SpeedRow("n14 initial", 14, tuple(range(1, 14)), "initial"),
        SpeedRow("n14 row-parent", 14, ladder(14, 7, 6), "row-parent"),
        SpeedRow("n14 gate", 14, ladder(14, 14, 6), "gate"),
        SpeedRow("n14 double-gate", 14, ladder(14, 28, 6), "double-gate"),
        SpeedRow("n18 initial", 18, tuple(range(1, 18)), "initial"),
        SpeedRow("n18 row-parent", 18, ladder(18, 9, 8), "row-parent"),
        SpeedRow("n18 gate", 18, ladder(18, 18, 8), "gate"),
        SpeedRow("n18 double-gate", 18, ladder(18, 36, 8), "double-gate"),
    )


def positions(speeds0: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in speeds0)


def adjacent_gap_data(pos: tuple[Fraction, ...]) -> tuple[list[Fraction], list[Fraction], list[Fraction]]:
    ordered = sorted(range(len(pos)), key=lambda i: (pos[i], i))
    index = {vertex: idx for idx, vertex in enumerate(ordered)}
    gaps: list[Fraction] = []
    for idx, vertex in enumerate(ordered):
        nxt = ordered[(idx + 1) % len(ordered)]
        gaps.append(clockwise(pos[vertex], pos[nxt]))

    left = [Fraction(0)] * len(pos)
    right = [Fraction(0)] * len(pos)
    moat = [Fraction(0)] * len(pos)
    for vertex in range(len(pos)):
        idx = index[vertex]
        left_gap = gaps[(idx - 1) % len(gaps)]
        right_gap = gaps[idx]
        left[vertex] = left_gap
        right[vertex] = right_gap
        moat[vertex] = min(left_gap, right_gap)
    return left, right, moat


def safe_vertex_count(pos: tuple[Fraction, ...], threshold: Fraction) -> int:
    left, right, _moat = adjacent_gap_data(pos)
    return sum(1 for lft, rgt in zip(left, right) if lft >= threshold and rgt >= threshold)


def origin_margin(speeds: tuple[int, ...], t: Fraction, n: int) -> Fraction:
    threshold = Fraction(1, n)
    return min(dist0(speed * t) for speed in speeds) / threshold


def nearest_distances(pos: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    out: list[Fraction] = []
    for i in range(len(pos)):
        out.append(min(circular_distance(pos[i], pos[j]) for j in range(len(pos)) if j != i))
    return tuple(out)


def nearest_without(pos: tuple[Fraction, ...], i: int, removed: int) -> Fraction:
    return min(
        circular_distance(pos[i], pos[j])
        for j in range(len(pos))
        if j != i and j != removed
    )


def distance_matrix(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(0) if i == j else circular_distance(pos[i], pos[j]) for j in range(len(pos)))
        for i in range(len(pos))
    )


def k_nearest_sum(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    k: int,
    removed: int | None = None,
) -> Fraction:
    values = sorted(
        dist[i][j]
        for j in range(len(dist))
        if j != i and j != removed
    )
    return sum(values[:k], Fraction(0))


def first_distances(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    k: int,
    removed: int | None = None,
) -> tuple[Fraction, ...]:
    values = sorted(
        dist[i][j]
        for j in range(len(dist))
        if j != i and j != removed
    )
    return tuple(values[:k])


def deficit_score(values: tuple[Fraction, ...], threshold: Fraction) -> Fraction:
    return sum((threshold - value for value in values if value < threshold), Fraction(0))


def orient_by_values(values: tuple[Fraction, ...], low_points_high: bool) -> list[list[bool]]:
    n = len(values)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if values[i] == values[j]:
            winner = min(i, j)
        elif low_points_high:
            winner = i if values[i] < values[j] else j
        else:
            winner = i if values[i] > values[j] else j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def phase_from_positions(pos: tuple[Fraction, ...]) -> list[list[bool]]:
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        if delta == 0 or delta == HALF:
            winner = min(i, j)
        elif delta < HALF:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def close_pair_switch(pos: tuple[Fraction, ...], threshold: Fraction) -> list[list[bool]]:
    """Flip the base path exactly on close pairs."""
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        close = circular_distance(pos[i], pos[j]) < threshold
        winner = j if close else i
        loser = i if close else j
        adj[winner][loser] = True
    return adj


def band_pair_switch(pos: tuple[Fraction, ...], threshold: Fraction) -> list[list[bool]]:
    """Flip base path on the two-neighbor danger annulus [1/n,2/n)."""
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        distance = circular_distance(pos[i], pos[j])
        in_band = threshold <= distance < 2 * threshold
        winner = j if in_band else i
        loser = i if in_band else j
        adj[winner][loser] = True
    return adj


def quarter_lead_switch(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Orient only decisive quarter-turn leads; otherwise use base path."""
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        if Fraction(0) < delta < Fraction(1, 4):
            winner = i
        elif Fraction(3, 4) < delta < ONE:
            winner = j
        else:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def blocker_relief(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Orient blocker -> blocked using nearest-neighbor deletion relief."""
    n = len(pos)
    nearest = nearest_distances(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        relief_i_if_j_removed = nearest_without(pos, i, j) - nearest[i]
        relief_j_if_i_removed = nearest_without(pos, j, i) - nearest[j]
        if relief_i_if_j_removed == relief_j_if_i_removed:
            winner = min(i, j)
        elif relief_i_if_j_removed > relief_j_if_i_removed:
            winner = j
        else:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def two_neighbor_relief(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Orient blocker -> blocked using sum-of-two-nearest deletion relief."""
    n = len(pos)
    dist = distance_matrix(pos)
    base = [k_nearest_sum(dist, i, 2) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        relief_i_if_j_removed = k_nearest_sum(dist, i, 2, removed=j) - base[i]
        relief_j_if_i_removed = k_nearest_sum(dist, j, 2, removed=i) - base[j]
        if relief_i_if_j_removed == relief_j_if_i_removed:
            winner = min(i, j)
        elif relief_i_if_j_removed > relief_j_if_i_removed:
            winner = j
        else:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def threshold_deficit_relief(pos: tuple[Fraction, ...], threshold: Fraction) -> list[list[bool]]:
    """Orient blocker -> blocked by how much deletion removes two-neighbor deficit."""
    n = len(pos)
    dist = distance_matrix(pos)
    base = [deficit_score(first_distances(dist, i, 2), threshold) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        after_i = deficit_score(first_distances(dist, i, 2, removed=j), threshold)
        after_j = deficit_score(first_distances(dist, j, 2, removed=i), threshold)
        relief_i_if_j_removed = base[i] - after_i
        relief_j_if_i_removed = base[j] - after_j
        if relief_i_if_j_removed == relief_j_if_i_removed:
            winner = min(i, j)
        elif relief_i_if_j_removed > relief_j_if_i_removed:
            winner = j
        else:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def open_arc_density_switch(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Orient toward the less crowded clockwise arc between each pair."""
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        if delta == 0:
            winner = min(i, j)
        else:
            inside_cw = sum(
                1
                for k in range(n)
                if k not in (i, j) and Fraction(0) < clockwise(pos[i], pos[k]) < delta
            )
            inside_ccw = n - 2 - inside_cw
            density_cw = Fraction(inside_cw, 1) / delta
            density_ccw = Fraction(inside_ccw, 1) / (ONE - delta)
            if density_cw == density_ccw:
                winner = i if delta < HALF else j
            elif density_cw < density_ccw:
                winner = i
            else:
                winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def open_arc_count_switch(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Combinatorial version of open-arc density: fewer intervening vertices wins."""
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        if delta == 0:
            winner = min(i, j)
        else:
            inside_cw = sum(
                1
                for k in range(n)
                if k not in (i, j) and Fraction(0) < clockwise(pos[i], pos[k]) < delta
            )
            inside_ccw = n - 2 - inside_cw
            if inside_cw == inside_ccw:
                winner = i if delta < HALF else j
            elif inside_cw < inside_ccw:
                winner = i
            else:
                winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def kinetic_escape(row: SpeedRow, pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Orient along the short arc that is currently opening under relative speed."""
    speeds0 = (0,) + row.speeds
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        dv = speeds0[j] - speeds0[i]
        if delta == 0 or delta == HALF or dv == 0:
            winner = min(i, j)
        elif delta < HALF:
            winner = i if dv > 0 else j
        else:
            winner = j if dv < 0 else i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def kinetic_close_switch(row: SpeedRow, pos: tuple[Fraction, ...], threshold: Fraction) -> list[list[bool]]:
    """Flip the base path on pairs that are within 2/n and still closing."""
    speeds0 = (0,) + row.speeds
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = clockwise(pos[i], pos[j])
        dv = speeds0[j] - speeds0[i]
        distance = min(delta, ONE - delta)
        rate = dv if delta < HALF else -dv
        active = distance < 2 * threshold and rate < 0
        winner = j if active else i
        loser = i if active else j
        adj[winner][loser] = True
    return adj


def observer_two_guard(pos: tuple[Fraction, ...], threshold: Fraction) -> list[list[bool]]:
    """Orient runners by two-sided observer relief; observer edges mark danger."""
    n = len(pos)
    dist = distance_matrix(pos)
    base = k_nearest_sum(dist, 0, 2)
    relief = [Fraction(0)] * n
    for runner in range(1, n):
        relief[runner] = k_nearest_sum(dist, 0, 2, removed=runner) - base

    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if i == 0 or j == 0:
            runner = j if i == 0 else i
            winner = runner if dist[0][runner] < threshold else 0
        elif relief[i] == relief[j]:
            winner = min(i, j)
        elif relief[i] > relief[j]:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def observer_relief_rank(pos: tuple[Fraction, ...], _threshold: Fraction) -> list[list[bool]]:
    """Rank runners by how much deleting them improves observer nearest distance."""
    current = min(circular_distance(pos[0], pos[j]) for j in range(1, len(pos)))
    values: list[Fraction] = [Fraction(-1)]
    for removed in range(1, len(pos)):
        after = min(
            circular_distance(pos[0], pos[j])
            for j in range(1, len(pos))
            if j != removed
        )
        values.append(after - current)
    return orient_by_values(tuple(values), low_points_high=False)


def source_sink_counts(adj: list[list[bool]]) -> tuple[int, int]:
    n = len(adj)
    indeg = [0] * n
    outdeg = [0] * n
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                outdeg[i] += 1
                indeg[j] += 1
    return sum(1 for d in indeg if d == 0), sum(1 for d in outdeg if d == 0)


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
    seen.clear()
    sizes: list[int] = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def directed_triangles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def score_hist(scores: tuple[int, ...]) -> str:
    counts = Counter(scores)
    return " ".join(f"{score}:{counts[score]}" for score in sorted(counts))


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
def hamiltonian_paths(n: int, masks: tuple[int, ...]) -> int:
    full = (1 << n) - 1
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
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
                dp[(mask | bit) * n + nxt] += value
                allowed -= bit
    return sum(dp[full * n : (full + 1) * n])


def metrics_for_adj(gauge: str, adj: list[list[bool]], observer: bool) -> GaugeMetrics:
    n = len(adj)
    scores = tuple(sorted(sum(row) for row in adj))
    hp: int | None = None
    h_ratio: Fraction | None = None
    if n <= 14:
        hp = hamiltonian_paths(n, row_masks(adj))
        h_ratio = Fraction(hp * (2 ** (n - 1)), factorial(n))
    sources, sinks = source_sink_counts(adj)
    return GaugeMetrics(
        gauge=gauge,
        vertices=n,
        hamiltonian_paths=hp,
        h_ratio=h_ratio,
        score_width=max(scores) - min(scores),
        score_hist=score_hist(scores),
        cyclic_triples=directed_triangles(adj),
        largest_scc=scc_sizes(adj)[0],
        source_count=sources,
        sink_count=sinks,
        observer_out=sum(adj[0]) if observer else None,
    )


TIME_GAUGES = (
    ("phase_half", lambda row, pos, th: phase_from_positions(pos)),
    ("origin_danger_rank", lambda row, pos, th: orient_by_values(tuple(dist0(x) for x in pos), True)),
    ("local_moat_rank", lambda row, pos, th: orient_by_values(tuple(adjacent_gap_data(pos)[2]), True)),
    ("close_pair_switch", lambda row, pos, th: close_pair_switch(pos, th)),
    ("danger_band_switch", lambda row, pos, th: band_pair_switch(pos, th)),
    ("quarter_lead_switch", lambda row, pos, th: quarter_lead_switch(pos, th)),
    ("open_arc_density", lambda row, pos, th: open_arc_density_switch(pos, th)),
    ("open_arc_count", lambda row, pos, th: open_arc_count_switch(pos, th)),
    ("blocker_relief", lambda row, pos, th: blocker_relief(pos, th)),
    ("two_neighbor_relief", lambda row, pos, th: two_neighbor_relief(pos, th)),
    ("threshold_deficit", lambda row, pos, th: threshold_deficit_relief(pos, th)),
    ("observer_two_guard", lambda row, pos, th: observer_two_guard(pos, th)),
    ("observer_relief_rank", lambda row, pos, th: observer_relief_rank(pos, th)),
    ("kinetic_escape", kinetic_escape),
    ("kinetic_close_switch", kinetic_close_switch),
)


def endpoint_tournament(
    speeds: tuple[int, ...], mode: str
) -> list[list[bool]]:
    n = len(speeds)
    endpoints = S360.endpoints(speeds)
    values = sorted({endpoint.value for endpoint in endpoints})
    endpoints_by_owner: dict[int, list[S360.Endpoint]] = defaultdict(list)
    for endpoint in endpoints:
        endpoints_by_owner[endpoint.speed].append(endpoint)

    protection = {
        (protector, value): S360.direct_protects(speeds, protector, value)
        for protector in speeds
        for value in values
    }
    speed_index = {speed: idx for idx, speed in enumerate(speeds)}
    adj = [[False] * n for _ in range(n)]

    if mode == "endpoint_debt_rank":
        debt_by_speed: dict[int, Fraction] = {}
        for speed in speeds:
            debt = 0
            for endpoint in endpoints_by_owner[speed]:
                if not any(protection[(protector, endpoint.value)] for protector in speeds):
                    debt += 1
            debt_by_speed[speed] = Fraction(debt)
        return orient_by_values(tuple(debt_by_speed[speed] for speed in speeds), False)

    for i, j in combinations(range(n), 2):
        a = speeds[i]
        b = speeds[j]
        a_protects_b = sum(
            1 for endpoint in endpoints_by_owner[b] if protection[(a, endpoint.value)]
        )
        b_protects_a = sum(
            1 for endpoint in endpoints_by_owner[a] if protection[(b, endpoint.value)]
        )
        if a_protects_b == b_protects_a:
            winner = min(speed_index[a], speed_index[b])
        elif a_protects_b > b_protects_a:
            winner = speed_index[a]
        else:
            winner = speed_index[b]
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def selected_time(row: SpeedRow) -> tuple[str, Fraction]:
    kind, t, _gap = S502.selected_time_and_gap(row.speeds)
    return kind, t


def selected_records() -> list[RowGaugeRecord]:
    out: list[RowGaugeRecord] = []
    for row in rows():
        kind, t = selected_time(row)
        speeds0 = (0,) + row.speeds
        pos = positions(speeds0, t)
        threshold = Fraction(1, row.n)
        report = S356.report(row.label, list(row.speeds))
        summary = S360.summarize(list(row.speeds))
        for name, builder in TIME_GAUGES:
            adj = builder(row, pos, threshold)
            out.append(
                RowGaugeRecord(
                    row=row,
                    t_kind=kind,
                    t=t,
                    origin_margin=origin_margin(row.speeds, t, row.n),
                    safe_vertices=safe_vertex_count(pos, threshold),
                    gap_ratio=report.max_gap / report.threshold,
                    endpoint_debt=summary.unprotected_count,
                    gauge=metrics_for_adj(name, adj, observer=True),
                )
            )

        for mode in ("endpoint_cross_protect", "endpoint_debt_rank"):
            adj = endpoint_tournament(row.speeds, mode)
            out.append(
                RowGaugeRecord(
                    row=row,
                    t_kind="endpoint",
                    t=t,
                    origin_margin=origin_margin(row.speeds, t, row.n),
                    safe_vertices=safe_vertex_count(pos, threshold),
                    gap_ratio=report.max_gap / report.threshold,
                    endpoint_debt=summary.unprotected_count,
                    gauge=metrics_for_adj(mode, adj, observer=False),
                )
            )
    return out


def sample_times(row: SpeedRow) -> tuple[Fraction, ...]:
    kind, t = selected_time(row)
    q = 1009
    base = [Fraction(((37 + 97 * k) % q), q) for k in range(24)]
    return tuple(sorted(set(base + [t])))


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def correlation_audit() -> list[tuple[str, int, float | None, float | None, float | None, float | None]]:
    sample_rows = [row for row in rows() if row.n == 14 and row.mode in {"row-parent", "gate"}]
    by_gauge: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in sample_rows:
        threshold = Fraction(1, row.n)
        for t in sample_times(row):
            pos = positions((0,) + row.speeds, t)
            margin = float(origin_margin(row.speeds, t, row.n))
            for name, builder in TIME_GAUGES:
                adj = builder(row, pos, threshold)
                metrics = metrics_for_adj(name, adj, observer=True)
                by_gauge[name]["margin"].append(margin)
                by_gauge[name]["score_width"].append(float(metrics.score_width))
                by_gauge[name]["cyclic_triples"].append(float(metrics.cyclic_triples))
                by_gauge[name]["observer_out"].append(float(metrics.observer_out or 0))
                if metrics.h_ratio is not None:
                    by_gauge[name]["h_ratio"].append(float(metrics.h_ratio))
    out = []
    for name, data in sorted(by_gauge.items()):
        n = len(data["margin"])
        out.append(
            (
                name,
                n,
                pearson(data["margin"], data.get("h_ratio", [])),
                pearson(data["margin"], data["score_width"]),
                pearson(data["margin"], data["cyclic_triples"]),
                pearson(data["margin"], data["observer_out"]),
            )
        )
    return out


def print_selected_tables(records: list[RowGaugeRecord]) -> None:
    print("LRC arc-gauge zoo (codex-2026-06-01 S506c)")
    print("=" * 120)
    print("Selected rows: gauges are evaluated at the LRC witness boundary or gap midpoint.")
    print("H is exact for gauges with <=14 vertices; n=18 rows still report score/c3/SCC.")
    print()

    for row in rows():
        print(f"[{row.label}]")
        print(
            f"  mode={row.mode} t={fmt(selected_time(row)[1])} "
            f"origin/th={fmt(origin_margin(row.speeds, selected_time(row)[1], row.n))}"
        )
        print(
            "  gauge                  V        H   Hratio scoreW c3  SCC src/sink obsOut score-hist"
        )
        for record in [item for item in records if item.row == row]:
            gauge = record.gauge
            print(
                f"  {gauge.gauge:<22} {gauge.vertices:>2} "
                f"{fmt(gauge.hamiltonian_paths):>8} {fmt_dec(gauge.h_ratio):>7} "
                f"{gauge.score_width:>6} {gauge.cyclic_triples:>3} "
                f"{gauge.largest_scc:>4} {gauge.source_count:>2}/{gauge.sink_count:<2} "
                f"{fmt(gauge.observer_out):>6} {gauge.score_hist}"
            )
        print()


def gauge_read(name: str) -> str:
    if name.endswith("rank"):
        return "scalar rank; readable but often H=1"
    if name in {"phase_half", "open_arc_density", "open_arc_count"}:
        return "global circular spread"
    if name in {"close_pair_switch", "danger_band_switch", "quarter_lead_switch"}:
        return "edge-local threshold switch"
    if name in {"blocker_relief", "two_neighbor_relief", "threshold_deficit"}:
        return "deletion relief / pressure"
    if name == "observer_two_guard":
        return "marked observer two-neighbor guard"
    if name.startswith("kinetic"):
        return "velocity-aware wall-crossing"
    return "endpoint/boundary summary"


def print_correlation_table() -> None:
    print("CORRELATION AUDIT ON N14 ROW-PARENT/GATE TIME SAMPLES")
    print("=" * 120)
    print("Target is actual observer loneliness margin min||v t||/(1/n).")
    print("Positive corr means the invariant rises when the observer is lonelier.")
    print()
    print(
        f"{'gauge':<22} {'samples':>7} {'Hratio':>9} {'scoreW':>9} "
        f"{'c3':>9} {'obsOut':>9} read"
    )
    for name, n, hcorr, wcorr, c3corr, obscorr in correlation_audit():
        read = gauge_read(name)
        print(
            f"{name:<22} {n:>7} {fmt_dec(hcorr):>9} {fmt_dec(wcorr):>9} "
            f"{fmt_dec(c3corr):>9} {fmt_dec(obscorr):>9} {read}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 120)
    print("1. Rank gauges are useful diagnostics but poor tournaments: H usually collapses to 1.")
    print("2. Phase/open-arc gauges are global circular-spread meters; phase_H is the HYP-1970 channel.")
    print("3. Close/danger/kinetic-close switches preserve local 1/n crowding and approach bits.")
    print("4. Blocker, two-neighbor, and threshold-deficit relief are pressure gauges; SCC > 1 is the core alarm.")
    print("5. Observer-two-guard is the marked loneliness channel: observer outdegree tracks the bracket directly.")
    print("6. Runner-level endpoint summaries were transitive here; endpoint data should stay labelled, not collapsed.")
    print("7. A useful LRC metric is therefore a bundle:")
    print("   phase_H + open_arc scores + threshold-switch c3 + pressure SCC/sinks + observer out + endpoint debt.")


def main() -> None:
    records = selected_records()
    print_selected_tables(records)
    print_correlation_table()
    print_synthesis()


if __name__ == "__main__":
    main()
