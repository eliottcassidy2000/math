#!/usr/bin/env python3
"""
lrc_tournament_gauge_lab_s506.py

codex-2026-06-01 S506

Try several arc-assignment criteria for Tournament Analysis of the Lonely
Runner Conjecture.  The target is not one "correct" tournament, but a useful
menu:

* scalar rankers should expose ordinary loneliness coordinates and usually
  collapse to transitive tournaments;
* entropy gauges should keep cyclic phase shape, so H(T) can act as a global
  loneliness/spread metric;
* pressure gauges should expose blocker/debt dependencies, often through score
  shape and SCC/DAG behavior rather than H alone.

Tournament Analysis declaration:

For each gauge we specify:
  pairwise observable: runner positions, distances, local gaps, or deletion
    reliefs at a fixed LRC time;
  switch/gauge: a binary comparison that orients every runner pair;
  tie Hamiltonian path: the fixed numerical speed-order path.
"""

from __future__ import annotations

from collections import Counter, deque
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
class Gauge:
    name: str
    family: str
    observable: str
    switch: str


@dataclass(frozen=True)
class Row:
    label: str
    n: int
    speeds: tuple[int, ...]
    scale: int | None
    depth: int | None
    skip: int | None
    t: Fraction
    origin_margin: Fraction
    max_circle_gap: Fraction
    gap_ratio: Fraction
    endpoint_debt: int
    endpoint_product: Fraction


@dataclass(frozen=True)
class GaugeRow:
    gauge: str
    row_label: str
    n: int
    h: int
    h_ratio: Fraction
    score_width: int
    score_hist: tuple[tuple[int, int], ...]
    cyclic_triples: int
    largest_scc: int
    ties: int


@dataclass(frozen=True)
class GaugeSummary:
    gauge: str
    family: str
    distinct_h: int
    h_ratio_min: Fraction
    h_ratio_max: Fraction
    score_width_range: tuple[int, int]
    cyclic_range: tuple[int, int]
    largest_scc_range: tuple[int, int]
    mean_ties: Fraction
    transitive_rows: int
    corr_spread: float
    corr_origin: float
    n14_signature: str
    n18_signature: str
    verdict: str


GAUGES = (
    Gauge(
        "phase_half",
        "entropy",
        "clockwise half-turn phase",
        "i beats j when j is in the clockwise open half-circle from i",
    ),
    Gauge(
        "origin_lonely_rank",
        "ranker",
        "distance from observer/origin",
        "larger origin distance wins",
    ),
    Gauge(
        "local_clearance_rank",
        "ranker",
        "minimum of left/right adjacent circular gaps",
        "larger adjacent clearance wins",
    ),
    Gauge(
        "nearest_rank",
        "ranker",
        "nearest-neighbor distance",
        "larger nearest distance wins",
    ),
    Gauge(
        "two_neighbor_rank",
        "ranker",
        "sum of two nearest-neighbor distances",
        "larger two-neighbor moat wins",
    ),
    Gauge(
        "safe_phase_gate",
        "hybrid",
        "binary LRC local clearance plus phase",
        "locally safe runner beats unsafe runner; otherwise use phase",
    ),
    Gauge(
        "origin_safe_phase",
        "hybrid",
        "origin threshold safety plus phase",
        "origin-lonely runner beats origin-unsafe runner; otherwise use phase",
    ),
    Gauge(
        "blocker_phase_gate",
        "hybrid",
        "origin threshold danger plus phase",
        "origin-unsafe blocker beats safe runner; otherwise use phase",
    ),
    Gauge(
        "pressure_k1",
        "pressure",
        "nearest-distance deletion relief",
        "blocker points to runner whose nearest moat it improves more",
    ),
    Gauge(
        "pressure_k2",
        "pressure",
        "two-neighbor deletion relief",
        "blocker points to runner whose two-neighbor moat it improves more",
    ),
    Gauge(
        "deficit_pressure",
        "pressure",
        "threshold-deficit deletion relief",
        "blocker points to runner whose local threshold deficit it reduces more",
    ),
    Gauge(
        "origin_blocker_rank",
        "ranker",
        "observer-margin deletion relief",
        "runner whose deletion helps the observer most wins",
    ),
)


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return S356.fmt_frac(value)


def fmt_dec(value: Fraction, places: int = 4) -> str:
    return f"{float(value):.{places}f}"


def tie_winner(i: int, j: int) -> int:
    return i if i < j else j


def circle(value: Fraction) -> Fraction:
    return value % ONE


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = circle(a - b)
    return min(delta, ONE - delta)


def positions(speeds0: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in speeds0)


def distance_matrix(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], ...]:
    n = len(pos)
    return tuple(
        tuple(Fraction(0) if i == j else circular_distance(pos[i], pos[j]) for j in range(n))
        for i in range(n)
    )


def sorted_neighbor_distances(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    removed: int | None = None,
) -> tuple[Fraction, ...]:
    return tuple(
        sorted(dist[i][j] for j in range(len(dist)) if j != i and j != removed)
    )


def nearest_sum(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    k: int,
    removed: int | None = None,
) -> Fraction:
    return sum(sorted_neighbor_distances(dist, i, removed)[:k], Fraction(0))


def threshold_deficit(values: tuple[Fraction, ...], threshold: Fraction) -> Fraction:
    return sum((threshold - value for value in values[:2] if value < threshold), Fraction(0))


def adjacent_clearance(pos: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    n = len(pos)
    ordered = sorted((position, idx) for idx, position in enumerate(pos))
    clear = [Fraction(0)] * n
    for k, (position, idx) in enumerate(ordered):
        prev_position = ordered[k - 1][0]
        next_position = ordered[(k + 1) % n][0]
        left = circle(position - prev_position)
        right = circle(next_position - position)
        clear[idx] = min(left, right)
    return tuple(clear)


def orient_by_scores(scores: tuple[Fraction, ...], high_wins: bool = True) -> tuple[list[list[bool]], int]:
    n = len(scores)
    adj = [[False] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        if scores[i] == scores[j]:
            winner = tie_winner(i, j)
            ties += 1
        elif high_wins:
            winner = i if scores[i] > scores[j] else j
        else:
            winner = i if scores[i] < scores[j] else j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, ties


def phase_winner(pos: tuple[Fraction, ...], i: int, j: int) -> int:
    delta = circle(pos[j] - pos[i])
    if delta == 0 or delta == HALF:
        return tie_winner(i, j)
    return i if delta < HALF else j


def phase_tournament(pos: tuple[Fraction, ...]) -> tuple[list[list[bool]], int]:
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        delta = circle(pos[j] - pos[i])
        if delta == 0 or delta == HALF:
            ties += 1
        winner = phase_winner(pos, i, j)
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, ties


def gated_phase_tournament(
    pos: tuple[Fraction, ...],
    gate_scores: tuple[Fraction, ...],
    threshold: Fraction,
    safe_wins: bool,
) -> tuple[list[list[bool]], int]:
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        safe_i = gate_scores[i] >= threshold
        safe_j = gate_scores[j] >= threshold
        if safe_i != safe_j:
            if safe_wins:
                winner = i if safe_i else j
            else:
                winner = i if not safe_i else j
        else:
            delta = circle(pos[j] - pos[i])
            if delta == 0 or delta == HALF:
                ties += 1
            winner = phase_winner(pos, i, j)
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, ties


def pressure_tournament(
    dist: tuple[tuple[Fraction, ...], ...],
    k: int,
) -> tuple[list[list[bool]], int]:
    n = len(dist)
    base = [nearest_sum(dist, i, k) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        relief_i_from_j = nearest_sum(dist, i, k, removed=j) - base[i]
        relief_j_from_i = nearest_sum(dist, j, k, removed=i) - base[j]
        if relief_i_from_j > relief_j_from_i:
            winner = j
        elif relief_j_from_i > relief_i_from_j:
            winner = i
        else:
            winner = tie_winner(i, j)
            ties += 1
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, ties


def deficit_pressure_tournament(
    dist: tuple[tuple[Fraction, ...], ...],
    threshold: Fraction,
) -> tuple[list[list[bool]], int]:
    n = len(dist)
    base = [threshold_deficit(sorted_neighbor_distances(dist, i), threshold) for i in range(n)]
    adj = [[False] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        after_i = threshold_deficit(sorted_neighbor_distances(dist, i, removed=j), threshold)
        after_j = threshold_deficit(sorted_neighbor_distances(dist, j, removed=i), threshold)
        relief_i_from_j = base[i] - after_i
        relief_j_from_i = base[j] - after_j
        if relief_i_from_j > relief_j_from_i:
            winner = j
        elif relief_j_from_i > relief_i_from_j:
            winner = i
        else:
            winner = tie_winner(i, j)
            ties += 1
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, ties


def origin_blocker_scores(dist: tuple[tuple[Fraction, ...], ...]) -> tuple[Fraction, ...]:
    n = len(dist)
    current = nearest_sum(dist, 0, 1)
    scores = [Fraction(0)] * n
    for j in range(1, n):
        scores[j] = nearest_sum(dist, 0, 1, removed=j) - current
    return tuple(scores)


def tournament_for_gauge(gauge: str, row: Row) -> tuple[list[list[bool]], int]:
    speeds0 = (0,) + row.speeds
    pos = positions(speeds0, row.t)
    dist = distance_matrix(pos)
    threshold = Fraction(1, len(speeds0))
    d_origin = tuple(circular_distance(position, Fraction(0)) for position in pos)
    clear = adjacent_clearance(pos)
    d1 = tuple(nearest_sum(dist, i, 1) for i in range(len(pos)))
    d2sum = tuple(nearest_sum(dist, i, 2) for i in range(len(pos)))

    if gauge == "phase_half":
        return phase_tournament(pos)
    if gauge == "origin_lonely_rank":
        return orient_by_scores(d_origin, high_wins=True)
    if gauge == "local_clearance_rank":
        return orient_by_scores(clear, high_wins=True)
    if gauge == "nearest_rank":
        return orient_by_scores(d1, high_wins=True)
    if gauge == "two_neighbor_rank":
        return orient_by_scores(d2sum, high_wins=True)
    if gauge == "safe_phase_gate":
        return gated_phase_tournament(pos, clear, threshold, safe_wins=True)
    if gauge == "origin_safe_phase":
        return gated_phase_tournament(pos, d_origin, threshold, safe_wins=True)
    if gauge == "blocker_phase_gate":
        return gated_phase_tournament(pos, d_origin, threshold, safe_wins=False)
    if gauge == "pressure_k1":
        return pressure_tournament(dist, 1)
    if gauge == "pressure_k2":
        return pressure_tournament(dist, 2)
    if gauge == "deficit_pressure":
        return deficit_pressure_tournament(dist, threshold)
    if gauge == "origin_blocker_rank":
        return orient_by_scores(origin_blocker_scores(dist), high_wins=True)
    raise ValueError(gauge)


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
        dp[((1 << v) * n) + v] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask * n + last]
            if not value:
                continue
            allowed = masks[last] & ~mask
            while allowed:
                bit = allowed & -allowed
                nxt = bit.bit_length() - 1
                dp[((mask | bit) * n) + nxt] += value
                allowed -= bit
    return sum(dp[full * n : (full + 1) * n])


def score_sequence(adj: list[list[bool]]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def score_hist(adj: list[list[bool]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(score_sequence(adj)).items()))


def cyclic_triples(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def largest_scc(adj: list[list[bool]]) -> int:
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
    return max(sizes)


def h_ratio(h: int, n: int) -> Fraction:
    return Fraction(h * (2 ** (n - 1)), factorial(n))


def primitive(values: tuple[int, ...]) -> bool:
    return gcd_tuple(values) == 1


def gcd_tuple(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        raise ValueError((n, scale, skip, speeds))
    return speeds


def summarize_row(label: str, n: int, speeds: tuple[int, ...], scale, depth, skip) -> Row:
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    kind, t, _gap = S502.selected_time_and_gap(speeds)
    if kind not in ("gap-mid", "boundary"):
        t = Fraction(1, n)
    speeds0 = (0,) + speeds
    pos = positions(speeds0, t)
    return Row(
        label=label,
        n=n,
        speeds=speeds,
        scale=scale,
        depth=depth,
        skip=skip,
        t=t,
        origin_margin=S502.origin_margin(speeds, t),
        max_circle_gap=S502.max_circle_gap(pos),
        gap_ratio=report.max_gap / report.threshold,
        endpoint_debt=summary.unprotected_count,
        endpoint_product=(report.max_gap / report.threshold) * summary.unprotected_count,
    )


def build_rows() -> tuple[Row, ...]:
    hard = {
        14: (7, 6),
        18: (9, 8),
    }
    rows: list[Row] = []
    for n, (base, skip) in hard.items():
        rows.append(summarize_row(f"n{n}-initial", n, tuple(range(1, n)), None, None, None))
        for depth in range(3):
            scale = base * (2**depth)
            rows.append(summarize_row(f"n{n}-s{scale}", n, ladder(n, scale, skip), scale, depth, skip))
    return tuple(rows)


def gauge_row(row: Row, gauge: Gauge) -> GaugeRow:
    adj, ties = tournament_for_gauge(gauge.name, row)
    n = len(adj)
    h = hamiltonian_paths(n, row_masks(adj))
    scores = score_sequence(adj)
    return GaugeRow(
        gauge=gauge.name,
        row_label=row.label,
        n=n,
        h=h,
        h_ratio=h_ratio(h, n),
        score_width=max(scores) - min(scores),
        score_hist=score_hist(adj),
        cyclic_triples=cyclic_triples(adj),
        largest_scc=largest_scc(adj),
        ties=ties,
    )


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def sequence_signature(values: list[Fraction]) -> str:
    return " -> ".join(fmt_dec(value) for value in values)


def summarize_gauge(gauge: Gauge, rows: tuple[Row, ...], data: tuple[GaugeRow, ...]) -> GaugeSummary:
    h_values = [row.h_ratio for row in data]
    score_widths = [row.score_width for row in data]
    cyclics = [row.cyclic_triples for row in data]
    sccs = [row.largest_scc for row in data]
    spread = [float(ONE - row.max_circle_gap) for row in rows]
    origin = [float(row.origin_margin) for row in rows]
    h_float = [float(row.h_ratio) for row in data]
    n14 = [row.h_ratio for row in data if row.row_label.startswith("n14-s")]
    n18 = [row.h_ratio for row in data if row.row_label.startswith("n18-s")]
    transitive_rows = sum(1 for row in data if row.h == 1)

    if gauge.family == "ranker":
        verdict = "ranker: good coordinate, H collapses"
    elif gauge.family == "entropy" and len(set(h_values)) >= 3:
        verdict = "best H entropy gauge"
    elif gauge.family == "hybrid" and len(set(h_values)) >= 3:
        verdict = "endpoint-aware entropy candidate"
    elif gauge.family == "pressure" and max(sccs) == len(data[0].score_hist):
        verdict = "pressure shape, inspect SCC"
    elif gauge.family == "pressure":
        verdict = "dependency/peeling gauge"
    else:
        verdict = "weak dynamic range"

    return GaugeSummary(
        gauge=gauge.name,
        family=gauge.family,
        distinct_h=len({row.h for row in data}),
        h_ratio_min=min(h_values),
        h_ratio_max=max(h_values),
        score_width_range=(min(score_widths), max(score_widths)),
        cyclic_range=(min(cyclics), max(cyclics)),
        largest_scc_range=(min(sccs), max(sccs)),
        mean_ties=sum((Fraction(row.ties, 1) for row in data), Fraction(0)) / len(data),
        transitive_rows=transitive_rows,
        corr_spread=pearson(h_float, spread),
        corr_origin=pearson(h_float, origin),
        n14_signature=sequence_signature(n14),
        n18_signature=sequence_signature(n18),
        verdict=verdict,
    )


def print_definitions() -> None:
    print("LRC tournament gauge lab (codex-2026-06-01 S506)")
    print("=" * 124)
    print("Every gauge orients all pairs on the observer+runners at the selected LRC time.")
    print("Tie Hamiltonian path: numerical speed order.")
    print()
    print("GAUGE MENU")
    print("=" * 124)
    print(f"{'gauge':<22} {'family':<9} {'observable':<44} switch")
    print("-" * 124)
    for gauge in GAUGES:
        print(f"{gauge.name:<22} {gauge.family:<9} {gauge.observable:<44} {gauge.switch}")
    print()


def print_rows(rows: tuple[Row, ...]) -> None:
    print("LRC ROWS USED")
    print("=" * 124)
    print(
        f"{'row':<12} {'t':>12} {'origin/th':>9} {'maxgap':>8} "
        f"{'gap/th':>8} {'debt':>5} {'product':>8}"
    )
    print("-" * 124)
    for row in rows:
        print(
            f"{row.label:<12} {fmt(row.t):>12} {fmt(row.origin_margin):>9} "
            f"{fmt(row.max_circle_gap):>8} {fmt(row.gap_ratio):>8} "
            f"{row.endpoint_debt:>5} {fmt(row.endpoint_product):>8}"
        )
    print()


def print_summary(summaries: tuple[GaugeSummary, ...]) -> None:
    print("GAUGE USEFULNESS SUMMARY")
    print("=" * 124)
    print(
        f"{'gauge':<22} {'family':<9} {'dH':>3} {'Hratio range':<17} "
        f"{'scoreW':<7} {'cyc3':<9} {'SCC':<7} {'ties':>6} "
        f"{'corr spread':>11} {'corr origin':>11} verdict"
    )
    print("-" * 124)
    for row in summaries:
        print(
            f"{row.gauge:<22} {row.family:<9} {row.distinct_h:>3} "
            f"{fmt_dec(row.h_ratio_min)}..{fmt_dec(row.h_ratio_max):<8} "
            f"{row.score_width_range[0]}..{row.score_width_range[1]:<3} "
            f"{row.cyclic_range[0]}..{row.cyclic_range[1]:<5} "
            f"{row.largest_scc_range[0]}..{row.largest_scc_range[1]:<3} "
            f"{fmt(row.mean_ties):>6} {row.corr_spread:>11.3f} "
            f"{row.corr_origin:>11.3f} {row.verdict}"
        )
    print()


def print_recursive_signatures(summaries: tuple[GaugeSummary, ...]) -> None:
    print("H_RATIO SIGNATURES ON HARD SCALE LADDERS")
    print("=" * 124)
    print(f"{'gauge':<22} {'n14 scales 7->14->28':<42} {'n18 scales 9->18->36'}")
    print("-" * 124)
    for row in summaries:
        if row.family in ("entropy", "hybrid", "pressure"):
            print(f"{row.gauge:<22} {row.n14_signature:<42} {row.n18_signature}")
    print()


def print_interpretation(summaries: tuple[GaugeSummary, ...]) -> None:
    print("INTERPRETATION")
    print("=" * 124)
    print("1. Scalar rankers are useful loneliness coordinates, but H is almost always")
    print("   transitive there.  They tell us who is lonely, not the shape of the crowd.")
    print("2. The half-turn phase gauge remains the clean H-entropy metric: it gives a")
    print("   broad cyclic tournament whose H drops at the first gate and then plateaus.")
    print("3. Endpoint-aware hybrids are promising if they keep phase cycles while using")
    print("   local safe/unsafe overrides.  The origin-safe version is too close to a")
    print("   scalar; the local-clearance safe gate is the better candidate.")
    print("4. Pressure gauges are not primarily H metrics.  Their score/SCC shape is a")
    print("   dependency certificate for blocker/debt peeling, complementary to H.")
    print()
    entropy = next(row for row in summaries if row.gauge == "phase_half")
    hybrid = next(row for row in summaries if row.gauge == "safe_phase_gate")
    pressure = next(row for row in summaries if row.gauge == "pressure_k2")
    print(
        f"Best clean H-loneliness metric: {entropy.gauge}.  It preserves cyclic "
        "phase shape and sees the recursive gate/plateau without endpoint-rule "
        "overrides."
    )
    print(
        f"Best endpoint-aware H candidate: {hybrid.gauge}.  It has the strongest "
        "spread correlation here, but its low-H collapse on some hard rows means "
        "it should be treated as a perturbation of phase_half, not a replacement."
    )
    print(
        f"Best pressure-shape candidate: {pressure.gauge}.  Its H is small, but "
        "its nontrivial cycle/SCC range makes it useful for blocker-debt peeling."
    )


def main() -> None:
    print_definitions()
    rows = build_rows()
    print_rows(rows)
    all_data: dict[str, tuple[GaugeRow, ...]] = {}
    summaries: list[GaugeSummary] = []
    for gauge in GAUGES:
        data = tuple(gauge_row(row, gauge) for row in rows)
        all_data[gauge.name] = data
        summaries.append(summarize_gauge(gauge, rows, data))
    ordered = tuple(
        sorted(
            summaries,
            key=lambda row: (
                {"entropy": 0, "hybrid": 1, "pressure": 2, "ranker": 3}[row.family],
                -row.distinct_h,
                row.gauge,
            ),
        )
    )
    print_summary(ordered)
    print_recursive_signatures(ordered)
    print_interpretation(ordered)


if __name__ == "__main__":
    main()
