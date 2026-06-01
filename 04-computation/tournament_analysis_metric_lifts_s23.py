#!/usr/bin/env python3
"""
tournament_analysis_metric_lifts_s23.py

oracle-2026-06-01 S23

Tournament Analysis, in this repo's operational sense:

    raw pairwise or geometric data
      -> comparator / tie-break rule
      -> tournament
      -> tournament invariants and wall-crossing as variables move.

This script tries several metric-to-tournament lifts:

1. directed pair data, modeled by basketball pass counts;
2. runners moving on a circle, using arc/chord/nearest-cell metrics;
3. cuboid/sphere/simplex lifts from continuous point clouds;
4. LRC witness rows from the active n=14,15,16 frontier.

The main test is whether a lift produces genuine cyclic tournament structure
or collapses to a scalar ranking/transitive tournament.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import acos, cos, gcd, log, pi, sin, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

EPS = 1.0e-12


@dataclass(frozen=True)
class TournamentStats:
    n: int
    bits: int
    score_hist: tuple[tuple[int, int], ...]
    cyclic_triples: int
    cycles_through_0: int
    hamiltonian_paths: int | None
    regularity_defect2: int
    scc_count: int


@dataclass(frozen=True)
class SeriesSummary:
    family: str
    lift: str
    kind: str
    n: int
    samples: int
    distinct_states: int
    total_flips: int
    max_step_flips: int
    transitive_samples: int
    cyclic_range: tuple[int, int]
    h_range: tuple[int | None, int | None]
    regularity_range2: tuple[int, int]
    scc_range: tuple[int, int]


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_range(values: tuple[int | None, int | None]) -> str:
    lo, hi = values
    if lo is None or hi is None:
        return "-"
    if lo == hi:
        return str(lo)
    return f"{lo}..{hi}"


def tie_beats(i: int, j: int) -> bool:
    """Fixed Hamiltonian path tie-break: 1 -> 2 -> ... -> n."""

    return i < j


def orient_by_score(scores: list[float]) -> list[list[int]]:
    n = len(scores)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if abs(scores[i] - scores[j]) <= EPS:
            winner = i if tie_beats(i, j) else j
        else:
            winner = i if scores[i] > scores[j] else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def orient_by_edge_rule(n: int, beats) -> list[list[int]]:
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        winner = i if beats(i, j) else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def bits_from_adj(adj: list[list[int]]) -> int:
    n = len(adj)
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                bits |= 1 << idx
            idx += 1
    return bits


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        row = dp[mask]
        for v, count in enumerate(row):
            if count == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += count
    return sum(dp[full])


def scc_count(adj: list[list[int]]) -> int:
    n = len(adj)
    reach = [[bool(adj[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if not reach[i][k]:
                continue
            for j in range(n):
                reach[i][j] = reach[i][j] or reach[k][j]
    seen: set[int] = set()
    count = 0
    for i in range(n):
        if i in seen:
            continue
        count += 1
        for j in range(n):
            if reach[i][j] and reach[j][i]:
                seen.add(j)
    return count


def stats(adj: list[list[int]], compute_h: bool = True) -> TournamentStats:
    n = len(adj)
    outdegrees = [sum(row) for row in adj]
    cyclic = 0
    cyclic0 = 0
    for triple in combinations(range(n), 3):
        local_out = Counter()
        for i, j in combinations(triple, 2):
            winner = i if adj[i][j] else j
            local_out[winner] += 1
        if sorted(local_out.values()) == [1, 1, 1]:
            cyclic += 1
            if 0 in triple:
                cyclic0 += 1
    h = hamiltonian_path_count(adj) if compute_h and n <= 10 else None
    return TournamentStats(
        n=n,
        bits=bits_from_adj(adj),
        score_hist=tuple(sorted(Counter(outdegrees).items())),
        cyclic_triples=cyclic,
        cycles_through_0=cyclic0,
        hamiltonian_paths=h,
        regularity_defect2=sum(abs(2 * d - (n - 1)) for d in outdegrees),
        scc_count=scc_count(adj),
    )


def summarize_series(
    family: str, lift: str, kind: str, adjs: list[list[list[int]]]
) -> SeriesSummary:
    rows = [stats(adj) for adj in adjs]
    bits = [row.bits for row in rows]
    flips = [
        (bits[idx] ^ bits[idx - 1]).bit_count()
        for idx in range(1, len(bits))
    ]
    h_values = [row.hamiltonian_paths for row in rows if row.hamiltonian_paths is not None]
    h_range = (min(h_values), max(h_values)) if h_values else (None, None)
    return SeriesSummary(
        family=family,
        lift=lift,
        kind=kind,
        n=rows[0].n,
        samples=len(rows),
        distinct_states=len(set(bits)),
        total_flips=sum(flips),
        max_step_flips=max(flips) if flips else 0,
        transitive_samples=sum(1 for row in rows if row.cyclic_triples == 0),
        cyclic_range=(min(row.cyclic_triples for row in rows), max(row.cyclic_triples for row in rows)),
        h_range=h_range,
        regularity_range2=(min(row.regularity_defect2 for row in rows), max(row.regularity_defect2 for row in rows)),
        scc_range=(min(row.scc_count for row in rows), max(row.scc_count for row in rows)),
    )


def print_series_table(title: str, rows: list[SeriesSummary]) -> None:
    print(title)
    print("=" * 118)
    print(
        "  family      lift                    kind        n samp states flips maxF "
        "trans cyc-range H-range reg2 scc"
    )
    for row in rows:
        cyc = fmt_range(row.cyclic_range)
        reg = fmt_range(row.regularity_range2)
        scc = fmt_range(row.scc_range)
        print(
            f"  {row.family:<11} {row.lift:<23} {row.kind:<10} "
            f"{row.n:>2} {row.samples:>4} {row.distinct_states:>6} "
            f"{row.total_flips:>5} {row.max_step_flips:>4} "
            f"{row.transitive_samples:>5} {cyc:<9} {fmt_range(row.h_range):<8} "
            f"{reg:<7} {scc}"
        )
    print()


def pass_matrix(n: int, pairs: dict[tuple[int, int], tuple[int, int]]) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for (i, j), (ij, ji) in pairs.items():
        matrix[i][j] = ij
        matrix[j][i] = ji
    return matrix


def pass_tournament(matrix: list[list[int]]) -> tuple[list[list[int]], int]:
    n = len(matrix)
    ties = 0

    def beats(i: int, j: int) -> bool:
        nonlocal ties
        if matrix[i][j] == matrix[j][i]:
            ties += 1
            return tie_beats(i, j)
        return matrix[i][j] > matrix[j][i]

    return orient_by_edge_rule(n, beats), ties


def basketball_examples() -> None:
    quarters = [
        (
            "Q1 guard hub",
            pass_matrix(
                5,
                {
                    (0, 1): (8, 8),
                    (0, 2): (7, 2),
                    (0, 3): (6, 1),
                    (0, 4): (9, 3),
                    (1, 2): (5, 4),
                    (1, 3): (4, 2),
                    (1, 4): (6, 6),
                    (2, 3): (3, 5),
                    (2, 4): (4, 7),
                    (3, 4): (6, 4),
                },
            ),
        ),
        (
            "Q2 triangle",
            pass_matrix(
                5,
                {
                    (0, 1): (9, 4),
                    (0, 2): (3, 8),
                    (0, 3): (6, 6),
                    (0, 4): (5, 7),
                    (1, 2): (8, 5),
                    (1, 3): (4, 7),
                    (1, 4): (8, 3),
                    (2, 3): (7, 4),
                    (2, 4): (5, 5),
                    (3, 4): (4, 9),
                },
            ),
        ),
        (
            "Q3 post split",
            pass_matrix(
                5,
                {
                    (0, 1): (4, 7),
                    (0, 2): (5, 5),
                    (0, 3): (3, 6),
                    (0, 4): (7, 4),
                    (1, 2): (6, 2),
                    (1, 3): (7, 3),
                    (1, 4): (3, 8),
                    (2, 3): (8, 6),
                    (2, 4): (4, 9),
                    (3, 4): (6, 6),
                },
            ),
        ),
        (
            "Q4 scramble",
            pass_matrix(
                5,
                {
                    (0, 1): (6, 9),
                    (0, 2): (8, 3),
                    (0, 3): (4, 7),
                    (0, 4): (5, 5),
                    (1, 2): (3, 8),
                    (1, 3): (7, 4),
                    (1, 4): (6, 2),
                    (2, 3): (5, 9),
                    (2, 4): (8, 6),
                    (3, 4): (9, 4),
                },
            ),
        ),
    ]
    print("A. DISCRETE DIRECTED DATA: BASKETBALL PASS TOURNAMENT")
    print("=" * 118)
    print("  tie-break is fixed path 1 -> 2 -> 3 -> 4 -> 5.")
    print("  quarter        ties score_hist      cycles H  scc reg2 flips")
    previous_bits: int | None = None
    for label, matrix in quarters:
        adj, ties = pass_tournament(matrix)
        row = stats(adj)
        flips = "-" if previous_bits is None else str((row.bits ^ previous_bits).bit_count())
        previous_bits = row.bits
        print(
            f"  {label:<14} {ties:>4} {row.score_hist!s:<15} "
            f"{row.cyclic_triples:>6} {row.hamiltonian_paths:>2} "
            f"{row.scc_count:>4} {row.regularity_defect2:>4} {flips:>5}"
        )
    print()


def circle_positions(speeds: list[int], t: float) -> list[float]:
    return [(speed * t) % 1.0 for speed in speeds]


def circ_delta(a: float, b: float) -> float:
    return (b - a) % 1.0


def circ_dist(a: float, b: float) -> float:
    delta = circ_delta(a, b)
    return min(delta, 1.0 - delta)


def chord_dist(a: float, b: float) -> float:
    return sqrt(max(0.0, 2.0 - 2.0 * cos(2.0 * pi * circ_delta(a, b))))


def entropy(values: list[float]) -> float:
    total = sum(values)
    if total <= EPS:
        return 0.0
    out = 0.0
    for value in values:
        if value <= EPS:
            continue
        p = value / total
        out -= p * log(p)
    return out


def phase_halfturn_tournament(positions: list[float]) -> list[list[int]]:
    n = len(positions)

    def beats(i: int, j: int) -> bool:
        delta = circ_delta(positions[i], positions[j])
        if abs(delta) <= EPS or abs(delta - 0.5) <= EPS:
            return tie_beats(i, j)
        return delta < 0.5

    return orient_by_edge_rule(n, beats)


def pivot_area_tournament(positions: list[float]) -> list[list[int]]:
    n = len(positions)

    def beats(i: int, j: int) -> bool:
        value = sin(2.0 * pi * (positions[j] - positions[i]))
        if abs(value) <= EPS:
            return tie_beats(i, j)
        return value > 0.0

    return orient_by_edge_rule(n, beats)


def row_sum_circle_tournament(positions: list[float], dist) -> list[list[int]]:
    scores = [
        sum(dist(positions[i], positions[j]) for j in range(len(positions)) if i != j)
        for i in range(len(positions))
    ]
    return orient_by_score(scores)


def anchor_distance_tournament(positions: list[float], dist) -> list[list[int]]:
    scores = [dist(positions[i], positions[0]) for i in range(len(positions))]
    return orient_by_score(scores)


def local_cell_tournament(positions: list[float]) -> list[list[int]]:
    n = len(positions)
    order = sorted(range(n), key=lambda idx: (positions[idx], idx))
    scores = [0.0] * n
    for order_idx, vertex in enumerate(order):
        prev_vertex = order[(order_idx - 1) % n]
        next_vertex = order[(order_idx + 1) % n]
        left = circ_delta(positions[prev_vertex], positions[vertex])
        right = circ_delta(positions[vertex], positions[next_vertex])
        scores[vertex] = min(left, right)
    return orient_by_score(scores)


def entropy_circle_tournament(positions: list[float], dist) -> list[list[int]]:
    scores = [
        entropy([dist(positions[i], positions[j]) for j in range(len(positions)) if i != j])
        for i in range(len(positions))
    ]
    return orient_by_score(scores)


def label_shift_circle_tournament(positions: list[float], dist, shift: int = 1) -> list[list[int]]:
    n = len(positions)

    def beats(i: int, j: int) -> bool:
        left = dist(positions[i], positions[(j + shift) % n])
        right = dist(positions[j], positions[(i + shift) % n])
        if abs(left - right) <= EPS:
            return tie_beats(i, j)
        return left < right

    return orient_by_edge_rule(n, beats)


def metric_switch_circle_tournament(positions: list[float], dist, resonance: int | None = None) -> list[list[int]]:
    """Use a symmetric pair metric to switch each edge against the base path."""

    n = len(positions)
    pair_distances = [
        dist(positions[i], positions[j])
        for i, j in combinations(range(n), 2)
    ]
    threshold = sorted(pair_distances)[len(pair_distances) // 2]

    def switch_value(i: int, j: int) -> float:
        value = dist(positions[i], positions[j])
        if resonance is None:
            return value - threshold
        return sin(resonance * pi * value)

    def beats(i: int, j: int) -> bool:
        value = switch_value(i, j)
        if abs(value) <= EPS:
            return tie_beats(i, j)
        follows_base_path = value > 0.0
        return tie_beats(i, j) if follows_base_path else tie_beats(j, i)

    return orient_by_edge_rule(n, beats)


def receding_from_anchor_tournament(speeds: list[int], t: float) -> list[list[int]]:
    here = circle_positions(speeds, t)
    there = circle_positions(speeds, t + 1.0e-5)
    scores = [
        chord_dist(there[i], there[0]) - chord_dist(here[i], here[0])
        for i in range(len(speeds))
    ]
    return orient_by_score(scores)


def circle_series() -> tuple[list[SeriesSummary], dict[str, tuple[int, ...]]]:
    speeds = [0, 1, 2, 3, 5, 8, 13, 21]
    times = [k / 29.0 for k in range(1, 29)]
    lifts = {
        "phase-halfturn": ("edge-local", lambda pos, _t: phase_halfturn_tournament(pos)),
        "pivot-area": ("edge-local", lambda pos, _t: pivot_area_tournament(pos)),
        "anchor-arc-far": ("score", lambda pos, _t: anchor_distance_tournament(pos, circ_dist)),
        "anchor-chord-far": ("score", lambda pos, _t: anchor_distance_tournament(pos, chord_dist)),
        "row-sum-arc": ("score", lambda pos, _t: row_sum_circle_tournament(pos, circ_dist)),
        "row-sum-chord": ("score", lambda pos, _t: row_sum_circle_tournament(pos, chord_dist)),
        "local-cell-radius": ("score", lambda pos, _t: local_cell_tournament(pos)),
        "profile-entropy": ("score", lambda pos, _t: entropy_circle_tournament(pos, chord_dist)),
        "label-shift-arc": ("edge-local", lambda pos, _t: label_shift_circle_tournament(pos, circ_dist)),
        "arc-median-switch": ("edge-switch", lambda pos, _t: metric_switch_circle_tournament(pos, circ_dist)),
        "chord-resonance": ("edge-switch", lambda pos, _t: metric_switch_circle_tournament(pos, chord_dist, resonance=5)),
        "receding-anchor": ("score", lambda _pos, t: receding_from_anchor_tournament(speeds, t)),
    }
    summaries: list[SeriesSummary] = []
    sequences: dict[str, tuple[int, ...]] = {}
    for lift_name, (kind, maker) in lifts.items():
        adjs = []
        for t in times:
            positions = circle_positions(speeds, t)
            adjs.append(maker(positions, t))
        summaries.append(summarize_series("circle8", lift_name, kind, adjs))
        sequences[lift_name] = tuple(bits_from_adj(adj) for adj in adjs)
    return summaries, sequences


def print_equivalence_clusters(family: str, sequences: dict[str, tuple[int, ...]]) -> None:
    buckets: defaultdict[tuple[int, ...], list[str]] = defaultdict(list)
    for name, seq in sequences.items():
        buckets[seq].append(name)
    print(f"{family} exact state-sequence clusters")
    print("=" * 118)
    for names in buckets.values():
        if len(names) > 1:
            print(f"  same tournament path: {', '.join(names)}")
    if all(len(names) == 1 for names in buckets.values()):
        print("  no exact coincidences")
    print()


def cube_points(speed_vectors: list[tuple[int, int, int]], t: float) -> list[tuple[float, float, float]]:
    return [tuple((component * t) % 1.0 for component in vector) for vector in speed_vectors]


def lp_dist(p: tuple[float, ...], q: tuple[float, ...], power: float) -> float:
    diffs = [abs(a - b) for a, b in zip(p, q)]
    if power == float("inf"):
        return max(diffs)
    return sum(value**power for value in diffs) ** (1.0 / power)


def row_sum_points_tournament(points: list[tuple[float, ...]], power: float) -> list[list[int]]:
    scores = [
        sum(lp_dist(points[i], points[j], power) for j in range(len(points)) if i != j)
        for i in range(len(points))
    ]
    return orient_by_score(scores)


def entropy_points_tournament(points: list[tuple[float, ...]], power: float) -> list[list[int]]:
    scores = [
        entropy([lp_dist(points[i], points[j], power) for j in range(len(points)) if i != j])
        for i in range(len(points))
    ]
    return orient_by_score(scores)


def label_shift_points_tournament(points: list[tuple[float, ...]], power: float, shift: int = 1) -> list[list[int]]:
    n = len(points)

    def beats(i: int, j: int) -> bool:
        left = lp_dist(points[i], points[(j + shift) % n], power)
        right = lp_dist(points[j], points[(i + shift) % n], power)
        if abs(left - right) <= EPS:
            return tie_beats(i, j)
        return left < right

    return orient_by_edge_rule(n, beats)


def metric_switch_points_tournament(
    points: list[tuple[float, ...]], power: float, resonance: int | None = None
) -> list[list[int]]:
    """Use an Lp pair distance to switch each edge against the base path."""

    n = len(points)
    pair_distances = [
        lp_dist(points[i], points[j], power)
        for i, j in combinations(range(n), 2)
    ]
    threshold = sorted(pair_distances)[len(pair_distances) // 2]

    def switch_value(i: int, j: int) -> float:
        value = lp_dist(points[i], points[j], power)
        if resonance is None:
            return value - threshold
        return sin(resonance * pi * value)

    def beats(i: int, j: int) -> bool:
        value = switch_value(i, j)
        if abs(value) <= EPS:
            return tie_beats(i, j)
        follows_base_path = value > 0.0
        return tie_beats(i, j) if follows_base_path else tie_beats(j, i)

    return orient_by_edge_rule(n, beats)


def det3(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> float:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def oriented_volume_tournament(points: list[tuple[float, float, float]]) -> list[list[int]]:
    n = len(points)
    center = (0.5, 0.5, 0.5)
    drift = (1.0, (1.0 + sqrt(5.0)) / 2.0, 2.0)

    def shifted(point: tuple[float, float, float]) -> tuple[float, float, float]:
        return (point[0] - center[0], point[1] - center[1], point[2] - center[2])

    def beats(i: int, j: int) -> bool:
        value = det3(shifted(points[i]), shifted(points[j]), drift)
        if abs(value) <= EPS:
            return tie_beats(i, j)
        return value > 0.0

    return orient_by_edge_rule(n, beats)


def xy_area_tournament(points: list[tuple[float, float, float]]) -> list[list[int]]:
    n = len(points)

    def beats(i: int, j: int) -> bool:
        ax, ay = points[i][0] - 0.5, points[i][1] - 0.5
        bx, by = points[j][0] - 0.5, points[j][1] - 0.5
        value = ax * by - ay * bx
        if abs(value) <= EPS:
            return tie_beats(i, j)
        return value > 0.0

    return orient_by_edge_rule(n, beats)


def cube_series() -> tuple[list[SeriesSummary], dict[str, tuple[int, ...]]]:
    speed_vectors = [
        (0, 0, 0),
        (1, 2, 3),
        (2, 3, 5),
        (3, 5, 8),
        (5, 8, 13),
        (8, 13, 21),
        (13, 21, 34),
        (21, 34, 55),
    ]
    times = [k / 31.0 for k in range(1, 31)]
    lifts = {
        "cuboid-L1-row": ("score", lambda pts: row_sum_points_tournament(pts, 1.0)),
        "cuboid-L2-row": ("score", lambda pts: row_sum_points_tournament(pts, 2.0)),
        "cuboid-Linf-row": ("score", lambda pts: row_sum_points_tournament(pts, float("inf"))),
        "simplex-entropy": ("score", lambda pts: entropy_points_tournament(pts, 2.0)),
        "label-shift-L2": ("edge-local", lambda pts: label_shift_points_tournament(pts, 2.0)),
        "L2-median-switch": ("edge-switch", lambda pts: metric_switch_points_tournament(pts, 2.0)),
        "Linf-resonance": ("edge-switch", lambda pts: metric_switch_points_tournament(pts, float("inf"), resonance=4)),
        "xy-area": ("edge-local", lambda pts: xy_area_tournament(pts)),
        "drift-volume": ("edge-local", lambda pts: oriented_volume_tournament(pts)),
    }
    summaries: list[SeriesSummary] = []
    sequences: dict[str, tuple[int, ...]] = {}
    for lift_name, (kind, maker) in lifts.items():
        adjs = []
        for t in times:
            adjs.append(maker(cube_points(speed_vectors, t)))
        summaries.append(summarize_series("cube8", lift_name, kind, adjs))
        sequences[lift_name] = tuple(bits_from_adj(adj) for adj in adjs)
    return summaries, sequences


def primitive_gcd(values: tuple[int, ...]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or primitive_gcd(speeds) != 1:
        raise ValueError(f"bad ladder n={n}, scale={scale}, skip={skip}")
    return speeds


def lrc_witness_time(speeds: tuple[int, ...]) -> Fraction:
    report = S356.report("lrc-metric-lift", list(speeds))
    if report.witness is not None:
        return report.witness
    if report.boundary_witness is not None:
        return report.boundary_witness
    raise ValueError(f"no witness for {speeds}")


def bracket_margin(positions: list[float], threshold: float) -> float:
    right = min(circ_delta(positions[0], positions[i]) for i in range(1, len(positions)))
    left = min(circ_delta(positions[i], positions[0]) for i in range(1, len(positions)))
    return min(left, right) - threshold


def lrc_witness_table() -> None:
    samples = [
        ("n14 d=7", ladder(14, 7, 6)),
        ("n14 d=14", ladder(14, 14, 6)),
        ("n15 d=3", ladder(15, 3, 6)),
        ("n15 d=5", ladder(15, 5, 6)),
        ("n16 d=2", ladder(16, 2, 14)),
        ("n16 d=8", ladder(16, 8, 14)),
    ]
    print("D. ACTIVE LRC WITNESS ROWS UNDER MULTIPLE TOURNAMENT LIFTS")
    print("=" * 118)
    print(
        "  sample     t          margin      phase cyc/H  anchor cyc/H "
        "cell cyc/H  shift cyc/H switch phase hist"
    )
    for label, moving_speeds in samples:
        n = len(moving_speeds) + 1
        speeds = (0,) + moving_speeds
        t_frac = lrc_witness_time(moving_speeds)
        t = float(t_frac)
        pos = circle_positions(list(speeds), t)
        margin = bracket_margin(pos, 1.0 / n)
        phase = stats(phase_halfturn_tournament(pos), compute_h=False)
        anchor = stats(anchor_distance_tournament(pos, circ_dist), compute_h=False)
        cell = stats(local_cell_tournament(pos), compute_h=False)
        shift = stats(label_shift_circle_tournament(pos, circ_dist), compute_h=False)
        switch = stats(metric_switch_circle_tournament(pos, circ_dist), compute_h=False)

        def cyc_h(row: TournamentStats) -> str:
            return f"{row.cyclic_triples}/-"

        print(
            f"  {label:<9} {fmt_frac(t_frac):<10} {margin:+.6f} "
            f"{cyc_h(phase):>11} {cyc_h(anchor):>13} "
            f"{cyc_h(cell):>10} {cyc_h(shift):>11} "
            f"{switch.cyclic_triples:>4}/- {phase.score_hist}"
        )
    print()


def print_method_synthesis(circle_rows: list[SeriesSummary], cube_rows: list[SeriesSummary]) -> None:
    all_rows = circle_rows + cube_rows
    score_rows = [row for row in all_rows if row.kind == "score"]
    edge_rows = [row for row in all_rows if row.kind in {"edge-local", "edge-switch"}]
    score_transitive = sum(row.transitive_samples for row in score_rows)
    score_samples = sum(row.samples for row in score_rows)
    edge_transitive = sum(row.transitive_samples for row in edge_rows)
    edge_samples = sum(row.samples for row in edge_rows)
    print("E. SYNTHESIS: WHAT TOURNAMENT ANALYSIS IS TESTING")
    print("=" * 118)
    print("  abstraction stack:")
    print("    raw data -> metric/sensor -> comparator -> tournament -> invariants -> wall-crossing path")
    print()
    print(
        "  Score lifts are rankers.  They orient i,j by scalar s_i - s_j, so "
        "distinct scores force transitive tournaments."
    )
    print(
        f"    observed score transitive samples: {score_transitive}/{score_samples}"
    )
    print(
        "  Edge-local and edge-switch lifts are analyzers.  They decide every "
        "pair from pair- or lens-specific data, so cycles can survive."
    )
    print(
        f"    observed edge-local transitive samples: {edge_transitive}/{edge_samples}"
    )
    print()
    print("  useful layer names:")
    print("    flux lift:       directed counts, e.g. basketball passes i->j vs j->i")
    print("    phase lift:      circular half-turn/chirality from moving runners")
    print("    rank lift:       any centrality, distance-to-anchor, entropy, or isolation score")
    print("    lens lift:       pair i,j compared through a third label, anchor, or shifted coordinate")
    print("    switch lift:     symmetric distance D_ij toggles a base Hamiltonian-path edge")
    print("    volume lift:     oriented area/volume after embedding in sphere/cuboid/simplex")
    print("    dynamic lift:    edge flips as the continuous variable crosses comparator walls")
    print()
    print(
        "  repo connection: the Hamiltonian tie-break path is not a detail.  It "
        "is the fixed labelled path used to complete ties, just as the tiling "
        "model uses a base path and S22 used lex completion for degenerate "
        "circle states."
    )
    print()


def main() -> None:
    print("TOURNAMENT ANALYSIS METRIC LIFTS (oracle-2026-06-01-S23)")
    print("=" * 118)
    print(
        "Tournament Analysis means choosing a meaningful binary comparator for "
        "each pair, then studying the resulting tournament and its changes."
    )
    print()

    basketball_examples()

    circle_rows, circle_sequences = circle_series()
    print_series_table("B. CONTINUOUS CIRCLE RUNNERS: MANY METRIC LIFTS", circle_rows)
    print_equivalence_clusters("circle8", circle_sequences)

    cube_rows, cube_sequences = cube_series()
    print_series_table("C. CUBOID / SPHERE / SIMPLEX-STYLE LIFTS", cube_rows)
    print_equivalence_clusters("cube8", cube_sequences)

    lrc_witness_table()
    print_method_synthesis(circle_rows, cube_rows)


if __name__ == "__main__":
    main()
