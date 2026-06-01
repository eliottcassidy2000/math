#!/usr/bin/env python3
"""
tournament_analysis_metric_lifts_s480.py

codex-2026-05-31 S480

Tournament Analysis:

    pairwise or geometric data + a binary decision gauge
        -> tournament-valued structure
        -> tournament invariants as the original variables move.

The user suggested basketball passing, runner distances on a circle, chord
metrics, and higher-dimensional geometric metrics.  This script makes those
examples explicit and compares several gauges that turn continuous or weighted
data into complete binary relations.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import cos, pi, sin, sqrt


EPS = 1e-10


@dataclass(frozen=True)
class Tournament:
    labels: tuple[str, ...]
    arcs: frozenset[tuple[int, int]]

    @property
    def n(self) -> int:
        return len(self.labels)

    def beats(self, i: int, j: int) -> bool:
        return (i, j) in self.arcs

    def score_sequence(self) -> tuple[int, ...]:
        scores = [0] * self.n
        for i, _j in self.arcs:
            scores[i] += 1
        return tuple(sorted(scores))

    def score_hist(self) -> tuple[tuple[int, int], ...]:
        return tuple(sorted(Counter(self.score_sequence()).items()))

    def cyclic_triples(self) -> int:
        out = 0
        for a, b, c in combinations(range(self.n), 3):
            s = sum(
                1
                for x, y in ((a, b), (a, c), (b, c))
                if self.beats(x, y)
            )
            if s in (1, 2):
                # The transitive triples have a vertex with outdegree 2 and a
                # vertex with outdegree 0.  Cyclic triples have all local
                # outdegrees equal to 1.
                local = Counter()
                for x, y in ((a, b), (a, c), (b, c)):
                    local[x if self.beats(x, y) else y] += 1
                if sorted(local.values()) == [1, 1, 1]:
                    out += 1
        return out

    def scc_count(self) -> int:
        adj = [[] for _ in range(self.n)]
        radj = [[] for _ in range(self.n)]
        for a, b in self.arcs:
            adj[a].append(b)
            radj[b].append(a)

        seen: set[int] = set()
        order: list[int] = []

        def dfs(v: int) -> None:
            seen.add(v)
            for w in adj[v]:
                if w not in seen:
                    dfs(w)
            order.append(v)

        for v in range(self.n):
            if v not in seen:
                dfs(v)

        seen.clear()
        components = 0

        def rdfs(v: int) -> None:
            seen.add(v)
            for w in radj[v]:
                if w not in seen:
                    rdfs(w)

        for v in reversed(order):
            if v not in seen:
                components += 1
                rdfs(v)
        return components

    def hamiltonian_paths(self) -> int:
        n = self.n
        dp = [[0] * n for _ in range(1 << n)]
        for i in range(n):
            dp[1 << i][i] = 1
        for mask in range(1 << n):
            for last in range(n):
                total = dp[mask][last]
                if not total:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    if self.beats(last, nxt):
                        dp[mask | (1 << nxt)][nxt] += total
        return sum(dp[(1 << n) - 1])

    def signature(self) -> str:
        bits = []
        for i, j in combinations(range(self.n), 2):
            bits.append("1" if self.beats(i, j) else "0")
        return "".join(bits)

    def summary(self) -> str:
        return (
            f"score={self.score_sequence()} "
            f"scc={self.scc_count()} "
            f"c3={self.cyclic_triples()} "
            f"H={self.hamiltonian_paths()}"
        )


def lower_label_tiebreak(i: int, j: int) -> int:
    return i if i < j else j


def tournament_from_winner(
    labels: tuple[str, ...], winner
) -> Tournament:
    arcs = set()
    for i, j in combinations(range(len(labels)), 2):
        w = winner(i, j)
        arcs.add((w, j if w == i else i))
    return Tournament(labels, frozenset(arcs))


def sign_tiebreak(value: float, i: int, j: int) -> int:
    if value > EPS:
        return i
    if value < -EPS:
        return j
    return lower_label_tiebreak(i, j)


def distance(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sum(x * y for x, y in zip(a, b))


def normalize(v: tuple[float, ...]) -> tuple[float, ...]:
    norm = sqrt(sum(x * x for x in v))
    return tuple(x / norm for x in v)


def circle_point(speed: int, t: Fraction) -> tuple[float, float]:
    angle = 2 * pi * float((speed * t) % 1)
    return (cos(angle), sin(angle))


def circle_velocity(speed: int, t: Fraction) -> tuple[float, float]:
    angle = 2 * pi * float((speed * t) % 1)
    return (-2 * pi * speed * sin(angle), 2 * pi * speed * cos(angle))


def clockwise_delta(a: Fraction, b: Fraction) -> Fraction:
    return (b - a) % 1


def circle_distance_fraction(a: Fraction, b: Fraction) -> Fraction:
    d = clockwise_delta(a, b)
    return min(d, 1 - d)


def positions_for_speeds(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple((v * t) % 1 for v in speeds)


def chord_metric(positions: tuple[tuple[float, float], ...]) -> list[list[float]]:
    return [[distance(a, b) for b in positions] for a in positions]


def tournament_from_vertex_scores(
    labels: tuple[str, ...], scores: tuple[float, ...], high_wins: bool = True
) -> Tournament:
    def winner(i: int, j: int) -> int:
        delta = scores[i] - scores[j]
        return sign_tiebreak(delta if high_wins else -delta, i, j)

    return tournament_from_winner(labels, winner)


def metric_row_sum_tournament(
    labels: tuple[str, ...], metric: list[list[float]]
) -> Tournament:
    scores = tuple(sum(row) for row in metric)
    return tournament_from_vertex_scores(labels, scores, high_wins=True)


def nearest_radius_tournament(
    labels: tuple[str, ...], metric: list[list[float]]
) -> Tournament:
    radii = []
    for i in range(len(labels)):
        radii.append(min(metric[i][j] for j in range(len(labels)) if i != j))
    return tournament_from_vertex_scores(labels, tuple(radii), high_wins=True)


def metric_threshold_switch_tournament(
    labels: tuple[str, ...], metric: list[list[float]], threshold: float
) -> Tournament:
    def winner(i: int, j: int) -> int:
        if metric[i][j] > threshold + EPS:
            return j
        return i

    return tournament_from_winner(labels, winner)


def metric_annulus_switch_tournament(
    labels: tuple[str, ...], metric: list[list[float]], low: float, high: float
) -> Tournament:
    def winner(i: int, j: int) -> int:
        if low + EPS < metric[i][j] < high - EPS:
            return j
        return i

    return tournament_from_winner(labels, winner)


def projection_tournament(
    labels: tuple[str, ...],
    points: tuple[tuple[float, ...], ...],
    axis: tuple[float, ...],
) -> Tournament:
    axis = normalize(axis)
    scores = tuple(dot(point, axis) for point in points)
    return tournament_from_vertex_scores(labels, scores, high_wins=True)


def circular_half_tournament(
    labels: tuple[str, ...], positions: tuple[Fraction, ...]
) -> Tournament:
    def winner(i: int, j: int) -> int:
        d = clockwise_delta(positions[i], positions[j])
        if d == 0 or d == Fraction(1, 2):
            return lower_label_tiebreak(i, j)
        return i if d < Fraction(1, 2) else j

    return tournament_from_winner(labels, winner)


def zero_distance_tournament(
    labels: tuple[str, ...], positions: tuple[Fraction, ...]
) -> Tournament:
    scores = tuple(float(circle_distance_fraction(Fraction(0), p)) for p in positions)
    return tournament_from_vertex_scores(labels, scores, high_wins=True)


def pull_apart_tournament(
    labels: tuple[str, ...], speeds: tuple[int, ...], t: Fraction
) -> Tournament:
    points = tuple(circle_point(v, t) for v in speeds)
    velocities = tuple(circle_velocity(v, t) for v in speeds)

    def winner(i: int, j: int) -> int:
        midpoint = tuple((points[i][k] + points[j][k]) / 2 for k in range(2))
        out_i = tuple(points[i][k] - midpoint[k] for k in range(2))
        out_j = tuple(points[j][k] - midpoint[k] for k in range(2))
        score_i = dot(velocities[i], out_i)
        score_j = dot(velocities[j], out_j)
        return sign_tiebreak(score_i - score_j, i, j)

    return tournament_from_winner(labels, winner)


def basketball_tournament() -> tuple[Tournament, int]:
    labels = ("1", "2", "3", "4", "5")
    # passes[i][j] = passes from player i+1 to player j+1.
    passes = [
        [0, 18, 9, 7, 12],
        [12, 0, 16, 6, 10],
        [9, 8, 0, 11, 4],
        [5, 6, 13, 0, 15],
        [14, 4, 10, 15, 0],
    ]
    ties = 0

    def winner(i: int, j: int) -> int:
        nonlocal ties
        delta = passes[i][j] - passes[j][i]
        if delta == 0:
            ties += 1
        return sign_tiebreak(float(delta), i, j)

    return tournament_from_winner(labels, winner), ties


@dataclass(frozen=True)
class MovieSummary:
    gauge: str
    distinct: int
    changes: int
    h_range: tuple[int, int]
    c3_range: tuple[int, int]
    scc_values: tuple[int, ...]
    score_patterns: int
    first_signature: str
    last_signature: str


def summarize_movie(
    gauge: str, tournaments: tuple[Tournament, ...]
) -> MovieSummary:
    signatures = [t.signature() for t in tournaments]
    h_values = [t.hamiltonian_paths() for t in tournaments]
    c3_values = [t.cyclic_triples() for t in tournaments]
    return MovieSummary(
        gauge=gauge,
        distinct=len(set(signatures)),
        changes=sum(a != b for a, b in zip(signatures, signatures[1:])),
        h_range=(min(h_values), max(h_values)),
        c3_range=(min(c3_values), max(c3_values)),
        scc_values=tuple(sorted({t.scc_count() for t in tournaments})),
        score_patterns=len({t.score_sequence() for t in tournaments}),
        first_signature=signatures[0],
        last_signature=signatures[-1],
    )


def runner_tournament(
    gauge: str, speeds: tuple[int, ...], t: Fraction
) -> Tournament:
    labels = tuple(str(v) for v in speeds)
    positions_f = positions_for_speeds(speeds, t)
    points = tuple(circle_point(v, t) for v in speeds)
    metric = chord_metric(points)

    if gauge == "round-half":
        return circular_half_tournament(labels, positions_f)
    if gauge == "zero-distance":
        return zero_distance_tournament(labels, positions_f)
    if gauge == "chord-row-sum":
        return metric_row_sum_tournament(labels, metric)
    if gauge == "nearest-radius":
        return nearest_radius_tournament(labels, metric)
    if gauge == "chord-threshold":
        return metric_threshold_switch_tournament(labels, metric, sqrt(2))
    if gauge == "chord-annulus":
        return metric_annulus_switch_tournament(labels, metric, 0.75, 1.55)
    if gauge == "x-projection":
        return projection_tournament(labels, points, (1.0, 0.0))
    if gauge == "pull-apart":
        return pull_apart_tournament(labels, speeds, t)
    raise ValueError(gauge)


def runner_movie() -> tuple[MovieSummary, ...]:
    speeds = (0, 1, 4, 5, 6, 7, 11)
    times = tuple(Fraction(k, 84) for k in range(1, 84))
    gauges = (
        "round-half",
        "zero-distance",
        "chord-row-sum",
        "nearest-radius",
        "chord-threshold",
        "chord-annulus",
        "x-projection",
        "pull-apart",
    )
    return tuple(
        summarize_movie(gauge, tuple(runner_tournament(gauge, speeds, t) for t in times))
        for gauge in gauges
    )


def geometric_cloud_tournaments() -> tuple[tuple[str, Tournament], ...]:
    labels = ("A", "B", "C", "D", "E", "F")
    raw = (
        (1.0, 0.0, 0.0),
        (0.1, 0.95, 0.2),
        (-0.7, 0.35, 0.6),
        (-0.8, -0.45, 0.4),
        (0.25, -0.85, -0.45),
        (0.65, 0.1, -0.75),
    )
    points = tuple(normalize(p) for p in raw)
    metric = [[distance(a, b) for b in points] for a in points]
    centroid = tuple(sum(p[k] for p in points) / len(points) for k in range(3))
    dist_to_centroid = tuple(distance(p, centroid) for p in points)
    simplex_vertex = (1.0, 0.0, 0.0)
    close_to_vertex = tuple(distance(p, simplex_vertex) for p in points)
    pair_distances = sorted(metric[i][j] for i, j in combinations(range(len(labels)), 2))
    median_distance = pair_distances[len(pair_distances) // 2]
    low_annulus = pair_distances[len(pair_distances) // 4]
    high_annulus = pair_distances[(3 * len(pair_distances)) // 4]

    return (
        ("3d projection axis", projection_tournament(labels, points, (1.0, 2.0, -1.0))),
        ("centroid outside", tournament_from_vertex_scores(labels, dist_to_centroid, True)),
        ("simplex vertex affinity", tournament_from_vertex_scores(labels, close_to_vertex, False)),
        ("metric row-sum isolation", metric_row_sum_tournament(labels, metric)),
        ("nearest-radius isolation", nearest_radius_tournament(labels, metric)),
        (
            "median-distance switch",
            metric_threshold_switch_tournament(labels, metric, median_distance),
        ),
        (
            "middle-annulus switch",
            metric_annulus_switch_tournament(labels, metric, low_annulus, high_annulus),
        ),
    )


def print_pipeline() -> None:
    print("TOURNAMENT ANALYSIS PIPELINE")
    print("=" * 96)
    rows = [
        (
            "raw pairwise direction",
            "pass_i_j - pass_j_i",
            "sports, votes, traffic, trades",
            "already antisymmetric; ties need a path",
        ),
        (
            "marked reference",
            "d(i, base) - d(j, base)",
            "LRC stationary runner, hub audits",
            "usually transitive but tracks threat order",
        ),
        (
            "circular orientation",
            "clockwise displacement < 1/2",
            "runners, rotations, phases",
            "round tournament with tie defects",
        ),
        (
            "global metric row",
            "sum_k d(i,k) - sum_k d(j,k)",
            "isolation/spread/spacing",
            "turns symmetric metric into vertex scores",
        ),
        (
            "local metric row",
            "nearest_radius(i)-nearest_radius(j)",
            "crowding, local density",
            "detects private leaves and packed clusters",
        ),
        (
            "metric switch",
            "1[d(i,j)>theta] flips base edge",
            "pairwise distances, chords, polls",
            "the N choose 2 metrics become a tiling coordinate",
        ),
        (
            "metric annulus",
            "1[lo<d(i,j)<hi] flips base edge",
            "middle-distance bonds, shell tests",
            "switches edges by belonging to a chosen shell",
        ),
        (
            "projection gauge",
            "<x_i,u>-<x_j,u>",
            "sphere/cube/simplex shadows",
            "sweeping u gives a tournament movie",
        ),
        (
            "kinetic gauge",
            "outward velocity from pair midpoint",
            "moving runners, fluids, teams in motion",
            "uses derivative of metric, not only metric",
        ),
        (
            "tie Hamiltonian path",
            "fixed label order 1->2->...->n",
            "basketball jerseys, base-path tilings",
            "same convenience as the repo tiling model",
        ),
    ]
    print(f"{'gauge':<24} {'binary function':<36} {'example'}")
    print("-" * 96)
    for gauge, fn, example, warning in rows:
        print(f"{gauge:<24} {fn:<36} {example}")
        print(f"{'':<24} {'':<36} {warning}")
    print()


def print_basketball() -> None:
    print("BASKETBALL PASS-FLOW TOURNAMENT")
    print("=" * 96)
    tournament, ties = basketball_tournament()
    print("players = 1,2,3,4,5; ties broken by Hamiltonian path 1->2->3->4->5")
    print(f"tied pair count={ties}")
    print(f"summary: {tournament.summary()}")
    print(f"signature={tournament.signature()}")
    print(
        "reading: the pass matrix is a weighted directed relation; Tournament "
        "Analysis keeps the binary skeleton first, then asks which weighted "
        "margins explain its cycles."
    )
    print()


def print_runner_movie() -> None:
    print("RUNNER MOVIE: ONE SPEED SET, MANY TOURNAMENT GAUGES")
    print("=" * 96)
    print("speeds=(0,1,4,5,6,7,11), sampled at t=k/84 for k=1..83")
    print(
        f"{'gauge':<18} {'distinct':>8} {'changes':>8} {'H-range':>16} "
        f"{'c3-range':>12} {'scc':>10} {'score-patterns':>15}"
    )
    print("-" * 96)
    for row in runner_movie():
        print(
            f"{row.gauge:<18} {row.distinct:>8} {row.changes:>8} "
            f"{str(row.h_range):>16} {str(row.c3_range):>12} "
            f"{str(row.scc_values):>10} {row.score_patterns:>15}"
        )
    print()

    print("selected time t=181/588 from the S431 near-tight k=6 example")
    speeds = (0, 1, 4, 5, 6, 7, 11)
    t = Fraction(181, 588)
    for gauge in (
        "round-half",
        "zero-distance",
        "chord-row-sum",
        "nearest-radius",
        "chord-threshold",
        "chord-annulus",
        "x-projection",
        "pull-apart",
    ):
        tour = runner_tournament(gauge, speeds, t)
        print(f"  {gauge:<18} {tour.summary()}")
    print()


def print_geometric_cloud() -> None:
    print("SPHERE/CUBE/SIMPLEX POINT CLOUD")
    print("=" * 96)
    print("same six points; seven gauges")
    for name, tour in geometric_cloud_tournaments():
        print(f"  {name:<25} {tour.summary()} signature={tour.signature()}")
    print()


def print_patterns() -> None:
    print("PATTERNS / QUESTIONS")
    print("=" * 96)
    patterns = [
        (
            "Gauge non-uniqueness is the point",
            "A symmetric metric does not become a tournament until a gauge chooses what a pairwise comparison means.",
        ),
        (
            "Ties are structure, not noise",
            "Basketball equal pass counts and runner antipodes are both incomplete relations completed by a base Hamiltonian path.",
        ),
        (
            "Movies matter more than frames",
            "For continuous variables, the flip events and invariant ranges are often more meaningful than a single tournament.",
        ),
        (
            "LRC is marked Tournament Analysis",
            "The stationary runner supplies the mark; S22 keeps its two bracket gaps, while S431 keeps the surrounding pairwise tournament shadow.",
        ),
        (
            "Higher geometry supplies gauges",
            "Spheres, cubes, and simplices naturally give projection, centroid, face-affinity, and isolation tournaments.",
        ),
        (
            "Repo bridge",
            "The fixed tie path is the same mental move as the tiling model's fixed Hamiltonian path: it turns ambiguity into a controlled coordinate system.",
        ),
    ]
    for title, body in patterns:
        print(f"- {title}: {body}")


def main() -> None:
    print("Tournament Analysis metric lifts (codex-2026-05-31 S480)")
    print("Pairwise data becomes tournament structure only after choosing a gauge.\n")
    print_pipeline()
    print_basketball()
    print_runner_movie()
    print_geometric_cloud()
    print_patterns()


if __name__ == "__main__":
    main()
