#!/usr/bin/env python3
"""
tournament_analysis_framework_s471.py

codex-2026-06-01 S471

Tournament Analysis: turn pairwise data into tournament-valued observables.

The general move is:

    objects + pairwise measurements + switch functional + tie Hamiltonian path
      -> a tournament
      -> tournament fingerprints
      -> trajectory as the measurement variable changes.

This script tries several concrete switch functionals:

* basketball pass majority with role-order tie-breaks;
* circle runner semicircle orientation;
* chord-threshold cuboid orientation;
* chord opening/closing orientation;
* circle and Fourier-lift isolation tournaments;
* simplex stress tournaments from an arbitrary symmetric metric;
* two-neighbor pressure tournaments.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import cos, exp, pi, sin, sqrt


EPS = 1e-10


@dataclass(frozen=True)
class TournamentSummary:
    label: str
    n: int
    hamiltonian_paths: int
    score_hist: str
    directed_triangles: int
    upsets: int
    top_scores: tuple[tuple[int, int], ...]
    bitstring: str


def base_labels(n: int) -> tuple[int, ...]:
    return tuple(range(n))


def tie_winner(i: int, j: int) -> int:
    """Base Hamiltonian path tie-break: lower label beats higher label."""

    return i if i < j else j


def tournament_from_pair_values(
    n: int, value
) -> list[list[bool]]:
    """Orient i->j when value(i,j)>0; ties use the base Hamiltonian path."""

    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = value(i, j)
        if delta > EPS:
            winner = i
        elif delta < -EPS:
            winner = j
        else:
            winner = tie_winner(i, j)
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def tournament_from_bits(
    n: int, edge_bit
) -> list[list[bool]]:
    """For i<j, bit 1 keeps base edge i->j; bit 0 flips it."""

    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if edge_bit(i, j):
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def score_hist(adj: list[list[bool]]) -> str:
    scores = [sum(row) for row in adj]
    pieces = []
    for value in sorted(set(scores)):
        count = scores.count(value)
        pieces.append(f"{value}^{count}" if count > 1 else str(value))
    return " ".join(pieces)


def top_scores(adj: list[list[bool]], limit: int = 5) -> tuple[tuple[int, int], ...]:
    scores = [(i, sum(adj[i])) for i in range(len(adj))]
    return tuple(sorted(scores, key=lambda item: (-item[1], item[0]))[:limit])


def directed_triangles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n + 1):
        for mask in range(1 << n):
            if mask.bit_count() != size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if prev & (1 << u) and adj[u][v]:
                        total += dp.get((prev, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def bitstring(adj: list[list[bool]]) -> str:
    bits = []
    for i, j in combinations(range(len(adj)), 2):
        bits.append("1" if adj[i][j] else "0")
    return "".join(bits)


def upsets_against_base(adj: list[list[bool]]) -> int:
    return sum(1 for i, j in combinations(range(len(adj)), 2) if adj[j][i])


def summarize(label: str, adj: list[list[bool]]) -> TournamentSummary:
    return TournamentSummary(
        label=label,
        n=len(adj),
        hamiltonian_paths=hamiltonian_paths(adj),
        score_hist=score_hist(adj),
        directed_triangles=directed_triangles(adj),
        upsets=upsets_against_base(adj),
        top_scores=top_scores(adj),
        bitstring=bitstring(adj),
    )


def print_summary_table(title: str, rows: list[TournamentSummary]) -> None:
    print(title)
    print("=" * 104)
    print(
        f"{'label':<34} {'n':>2} {'H':>8} {'3cyc':>5} {'up':>4} "
        f"{'score_hist':<16} top_scores"
    )
    print("-" * 104)
    for row in rows:
        print(
            f"{row.label:<34} {row.n:>2} {row.hamiltonian_paths:>8} "
            f"{row.directed_triangles:>5} {row.upsets:>4} "
            f"{row.score_hist:<16} {row.top_scores}"
        )
    print()


def basketball_pass_examples() -> list[TournamentSummary]:
    """Synthetic 5-starter passing matrices.  Labels 0..4 correspond to 1..5."""

    matrices: list[tuple[str, list[list[int]]]] = []

    matrices.append(
        (
            "passes: point-hub",
            [
                [0, 18, 15, 10, 8],
                [10, 0, 9, 7, 4],
                [8, 12, 0, 11, 6],
                [6, 9, 7, 0, 13],
                [5, 5, 5, 9, 0],
            ],
        )
    )
    matrices.append(
        (
            "passes: motion-ties",
            [
                [0, 12, 9, 8, 7],
                [12, 0, 11, 8, 6],
                [7, 11, 0, 10, 10],
                [8, 9, 10, 0, 12],
                [7, 6, 10, 12, 0],
            ],
        )
    )
    matrices.append(
        (
            "passes: inverted-big",
            [
                [0, 9, 7, 5, 11],
                [13, 0, 8, 6, 10],
                [10, 9, 0, 12, 14],
                [8, 7, 10, 0, 16],
                [17, 12, 13, 15, 0],
            ],
        )
    )

    out = []
    for label, mat in matrices:
        out.append(
            summarize(
                label,
                tournament_from_pair_values(
                    5, lambda i, j, mat=mat: mat[i][j] - mat[j][i]
                ),
            )
        )
    return out


def mod1(x: float) -> float:
    return x % 1.0


def circ_delta(a: float, b: float) -> float:
    return mod1(b - a)


def circ_dist(a: float, b: float) -> float:
    d = circ_delta(a, b)
    return min(d, 1.0 - d)


def chord(a: float, b: float) -> float:
    return 2.0 * sin(pi * circ_dist(a, b))


def positions(speeds: tuple[int, ...], t: Fraction) -> tuple[float, ...]:
    return tuple(mod1(float(t) * speed) for speed in speeds)


def chord_matrix(pos: tuple[float, ...]) -> list[list[float]]:
    n = len(pos)
    return [[0.0 if i == j else chord(pos[i], pos[j]) for j in range(n)] for i in range(n)]


def fourier_point(x: float, depth: int = 2) -> tuple[float, ...]:
    coords: list[float] = []
    for k in range(1, depth + 1):
        coords.append(cos(2 * pi * k * x) / sqrt(depth))
        coords.append(sin(2 * pi * k * x) / sqrt(depth))
    return tuple(coords)


def fourier_velocity(speed: int, x: float, depth: int = 2) -> tuple[float, ...]:
    coords: list[float] = []
    for k in range(1, depth + 1):
        scale = 2 * pi * k * speed / sqrt(depth)
        coords.append(-scale * sin(2 * pi * k * x))
        coords.append(scale * cos(2 * pi * k * x))
    return tuple(coords)


def euclidean(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def runner_tournaments(
    label: str, speeds: tuple[int, ...], t: Fraction
) -> list[TournamentSummary]:
    pos = positions(speeds, t)
    n = len(speeds)
    chords = chord_matrix(pos)
    threshold = 2.0 * sin(pi / n)

    rows: list[TournamentSummary] = []

    rows.append(
        summarize(
            f"{label}: semicircle",
            tournament_from_pair_values(
                n,
                lambda i, j: (
                    1.0
                    if 0.0 < circ_delta(pos[i], pos[j]) < 0.5
                    else -1.0
                    if 0.5 < circ_delta(pos[i], pos[j]) < 1.0
                    else 0.0
                ),
            ),
        )
    )

    rows.append(
        summarize(
            f"{label}: chord-threshold",
            tournament_from_bits(n, lambda i, j: chords[i][j] >= threshold),
        )
    )

    def opening_value(i: int, j: int) -> float:
        delta = circ_delta(pos[i], pos[j])
        raw = (speeds[j] - speeds[i]) * sin(2 * pi * delta)
        return (speeds[i] - speeds[j]) * raw

    rows.append(
        summarize(
            f"{label}: chord-opening",
            tournament_from_pair_values(n, opening_value),
        )
    )

    isolation = [sum(chords[i]) for i in range(n)]
    rows.append(
        summarize(
            f"{label}: chord-isolation",
            tournament_from_pair_values(n, lambda i, j: isolation[i] - isolation[j]),
        )
    )

    heat = [sum(exp(-3.0 * chords[i][j] ** 2) for j in range(n) if j != i) for i in range(n)]
    rows.append(
        summarize(
            f"{label}: heat-escape",
            tournament_from_pair_values(n, lambda i, j: heat[j] - heat[i]),
        )
    )

    fpoints = [fourier_point(x, 2) for x in pos]
    fvels = [fourier_velocity(speeds[i], pos[i], 2) for i in range(n)]
    fdist = [[0.0 if i == j else euclidean(fpoints[i], fpoints[j]) for j in range(n)] for i in range(n)]
    fstress = [sum(fdist[i][k] for k in range(n) if k != i) for i in range(n)]
    rows.append(
        summarize(
            f"{label}: fourier-isolation",
            tournament_from_pair_values(n, lambda i, j: fstress[i] - fstress[j]),
        )
    )

    def fourier_opening(i: int, j: int) -> float:
        diff = [fpoints[i][k] - fpoints[j][k] for k in range(len(fpoints[i]))]
        vdiff = [fvels[i][k] - fvels[j][k] for k in range(len(fvels[i]))]
        raw = 2.0 * sum(a * b for a, b in zip(diff, vdiff))
        return (speeds[i] - speeds[j]) * raw

    rows.append(
        summarize(
            f"{label}: fourier-opening",
            tournament_from_pair_values(n, fourier_opening),
        )
    )

    def simplex_stress(i: int, j: int) -> float:
        left = sum(chords[i][k] for k in range(n) if k not in (i, j))
        right = sum(chords[j][k] for k in range(n) if k not in (i, j))
        return left - right

    rows.append(
        summarize(
            f"{label}: simplex-stress",
            tournament_from_pair_values(n, simplex_stress),
        )
    )

    def nearest_without(i: int, removed: int) -> float:
        vals = [circ_dist(pos[i], pos[k]) for k in range(n) if k not in (i, removed)]
        return min(vals) if vals else 0.0

    nearest = [
        min(circ_dist(pos[i], pos[k]) for k in range(n) if k != i)
        for i in range(n)
    ]

    def pressure_value(i: int, j: int) -> float:
        relief_i_from_j = nearest_without(i, j) - nearest[i]
        relief_j_from_i = nearest_without(j, i) - nearest[j]
        # j blocks i more -> j beats i, hence i loses.
        return relief_j_from_i - relief_i_from_j

    rows.append(
        summarize(
            f"{label}: pressure",
            tournament_from_pair_values(n, pressure_value),
        )
    )

    return rows


def trajectory_report() -> None:
    speeds = tuple(range(5))
    n = len(speeds)
    threshold = 2.0 * sin(pi / n)
    seen: Counter[str] = Counter()
    h_seen: Counter[int] = Counter()
    flips: Counter[tuple[int, int]] = Counter()
    previous: str | None = None
    previous_adj: list[list[bool]] | None = None
    for q in range(1, 61):
        t = Fraction(q, 60)
        pos = positions(speeds, t)
        chords = chord_matrix(pos)
        adj = tournament_from_bits(n, lambda i, j, chords=chords: chords[i][j] >= threshold)
        bits = bitstring(adj)
        seen[bits] += 1
        h_seen[hamiltonian_paths(adj)] += 1
        if previous is not None and previous_adj is not None:
            for i, j in combinations(range(n), 2):
                if adj[i][j] != previous_adj[i][j]:
                    flips[(i, j)] += 1
        previous = bits
        previous_adj = adj

    print("Chord-threshold cuboid trajectory: initial five runners, t=q/60")
    print("=" * 104)
    print(f"unique tournaments={len(seen)} over 60 samples")
    print(f"H spectrum={dict(sorted(h_seen.items()))}")
    print(f"most common bitstrings={seen.most_common(5)}")
    print(f"edge flip counts={dict(sorted(flips.items()))}")
    print()


def print_rule_atlas() -> None:
    print("Tournament Analysis schema")
    print("=" * 104)
    print("objects + pairwise observable + switch functional + tie path -> tournament trajectory")
    print("Tie path used here is the base Hamiltonian path 0->1->...->n-1.")
    print()
    print("Switch families tried")
    print("-" * 104)
    families = [
        ("pass majority", "asymmetric count", "i beats j when passes i->j exceed j->i"),
        ("semicircle", "positions on circle/sphere", "i beats j when j lies in i's forward open half"),
        ("chord threshold", "symmetric chord metric", "edge bit switches base path when chord crosses threshold"),
        ("chord opening", "metric derivative", "faster endpoint wins when chord is opening, slower when closing"),
        ("isolation/heat", "metric centrality", "less crowded or more isolated vertex wins"),
        ("simplex stress", "edge-weighted simplex", "compare incident stress after deleting the pair edge"),
        ("pressure", "first/second nearest", "irreplaceable blocker beats the runner it blocks"),
    ]
    for name, data, rule in families:
        print(f"{name:<18} {data:<28} {rule}")
    print()


def main() -> None:
    print("codex-2026-06-01-S471")
    print("Tournament Analysis framework: pairwise data -> binary relation -> fingerprints")
    print()
    print_rule_atlas()
    print_summary_table("Basketball pass tournaments with role-order tie-breaks", basketball_pass_examples())

    examples = [
        ("initial5 boundary", tuple(range(5)), Fraction(1, 5)),
        ("initial6 boundary", tuple(range(6)), Fraction(1, 6)),
        ("initial14 boundary", tuple(range(14)), Fraction(1, 14)),
        (
            "n14 seven boundary",
            (0, 1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
            Fraction(9, 98),
        ),
    ]
    for label, speeds, t in examples:
        print_summary_table(
            f"Runner Tournament Analysis at {label}, t={t}",
            runner_tournaments(label, speeds, t),
        )

    trajectory_report()

    print("Synthesis")
    print("=" * 104)
    print(
        "A metric alone gives an edge-weighted complete graph.  Tournament Analysis "
        "starts when we choose a switch functional for every edge."
    )
    print(
        "Symmetric metrics are not a problem: threshold them against a base "
        "Hamiltonian path and the edge set becomes a cuboid of tournament bits."
    )
    print(
        "As variables move, we should study the trajectory in tournament space: "
        "which edges flip, which H-values appear, where 3-cycle curvature forms, "
        "and whether pressure or handoff SCCs survive peeling."
    )
    print(
        "The basketball role-order tie-break and the LRC base-path/staircase "
        "encoding are the same abstraction: a convenient labelled Hamiltonian "
        "path turns ties or thresholds into honest tournaments."
    )


if __name__ == "__main__":
    main()
